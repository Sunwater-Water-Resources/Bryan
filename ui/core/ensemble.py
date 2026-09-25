"""Ensemble results - the PMF - and a notional AEP for it.

**Which event.** An ensemble run routes every temporal pattern of every duration
at one AEP, so a result has to be picked from the spread. Two conventions, and a
study uses both:

``highest``  The largest level of any event, with *that event's* inflow and
             outflow. What the PMF is taken as.
``median``   For each duration, the median pattern - ``lib/EnbAnalysis.py``'s own
             pick, the value at position ``round(n / 2)`` of the ascending sort,
             with numpy's rounding (the sixth of ten) - then the duration whose
             median is largest. What everything else is taken as, and what Bryan
             writes to ``csv/<name>_critical.csv``. Here the inflow and outflow
             are the same event's, where EnbAnalysis takes each result's median
             separately; the level is the same either way.

At Callide the two differ by 0.14-0.28 m of PMF level (E012, RFSL).

**A notional AEP.** The PMF has no AEP of its own, but a risk assessment needs
one, and the Monte Carlo realisations can give it: they extend past the AEP of
the PMP as far as the storm config extrapolates the rainfall (a GEV at
Callide), so the PMF level usually sits *inside* the sample rather than beyond
it - at Callide E012 the realisations reach 1 in 38.8 million and 221.17 m,
above the 221.05 m PMF. The estimate reads the PMF level against the realisations
of one duration (by default the level's critical duration at the AEP of the
PMP) with a polynomial in log10(level) fitted to the standard normal variate over
an AEP window at the top of the sample - ``AEPofPMF.py``'s method, made
explicit. A straight read off the realisations without the fit is too noisy up
there (a dozen or so events within 50 mm of the PMF level).

The realisations are those of **the PMF event's own duration** by default (9 h
at Callide GWL 0, where the level's Monte Carlo critical duration at the AEP of
the PMP is 12 h) - the choice ``AEPofPMF.py`` made by hand. The window's upper
bound defaults to the top of the sample: the script's 1 in 4,000,000 left the
PMF above every realisation in its window, so it was extrapolating off the end
of its own fit (1 in 7.8 million at degree 1, against 9.0 million with the
whole top of the sample).

The window's lower bound and the degree change the answer, so the page shows a
grid of both. That grid is the sensitivity of the **fit**, not the whole
uncertainty: the shape of the tail is set by the rainfall extrapolation and by
the rating and storage curves, which no fit setting touches.

numpy, pandas and ``statistics.NormalDist`` only - no scipy - so it runs in the
UI's own environment. Nothing here imports nicegui.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

from . import results
from .outputs import find_database
from .paths import cell_text
from .study import Study, StudyError

HIGHEST = "highest"
MEDIAN = "median"
RESULT_TYPES = ("level", "inflow", "outflow")

DEFAULT_PMP_AEP = 1_900_000
DEFAULT_LOWER_AEP = 500_000
DEGREES = (1, 2, 3)
# Fewer realisations than this in the window and there is no fit worth quoting.
MIN_POINTS = 30
THIN_POINTS = 100

_NORMAL = NormalDist()


# -- the ensemble database -----------------------------------------------------

def load(study: Study, run_name: str, group: str) -> pd.DataFrame:
    """A group's ensemble database: one row per event, rows of every member."""
    run = study.open_run(run_name)
    rows = run.rows_in(group)
    if not rows:
        raise StudyError(f"{run_name} has no group {group!r}")
    frames = []
    for index in rows:
        path = find_database(run.frame.loc[index], run.project_folder)
        if path is None:
            continue
        try:
            frame = (pd.read_parquet(path) if path.suffix == ".parquet"
                     else pd.read_csv(path, index_col=0))
        except Exception as exc:                 # noqa: BLE001 - report, never crash
            raise StudyError(f"{path.name} could not be read ({exc})") from exc
        frames.append(frame)
    if not frames:
        raise StudyError(f"{group}: no results database on disk")
    frame = pd.concat(frames, ignore_index=True)
    missing = [name for name in ("duration", "level", "inflow", "outflow")
               if name not in frame.columns]
    if missing:
        raise StudyError(f"{group}: the database has no {', '.join(missing)} column - "
                         f"is this an ensemble run?")
    for name in ("duration", "level", "inflow", "outflow"):
        frame[name] = pd.to_numeric(frame[name], errors="coerce")
    return frame


def pattern_label(row) -> str:
    """'GSDM: 3', as EnbAnalysis labels a pattern."""
    method = cell_text(row.get("storm_method")) if hasattr(row, "get") else ""
    tp = row.get("tp") if hasattr(row, "get") else None
    tp_text = "" if tp is None or (isinstance(tp, float) and math.isnan(tp)) else (
        f"{int(tp)}" if float(tp).is_integer() else f"{tp}")
    return f"{method}: {tp_text}" if method else tp_text


def median_position(count: int) -> int:
    """lib/EnbAnalysis.py: ``int(np.around(n / 2, 0))`` of the ascending sort."""
    return int(np.around(count / 2, 0))


@dataclass
class Pick:
    convention: str
    event: object = None           # the database index of the event
    duration: float = math.nan
    level: float = math.nan
    inflow: float = math.nan
    outflow: float = math.nan
    pattern: str = ""

    @property
    def found(self) -> bool:
        return self.event is not None


def _pick_event(frame, event, convention) -> Pick:
    row = frame.loc[event]
    return Pick(convention=convention, event=event, duration=float(row["duration"]),
                level=float(row["level"]), inflow=float(row["inflow"]),
                outflow=float(row["outflow"]), pattern=pattern_label(row))


def median_events(frame: pd.DataFrame, result: str = "level") -> pd.Series:
    """Each duration's median-pattern event, by EnbAnalysis's position rule."""
    picked = {}
    for duration, local in frame.dropna(subset=[result]).groupby("duration"):
        ordered = local.sort_values(result, ascending=True, kind="mergesort")
        picked[duration] = ordered.index[median_position(len(ordered))]
    return pd.Series(picked, dtype=object)


def pick(frame: pd.DataFrame, convention: str = HIGHEST, result: str = "level") -> Pick:
    """The event a table quotes, with its own inflow and outflow.

    Picked on ``result`` - the level for the PMF; the Ensemble page picks on
    whichever result it is showing.
    """
    values = frame[result].dropna()
    if values.empty:
        return Pick(convention=convention)
    if convention == MEDIAN:
        events = median_events(frame, result)
        best = max(events.index, key=lambda duration: frame.loc[events[duration], result])
        return _pick_event(frame, events[best], convention)
    return _pick_event(frame, values.idxmax(), convention)


def by_duration(frame: pd.DataFrame, result: str = "level") -> pd.DataFrame:
    """Per duration: the spread of the patterns, the median pick and the highest."""
    rows = []
    for duration, local in frame.dropna(subset=[result]).groupby("duration"):
        values = local[result].sort_values(kind="mergesort")
        median_event = values.index[median_position(len(values))]
        top = values.idxmax()
        rows.append({"duration": float(duration), "n": len(values),
                     "min": float(values.min()),
                     "q1": float(values.quantile(0.25)),
                     "median": float(values.loc[median_event]),
                     "median_pattern": pattern_label(frame.loc[median_event]),
                     "q3": float(values.quantile(0.75)),
                     "max": float(values.max()),
                     "max_pattern": pattern_label(frame.loc[top])})
    return pd.DataFrame(rows).set_index("duration") if rows else pd.DataFrame()


def pattern_points(frame: pd.DataFrame, result: str = "level") -> list:
    """(duration, value, pattern) for every event - the dots over the boxes."""
    return [(float(row["duration"]), float(row[result]), pattern_label(row))
            for _, row in frame.dropna(subset=[result]).iterrows()]


# -- the Monte Carlo realisations ----------------------------------------------

def variate(probability) -> float:
    """z for an exceedance probability, from the lower tail so it keeps its digits.

    ``ndtri(1 - p)`` loses the precision of p near 1e-8 in the subtraction;
    ``-inv_cdf(p)`` does not.
    """
    try:
        p = float(probability)
    except (TypeError, ValueError):
        return math.nan
    if not 0.0 < p < 1.0:
        return math.nan
    return -_NORMAL.inv_cdf(p)


def aep_of_variate(z: float) -> float:
    """'1 in X' for a standard normal variate, from the upper tail directly."""
    if not math.isfinite(z):
        return math.nan
    tail = 0.5 * math.erfc(z / math.sqrt(2.0))
    return 1.0 / tail if tail > 0 else math.inf


@dataclass
class Realisations:
    """One duration's Monte Carlo events, as (z, value) of one result."""

    path: Path
    result: str
    z: np.ndarray
    value: np.ndarray

    @property
    def top_aep(self) -> float:
        return aep_of_variate(float(self.z.max())) if len(self.z) else math.nan

    @property
    def top_value(self) -> float:
        return float(self.value.max()) if len(self.value) else math.nan


_REALISATION_CACHE: dict = {}


def read_realisations(path: Path, result: str = "level") -> Realisations:
    """The ``<result>`` and ``<result>_aep`` columns of an mcdf.

    The mcdf stores the AEP as a **probability** (``TotalProbTheorem.assign_aep``),
    not a '1 in X' - see lib/RepresentativeEvents.py.
    """
    path = Path(path)
    stamp = (str(path), path.stat().st_mtime_ns, result)
    if stamp in _REALISATION_CACHE:
        return _REALISATION_CACHE[stamp]
    column = f"{result}_aep"
    try:
        frame = pd.read_csv(path, usecols=[result, column])
    except ValueError as exc:
        raise StudyError(f"{path.name} has no {result!r} and {column!r} columns - "
                         f"was it analysed?") from exc
    z = np.array([variate(p) for p in frame[column]], dtype=float)
    value = pd.to_numeric(frame[result], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(z) & np.isfinite(value) & (value > 0)
    order = np.argsort(z[keep], kind="mergesort")
    out = Realisations(path=path, result=result, z=z[keep][order], value=value[keep][order])
    if len(_REALISATION_CACHE) > 12:        # an mcdf is megabytes; keep a handful
        _REALISATION_CACHE.clear()
    _REALISATION_CACHE[stamp] = out
    return out


# -- the fit -------------------------------------------------------------------

@dataclass
class Estimate:
    target: float
    lower_aep: float
    upper_aep: float
    degree: int
    aep: float = math.nan
    z: float = math.nan
    n: int = 0
    coefficients: list = field(default_factory=list)   # z = poly(log10 value)
    warnings: list = field(default_factory=list)
    refused: str = ""

    @property
    def ok(self) -> bool:
        return not self.refused and math.isfinite(self.aep)


def fit(real: Realisations, target: float, lower_aep: float = DEFAULT_LOWER_AEP,
        upper_aep: float | None = None, degree: int = 1) -> Estimate:
    """The AEP at which the fitted top of the sample reaches ``target``."""
    upper = float(upper_aep) if upper_aep else real.top_aep
    estimate = Estimate(target=float(target), lower_aep=float(lower_aep),
                        upper_aep=upper, degree=int(degree))
    if not (target and target > 0):
        estimate.refused = "no PMF value to place"
        return estimate
    z_low = variate(1.0 / float(lower_aep))
    z_high = variate(1.0 / upper) if math.isfinite(upper) else math.inf
    inside = (real.z >= z_low) & (real.z <= z_high)
    estimate.n = int(inside.sum())
    if estimate.n < MIN_POINTS or estimate.n <= degree + 1:
        estimate.refused = (f"only {estimate.n} realisations between 1 in "
                            f"{results.format_aep(lower_aep)} and 1 in "
                            f"{results.format_aep(upper)} - widen the window")
        return estimate
    x = np.log10(real.value[inside])
    y = real.z[inside]
    coefficients = np.polyfit(x, y, int(degree))
    estimate.coefficients = [float(c) for c in coefficients]
    x_target = math.log10(float(target))
    estimate.z = float(np.polyval(coefficients, x_target))
    estimate.aep = aep_of_variate(estimate.z)

    if estimate.n < THIN_POINTS:
        estimate.warnings.append(f"{estimate.n} realisations in the window - thin")
    if target > real.top_value:
        estimate.warnings.append(
            f"the PMF ({target:.2f}) is above every realisation (top "
            f"{real.top_value:.2f}) - this is an extrapolation, and the degree decides it")
    elif target > float(real.value[inside].max()):
        estimate.warnings.append("the PMF is above the window - raise its upper bound")
    if target < float(real.value[inside].min()):
        estimate.warnings.append("the PMF is below the window - lower its lower bound")
    span = np.linspace(min(float(x.min()), x_target), max(float(x.max()), x_target), 200)
    slope = np.polyval(np.polyder(coefficients), span)
    if np.any(slope <= 0):
        estimate.warnings.append("the fitted curve turns over inside the range it is "
                                 "read on - use a lower degree")
    return estimate


def sensitivity(real: Realisations, target: float, lower_aep: float,
                upper_aep: float | None = None) -> pd.DataFrame:
    """AEP of the target for window lower bounds around the chosen one, by degree."""
    lowers = sorted({float(lower_aep) / 2.5, float(lower_aep), float(lower_aep) * 2.0})
    grid = pd.DataFrame(index=lowers, columns=list(DEGREES), dtype=float)
    for lower in lowers:
        for degree in DEGREES:
            estimate = fit(real, target, lower, upper_aep, degree)
            grid.loc[lower, degree] = estimate.aep if estimate.ok else math.nan
    grid.index.name = "window from 1 in"
    return grid


def fitted_curve(estimate: Estimate, real: Realisations, points: int = 60) -> list:
    """(z, value) along the fit, over the window and out to the target."""
    if not estimate.coefficients:
        return []
    z_low = variate(1.0 / estimate.lower_aep)
    inside = real.value[real.z >= z_low]
    if not len(inside):
        return []
    low = math.log10(float(inside.min()))
    high = max(math.log10(float(inside.max())), math.log10(estimate.target))
    out = []
    for x in np.linspace(low, high, points):
        out.append((float(np.polyval(estimate.coefficients, x)), float(10 ** x)))
    return out


# -- which realisations: the Monte Carlo group's rows ----------------------------

def mc_databases(study: Study, run_name: str, group: str) -> dict:
    """Duration label ('9h') -> the row's mcdf path, for one Monte Carlo group."""
    run = study.open_run(run_name)
    found = {}
    for index in run.rows_in(group):
        row = run.frame.loc[index]
        path = find_database(row, run.project_folder)
        if path is None:
            continue
        duration = results.duration_of(row, cell_text(row.get("Output file")))
        label = f"{duration:g}h" if duration is not None else f"row {index + 2}"
        found[label] = path
    return dict(sorted(found.items(), key=lambda item: _hours(item[0])))


def _hours(label: str) -> float:
    try:
        return float(str(label).rstrip("h"))
    except ValueError:
        return math.inf


def nearest_duration(labels, hours) -> str | None:
    """The label of the duration closest to ``hours`` ('9h' for 9.0)."""
    labels = [label for label in labels if math.isfinite(_hours(label))]
    if not labels or not math.isfinite(float(hours)):
        return None
    return min(labels, key=lambda label: abs(_hours(label) - float(hours)))


def critical_at(study: Study, run_name: str, group: str, aep: float) -> str | None:
    """The level's critical duration at ``aep`` (nearest standard AEP)."""
    from .reporttables import group_curves          # local: reporttables imports this
    curves = group_curves(study, run_name, group)
    critical = curves.critical("level").dropna()
    if critical.empty:
        return None
    z = results.normal_variate(aep)
    nearest = min(critical.index, key=lambda each: abs(results.normal_variate(each) - z))
    return str(critical.loc[nearest])


# -- what the study keeps ------------------------------------------------------

PMF_KEY = "pmf"


def settings(study: Study) -> dict:
    """The study's PMF section, completed with defaults."""
    stored = study.extra.get(PMF_KEY) or {}
    out = {"adopted_aep": stored.get("adopted_aep"),
           "adopted_note": stored.get("adopted_note", ""),
           "pmp_aep": stored.get("pmp_aep", DEFAULT_PMP_AEP),
           "reference_levels": list(stored.get("reference_levels") or []),
           "groups": [dict(entry) for entry in stored.get("groups") or []]}
    for entry in out["groups"]:
        entry.setdefault("label", entry.get("ensemble", {}).get("group", ""))
        entry.setdefault("ensemble", {"run": "", "group": ""})
        entry.setdefault("mc", {"run": "", "group": ""})
        entry.setdefault("duration", "")          # '' = the PMF event's own duration
        entry.setdefault("lower_aep", DEFAULT_LOWER_AEP)
        entry.setdefault("upper_aep", None)       # None = the top of the sample
        entry.setdefault("degree", 1)
    return out


def store(study: Study, section: dict) -> None:
    study.extra[PMF_KEY] = section


def aep_at_value(real: Realisations, value: float) -> float:
    """The AEP at which one duration's realisations reach ``value``, unfitted.

    What ``AEPofDCF_v2.py`` does for the dam crest flood: sort the events by the
    result, interpolate the standard normal variate linearly in log(result).
    Dense enough at a DCF (thousands of events); too sparse at the PMF, which is
    why that one is fitted instead. NaN outside the range the run reached.
    """
    if not len(real.value) or not (value and value > 0):
        return math.nan
    order = np.argsort(real.value, kind="mergesort")
    values, variates = real.value[order], real.z[order]
    if value < values[0] or value > values[-1]:
        return math.nan
    return aep_of_variate(float(np.interp(math.log(value), np.log(values), variates)))


def value_at_aep(real: Realisations, aep: float) -> float:
    """One duration's result at a '1 in X' AEP, read straight off its realisations.

    Linear in the AEP, as ``List_1EY_results.py`` reads the 1 EY level
    (``np.interp(aep, 1 / level_aep, level)``) - an AEP more frequent than the
    quantile tables go (they start at 1 in 2), so the database is the only place
    it can come from. NaN outside the range the run reached.
    """
    if not len(real.z) or not (aep and aep > 1):
        return math.nan
    aeps = np.array([aep_of_variate(float(z)) for z in real.z])
    keep = np.isfinite(aeps)
    aeps, values = aeps[keep], real.value[keep]
    if not len(aeps) or aep < aeps[0] or aep > aeps[-1]:
        return math.nan
    return float(np.interp(float(aep), aeps, values))
