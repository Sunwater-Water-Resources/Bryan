"""Ensemble results: the median pattern by AEP and duration, and the critical duration.

An ensemble run routes every temporal pattern of every storm duration at each
standard AEP. Bryan's own analysis (``lib/EnbAnalysis.py``) takes, for each AEP
and duration, the **median pattern** - position ``int(np.around(n / 2))`` of the
ascending sort, the sixth of ten - and calls the duration whose median is largest
**critical**. It writes ``csv/<name>_critical.csv`` and one box-plot image per AEP.

This reads the database rather than those files. The medians, AEP against
duration, are the same shape as the Monte Carlo quantile tables, so they go
through ``core/results.analyse`` unchanged: the same envelope, the same crossover
checks, and the margin in metres for lake level. Where Bryan's
``_critical.csv`` is on disk the answer is checked against it (``check``).

Each result type is analysed on its own, as EnbAnalysis does: the critical
duration for inflow is the inflow's own, not the level's.

Two small differences from EnbAnalysis, neither of which moves a median's value:

- ties are broken by a stable sort, so where several patterns give the same
  value (a lake held at full supply) the pattern named may differ from Bryan's;
- an event with no result (a failed run) is left out of the count, where
  EnbAnalysis sorts it last and counts it. ``load`` reports any.

pandas and the standard library only; no nicegui.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import ensemble, grouping, results
from .columns import ENSEMBLE, RESERVOIR_ROUTING, normalise_method
from .outputs import find_database

RESULT_TYPES = ("level", "inflow", "outflow")
REQUIRED = ("rain_aep", "duration", "tp") + RESULT_TYPES


class EnsembleError(Exception):
    """A database that cannot be read as an ensemble."""


@dataclass(frozen=True)
class Source:
    """One sims-list row's ensemble database."""

    label: str
    path: Path


def duration_label(hours) -> str:
    return f"{float(hours):g}h"


# -- finding the databases -------------------------------------------------------

def _columns(path: Path) -> set:
    try:
        if path.suffix.lower() == ".parquet":
            return set(pd.read_parquet(path).columns)
        return set(pd.read_csv(path, nrows=0).columns)
    except Exception:                            # noqa: BLE001 - unreadable is "no"
        return set()


def looks_like_ensemble(path: Path) -> bool:
    """The test ``ReservoirRouting._detect_scheme`` makes: duration and pattern,
    and no Monte Carlo sample position."""
    columns = _columns(Path(path))
    return {"duration", "tp"} <= columns and not {"m", "n"} <= columns


def sources_by_group(project) -> dict:
    """Group key -> the ensemble databases on disk, in sims-list order.

    Ensemble rows, and reservoir-routing rows that re-routed an ensemble - told
    apart from re-routed Monte Carlo runs by the database's own columns.
    """
    frame = project.frame
    folder = project.config.project_folder
    keys = project.group_keys()
    found: dict = {}
    for key in grouping.groups_in_order(frame):
        sources = []
        for index in frame.index:
            if keys[index] != key:
                continue
            row = frame.loc[index]
            method = normalise_method(row.get("Method"))
            if method not in (ENSEMBLE, RESERVOIR_ROUTING):
                continue
            path = find_database(row, folder)
            if path is None:
                continue
            if method == RESERVOIR_ROUTING and not looks_like_ensemble(path):
                continue
            sources.append(Source(label=path.stem, path=path))
        if sources:
            found[key] = sources
    return found


# -- reading ---------------------------------------------------------------------

@dataclass
class Database:
    frame: pd.DataFrame
    sources: list
    notes: list = field(default_factory=list)

    @property
    def aeps(self) -> list:
        return sorted(float(aep) for aep in self.frame["rain_aep"].dropna().unique())


def load(sources) -> Database:
    """Every event of a group's databases, with the numbers made numeric."""
    frames = []
    for source in sources:
        try:
            frame = (pd.read_parquet(source.path) if source.path.suffix.lower() == ".parquet"
                     else pd.read_csv(source.path, index_col=0))
        except Exception as exc:                 # noqa: BLE001 - report, never crash
            raise EnsembleError(f"{source.path.name} could not be read ({exc})") from exc
        frames.append(frame.assign(source=source.label))
    if not frames:
        raise EnsembleError("no ensemble database on disk")
    frame = pd.concat(frames, ignore_index=True)
    missing = [name for name in REQUIRED if name not in frame.columns]
    if missing:
        raise EnsembleError(f"the database has no {', '.join(missing)} column - "
                            f"is this an ensemble run?")
    for name in ("rain_aep", "duration") + RESULT_TYPES:
        frame[name] = pd.to_numeric(frame[name], errors="coerce")
    if "storm_method" not in frame.columns:
        frame["storm_method"] = ""

    notes = []
    failed = int(frame[list(RESULT_TYPES)].isna().any(axis=1).sum())
    if failed:
        notes.append(f"{failed} event(s) have no result and are left out of the medians. "
                     f"Bryan's own analysis counts them, so its median can differ.")
    repeated = frame.duplicated(subset=["rain_aep", "duration", "storm_method", "tp"])
    if repeated.any():
        notes.append(f"{int(repeated.sum())} event(s) appear in more than one of the "
                     f"group's databases - the medians mix those runs.")
    return Database(frame=frame, sources=list(sources), notes=notes)


# -- the medians -----------------------------------------------------------------

@dataclass
class Medians:
    """Each AEP and duration's median-pattern value, and which pattern gave it."""

    comparison: results.Comparison       # AEP x duration label, key = the result
    patterns: pd.DataFrame               # the same cells: the median pattern's label
    counts: pd.DataFrame                 # the same cells: how many patterns


def medians(frame: pd.DataFrame, result: str = "level") -> Medians:
    """The median pattern of every AEP and duration, by EnbAnalysis's rule.

    Durations keep the order the database first lists them in, as EnbAnalysis
    takes them, so a tie for critical goes to the same duration it gives it to.
    """
    order = [float(d) for d in pd.unique(frame["duration"].dropna())]
    labels = [duration_label(d) for d in order]
    values, patterns, counts = {}, {}, {}
    for aep, local in frame.dropna(subset=["rain_aep"]).groupby("rain_aep", sort=True):
        picked = ensemble.median_events(local, result)
        sizes = local.dropna(subset=[result]).groupby("duration").size()
        for duration, event in picked.items():
            label = duration_label(duration)
            values[(float(aep), label)] = float(local.loc[event, result])
            patterns[(float(aep), label)] = ensemble.pattern_label(local.loc[event])
            counts[(float(aep), label)] = int(sizes.get(duration, 0))

    def table(cells, fill):
        out = pd.Series(cells, dtype=object).unstack() if cells else pd.DataFrame()
        out = out.reindex(columns=[label for label in labels if label in out.columns])
        out.index.name = results.AEP_COLUMN
        return out if fill is None else out.where(out.notna(), fill)

    frame_values = table(values, None).astype(float)
    comparison = results.Comparison(
        frame=frame_values.dropna(axis=0, how="all"),
        durations={label: d for label, d in zip(labels, order)
                   if label in frame_values.columns},
        key=result)
    return Medians(comparison=comparison, patterns=table(patterns, ""),
                   counts=table(counts, 0))


def highest(frame: pd.DataFrame, result: str = "level") -> pd.DataFrame:
    """The single highest event at each AEP, over every pattern and duration."""
    rows = {}
    for aep, local in frame.dropna(subset=["rain_aep"]).groupby("rain_aep", sort=True):
        picked = ensemble.pick(local, ensemble.HIGHEST, result)
        if picked.found:
            rows[float(aep)] = {"highest": getattr(picked, result),
                                "highest duration": duration_label(picked.duration),
                                "highest pattern": picked.pattern}
    out = pd.DataFrame.from_dict(rows, orient="index")
    out.index.name = results.AEP_COLUMN
    return out


def critical_table(found: Medians, analysis: results.Analysis,
                   top: pd.DataFrame) -> pd.DataFrame:
    """The Results page's table, with the median pattern and the highest event."""
    frame = results.table(found.comparison, analysis)
    if frame.empty:
        return frame
    frame["median pattern"] = [
        found.patterns.loc[aep, owner] if isinstance(owner, str) else ""
        for aep, owner in analysis.critical.items()]
    return frame.join(top, how="left")


# -- one AEP ---------------------------------------------------------------------

def at_aep(frame: pd.DataFrame, aep: float) -> pd.DataFrame:
    return frame[frame["rain_aep"] == float(aep)]


# -- the check against Bryan -----------------------------------------------------

@dataclass
class Check:
    path: Path
    compared: int = 0
    differences: list = field(default_factory=list)
    older: bool = False            # the csv predates the database

    @property
    def agrees(self) -> bool:
        return self.compared > 0 and not self.differences


def critical_csv(source: Source) -> Path:
    """Where lib/EnbAnalysis.py writes the critical durations for this database."""
    return source.path.parent / "csv" / f"{source.path.stem}_critical.csv"


def check(sources, analysis: results.Analysis, result: str) -> Check | None:
    """This page's critical durations against the ``_critical.csv`` Bryan wrote.

    Only for a group of one database - Bryan analyses each run on its own - and
    only where the file is on disk. Compares the critical duration and its median
    at every AEP both have; the pattern label is not compared, see the module note.
    """
    if len(sources) != 1:
        return None
    path = critical_csv(sources[0])
    if not path.is_file():
        return None
    try:
        bryan = pd.read_csv(path, index_col=0)
    except Exception:                            # noqa: BLE001 - no check rather than a crash
        return None
    if result not in bryan.columns or f"{result}_duration" not in bryan.columns:
        return None
    out = Check(path=path, older=path.stat().st_mtime < sources[0].path.stat().st_mtime)
    for aep, owner in analysis.critical.items():
        if not isinstance(owner, str) or float(aep) not in bryan.index.astype(float):
            continue
        row = bryan.loc[bryan.index.astype(float) == float(aep)].iloc[0]
        theirs_duration = duration_label(row[f"{result}_duration"])
        theirs = float(row[result])
        ours = float(analysis.envelope.loc[aep])
        out.compared += 1
        if theirs_duration != owner or not math.isclose(theirs, ours, rel_tol=1e-9,
                                                        abs_tol=1e-9):
            out.differences.append(
                f"1 in {results.format_aep(aep)}: Bryan {theirs:,.6g} at "
                f"{theirs_duration}, here {ours:,.6g} at {owner}")
    return out
