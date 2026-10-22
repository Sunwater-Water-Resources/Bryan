"""Representative events from an ensemble run: the PMF, and design AEPs by the median pattern.

An ensemble routes every temporal pattern of every duration at each AEP, so there is
nothing to rank: the event is decided by a rule, and the rules are the PMF page's own
(``core/ensemble.py``), so the two pages cannot disagree.

``pmf``  The highest event of any duration and pattern (``ensemble.HIGHEST``).
``aep``  At that AEP, the median pattern of each duration by EnbAnalysis's position rule,
         then the duration whose median is largest (``ensemble.MEDIAN``) - the critical
         duration. A duration given in ``source`` overrides the critical one. The PMPF is
         this rule at the AEP of the PMP.

The event is identified by its row in its own database file - the index Bryan wrote, and
the ``sim_<id>`` column of the stored hydrographs - never by its position in a frame that
joined several files.

pandas and the standard library only; no nicegui.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import ensemble, grouping, results
from .bryan import representative_events
from .columns import ENSEMBLE, RESERVOIR_ROUTING, normalise_method
from .enbresults import looks_like_ensemble
from .outputs import find_database
from .paths import cell_text

EVENTS = representative_events()

PMF = "pmf"
AEP = "aep"
RULES = {PMF: ensemble.HIGHEST, AEP: ensemble.MEDIAN}
RULE_TEXT = {ensemble.HIGHEST: "the highest event of any duration and pattern",
             ensemble.MEDIAN: "the median pattern at the critical duration"}


class EnsembleEventError(Exception):
    """A database that cannot be read as an ensemble."""


@dataclass(frozen=True)
class Source:
    """One sims-list row's ensemble database."""

    output_name: str
    path: Path

    @property
    def label(self) -> str:
        return self.path.stem


def sources_by_group(project) -> dict:
    """Group key -> the ensemble databases on disk, in sims-list order: ensemble rows, and
    reservoir routing rows that re-routed an ensemble (told apart by the database's
    columns, as ``enbresults.sources_by_group`` does)."""
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
            if path is None or (method == RESERVOIR_ROUTING and not looks_like_ensemble(path)):
                continue
            sources.append(Source(output_name=cell_text(row.get("Output file")), path=path))
        if sources:
            found[key] = sources
    return found


_CACHE: dict = {}


def load(sources) -> pd.DataFrame:
    """Every event of a group's databases, one row each, with ``event`` (the row in its own
    file) and ``source`` (which file) beside the results. Cached on each file's mtime."""
    frames = []
    for number, source in enumerate(sources):
        try:
            stat = source.path.stat()
        except OSError as exc:
            raise EnsembleEventError(f"{source.path.name} is not on disk") from exc
        key = (str(source.path), stat.st_mtime, stat.st_size)
        if key not in _CACHE:
            try:
                frame = (pd.read_parquet(source.path)
                         if source.path.suffix.lower() == ".parquet"
                         else pd.read_csv(source.path, index_col=0))
            except Exception as exc:                 # noqa: BLE001 - report, never crash
                raise EnsembleEventError(f"{source.path.name} could not be read ({exc})") from exc
            if "Unnamed: 0" in frame.columns:        # parquet written with index=False
                frame = frame.set_index("Unnamed: 0")
            _CACHE[key] = frame
        frame = _CACHE[key].copy()
        frame["event"] = frame.index
        frame["source"] = number
        frames.append(frame)
    if not frames:
        raise EnsembleEventError("no ensemble database on disk")
    frame = pd.concat(frames, ignore_index=True)
    missing = [name for name in ("rain_aep", "duration", "tp", "level") if name not in frame]
    if missing:
        raise EnsembleEventError(f"the database has no {', '.join(missing)} column - "
                                 f"is this an ensemble run?")
    for name in ("rain_aep", "duration", "level", "inflow", "outflow", "ADV"):
        if name in frame:
            frame[name] = pd.to_numeric(frame[name], errors="coerce")
    return frame


def forget_cached() -> None:
    _CACHE.clear()


def aeps(frame) -> list:
    """The AEPs the ensemble ran, 1 in X, rarest last."""
    return sorted(float(value) for value in frame["rain_aep"].dropna().unique())


def durations(frame) -> list:
    return sorted(float(value) for value in frame["duration"].dropna().unique())


def default_targets(frame, result_type="level") -> list:
    """What a new ensemble group starts with. A run at one AEP is a PMF run - the PMF and
    the PMPF; otherwise the three rarest AEPs it ran."""
    ran = aeps(frame)
    if len(ran) == 1:
        return [EVENTS.Target(kind=PMF, result_type=result_type, name="PMF"),
                EVENTS.Target(kind=AEP, value=ran[0], result_type=result_type, name="PMPF")]
    return [EVENTS.Target(kind=AEP, value=aep, result_type=result_type) for aep in ran[-3:]]


@dataclass
class Choice:
    """The event a loading takes, and what it was chosen from."""

    target: object
    rule: str = ""
    event: int | None = None           # the row in its own file
    source: Source | None = None
    duration: float = math.nan
    critical: float = math.nan          # the duration the rule settled on
    pattern: str = ""
    values: dict = field(default_factory=dict)
    table: pd.DataFrame | None = None   # ensemble.by_duration over what was chosen from
    problem: str = ""

    @property
    def found(self) -> bool:
        return self.event is not None


def _duration_override(target) -> float | None:
    text = str(target.source or "").strip().lower().rstrip("h")
    try:
        return float(text) if text else None
    except ValueError:
        return None


def choose(frame, target, sources) -> Choice:
    """The event for one loading, by its rule."""
    rule = RULES.get(target.kind)
    choice = Choice(target=target, rule=rule or "")
    if rule is None:
        choice.problem = "an ensemble loading is the PMF or a design AEP"
        return choice
    result = target.result_type if target.result_type in frame else "level"
    pool = frame
    if target.kind == AEP:
        pool = frame[frame["rain_aep"] == float(target.value)]
        if pool.empty:
            ran = ", ".join(f"1 in {results.format_aep(aep)}" for aep in aeps(frame))
            choice.problem = (f"this ensemble did not run 1 in {results.format_aep(target.value)}"
                              f" - it ran {ran}")
            return choice
    override = _duration_override(target)
    choice.table = ensemble.by_duration(pool, result)
    if override is not None:
        pool = pool[pool["duration"] == override]
        if pool.empty:
            choice.problem = f"no {override:g} h events at this loading"
            return choice
    picked = ensemble.pick(pool, rule, result)
    if not picked.found:
        choice.problem = f"no event has a {result}"
        return choice
    row = frame.loc[picked.event]
    choice.event = int(row["event"])
    choice.source = sources[int(row["source"])]
    choice.duration = choice.critical = float(row["duration"])
    choice.pattern = picked.pattern
    choice.values = {name: float(row[name]) for name in
                     ("level", "inflow", "outflow", "rain_aep", "ADV")
                     if name in row and pd.notna(row[name])}
    return choice


def stored_target(choice, project_relative) -> object:
    """The loading as the selection file keeps it: the event, its duration and its rule,
    and the database relative to the project. ``project_relative`` is events.project_relative
    bound to the project."""
    target = EVENTS.Target.from_dict(choice.target.to_dict())
    target.rule = choice.rule
    if choice.found:
        target.picked = choice.event
        target.duration = choice.duration
        target.output_file = choice.source.output_name
        target.database = project_relative(choice.source.path)
    else:
        target.picked = None
    return target


def summary_rows(choices) -> list:
    """One row per loading with an event - what the Events summary table and its csv hold."""
    rows = []
    for choice in choices:
        if not choice.found:
            continue
        target, values = choice.target, choice.values
        rows.append({
            "loading": target.title,
            "rule": choice.rule,
            "event": choice.event,
            "hydrograph": f"sim_{choice.event:05d}",
            "duration (h)": f"{choice.duration:g}",
            "pattern": choice.pattern,
            "level": _fmt(values.get("level"), 3),
            "inflow": _fmt(values.get("inflow"), 1),
            "outflow": _fmt(values.get("outflow"), 1),
            "rain AEP": (f"1 in {results.format_aep(values['rain_aep'])}"
                         if "rain_aep" in values else ""),
            "notional AEP": (f"1 in {results.format_aep(target.notional_aep)}"
                             if target.notional_aep else ""),
            "ADV (ML)": _fmt(values.get("ADV"), 0),
            "database": choice.source.label,
        })
    return rows


def duration_rows(choice) -> list:
    """The per-duration table under a loading's card, the duration the rule took marked."""
    table = choice.table
    if table is None or table.empty:
        return []
    rows = []
    for duration, row in table.iterrows():
        rows.append({"duration": f"{duration:g} h", "n": int(row["n"]),
                     "median": f"{row['median']:.3f}", "median pattern": row["median_pattern"],
                     "highest": f"{row['max']:.3f}", "highest pattern": row["max_pattern"],
                     "chosen": bool(choice.found and duration == choice.duration)})
    return rows


def _fmt(value, digits) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{value:,.{digits}f}"
