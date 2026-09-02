"""Finding the Monte Carlo databases, and asking lib/RepresentativeEvents about them.

The analysis itself is Bryan's (``lib/RepresentativeEvents.py``, reached through
``core/bryan``), so that the util script which plots the chosen events cannot
rank them differently from the page that chose them. What lives here is
everything the analysis should not have to know: which sims-list row holds which
database, which duration is critical at a given AEP, and where the chosen list
is kept.

Two things worth stating.

**This is the first part of the UI that reads an mcdf.** The Results page
deliberately reads only the quantile tables - ten rows each - because that is
all a frequency curve needs. A representative event is a *realisation*, so
there is no version of this that avoids the database, and an mcdf is m x n rows
(10,000 is ordinary). Hence the cache: the expensive half of the analysis
depends only on the file and the result type, not on the loading being matched,
so it is computed once and every re-rank after that is arithmetic on a frame
already in memory.

**The source of an event should be the critical duration.** A design quantile
at a given AEP comes from whichever duration governs there, and for lake level
that duration changes across the frequency range - long while the storage fills,
short on the rare tail. Picking a representative event out of whatever run
happened to be selected would quietly answer a different question, so
``critical_source`` reuses the Results page's own envelope to suggest the row.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import grouping, results
from .bryan import representative_events
from .columns import MONTE_CARLO, RESERVOIR_ROUTING, normalise_method
from .outputs import find_database
from .paths import atomic_write_json, cell_text, normalise_sep

# Bryan's own analysis. Imported through the seam so the allow-list test sees it.
EVENTS = representative_events()

Filters = EVENTS.Filters
Target = EVENTS.Target

# Only these two methods write a per-realisation database. An ensemble run has
# no realisations to choose between - every combination it ran is in the file
# by construction.
METHODS_WITH_A_DATABASE = (MONTE_CARLO, RESERVOIR_ROUTING)

# How many prepared frames to keep. Each is an mcdf plus a dozen columns; a
# handful covers the durations of one group, which is the working set.
CACHE_LIMIT = 6

_CACHE: dict = {}


@dataclass(frozen=True)
class EventSource:
    """One row's Monte Carlo database."""

    label: str                 # '48h', or the output name when there is no duration
    row_index: object
    output_name: str
    duration: float | None
    path: Path
    group: str

    @property
    def sort_key(self) -> tuple:
        return (0, self.duration) if self.duration is not None else (1, 0.0)


# -- what there is to choose from --------------------------------------------

def sources_for_rows(project, rows=None) -> list:
    """Every row with a results database on disk.

    ``find_database`` is the same probe the Select page's completion check
    uses, so a reservoir routing row resolves to its suffixed mcdf and, failing
    that, to the legacy unsuffixed one.
    """
    frame = project.frame
    folder = project.config.project_folder
    keys = project.group_keys()
    indices = list(frame.index) if rows is None else [i for i in rows if i in frame.index]

    found = []
    for index in indices:
        row = frame.loc[index]
        if normalise_method(row.get("Method")) not in METHODS_WITH_A_DATABASE:
            continue
        path = find_database(row, folder)
        if path is None:
            continue
        output_name = cell_text(row.get("Output file"))
        display = Path(normalise_sep(output_name)).name if output_name else ""
        duration = results.duration_of(row, output_name)
        found.append(EventSource(
            label=f"{duration:g}h" if duration is not None else (
                display or f"row {index + 2}"),
            row_index=index, output_name=display, duration=duration,
            path=path, group=str(keys[index]),
        ))
    return _uniquely_labelled(sorted(found, key=lambda source: source.sort_key))


def _uniquely_labelled(sources) -> list:
    """Two rows of one duration must not both be called '24h'.

    Same reason as ``results._uniquely_labelled``: a duplicate label silently
    picks the wrong database, and one duration run twice under different options
    is exactly the case this page is used for.
    """
    counts: dict = {}
    for source in sources:
        counts[source.label] = counts.get(source.label, 0) + 1

    resolved, used = [], set()
    for source in sources:
        label = source.label
        if counts[label] > 1:
            label = f"{label} ({source.output_name or source.path.stem})"
        while label in used:
            label = f"{label}'"
        used.add(label)
        resolved.append(EventSource(
            label=label, row_index=source.row_index,
            output_name=source.output_name, duration=source.duration,
            path=source.path, group=source.group))
    return resolved


def sources_by_group(project) -> dict:
    """Databases grouped as the Select and Results pages group runs."""
    out: dict = {}
    for source in sources_for_rows(project):
        out.setdefault(source.group, []).append(source)
    return {key: out[key] for key in grouping.groups_in_order(project.frame)
            if key in out}


# -- loading, with the expensive half cached ---------------------------------

def prepared(source: EventSource, result_type: str) -> pd.DataFrame:
    """The database with its target-independent metrics, ready to score.

    Keyed on the file's mtime and size as well as its path, so a re-run of the
    row is picked up without a restart.
    """
    try:
        stat = source.path.stat()
        key = (str(source.path), stat.st_mtime_ns, stat.st_size, result_type)
    except OSError as error:
        raise ValueError(f"{source.path.name} could not be read ({error})") from error

    if key in _CACHE:
        return _CACHE[key]

    try:
        frame = EVENTS.load_mcdf(source.path)
    except Exception as error:               # noqa: BLE001 - report, never crash
        raise ValueError(f"{source.path.name} could not be read ({error})") from error

    value = EVENTS.prepare(frame, result_type)
    if len(_CACHE) >= CACHE_LIMIT:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = value
    return value


def forget_cached() -> None:
    """Drop the cache - for the page's Reload button."""
    _CACHE.clear()


# -- which duration the event should come from -------------------------------

def _curve_sources(project, sources, result_type):
    """The Results page's curve sources for the same rows, keyed by row index."""
    rows = [source.row_index for source in sources]
    found = results.sources_for_rows(project.frame, project.config.project_folder, rows)
    return found.get(result_type, [])


def level_curve(project, sources) -> pd.Series:
    """The level frequency curve a level target is read off.

    The **envelope** over the durations, not one duration's curve, because the
    envelope is the design quantile - the number a loading of '220.5 m AHD' was
    quoted from in the first place.
    """
    curves = _curve_sources(project, sources, "level")
    if not curves:
        return pd.Series(dtype=float)
    return results.analyse(results.compare(curves)).envelope


def critical_source(project, sources, result_type, target_aep):
    """The row whose duration governs at this AEP, or None.

    Returns ``(source, note)``. The note says which AEP on the curve the answer
    was read at, because the critical duration is only evaluated at the standard
    AEPs and a target between two of them takes the nearer one.
    """
    if not sources or not target_aep:
        return None, ""
    curves = _curve_sources(project, sources, result_type)
    if not curves:
        return None, ""

    analysis = results.analyse(results.compare(curves))
    critical = analysis.critical.dropna()
    if critical.empty:
        return None, ""

    z_target = EVENTS.normal_variate(target_aep)
    nearest = min(critical.index,
                  key=lambda aep: abs(EVENTS.normal_variate(aep) - z_target))
    label = critical.loc[nearest]

    by_row = {source.row_index: source for source in sources}
    for curve in curves:
        if curve.label == label and curve.row_index in by_row:
            source = by_row[curve.row_index]
            note = (f"{source.label} is critical for {result_type} at "
                    f"1 in {results.format_aep(nearest)}")
            return source, note
    return None, ""


# -- evaluating one loading --------------------------------------------------

@dataclass
class Outcome:
    """One target's answer, ready to render."""

    target: Target
    aep: float | None = None
    source: EventSource | None = None
    ranking: object = None                   # EVENTS.Ranking
    notes: list = field(default_factory=list)
    problem: str = ""

    @property
    def candidates(self) -> pd.DataFrame:
        if self.ranking is None:
            return pd.DataFrame()
        return self.ranking.candidates

    @property
    def picked(self):
        """The chosen row, defaulting to the best candidate."""
        frame = self.candidates
        if frame.empty:
            return None
        if self.target.picked is not None and self.target.picked in frame.index:
            return frame.loc[self.target.picked]
        return frame.iloc[0]

    @property
    def picked_id(self):
        row = self.picked
        return None if row is None else row.name


def evaluate(project, sources, target: Target, filters: Filters,
             curve=None) -> Outcome:
    """Rank the realisations of one database against one loading."""
    outcome = Outcome(target=target)
    if not sources:
        outcome.problem = "No Monte Carlo database has been run for this group."
        return outcome

    aep, above_curve = target.value, False
    if target.kind == "level":
        curve = level_curve(project, sources) if curve is None else curve
        if curve is None or len(curve) == 0:
            outcome.problem = ("There is no level frequency curve to read this "
                               "loading off - the rows have not been analysed.")
            return outcome
        lookup = EVENTS.aep_for_level(curve, target.value)
        above_curve = lookup.above_curve
        aep = lookup.aep
        if lookup.note:
            outcome.notes.append(lookup.note)
        if aep is None and not above_curve:
            outcome.problem = lookup.note or "The level is off the frequency curve."
            return outcome
    outcome.aep = aep

    source = _source_named(sources, target.source)
    if source is None:
        source, note = critical_source(project, sources, target.result_type,
                                       aep or target.value)
        if source is None:
            source = sources[0]
        elif note:
            outcome.notes.append(note)
    outcome.source = source

    try:
        frame = prepared(source, target.result_type)
    except ValueError as error:
        outcome.problem = str(error)
        return outcome

    # Above the curve there is no target to be near; rank() sorts by level
    # instead, and the scoring AEP is only there to fill the distance columns.
    scored = EVENTS.score(frame, aep or _top_aep(frame), target.rain_aep)
    outcome.ranking = EVENTS.rank(scored, filters, target.count,
                                  above_curve=above_curve)
    return outcome


def _top_aep(frame) -> float:
    """A stand-in target for the above-the-curve case: the rarest result there.

    Finite values only - a realisation whose TPT probability underflowed has an
    infinite AEP, and that is not a target anything can be scored against.
    """
    values = pd.to_numeric(frame.get("result_aep"), errors="coerce")
    values = values[values.notna() & (values != float("inf"))]
    return float(values.max()) if len(values) else 100.0


def _source_named(sources, name):
    if not name:
        return None
    for source in sources:
        if source.label == name or source.output_name == name:
            return source
    return None


# -- the chosen list ---------------------------------------------------------

def summary_rows(outcomes) -> list:
    """The list this whole page exists to produce: one row per loading."""
    rows = []
    for outcome in outcomes:
        row = {
            "loading": outcome.target.label,
            "result": outcome.target.result_type,
            "target aep (1 in x)": (round(outcome.aep, 1)
                                    if outcome.aep else ""),
            "source": outcome.source.label if outcome.source else "",
            "output file": outcome.source.output_name if outcome.source else "",
            "sim": "", "hydrograph": "",
            "rain aep (1 in x)": "", "result aep (1 in x)": "",
            "delta z": "", "level": "", "inflow": "", "outflow": "",
            "sub-burst ratio": "", "preburst percentile": "",
            "lake z": "", "adv (ml)": "", "flags": outcome.problem,
        }
        picked = outcome.picked
        if picked is not None:
            row.update({
                "sim": int(picked.name),
                "hydrograph": EVENTS.sim_label(picked.name),
                "rain aep (1 in x)": _round(picked.get("rain_aep"), 1),
                "result aep (1 in x)": _round(picked.get("result_aep"), 1),
                "delta z": _round(picked.get("delta_z"), 3),
                "level": _round(picked.get("level"), 2),
                "inflow": _round(picked.get("inflow"), 1),
                "outflow": _round(picked.get("outflow"), 1),
                "sub-burst ratio": _round(picked.get("subburst_ratio"), 2),
                "preburst percentile": _round(picked.get("preburst_p"), 2),
                "lake z": _round(picked.get("lake_z"), 2),
                "adv (ml)": _round(picked.get("ADV"), 0),
                "flags": "; ".join(picked.get("flags") or ()),
            })
        rows.append(row)
    return rows


def _round(value, places):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if number != number:                     # NaN
        return ""
    return round(number, places)


# The candidate table, in the order the columns are worth reading: how close
# the event is, then what would make it indefensible, then what it produced.
CANDIDATE_COLUMNS = (
    ("sim", "Sim"),
    ("rain_aep", "Rain AEP"),
    ("result_aep", "Result AEP"),
    ("delta_z", "\u0394z"),
    ("subburst_ratio", "Sub-burst"),
    ("preburst_p", "Pre-burst p"),
    ("lake_z", "Lake z"),
    ("level", "Level"),
    ("inflow", "Inflow"),
    ("outflow", "Outflow"),
    ("flags", "Flags"),
)

_PLACES = {"rain_aep": 0, "result_aep": 0, "delta_z": 3, "subburst_ratio": 2,
           "preburst_p": 2, "lake_z": 2, "level": 2, "inflow": 0, "outflow": 0}


def candidate_rows(outcome) -> list:
    """The ranked candidates as table rows - plain types, ready for ui.table."""
    frame = outcome.candidates
    if frame is None or frame.empty:
        return []
    picked = outcome.picked_id
    rows = []
    for sim_id, row in frame.iterrows():
        # Plain bool, not numpy's: these rows are serialised to the browser,
        # and json.dumps refuses a numpy bool_ - which is a blank table, not an
        # error anyone would see in a test that only builds the rows.
        entry = {"sim": int(sim_id), "picked": bool(sim_id == picked),
                 "flags": "; ".join(row.get("flags") or ()) or "-"}
        for column, places in _PLACES.items():
            value = _round(row.get(column), places)
            entry[column] = value if value != "" else "-"
        rows.append(entry)
    return rows


def selection_path(folder, group: str) -> Path:
    """Where a group's chosen events are kept.

    Named for the group because several groups commonly share one results
    folder - a GWL series usually does - and one file per scenario is what
    stops them overwriting each other.
    """
    safe = "".join(character if character.isalnum() or character in "-_." else "_"
                   for character in str(group)).strip("_")
    stem = f"{safe}_{EVENTS.SELECTION_FILE}" if safe else EVENTS.SELECTION_FILE
    return Path(folder) / stem


def default_folder(sources) -> Path | None:
    """The results folder the databases are in."""
    return sources[0].path.parent if sources else None


def load_targets(path) -> tuple:
    return EVENTS.read_selection(path)


def save_targets(path, targets, settings=None) -> None:
    atomic_write_json(Path(path), EVENTS.selection_payload(targets, settings))
