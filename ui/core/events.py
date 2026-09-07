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


# How close two results count as the same, and the units the page asks for it
# in. A lake level is quoted in millimetres because that is the precision it is
# argued about at; a flow in m3/s, which is its own unit.
BAND_UNITS = {"level": ("mm", 1000.0),
              "inflow": ("m\u00b3/s", 1.0),
              "outflow": ("m\u00b3/s", 1.0)}

# 20 mm on a level: inside that the rating curve and the routing timestep are
# deciding the order, not the hydrology.
DEFAULT_BANDS = {"level": 20.0, "inflow": 10.0, "outflow": 10.0}


def band_units(result_type) -> tuple:
    """(label, how many of them make one of the result's own units)."""
    return BAND_UNITS.get(result_type, ("", 1.0))


def band_in_result_units(result_type, shown) -> float:
    """The page's number, in the units the ranking measures in."""
    try:
        shown = float(shown)
    except (TypeError, ValueError):
        return 0.0
    scale = band_units(result_type)[1]
    return max(shown, 0.0) / scale


def envelope_curve(project, sources, result_type="level") -> pd.Series:
    """The design frequency curve a loading is read against.

    The **envelope** over the durations, not one duration's curve, because the
    envelope is the design quantile - the number a loading of '220.5 m AHD' was
    quoted from in the first place. The same curve read the other way gives the
    design value at a design AEP, which is what ranking on the result needs.
    """
    curves = _curve_sources(project, sources, result_type)
    if not curves:
        return pd.Series(dtype=float)
    return results.analyse(results.compare(curves)).envelope


def level_curve(project, sources) -> pd.Series:
    """The level envelope - the curve a level loading is read off."""
    return envelope_curve(project, sources, "level")


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
    target_value: float | None = None        # the loading in the result's units
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
             curve=None, order=EVENTS.DELTA_Z, curves=None, band=0.0) -> Outcome:
    """Rank the realisations of one database against one loading.

    ``order`` is what closeness means - see ``EVENTS.rank``. Ranking on the
    result needs the loading in the result's own units: a level loading names
    it outright, and a design AEP has it read off the design curve for the
    result type. ``curve`` is the level envelope (a level loading is always
    read off the level curve, whatever result is being ranked); ``curves`` is
    the caller's cache of envelopes by result type, so the page does not
    rebuild one per loading.
    """
    outcome = Outcome(target=target)
    if not sources:
        outcome.problem = "No Monte Carlo database has been run for this group."
        return outcome

    def envelope(result_type):
        if result_type == "level" and curve is not None and len(curve):
            return curve
        if curves is not None and result_type in curves:
            return curves[result_type]
        found = envelope_curve(project, sources, result_type)
        if curves is not None:
            curves[result_type] = found
        return found

    aep, above_curve, target_value = target.value, False, None
    if target.kind == "level":
        # The loading is a lake level, so it is read off the level curve even
        # when the ranking is on inflow.
        target_value = target.value if target.result_type == "level" else None
        level = envelope("level")
        if level is None or len(level) == 0:
            outcome.problem = ("There is no level frequency curve to read this "
                               "loading off - the rows have not been analysed.")
            return outcome
        lookup = EVENTS.aep_for_level(level, target.value)
        above_curve = lookup.above_curve
        aep = lookup.aep
        if lookup.note:
            outcome.notes.append(lookup.note)
        if aep is None and not above_curve:
            outcome.problem = lookup.note or "The level is off the frequency curve."
            return outcome
    outcome.aep = aep

    if order == EVENTS.RESULT and target_value is None and aep and not above_curve:
        # The loading is named in the wrong units for this ranking - a design
        # AEP, or a level while the inflow is being ranked - so the design
        # value at that AEP is read off the curve for the result type.
        design = envelope(target.result_type)
        value = (EVENTS.value_for_aep(design, aep)
                 if design is not None and len(design) else float("nan"))
        target_value = None if value != value else float(value)
        if target_value is None:
            outcome.notes.append(
                f"No design {target.result_type} at 1 in "
                f"{results.format_aep(aep)} on the envelope, so the ranking "
                f"falls back to the result AEP.")
        else:
            outcome.notes.append(
                f"Ranking on {target.result_type}: the design value at 1 in "
                f"{results.format_aep(aep)} is {target_value:,.2f}.")
    if order == EVENTS.RESULT and band:
        unit, scale = band_units(target.result_type)
        outcome.notes.append(
            f"Events within {band * scale:g} {unit} of the loading count as "
            f"reaching it and are ordered by AEP neutrality.")
    outcome.target_value = target_value

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
    scored = EVENTS.score(frame, aep or _top_aep(frame), target.rain_aep,
                          target_value=target_value)
    outcome.ranking = EVENTS.rank(scored, filters, target.count,
                                  above_curve=above_curve, order=order,
                                  band=band)
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
            "target value": _round(outcome.target_value, 2),
            "output file": outcome.source.output_name if outcome.source else "",
            "sim": "", "hydrograph": "",
            "rain aep (1 in x)": "", "result aep (1 in x)": "",
            "delta z": "", "delta target": "",
            "level": "", "inflow": "", "outflow": "",
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
                "delta target": _round(picked.get("delta_value"), 2),
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
    ("delta_value", "\u0394 target"),
    ("subburst_ratio", "Sub-burst"),
    ("preburst_p", "Pre-burst p"),
    ("lake_z", "Lake z"),
    ("level", "Level"),
    ("inflow", "Inflow"),
    ("outflow", "Outflow"),
    ("flags", "Flags"),
)

_PLACES = {"rain_aep": 0, "result_aep": 0, "delta_z": 3, "delta_value": 2,
           "subburst_ratio": 2,
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


def project_relative(project, path) -> str:
    """A path as the sims list would carry it - relative to the project folder.

    The CLI resolves what it finds in the selection file against that folder,
    the way Bryan resolves a sims-list value, so a saved selection survives the
    project being moved or opened from another machine.
    """
    try:
        return str(Path(path).resolve().relative_to(
            Path(project.config.project_folder).resolve()))
    except (ValueError, OSError):
        return str(path)


def extract_command(project, selection_path) -> str:
    """The util command that turns a saved selection into hydrographs and plots."""
    return ("python util/RepresentativeEvents.py --config "
            f"{project.config.config_path} --selection {selection_path}")


def load_targets(path) -> tuple:
    return EVENTS.read_selection(path)


def save_targets(path, targets, settings=None) -> None:
    atomic_write_json(Path(path), EVENTS.selection_payload(targets, settings))
