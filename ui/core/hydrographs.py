"""One realisation's stored hydrographs, for previewing a candidate event.

Closeness to a loading is not the only thing that makes an event usable. A
representative event is usually wanted **simple**: one rise, one peak, one
recession. A double-peaked outflow is defensible as a flood but awkward as a
design event - it makes a gate operation ambiguous and a dambreak run arguable
- and nothing in the mcdf says which candidates have one, because the database
records the peaks and not the shape between them.

So this reads the stored hydrographs, which exist only where the run was made
with ``Store hydrographs`` set to ``yes``. That is the whole caveat: no files,
no preview, and the page says so rather than pretending the events are all
alike.

**Reading is deliberately explicit and cached.** A stored hydrograph file is
one column per simulation - 10,000 of them is ordinary - so it is tens of
megabytes and takes seconds to parse. Nothing here happens on a redraw: the
page asks for a run's hydrographs when the user presses the button, the frame
is kept per file (mtime and size, so a re-run is picked up), and every preview
after that is a lookup. Two files are kept, which covers looking at two
durations of one group in turn.

Peak counting is done here rather than in ``lib/RepresentativeEvents.py``
because it is about a hydrograph, not about a realisation - and this module,
like the rest of ``ui/core``, is pandas and the standard library only, so there
is no ``scipy.signal.find_peaks`` to reach for.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .bryan import representative_events
from .columns import RESERVOIR_ROUTING, normalise_method
from .paths import cell_text, resolve_value

EVENTS = representative_events()

# The order they are drawn in: the outflow is what the question is usually
# about, the inflow says whether the shape came in with the flood or was made
# by the dam, and the level is what a level loading was quoted from.
KINDS = ("outflows", "inflows", "levels")

# Two runs' hydrographs, which is a duration and the one being compared with it.
CACHE_LIMIT = 2

# A local maximum counts as a peak of its own when it stands this far above the
# deepest trough between it and a higher one, as a fraction of the series peak.
# Below that it is a shoulder or a wobble on the recession, not a second flood.
PROMINENCE = 0.10

# And anything below this fraction of the peak is not worth calling a peak at
# all - a small early rise before the main event is not what makes a hydrograph
# awkward to defend.
FLOOR = 0.20

_CACHE: dict = {}


@dataclass(frozen=True)
class Shape:
    """What a hydrograph looks like, in the few terms a choice is made on."""

    peaks: int = 0
    peak_time: float | None = None
    summary: str = ""

    @property
    def single(self) -> bool:
        return self.peaks == 1


def paths_for(project, source) -> dict:
    """The stored hydrograph files of one run, by kind - existing ones only.

    ``EVENTS.hydrograph_paths`` is the layout Bryan's own writers use, and it
    knows the two of them: a monte carlo row writes into the model's results
    folder beside the run, while a reservoir routing row writes the levels and
    outflows into its ``Hydrographs folder`` and never writes the inflows at
    all - those are its input, named by the sims-list ``Inflow``.
    """
    row = project.frame.loc[source.row_index]
    folder = project.config.project_folder
    routed = normalise_method(row.get("Method")) == RESERVOIR_ROUTING

    found = EVENTS.hydrograph_paths(
        cell_text(row.get("Output file")),
        model=cell_text(row.get("Model type")) or "urbs",
        hydrographs_folder=(resolve_value(folder, row.get("Hydrographs folder"))
                            if routed else None),
        suffix=cell_text(row.get("Output suffix")),
        inflow=resolve_value(folder, row.get("Inflow")),
    )
    out = {}
    for kind, path in found.items():
        candidate = resolve_value(folder, path)
        if candidate is not None and Path(candidate).is_file():
            out[kind] = Path(candidate)
    return out


def load(path) -> pd.DataFrame:
    """A stored hydrograph file, cached on its path, mtime and size.

    Raises ValueError with something worth showing rather than a traceback: a
    parquet inflow file needs pyarrow, which the UI environment is not obliged
    to have.
    """
    path = Path(path)
    try:
        stat = path.stat()
    except OSError as error:
        raise ValueError(f"{path.name} could not be read ({error})") from error

    key = (str(path), stat.st_mtime_ns, stat.st_size)
    if key in _CACHE:
        return _CACHE[key]

    try:
        frame = EVENTS.read_hydrographs(path)
    except Exception as error:               # noqa: BLE001 - report, never crash
        raise ValueError(f"{path.name} could not be read ({error})") from error

    frame = frame.apply(pd.to_numeric, errors="coerce")
    frame.index = pd.to_numeric(frame.index, errors="coerce")
    frame = frame[frame.index.notna()]

    if len(_CACHE) >= CACHE_LIMIT:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = frame
    return frame


def forget_cached() -> None:
    """Drop the frames - for the page's Reload button."""
    _CACHE.clear()


def series_for(project, source, sim_id, kinds=KINDS) -> tuple:
    """``({kind: series}, [notes])`` for one simulation.

    A kind that has no file, or no column for this simulation, is a note and
    not an error: a RORB run stores no levels, and a routed row's inflows are
    only there if the sims list named them.
    """
    column = EVENTS.sim_label(sim_id)
    paths = paths_for(project, source)
    series, notes = {}, []
    if not paths:
        return series, ["No stored hydrographs for this run - it was run with "
                        "'Store hydrographs' off."]

    for kind in kinds:
        path = paths.get(kind)
        if path is None:
            continue
        try:
            frame = load(path)
        except ValueError as error:
            notes.append(str(error))
            continue
        if column not in frame.columns:
            notes.append(f"{path.name} has no {column}")
            continue
        found = frame[column].dropna()
        if len(found):
            series[kind] = found
    if not series and not notes:
        notes.append(f"Nothing stored for {column}.")
    return series, notes


def shapes_for(project, source, sim_ids, kind="outflows") -> tuple:
    """``({sim: Shape}, [notes])`` for a whole table of candidates at once.

    One read of one file answers the question for every candidate, which is the
    point of doing it per loading rather than per hover.
    """
    paths = paths_for(project, source)
    path = paths.get(kind) or paths.get("inflows")
    if path is None:
        return {}, ["No stored hydrographs for this run - it was run with "
                    "'Store hydrographs' off."]
    try:
        frame = load(path)
    except ValueError as error:
        return {}, [str(error)]

    shapes, missing = {}, 0
    for sim_id in sim_ids:
        column = EVENTS.sim_label(sim_id)
        if column not in frame.columns:
            missing += 1
            continue
        shapes[int(sim_id)] = shape_of(frame[column].dropna())
    notes = [f"{missing} of the candidates are not in {path.name}"] if missing else []
    return shapes, notes


def shape_of(series) -> Shape:
    """How many peaks a hydrograph has, and when the biggest one is.

    A peak is a local maximum that is worth arguing about: at least ``FLOOR``
    of the series peak, and standing at least ``PROMINENCE`` of the peak above
    the deepest trough between it and the next higher one. That second test is
    what separates a genuine second flood from the steps and wobbles every
    routed outflow has on its recession - counting bare local maxima would call
    almost everything multi-peaked and the column would say nothing.
    """
    values = [float(value) for value in pd.to_numeric(series, errors="coerce").dropna()]
    times = [float(time) for time in
             pd.to_numeric(pd.Series(series.index), errors="coerce").dropna()]
    if not values:
        return Shape(summary="nothing stored")

    highest = max(values)
    if highest <= 0:
        return Shape(summary="no flow")

    tops = _local_maxima(values)
    kept = [position for position in tops if values[position] >= FLOOR * highest]
    kept = _prominent(values, kept, highest)

    peak_at = values.index(highest)
    peak_time = times[peak_at] if peak_at < len(times) else None
    if len(kept) <= 1:
        summary = "single peaked"
    else:
        summary = f"{len(kept)} peaks"
    return Shape(peaks=max(len(kept), 1), peak_time=peak_time, summary=summary)


def _local_maxima(values) -> list:
    """Positions of the local maxima, flat tops counted once, at their middle."""
    tops = []
    position = 1
    while position < len(values) - 1:
        if values[position] <= values[position - 1]:
            position += 1
            continue
        end = position
        while end + 1 < len(values) and values[end + 1] == values[position]:
            end += 1
        if end + 1 < len(values) and values[end + 1] < values[position]:
            tops.append((position + end) // 2)
        position = end + 1
    return tops


def _prominent(values, tops, highest) -> list:
    """Drop the maxima that never come back down before a bigger one.

    Walking left to right and keeping a peak only when the trough since the
    last kept peak is deep enough: a shoulder on the rising limb and a step on
    the recession both fail that, and a genuine second flood passes it.
    """
    kept = []
    for position in tops:
        if not kept:
            kept.append(position)
            continue
        previous = kept[-1]
        trough = min(values[previous:position + 1])
        lower = min(values[previous], values[position])
        if lower - trough >= PROMINENCE * highest:
            kept.append(position)
        elif values[position] > values[previous]:
            kept[-1] = position              # the same peak, seen better
    return kept


def for_plot(series, limit=1200) -> tuple:
    """``(times, values)`` thinned enough to draw without losing the shape."""
    step = max(len(series) // limit, 1)
    thinned = series.iloc[::step]
    return ([round(float(time), 4) for time in thinned.index],
            [round(float(value), 4) for value in thinned])
