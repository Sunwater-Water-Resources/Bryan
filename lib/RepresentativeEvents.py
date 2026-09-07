"""Choosing representative events out of a Monte Carlo database.

A design flood quantile is a statistic over thousands of realisations, but the
things that get built - a spillway gate operation, a dambreak run, an emergency
action plan trigger - need a single *event*: one hydrograph, one storm, one
starting lake level. Picking it is a judgement, and this module lays out what
the judgement is made on.

What makes an event representative of a loading:

- **AEP neutrality.** The realisation's rainfall should be about as rare as the
  flood it produced. An event that reaches the 1 in 2,000 lake level off 1 in
  200 rainfall got there through an unusual coincidence - a full lake, a huge
  pre-burst, a pattern with a spike in it - and it will not behave like a 1 in
  2,000 event when anything downstream of it is changed.
- **No borrowed rarity.** The same applies within the storm: an embedded burst
  means part of the pattern is rarer than the storm it sits in, and a pre-burst
  or an antecedent storage far off the median is the same coincidence wearing a
  different hat.

So the ranking is a distance from the target, and everything else is reported
alongside it rather than folded into a single number - the trade-off between a
close AEP match and a clean storm is the user's to make, not a weighting's.

**Distance is measured in standard normal variate space.** ``1 in X`` is not a
linear scale: at a 1 in 2,000 target, being 200 out is nothing, and at 1 in 100
it is everything. z spaces the frequency range evenly, which is also how every
frequency curve in Bryan is plotted. The old ``1 in X`` distance is kept as
``delta_aep`` so a selection made with util/GetRepresentativeEvents.py can still
be checked against its own measure.

**Units, which are mixed and easy to get wrong** (lib/MCScheme.py:26-31 and
lib/Simulator.py:1012):

- ``rain_aep`` is **1 in X**, and ``rain_z`` is its standard normal variate.
- ``level_aep``, ``inflow_aep``, ``outflow_aep`` are **probabilities** - what
  ``TotalProbTheorem.assign_aep`` returns. They need ``1/p`` before they can be
  compared with ``rain_aep``.

**Dependencies are deliberately pandas and the standard library only.** The
run launcher imports this module directly (``ui/core/bryan.py``), and the whole
point of that allow-list is that the UI environment carries neither scipy nor
matplotlib. ``statistics.NormalDist`` covers what ``ndtri``/``ndtr`` were doing.
Keep it that way - anything needing scipy belongs in the util script that plots
the events, not here.
"""

from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass, field
from statistics import NormalDist

import pandas as pd

RESULT_TYPES = ('level', 'inflow', 'outflow')

# What 'closest to the loading' means when the candidates are ranked. See rank().
DELTA_Z = 'delta_z'          # both axes: reach the loading, and be neutral about it
RESULT = 'result'            # the result alone: hit the level, report the neutrality
ORDERS = (DELTA_Z, RESULT)

# lib/Simulator.py:1115-1116 writes these in pairs, one per sub-duration
# shorter than the storm. Their ratio is the quantitative embedded burst
# measure; the 'embedded_bursts' comment is the same finding in words.
SUBBURST_RE = re.compile(r'^subburst_(\d+(?:\.\d+)?)h$')

# StormGenerator.embedded_burst_comment_from_measurement's clean answer.
NO_EMBEDDED_BURSTS = 'No embedded bursts'

_NORMAL = NormalDist()


# -- the standard normal variate ---------------------------------------------

def normal_variate(aep) -> float:
    """z for an AEP expressed as '1 in X'.

    NaN rather than an exception for anything that is not a real AEP, infinity
    included: ``aep_of_variate`` genuinely returns infinity once the upper tail
    underflows (around z = 8), and an AEP read off the top of a steep level
    curve can get there. Left to ``inv_cdf`` that is a StatisticsError out of
    the middle of a page redraw.
    """
    try:
        aep = float(aep)
    except (TypeError, ValueError):
        return math.nan
    if not aep > 1 or not math.isfinite(aep):
        return math.nan
    return _NORMAL.inv_cdf(1.0 - 1.0 / aep)


def variate_of_probability(probability) -> float:
    """z for an exceedance probability - what the TPT columns hold."""
    try:
        probability = float(probability)
    except (TypeError, ValueError):
        return math.nan
    if not 0.0 < probability < 1.0:
        return math.nan
    return _NORMAL.inv_cdf(1.0 - probability)


def aep_of_variate(z) -> float:
    """'1 in X' for a standard normal variate."""
    try:
        z = float(z)
    except (TypeError, ValueError):
        return math.nan
    tail = 1.0 - _NORMAL.cdf(z)
    return math.inf if tail <= 0.0 else 1.0 / tail


# -- reading the database ----------------------------------------------------

def load_mcdf(path) -> pd.DataFrame:
    """A Monte Carlo database, indexed by simulation id.

    The index is the integer sim id (``MCScheme.py:26``), and it is what ties a
    row to its hydrograph: the stored columns are ``sim_00000`` and up
    (``URBSmodel.py:653``), which is what ``sim_label`` below rebuilds.

    Parquet is accepted because reservoir routing already accepts it for its
    inputs, and it carries the index as an ordinary column when it was written
    with ``index=False`` - the same promotion ``ReservoirRouting._read_indexed``
    does.
    """
    path = str(path)
    if path.lower().endswith('.parquet'):
        frame = pd.read_parquet(path)
        if frame.index.name is None and frame.columns.size:
            first = frame.columns[0]
            if str(first).startswith('Unnamed') or str(first) in ('', 'index'):
                frame = frame.set_index(first)
                frame.index.name = None
        return frame
    return pd.read_csv(path, index_col=0)


def read_hydrographs(path) -> pd.DataFrame:
    """A stored hydrograph file, indexed by time, one column per simulation.

    Bryan only ever *writes* these as csv, but a reservoir routing row's
    inflows are not written by the run at all - they are its input, taken
    verbatim from the sims-list ``Inflow`` column, and that file is routinely
    ``.parquet`` (``ReservoirRouting._read_inflows``). Reading it with
    ``read_csv`` gets a UnicodeDecodeError on the parquet magic, which is what
    used to end the whole plotting run.

    The index promotion is unconditional here, unlike ``load_mcdf`` above:
    parquet written with ``index=False`` carries the time axis as an ordinary
    first column, and a hydrograph file's first column is always the time.
    """
    path = str(path)
    if path.lower().endswith('.parquet'):
        frame = pd.read_parquet(path)
        if frame.index.name is None and frame.columns.size:
            frame = frame.set_index(frame.columns[0])
            frame.index.name = None
        return frame
    return pd.read_csv(path, index_col=0)


def sim_label(sim_id) -> str:
    """The hydrograph column for a simulation id: 'sim_00042'."""
    return 'sim_{}'.format(str(int(sim_id)).zfill(5))


# Where a run leaves the series behind, by the code that writes them:
#   URBS  UrbsModel.store_hydrographs       <run folder>/../urbs_results/<name>_<kind>.csv
#   RORB  RorbModel.store_hydrographs       <run folder>/../rorb_results/<name>_<kind>.csv
#   routing ReservoirRouting._write_hydrographs
#                                <Hydrographs folder>/<name>_<kind><suffix>.csv
# Reservoir routing writes no inflows - they are its input, so the sims-list
# 'Inflow' column is where those come from. RORB writes no levels.
HYDROGRAPH_KINDS = ('inflows', 'levels', 'outflows')

MODEL_FOLDERS = {'urbs': 'urbs_results', 'rorb': 'rorb_results'}


def hydrograph_paths(output_file, model='urbs', hydrographs_folder=None,
                     suffix='', inflow=None) -> dict:
    """Where this run's stored hydrographs are, by kind.

    Layout only - nothing is opened, and a run with ``Store hydrographs`` off
    has none of them. ``output_file`` is the sims-list value, so the folders
    resolve the way Bryan's own writers resolve them.
    """
    output_file = str(output_file).replace('\\', os.sep).replace('/', os.sep)
    name = os.path.basename(output_file)

    if hydrographs_folder:
        # Reservoir routing: everything but the inflows, which it did not write.
        tag = f'_{suffix}' if suffix and not str(suffix).startswith('_') else (suffix or '')
        found = {kind: os.path.join(str(hydrographs_folder), f'{name}_{kind}{tag}.csv')
                 for kind in ('levels', 'outflows')}
        if inflow:
            found['inflows'] = str(inflow)
        return found

    folder = os.path.join(os.path.dirname(output_file),
                          '..', MODEL_FOLDERS.get(str(model).lower(), 'urbs_results'))
    folder = os.path.normpath(folder)
    kinds = HYDROGRAPH_KINDS if str(model).lower() != 'rorb' else ('inflows', 'outflows')
    return {kind: os.path.join(folder, f'{name}_{kind}.csv') for kind in kinds}


def subburst_durations(frame) -> list:
    """The sub-durations this mcdf recorded, shortest first.

    Empty for a database written before the sub-burst columns existed, which is
    not an error - the embedded burst comment still says whether the pattern
    carried one, just not by how much.
    """
    found = []
    for name in frame.columns:
        match = SUBBURST_RE.match(str(name))
        if match and f'ifd_{match.group(1)}h' in frame.columns:
            found.append(float(match.group(1)))
    return sorted(found)


# -- turning a lake level into an AEP ----------------------------------------

@dataclass(frozen=True)
class LevelLookup:
    """A target lake level, and the AEP the frequency curve puts it at."""

    level: float
    aep: float | None = None
    note: str = ''
    above_curve: bool = False

    @property
    def found(self) -> bool:
        return self.aep is not None


def aep_for_level(curve, level) -> LevelLookup:
    """Read a lake level off a level frequency curve.

    ``curve`` is a series of level indexed by AEP ('1 in X') - a quantile table
    as ``ui/core/results.read_curve`` returns it, or the envelope over several
    durations, which is the design curve proper.

    Interpolation is (log level, z), as util/GetRepresentativeEvents.py does it:
    both axes are then close to straight, so a linear read between the standard
    AEPs is a fair one.

    A level **above** the curve is not an error and not something to invent an
    AEP for - it is the case where the dam has been asked about a loading rarer
    than anything simulated. It comes back with ``above_curve`` set, and
    ``rank`` then offers the highest-level events instead of the closest ones.
    """
    try:
        level = float(level)
    except (TypeError, ValueError):
        return LevelLookup(level=math.nan, note='not a number')

    points = []
    for aep, value in dict(curve).items():
        z = normal_variate(aep)
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if math.isnan(z) or math.isnan(value) or value <= 0:
            continue
        points.append((math.log10(value), z))
    points.sort()

    if len(points) < 2:
        return LevelLookup(level=level, note='the level curve has fewer than two points')

    target = math.log10(level) if level > 0 else math.nan
    if math.isnan(target):
        return LevelLookup(level=level, note='a level must be above zero')

    lowest, highest = points[0][0], points[-1][0]
    if target > highest:
        return LevelLookup(
            level=level, above_curve=True,
            note=f'{level:g} m is above the top of the curve '
                 f'({10 ** highest:g} m); ranking by the highest levels reached')
    if target < lowest:
        return LevelLookup(
            level=level,
            note=f'{level:g} m is below the bottom of the curve ({10 ** lowest:g} m)')

    for (x0, z0), (x1, z1) in zip(points, points[1:]):
        if x0 <= target <= x1:
            z = z0 if x1 == x0 else z0 + (z1 - z0) * (target - x0) / (x1 - x0)
            aep = aep_of_variate(z)
            if not math.isfinite(aep):
                # The upper tail has underflowed - the loading is rarer than
                # the scheme can put a number on, which is the above-the-curve
                # case however the interpolation got there.
                return LevelLookup(
                    level=level, above_curve=True,
                    note=f'{level:g} m is rarer than an AEP can express here; '
                         f'ranking by the highest levels reached')
            return LevelLookup(level=level, aep=aep)
    return LevelLookup(level=level, note='could not interpolate the AEP')


def value_for_aep(curve, aep) -> float:
    """The design value a frequency curve gives at an AEP - the inverse of
    ``aep_for_level``, on the same (log value, z) interpolation.

    This is what lets a **design AEP** loading be ranked on the result itself:
    the loading is quoted as a 1 in X, but what the event has to hit is the
    lake level that goes with it. NaN where the curve cannot answer, which the
    caller reports rather than silently ranking on something else.
    """
    z_target = normal_variate(aep)
    if math.isnan(z_target):
        return math.nan

    points = []
    for curve_aep, value in dict(curve).items():
        z = normal_variate(curve_aep)
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if math.isnan(z) or math.isnan(value) or value <= 0:
            continue
        points.append((z, math.log10(value)))
    points.sort()
    if len(points) < 2:
        return math.nan

    if z_target <= points[0][0]:
        return 10 ** points[0][1] if z_target == points[0][0] else math.nan
    if z_target >= points[-1][0]:
        return 10 ** points[-1][1] if z_target == points[-1][0] else math.nan
    for (z0, v0), (z1, v1) in zip(points, points[1:]):
        if z0 <= z_target <= z1:
            if z1 == z0:
                return 10 ** v0
            return 10 ** (v0 + (v1 - v0) * (z_target - z0) / (z1 - z0))
    return math.nan


# -- the metrics -------------------------------------------------------------

def prepare(frame, result_type: str) -> pd.DataFrame:
    """Everything about a realisation that does not depend on the target.

    Split out from ``score`` because it is the expensive half - a z per row -
    and it is the same whatever loading is being matched, so the caller can
    compute it once per (database, result type) and re-rank as often as it
    likes.
    """
    if result_type not in RESULT_TYPES:
        raise ValueError(f'{result_type!r} is not one of {RESULT_TYPES}')
    column = f'{result_type}_aep'
    if column not in frame.columns:
        raise ValueError(
            f"the database has no {column!r} column - it holds the sampling but "
            f"not the analysis, so run the row with 'Analyse results' set")

    out = frame.copy()

    # The TPT columns are probabilities; rain_aep is 1 in X. Never mix them.
    probability = pd.to_numeric(out[column], errors='coerce')
    out['result_aep'] = probability.where(probability > 0).rdiv(1.0)
    out['z_result'] = probability.map(variate_of_probability)

    if 'rain_z' in out.columns:
        out['z_rain'] = pd.to_numeric(out['rain_z'], errors='coerce')
    else:
        out['z_rain'] = pd.to_numeric(
            out.get('rain_aep'), errors='coerce').map(normal_variate)

    # The result in its own units - metres for a level, m3/s for a flow - so a
    # loading can be ranked on what it actually has to hit as well as on its AEP.
    out['result_value'] = pd.to_numeric(out.get(result_type), errors='coerce')

    out['subburst_ratio'] = _subburst_ratio(out)
    out['preburst_offset'] = _offset(out, 'preburst_p')
    out['il_offset'] = _offset(out, 'il_p')
    out['cl_offset'] = _offset(out, 'cl_p')
    if 'lake_z' in out.columns:
        out['lake_z'] = pd.to_numeric(out['lake_z'], errors='coerce')
    return out


def _offset(frame, column) -> pd.Series:
    """How far a sampled percentile sits from the median."""
    if column not in frame.columns:
        return pd.Series(math.nan, index=frame.index)
    return pd.to_numeric(frame[column], errors='coerce') - 0.5


def _subburst_ratio(frame) -> pd.Series:
    """The worst sub-burst in each realisation, as a ratio of the IFD depth.

    Above 1.0 the pattern contains a window rarer than the storm around it -
    the embedded burst, measured rather than described. NaN where the database
    predates the sub-burst columns; that is a gap in what is known, not a clean
    storm, and the flags say so.
    """
    ratios = []
    for duration in subburst_durations(frame):
        measured = pd.to_numeric(frame[f'subburst_{duration:g}h'], errors='coerce')
        reference = pd.to_numeric(frame[f'ifd_{duration:g}h'], errors='coerce')
        ratios.append(measured / reference.where(reference > 0))
    if not ratios:
        return pd.Series(math.nan, index=frame.index)
    return pd.concat(ratios, axis=1).max(axis=1, skipna=True)


def score(prepared, target_aep, rain_aep=None, target_value=None) -> pd.DataFrame:
    """Distance from one target. Cheap - call it per loading.

    ``rain_aep`` overrides the rainfall AEP the event is judged against. It
    defaults to the target, which is the AEP-neutral case; setting it asks for
    events whose rainfall is deliberately more or less rare than their flood.

    ``target_value`` is the loading in the result's own units - the lake level
    a loading of '220.5 m AHD' names outright, or the level the design curve
    gives at a '1 in 2,000'. It fills ``d_value``/``delta_value``, which
    ``rank(order='result')`` sorts on: hitting the level is often what the
    event is for, and AEP neutrality the thing traded against it.
    """
    target_aep = float(target_aep)
    rain_target = float(rain_aep) if rain_aep else target_aep
    z_target = normal_variate(target_aep)
    z_rain_target = normal_variate(rain_target)

    out = prepared.copy()
    out['d_z_result'] = out['z_result'] - z_target
    out['d_z_rain'] = out['z_rain'] - z_rain_target
    out['delta_z'] = (out['d_z_result'] ** 2 + out['d_z_rain'] ** 2) ** 0.5

    # The measure util/GetRepresentativeEvents.py sorts on, kept so a selection
    # made with it can be compared. Dominated by the rare end - see the module
    # docstring for why it is not the rank key.
    d_result = out['result_aep'] - target_aep
    d_rain = pd.to_numeric(out.get('rain_aep'), errors='coerce') - rain_target
    out['delta_aep'] = (d_result ** 2 + d_rain ** 2) ** 0.5

    values = pd.to_numeric(out.get('result_value'), errors='coerce')
    if target_value is None or values is None or _isnan(target_value):
        out['d_value'] = pd.Series(math.nan, index=out.index)
    else:
        out['d_value'] = values - float(target_value)
    out['delta_value'] = out['d_value'].abs()
    return out


# -- filtering and ranking ---------------------------------------------------

@dataclass(frozen=True)
class Filters:
    """What disqualifies an event, and what merely gets said about it.

    Everything except the PMP cap and the distance limit **flags** by default
    and excludes only when asked. A tool that silently drops the event you were
    looking for is worse than no tool: the flags are the point, and the choice
    between a close AEP match and a clean storm belongs to whoever has to
    defend the event.
    """

    aep_of_pmp: float | None = None
    pmp_factor: float = 1.1          # as util/GetRepresentativeEvents.py caps it
    max_delta_z: float | None = None
    subburst_limit: float = 1.0      # a window rarer than the storm around it
    preburst_band: float = 0.25      # |preburst_p - 0.5|
    lake_band: float = 1.0           # |lake_z|, so about 1 in 6 either way
    loss_band: float = 0.35          # |il_p - 0.5| and |cl_p - 0.5|
    exclude_embedded: bool = False
    exclude_flagged: bool = False


@dataclass
class Ranking:
    """One target's answer: the candidates, and what was left out on the way."""

    candidates: pd.DataFrame = field(default_factory=pd.DataFrame)
    excluded: dict = field(default_factory=dict)     # reason -> how many
    notes: list = field(default_factory=list)
    considered: int = 0

    @property
    def is_empty(self) -> bool:
        return self.candidates.empty


def _by_result(frame, band: float = 0.0, rounding: float = 0.0) -> tuple:
    """Sorted by how close the result is to the loading, ties on ``delta_z``.

    ``band`` and ``rounding`` are both in the result's own units, and both
    exist because a lake level read to the millimetre is false precision - the
    rating curve, the routing timestep and the model itself are nowhere near
    that, so an unbanded sort is decided by noise and neutrality never gets a
    look in. ``band`` is how far off still counts as reaching the loading;
    ``rounding`` is the grid everything further out is measured on. Within
    either group the order is ``delta_z``. See ``banded``.

    The loading's own units where they are known. Where they are not - no
    design value for that AEP, or a database with no such column - the result
    axis in z is the same order for the same reason, and the note says so
    rather than leaving the table looking like it sorted on the level.
    """
    distance = pd.to_numeric(frame.get('delta_value'), errors='coerce')
    if distance is not None and distance.notna().any():
        return frame.assign(_distance=banded(distance, band, rounding)).sort_values(
            by=['_distance', 'delta_z']).drop(columns='_distance'), ''
    return (frame.assign(_distance=frame['d_z_result'].abs())
            .sort_values(by=['_distance', 'delta_z']).drop(columns='_distance'),
            'Ranked on the result AEP: there is no design value at this '
            'loading to measure against.')


def banded(distance, band: float = 0.0, rounding: float = 0.0):
    """The sort key that decides which results count as the same.

    Two separate things, because they answer two questions:

    - ``band`` is how far from the loading still counts as **reaching** it.
      Everything inside it is one group - key 0 - and the order within the
      group is left to whatever the caller sorts on next, which is neutrality.
    - ``rounding`` is the grid the rest are measured on. Two events that round
      to the same difference are the same distance away as far as anyone can
      defend, so they tie and neutrality separates them too.

    Either can be zero on its own: no band and the closest event still leads;
    no rounding and the events outside the band keep their exact order. NaN
    stays NaN - an unknown result is not a match, banded or otherwise.
    """
    band, rounding = _positive(band), _positive(rounding)
    key = distance
    if rounding > 0:
        key = (distance / rounding).apply(
            lambda value: value if _isnan(value) else math.floor(value + 0.5))
    if band > 0:
        inside = distance <= band
        if rounding > 0:
            # Outside the band is never key 0, however coarse the rounding.
            key = key.mask(key < 1, 1)
        key = key.mask(inside.fillna(False), 0)
    return key


def _positive(value) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 0.0
    return value if value > 0 and math.isfinite(value) else 0.0


def flags_for(row, filters: Filters = Filters()) -> tuple:
    """Everything worth saying about one realisation, in words."""
    out = []

    comment = str(row.get('embedded_bursts') or '').strip()
    if comment and not comment.startswith(NO_EMBEDDED_BURSTS):
        out.append(comment)

    ratio = row.get('subburst_ratio')
    if ratio is not None and not _isnan(ratio) and ratio > filters.subburst_limit:
        out.append(f'sub-burst {ratio:.2f}x the IFD depth')

    preburst = row.get('preburst_offset')
    if preburst is not None and not _isnan(preburst) and abs(preburst) > filters.preburst_band:
        out.append(f'pre-burst percentile {preburst + 0.5:.2f}')

    lake = row.get('lake_z')
    if lake is not None and not _isnan(lake) and abs(lake) > filters.lake_band:
        which = 'high' if lake > 0 else 'low'
        out.append(f'{which} antecedent storage (z = {lake:+.2f})')

    for key, label in (('il_offset', 'initial loss'), ('cl_offset', 'continuing loss')):
        value = row.get(key)
        if value is not None and not _isnan(value) and abs(value) > filters.loss_band:
            out.append(f'{label} percentile {value + 0.5:.2f}')

    filtered = str(row.get('preburst_filter') or '').strip()
    if filtered.startswith('ERROR'):
        out.append(filtered)
    return tuple(out)


def has_embedded_burst(row, filters: Filters = Filters()) -> bool:
    """Whether this realisation's main burst carries an embedded burst.

    Asked of the data, never of the flag text: 'pre-burst percentile 0.88'
    contains the word 'burst' and has nothing to do with an embedded one, which
    is exactly the kind of thing a substring match gets wrong and a reviewer
    never notices.

    Two independent records, either of which counts: the comment
    ``StormGenerator`` wrote, and the measured sub-burst ratio.
    """
    comment = str(row.get('embedded_bursts') or '').strip()
    if comment and not comment.startswith(NO_EMBEDDED_BURSTS):
        return True
    ratio = row.get('subburst_ratio')
    return not _isnan(ratio) and float(ratio) > filters.subburst_limit


def embedded_mask(frame, filters: Filters = Filters()) -> pd.Series:
    """``has_embedded_burst`` over a whole frame, vectorised.

    ``rank`` needs this over the whole database, and a row-wise ``apply`` on an
    mcdf is a second or two - on a page that re-ranks whenever a control moves.
    Pinned against the scalar version by a test, since two implementations of
    one rule is what this module exists to avoid elsewhere.
    """
    comment = frame.get('embedded_bursts')
    if comment is None:
        flagged = pd.Series(False, index=frame.index)
    else:
        text = comment.astype('string').fillna('').str.strip()
        flagged = (text != '') & ~text.str.startswith(NO_EMBEDDED_BURSTS)
    ratio = pd.to_numeric(frame.get('subburst_ratio'), errors='coerce')
    if ratio is not None:
        flagged = flagged | (ratio > filters.subburst_limit).fillna(False)
    return flagged.astype(bool)


def flag_mask(frame, filters: Filters = Filters()) -> pd.Series:
    """Whether each realisation has anything flagged at all, vectorised.

    The counterpart of ``flags_for`` for the same reason as ``embedded_mask``,
    and pinned against it the same way.
    """
    flagged = embedded_mask(frame, filters)
    for column, band in (('preburst_offset', filters.preburst_band),
                         ('lake_z', filters.lake_band),
                         ('il_offset', filters.loss_band),
                         ('cl_offset', filters.loss_band)):
        values = pd.to_numeric(frame.get(column), errors='coerce')
        if values is not None:
            flagged = flagged | (values.abs() > band).fillna(False)

    errors = frame.get('preburst_filter')
    if errors is not None:
        text = errors.astype('string').fillna('').str.strip()
        flagged = flagged | text.str.startswith('ERROR')
    return flagged.astype(bool)


def _isnan(value) -> bool:
    try:
        return math.isnan(float(value))
    except (TypeError, ValueError):
        return True


def rank(scored, filters: Filters = Filters(), count: int = 10,
         above_curve: bool = False, order: str = DELTA_Z,
         band: float = 0.0, rounding: float = 0.0) -> Ranking:
    """The best candidates for one target, and an account of the rest.

    ``order`` is what closeness means here:

    - ``'delta_z'`` (the default) is the distance on both axes at once - the
      event should reach the loading *and* be AEP neutral about it.
    - ``'result'`` sorts on the result alone: the difference from the loading
      in its own units (``delta_value``) where the target value is known, and
      failing that the distance on the result axis in z. Neutrality is then
      reported rather than ranked on, which is the right way round when the
      event exists to reach a particular lake level and the rainfall's rarity
      is a thing to be checked afterwards. ``delta_z`` breaks its ties, so of
      two events at the same level the neutral one still comes first, and
      ``band``/``rounding`` are what 'the same level' means - see ``banded``.

    ``above_curve`` comes from ``aep_for_level``: the loading is rarer than
    anything simulated, so there is no distance to minimise and the honest
    answer is the events that got highest, ordered by how rare their rainfall
    was.
    """
    result = Ranking(considered=len(scored))
    frame = scored

    if filters.aep_of_pmp:
        cap = float(filters.aep_of_pmp) * filters.pmp_factor
        keep = pd.to_numeric(frame.get('rain_aep'), errors='coerce') < cap
        dropped = int((~keep).sum())
        if dropped:
            result.excluded[f'rainfall rarer than 1 in {cap:,.0f}'] = dropped
        frame = frame[keep.fillna(False)]

    if filters.max_delta_z is not None and not above_curve:
        keep = frame['delta_z'] <= float(filters.max_delta_z)
        dropped = int((~keep).sum())
        if dropped:
            result.excluded[f'further than {filters.max_delta_z:g} in z'] = dropped
        frame = frame[keep.fillna(False)]

    if filters.exclude_flagged or filters.exclude_embedded:
        if filters.exclude_flagged:
            keep = ~flag_mask(frame, filters)
            reason = 'flagged'
        else:
            keep = ~embedded_mask(frame, filters)
            reason = 'carrying an embedded burst'
        dropped = int((~keep).sum())
        if dropped:
            result.excluded[reason] = dropped
        frame = frame[keep]

    if frame.empty:
        result.notes.append('Nothing is left once the filters are applied.')
        result.candidates = frame.assign(flags=pd.Series(dtype=object))
        return result

    if above_curve:
        frame = frame.sort_values(by=['level', 'rain_aep'], ascending=False)
    elif order == RESULT:
        frame, note = _by_result(frame, band, rounding)
        if note:
            result.notes.append(note)
    else:
        frame = frame.sort_values(by=['delta_z', 'd_z_result'])

    # Only now, on the handful being offered: writing the flags out in words is
    # string formatting per row, and the database it started from is m x n.
    candidates = frame.head(max(int(count), 1)).copy()
    candidates['flags'] = [flags_for(row, filters)
                           for _, row in candidates.iterrows()]
    result.candidates = candidates
    return result


# -- targets, and the file they live in --------------------------------------

@dataclass
class Target:
    """One loading to find an event for.

    ``kind`` is 'aep' - a design AEP - or 'level', a lake level read off the
    level frequency curve first. ``source`` names the sims-list row whose
    database the event comes from; the UI fills it in with the row whose
    duration is critical at that AEP.
    """

    kind: str = 'aep'
    value: float = 0.0
    result_type: str = 'level'
    rain_aep: float | None = None
    source: str = ''
    output_file: str = ''      # the sims-list value, so the CLI can find the row
    database: str = ''         # the mcdf the event was chosen from
    count: int = 10
    picked: int | None = None
    comment: str = ''

    @property
    def label(self) -> str:
        if self.kind == 'level':
            return f'{self.value:g} m AHD'
        return f'1 in {self.value:,.0f}'

    def to_dict(self) -> dict:
        return {
            'kind': self.kind, 'value': self.value,
            'result_type': self.result_type, 'rain_aep': self.rain_aep,
            'source': self.source, 'output_file': self.output_file,
            'database': self.database, 'count': self.count,
            'picked': self.picked, 'comment': self.comment,
        }

    @classmethod
    def from_dict(cls, data) -> 'Target':
        known = {key: data.get(key, getattr(cls, key, None))
                 for key in cls.__dataclass_fields__}
        known['value'] = float(known.get('value') or 0.0)
        known['count'] = int(known.get('count') or 10)
        picked = known.get('picked')
        known['picked'] = None if picked in (None, '') else int(picked)
        return cls(**known)


SELECTION_FILE = 'representative_events.json'


def read_selection(path) -> tuple:
    """(targets, settings) from a saved selection file. Missing is empty."""
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return [], {}
    targets = [Target.from_dict(item) for item in data.get('targets', [])]
    return targets, dict(data.get('settings', {}))


def selection_payload(targets, settings=None) -> dict:
    return {
        'targets': [target.to_dict() for target in targets],
        'settings': dict(settings or {}),
    }


# -- the AEP of the PMP ------------------------------------------------------

def pmp_aep_from_storm_config(storm_config):
    """The AEP of the PMP, from the storm config's IFD files config.

    Where a Monte Carlo run gets it: ``StormBurst.get_aep_of_pmp`` ->
    ``ifdCurves.config_data['AEP_of_PMP']``, so the number lives in the IFD
    files config that the storm config's ``rare_ifds`` points at, not in the
    method config file. ``ReservoirRouting._aep_of_pmp_from_storm_config``
    walks the same chain for the same reason.

    Best effort, like that one: a missing or broken chain returns None and the
    caller carries on without the cap rather than failing.
    """
    if not storm_config:
        return None
    try:
        with open(storm_config) as handle:
            storm_data = json.load(handle)
        rare_ifds = storm_data['file_paths']['rare_ifds']

        # StormBurst.set_filepaths resolves against the storm config's folder,
        # and configs are written with either separator.
        folder = os.path.dirname(str(storm_config))
        ifd_path = os.path.normpath(os.path.join(folder, rare_ifds))
        if not os.path.isfile(ifd_path):
            swapped = rare_ifds.replace('\\', os.sep).replace('/', os.sep)
            ifd_path = os.path.normpath(os.path.join(folder, swapped))

        with open(ifd_path) as handle:
            return float(json.load(handle)['AEP_of_PMP'])
    except (OSError, KeyError, ValueError, TypeError):
        return None
