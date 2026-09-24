"""Read the raw Callide Dam headwater level record.

Two structures gauge the headwater at Callide Dam, at different points::

    130314C   "HW (Dam)"      at the dam wall,   -24.366719, 150.615713
    130314B   "HW (Intake)"   at the intake,     -24.366451, 150.617484  (182 m)

and a third site number carries the pre-1989 record::

    130314A   1970-10-08 -> 1989-09-28   daily, with event bursts at 5-60 min

Both C and B run from 1989-09-28.  For the level record that drives the model
they behave as one measurement almost the whole time: their daily maxima are
bit-identical to 2013 and agree to 0.009 m rmse (max 0.06 m) through mid-2024 --
one connected pool, gauged twice.  The default chain is therefore A then C: C is
the level *at the wall*, which is what governs storage and spillway discharge.

From September 2024 the two separate.  As the lake falls to ~200.9 m the intake
gauge draws down -- 200.88 (Sep) -> 200.16 (Dec 2024) -- while the wall gauge
stays pinned at ~200.94.  A sustained 0.78 m difference over 182 m of open water
is a 0.4% surface slope, impossible in a connected still pool, so below ~200.9 m
the headwater is no longer one pool: a sediment bar (2013/2015 flood alluvium)
between the intake and the wall has emerged, and the wall pocket is held at its
crest while the intake pocket keeps draining under extraction.  Below that crest
the wall gauge no longer reads the working storage; the intake gauge does.  See
:func:`splice_intake`.

The earlier State-portal export made B look like a duplicate that ceased in 2002
-- it dropped B's post-2002 data (a gauge-ownership change) and reproduced a
derived "130314B" column in the daily-maximum export that drifts from the wall.
The full Hydstra record (``130314B_hydstra.csv``) is the real intake gauge and
supersedes both.

The daily-maximum export (``Callide_Historical.csv``) is stamped a day late: the
row dated D holds the maximum of day D-1.  Against the raw record the daily
maxima agree to 0.006 m rmse at that offset and 0.124 m without it.  The raw
exports read here carry true timestamps and need no correction.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("bryan.homogenise.gauges")

DEFAULT_DIRECTORY = Path("HeadWaterLevels")

# In the order they operated.  The chain is the wall record (A then C); the
# intake gauge B is not a time segment of it but an overlay that replaces the
# wall level once the pool partitions below the bar -- see splice_intake.
DEFAULT_CHAIN = ("130314A", "130314C")

# The intake gauge, read from the full Hydstra export.
INTAKE_GAUGE = "130314B_hydstra"

# Level at the dam wall below which the wall gauge no longer reads the working
# storage: the crest of the sediment bar between the intake and the wall, at
# which the wall pocket is held while the intake pocket drains.  Empirical (the
# level at which the two gauges separate); a bathymetric survey would refine it.
DEAD_STORAGE_LEVEL_WALL = 200.90

# Once the intake record ends (2024-12-19) its last value is held while the wall
# stays within this margin of the bar crest; the wall record resumes only when
# an inflow event lifts it clear, reconnecting the pools.
INTAKE_RECONNECT_MARGIN = 0.10

# WMIP quality codes that carry no usable value.  The Callide exports contain
# only 9 (CITEC normal) and 59 (CITEC derived height), so nothing is dropped in
# practice -- these are here so that a future export cannot smuggle a gap in.
UNUSABLE_QUALITY = (151, 255)

# The two date layouts WMIP emits.  A and B use the first, C the second.
_TIMESTAMP_FORMATS = ("%H:%M:%S %d/%m/%Y", "%d/%m/%Y %H:%M")

_HEADER_ROWS = 4


def read_wmip_export(path):
    """Read one WMIP gauge export into a timestamp-indexed frame.

    The export has four header rows, then data, then a two-line licence footer.
    Both the footer and any unparseable row fall out when the timestamp fails to
    parse, but the fraction that fails is checked rather than assumed: a layout
    change that silently dropped half the record would otherwise look like a
    short record.
    """
    path = Path(path)
    header = pd.read_csv(path, nrows=_HEADER_ROWS, header=None, dtype=str)
    labels = [str(v).strip().lower() for v in header.iloc[_HEADER_ROWS - 1]]
    quality_column = labels.index("quality") if "quality" in labels else None

    # More names than columns is harmless (the extras come back all-NaN); fewer
    # would raise.  Eight covers both layouts with room to spare.
    frame = pd.read_csv(
        path, skiprows=_HEADER_ROWS, header=None, names=range(8), dtype=str,
        engine="python",
    )
    stamps = frame[0].astype(str).str.strip()

    parsed, used = None, None
    for fmt in _TIMESTAMP_FORMATS:
        candidate = pd.to_datetime(stamps, format=fmt, errors="coerce")
        if parsed is None or candidate.notna().sum() > parsed.notna().sum():
            parsed, used = candidate, fmt
    # The footer is two lines plus a blank; anything beyond that is a surprise.
    unparsed = int(parsed.isna().sum())
    if unparsed > 5:
        raise ValueError(
            f"{path}: {unparsed} rows have no recognisable timestamp under "
            f"{used!r}; the export layout may have changed"
        )

    data = pd.DataFrame({"Level": pd.to_numeric(frame[1], errors="coerce")})
    if quality_column is not None:
        data["Quality"] = pd.to_numeric(frame[quality_column], errors="coerce")
    data.index = parsed
    data = data[data.index.notna() & data["Level"].notna()].sort_index()
    data.index.name = "Timestamp"

    if "Quality" in data:
        unusable = data["Quality"].isin(UNUSABLE_QUALITY)
        if unusable.any():
            log.warning("%s: dropping %d values with unusable quality codes %s",
                        path.name, int(unusable.sum()),
                        sorted(data.loc[unusable, "Quality"].unique().tolist()))
            data = data[~unusable]

    duplicates = data.index.duplicated(keep="last")
    if duplicates.any():
        log.warning("%s: %d duplicate timestamps, keeping the last of each",
                    path.name, int(duplicates.sum()))
        data = data[~data.index.duplicated(keep="last")]

    gaps = data.index.to_series().diff().dt.total_seconds() / 60.0
    log.info("%s: %d values, %s to %s, median step %.0f min, longest gap %.1f d",
             path.name, len(data), data.index.min(), data.index.max(),
             gaps.median(), gaps.max() / 1440.0)
    return data


def read_gauge_chain(directory=DEFAULT_DIRECTORY, chain=DEFAULT_CHAIN):
    """Splice the gauges in ``chain`` into one continuous headwater record.

    Each gauge owns the half-open interval from its own first reading to the
    first reading of the next gauge in the chain, so an overlap is resolved in
    favour of the earlier gauge and no timestamp is claimed twice.  A and C hand
    over at 1989-09-28 12:00, where both record 203.441 m.
    """
    directory = Path(directory)
    return read_chain([directory / f"{gauge}.csv" for gauge in chain])


def read_chain(paths):
    """Splice gauge exports, given as files in the order the gauges operated.

    The general form of :func:`read_gauge_chain`: each file names its gauge by
    its stem, and owns the record from its first reading to the next file's.
    """
    records, chain = {}, []
    for path in map(Path, paths):
        if not path.exists():
            raise FileNotFoundError(f"gauge export not found: {path}")
        if path.stem in records:
            raise ValueError(f"{path.name} is in the chain twice")
        records[path.stem] = read_wmip_export(path)
        chain.append(path.stem)
    if not chain:
        raise ValueError("no gauge exports given")

    starts = [records[g].index.min() for g in chain]
    if list(starts) != sorted(starts):
        raise ValueError(f"chain {chain} is not in chronological order")

    pieces = []
    for i, gauge in enumerate(chain):
        data = records[gauge]
        end = starts[i + 1] if i + 1 < len(chain) else None
        piece = data.loc[starts[i]:] if end is None else data.loc[starts[i]:end - pd.Timedelta(1, "ns")]
        piece = piece.assign(Gauge=gauge)
        pieces.append(piece)
        log.info("Chain: %s contributes %d values, %s to %s", gauge, len(piece),
                 piece.index.min(), piece.index.max())

    record = pd.concat(pieces).sort_index()
    if record.index.has_duplicates:
        raise ValueError("splicing produced duplicate timestamps")

    handovers = starts[1:]
    for boundary in handovers:
        before = record.loc[:boundary - pd.Timedelta(1, "ns"), "Level"]
        after = record.loc[boundary:, "Level"]
        if len(before) and len(after):
            step = float(after.iloc[0] - before.iloc[-1])
            log.info("Handover at %s: level step %+.3f m", boundary, step)
            if abs(step) > 0.05:
                log.warning("Handover at %s is a %+.3f m jump; the gauges may "
                            "not share a datum", boundary, step)
    return record


def splice_intake(record, intake, dead_storage_level=DEAD_STORAGE_LEVEL_WALL,
                  reconnect_margin=INTAKE_RECONNECT_MARGIN,
                  intake_name=INTAKE_GAUGE):
    """Overlay the intake gauge on the wall record below dead storage.

    The wall record (``record``) reads the working storage only while the
    headwater is one pool.  Below ``dead_storage_level`` -- the crest of the bar
    between the intake and the wall -- the wall pocket is held at the crest and
    the intake gauge (``intake``) reads the level that still matters.  So the
    wall level is kept as primary and replaced by the intake level wherever the
    intake reads below ``dead_storage_level``.

    The intake record ends before the wall's does.  Reverting to the held wall
    level there would inject a step of the full pool difference, so once the
    intake ends its last value is carried forward while the wall stays within
    ``reconnect_margin`` of the crest, and the wall record resumes only when an
    inflow event lifts it clear -- the point at which the pools reconnect and
    the wall reads the working storage again.

    Above ``dead_storage_level`` this is a no-op: over 1989-2024 the two gauges
    agree to 0.06 m, so nothing changes until they separate in late 2024.
    """
    idx = record.index
    ns = idx.astype("int64").to_numpy(dtype=float)
    known = intake.index.astype("int64").to_numpy(dtype=float)
    b = np.interp(ns, known, intake["Level"].to_numpy(), left=np.nan, right=np.nan)

    level = record["Level"].to_numpy().copy()
    gauge = record["Gauge"].astype(object).to_numpy().copy()

    use_b = ~np.isnan(b) & (b < dead_storage_level)
    level[use_b] = b[use_b]
    gauge[use_b] = intake_name

    replaced = int(use_b.sum())
    held = 0
    if use_b.any():
        # After the intake coverage ends, hold its last value until the wall
        # climbs clear of the crest band (the pools reconnect).
        hold_value = float(intake["Level"].iloc[-1])
        past = np.flatnonzero(ns > known[-1])
        if past.size:
            wall = record["Level"].to_numpy()[past]
            recovered = wall > dead_storage_level + reconnect_margin
            stop = int(np.argmax(recovered)) if recovered.any() else past.size
            hold_positions = past[:stop]
            level[hold_positions] = hold_value
            gauge[hold_positions] = f"{intake_name} (held)"
            held = hold_positions.size

    out = record.copy()
    out["Level"] = level
    out["Gauge"] = gauge
    log.info("Intake splice: %d values replaced below %.2f m, %d held after "
             "intake record ends %s", replaced, dead_storage_level, held,
             intake.index.max().date())
    return out


def cap_timestep(record, max_step, level_column="Level"):
    """Insert linearly interpolated points so that no step exceeds ``max_step``.

    Every observed timestamp is kept, so observed peaks survive exactly -- which
    resampling onto a fixed grid would not guarantee, since the February 2015
    peak falls at 21:45 and an hourly grid would straddle and under-read it.
    Only gaps longer than ``max_step`` are subdivided, and the inserted points
    are flagged.

    Interpolating the level linearly (rather than the volume) means the derived
    inflow is not quite constant across a filled gap.  The distinction only
    affects how the gap's total inflow is distributed within it, since the
    storage differences telescope; the level is what was measured, so the level
    is what is interpolated.
    """
    max_step = pd.Timedelta(max_step)
    if max_step <= pd.Timedelta(0):
        raise ValueError("max_step must be positive")

    index = record.index
    gaps = index.to_series().diff()
    inserted = []
    for position in np.flatnonzero((gaps > max_step).to_numpy()):
        start, end = index[position - 1], index[position]
        pieces = int(np.ceil((end - start) / max_step))
        step = (end - start) / pieces
        inserted.extend(start + step * k for k in range(1, pieces))

    if not inserted:
        out = record.copy()
        out["Interpolated"] = False
        return out

    grid = index.union(pd.DatetimeIndex(inserted))
    out = pd.DataFrame(index=grid)
    out.index.name = index.name
    seconds = grid.astype("int64").to_numpy(dtype=float)
    known = index.astype("int64").to_numpy(dtype=float)
    out[level_column] = np.interp(seconds, known, record[level_column].to_numpy())
    for column in record.columns:
        if column != level_column:
            out[column] = record[column].reindex(grid).ffill()
    out["Interpolated"] = ~grid.isin(index)

    log.info("Timestep capped at %s: %d observed + %d interpolated = %d steps",
             max_step, len(index), len(inserted), len(out))
    return out


def long_gaps(record, threshold="2D"):
    """Observed gaps longer than ``threshold``, for reporting."""
    observed = record.index[~record["Interpolated"]] if "Interpolated" in record else record.index
    gaps = observed.to_series().diff()
    long = gaps[gaps > pd.Timedelta(threshold)]
    return pd.DataFrame({"Gap_days": long.dt.total_seconds() / 86400.0,
                         "From": long.index - long}).set_index("From")


def daily_maximum(series):
    """Daily maximum on true calendar days, for comparison with the WMIP export."""
    return series.resample("D").max().dropna()
