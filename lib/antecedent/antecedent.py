"""Rainfall-driven antecedent storage for Callide Dam.

This is the method of section 8.3 of the design hydrology report ("Antecedent
Storage Conditions"), after Criste-Jones et al. (2025), driven by the Python
homogenisation model rather than the GoldSim water balance the earlier notebook
used.

For each water year:

1. Take the peak lake volume and the instant it occurs, from the homogenised
   record's annual maximum series (sub-daily -- the peak is placed to the hour).
2. Look back up to ``WINDOW_DAYS`` before the peak for the main rainfall burst.
   For each storm duration of 1..5 days, form the largest forward accumulation
   of catchment-average AWAP rainfall, apply the AWAP restricted->unrestricted
   scaling, and read an AEP off the Sunwater IFD depth-AEP curve for that
   duration.  The critical duration is the one whose accumulation is rarest.
3. A year yields an antecedent-storage sample only if at least one duration's
   scaled accumulation exceeds ``THRESHOLD_FRACTION`` of the 1-in-2 AEP depth
   for that duration -- otherwise the annual maximum was not produced by a
   significant storm and the year is dropped.
4. Antecedent storage is the lake volume on the day the main burst started
   (``adv_burst_vol``).  A second sample, ``adv_preburst_vol``, is the volume at
   the start of the pre-burst rainfall (walking back while daily rain exceeds
   ``RAINFALL_THRESHOLD_MM``), for the case where excess pre-burst is prepended
   to design storms instead of excluded.  The report adopts the burst-start
   volume; both are reported.
5. Cunnane plotting positions rank the retained samples on each volume.

The lake volume comes from the homogenised model at the routing resolution,
snapshotted at 09:00 each day to match AWAP's 9am rainfall-day boundary.
"""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

log = logging.getLogger("bryan.antecedent")

# Storm durations carried, in days.  AWAP rainfall is daily and the Sunwater IFD
# stops at 120 h, so 1..5 days.
DURATIONS_D = (1, 2, 3, 4, 5)

# AWAP fixed 9am-9am day-window correction (restricted -> unrestricted), carried
# from the earlier ADV notebook.  The bias is largest at 1 day and negligible by
# 5 days.  Applied only when reading an AEP, not when reporting the raw depth.
RESTRICTION = {1: 1.15, 2: 1.11, 3: 1.07, 4: 1.05, 5: 1.04}

WINDOW_DAYS = 30            # days before the annual peak searched for the burst
THRESHOLD_FRACTION = 0.8    # of the 1-in-2 depth: the significance test
RAINFALL_THRESHOLD_MM = 10  # daily rain marking the edge of the pre-burst
CUNNANE_A = 0.4             # Cunnane plotting-position parameter


def load_rainfall(path):
    """CLD catchment-average AWAP daily rainfall as a date-indexed mm series."""
    frame = pd.read_csv(path, parse_dates=["date"])
    return frame.set_index("date")["rain_mm"].sort_index()


def load_ifd(path):
    """Sunwater IFD depths: index = duration (h), columns = AEP (1-in-X floats)."""
    ifd = pd.read_csv(path, index_col="duration_h")
    ifd.columns = [float(re.sub(r"^1\s*in\s*", "", c)) for c in ifd.columns]
    return ifd


def daily_volume_at_9am(result):
    """Homogenised lake volume (ML) snapshotted at 09:00 each calendar day.

    The sub-daily trace is interpolated linearly onto a 09:00 grid, and each
    sample is labelled with the AWAP rain-day it *begins*, so the value for date
    *D* is the lake state entering the rain-day that ends at 09:00 on *D* -- the
    lake as the rain labelled *D* starts falling.  Indexing the 09:00 sample by
    its own date instead would return the state *leaving* that rain-day, a day
    of the event's own rise later, which for a catchment that responds inside a
    day can put the reading after the flood peak: February 2015 peaks at 21:45
    on the 20th under rain labelled the 21st, and reads 216.33 m rather than
    215.46 m if the burst day itself is used.

    Antecedent storage is otherwise slow-moving, so away from the sharpest
    events the one-day choice moves it by a few hundred ML.
    """
    volume = result["New_Volume_ML"]
    days = pd.date_range(volume.index.min().normalize(),
                         volume.index.max().normalize(), freq="D")
    grid = days + pd.Timedelta(hours=9)
    grid = grid[(grid >= volume.index.min()) & (grid <= volume.index.max())]
    interpolated = np.interp(grid.view("int64"),
                             volume.index.view("int64"), volume.to_numpy())
    series = pd.Series(interpolated,
                       index=(grid + pd.Timedelta(days=1)).normalize())
    series.index.name = "date"
    return series


def _aep_of_depths(depths, aep_x, values):
    """Read 1-in-X AEPs off one duration's depth-AEP curve, in log-log space.

    ``depths`` and ``aep_x`` are the IFD depths and their 1-in-X labels for a
    single duration.  Linear in log-log with extrapolation, matching the notebook.
    A non-positive depth (a dry window) has no meaningful AEP; it is reported as
    the most frequent (1 in 0), never critical.

    Vectorised because the whole-record scan reads a quarter of a million depths
    off these curves and building the interpolator once per depth dominated it.
    """
    values = np.asarray(values, dtype=float)
    out = np.zeros(values.shape, dtype=float)
    positive = values > 0
    if positive.any():
        curve = interp1d(np.log(depths), np.log(aep_x), kind="linear",
                         fill_value="extrapolate")
        out[positive] = np.exp(curve(np.log(values[positive])))
    return out


def _aep_of_depth(depths, aep_x, depth):
    """One depth's 1-in-X AEP; see :func:`_aep_of_depths`."""
    return float(_aep_of_depths(depths, aep_x, [depth])[0])


def forward_accumulation(series, d):
    """Rolling ``d``-day forward accumulation: the value at T sums [T, T+d-1].

    Its argmax is therefore the day a d-day burst *starts*, which is what an
    antecedent volume has to be read on.
    """
    return series.rolling(d).sum().shift(-(d - 1))


def burst_frame(series, ifd, durations=None, threshold_fraction=None):
    """Every start date's burst depth, AEP and significance, for each duration.

    ``series`` is any stretch of the daily rainfall: the days around a flood
    peak, or the whole record.  Scoring the stretch handed in, rather than
    scoring the whole record once and slicing afterwards, is what keeps the two
    callers consistent -- an accumulation is only formed where all ``d`` days lie
    inside the stretch.  The peak-conditioned caller therefore has to hand in a
    stretch that runs past the peak; see :func:`critical_burst`.

    Returns a long frame indexed by burst start date with one row per (date,
    duration): ``duration_d``, ``depth_mm`` (raw AWAP), ``scaled_mm`` (after the
    restricted->unrestricted correction), ``aep_1inx``, ``threshold_mm`` and
    ``over``.  Empty if no accumulation can be formed at any duration.
    """
    # Resolved when called, not when defined, so settings.apply() takes effect
    # (Bryan's one change to the callide-fsl-reinstate copy).
    durations = DURATIONS_D if durations is None else durations
    threshold_fraction = THRESHOLD_FRACTION if threshold_fraction is None else threshold_fraction
    aep_x = ifd.columns.to_numpy(dtype=float)
    blocks = []
    for d in durations:
        hours = d * 24
        depths = ifd.loc[hours].to_numpy(dtype=float)
        accum = forward_accumulation(series, d).dropna()
        if accum.empty:
            continue
        scaled = accum.to_numpy(dtype=float) * RESTRICTION[d]
        threshold = threshold_fraction * float(ifd.loc[hours, 2.0])
        blocks.append(pd.DataFrame({
            "duration_d": d,
            "depth_mm": accum.to_numpy(dtype=float),
            "scaled_mm": scaled,
            "aep_1inx": _aep_of_depths(depths, aep_x, scaled),
            "threshold_mm": threshold,
            "over": scaled > threshold,
        }, index=accum.index))
    if not blocks:
        return pd.DataFrame(columns=["duration_d", "depth_mm", "scaled_mm",
                                     "aep_1inx", "threshold_mm", "over"])
    frame = pd.concat(blocks)
    frame.index.name = "burst_date"
    return frame


def critical_burst(rain, peak_date, ifd):
    """Identify the main burst in the 30 days before ``peak_date``.

    The searched stretch runs to ``peak_date`` plus the longest duration, not to
    ``peak_date``.  Two reasons, and the record flood needs both: an
    accumulation is only formed where every day lies inside the stretch, so a
    burst still running when the lake peaks is otherwise invisible; and AWAP
    stamps a rain-day at the 09:00 that *ends* it, so rain falling on the day of
    the peak carries the following day's label.  Callide responds inside a day,
    so for the sharpest events the causative rain is labelled after the peak it
    caused.  February 2015 is the case in point -- 216.9 mm labelled the 21st
    against a peak at 21:45 on the 20th, which scored 0.52 of the significance
    threshold and dropped the largest flood in the record from the sample, where
    the stretch below scores it at 2.15.

    The burst start is still bounded, since a burst beginning after the lake has
    already peaked cannot be the one that filled it; the widened stretch lets
    such a burst be *seen*, not selected.  The bound is the day after the peak
    rather than the peak, because of that same label lag: a burst labelled ``D``
    starts falling at 09:00 on ``D-1``, so ``D <= peak + 1`` is the statement
    that it began before the lake peaked.

    Returns a dict of the per-duration maxima and the critical burst, plus a
    ``qualified`` flag, or ``None`` if no rainfall is available in the window.
    """
    peak_date = pd.Timestamp(peak_date).normalize()
    window = rain.loc[peak_date - pd.Timedelta(days=WINDOW_DAYS):
                      peak_date + pd.Timedelta(days=max(DURATIONS_D))]
    if window.empty:
        return None

    result = {"qualified": False}
    per_duration = {}
    any_over_threshold = False

    frame = burst_frame(window, ifd)
    # Seen but not selected: the stretch runs past the peak so that a burst
    # still in progress forms an accumulation, but its start has to precede the
    # lake's peak to have caused it. See the label-lag note above for the +1.
    frame = frame[frame.index <= peak_date + pd.Timedelta(days=1)]
    for d, block in frame.groupby("duration_d", sort=True):
        at = block["depth_mm"].idxmax()
        row = block.loc[block.index == at].iloc[0]
        any_over_threshold = any_over_threshold or bool(row["over"])
        per_duration[int(d)] = {
            "depth_mm": float(row["depth_mm"]), "start": at,
            "aep_1inx": float(row["aep_1inx"]),
            "scaled_mm": float(row["scaled_mm"]),
            "threshold_mm": float(row["threshold_mm"]),
            "over": bool(row["over"]),
        }

    result["per_duration"] = per_duration
    if not per_duration:
        return result

    # Critical duration = rarest AEP (largest 1-in-X) among all durations, chosen
    # only when at least one duration cleared its significance threshold.
    result["qualified"] = any_over_threshold
    if any_over_threshold:
        crit_d = max(per_duration, key=lambda d: per_duration[d]["aep_1inx"])
        info = per_duration[crit_d]
        result.update({
            "burst_duration_d": crit_d,
            "burst_date": info["start"],
            "burst_depth_mm": info["depth_mm"],
            "burst_aep_1inx": info["aep_1inx"],
        })
    return result


def preburst_start(rain, burst_date):
    """Walk back from the burst start while daily rain stays above the threshold.

    Returns the first dry day before the pre-burst rainfall -- the volume there
    excludes any rise the pre-burst caused.
    """
    date = pd.Timestamp(burst_date).normalize()
    while True:
        date -= pd.Timedelta(days=1)
        if date not in rain.index:
            return date + pd.Timedelta(days=1)
        if rain.loc[date] <= RAINFALL_THRESHOLD_MM:
            return date


def _cunnane(volumes):
    """Cunnane plotting position on a volume series, largest volume ranked 1."""
    rank = volumes.rank(ascending=False)
    n = volumes.notna().sum()
    return rank, (rank - CUNNANE_A) / (n + 1 - 2 * CUNNANE_A)


def antecedent_series(ams, rain, ifd, daily_volume):
    """Build the antecedent-storage table from the AMS and the rainfall/IFD.

    ``ams`` is ``callide.peaks.annual_maxima`` output (indexed by water year with
    ``Period``, ``New_Volume_max_ML`` and ``New_Level_max_at``).  Returns a frame
    with one row per water year, the critical-burst columns, both antecedent
    volumes, and Cunnane plotting positions over the retained (qualified) years.
    """
    rows = []
    for year, peak in ams.iterrows():
        peak_at = pd.Timestamp(peak["New_Level_max_at"])
        row = {
            "WaterYear": year,
            "Period": peak["Period"],
            "Peak_date": peak_at.normalize(),
            "Peak_volume_ML": peak["New_Volume_max_ML"],
            "Peak_level_m": peak["New_Level_max"],
        }
        burst = critical_burst(rain, peak_at, ifd)
        if burst and burst["qualified"]:
            burst_date = burst["burst_date"]
            pre_date = preburst_start(rain, burst_date)
            row.update({
                "qualified": True,
                "burst_duration_d": burst["burst_duration_d"],
                "burst_date": burst_date,
                "burst_depth_mm": round(burst["burst_depth_mm"], 1),
                "burst_aep_1inx": round(burst["burst_aep_1inx"], 1),
                "preburst_start_date": pre_date,
                "adv_burst_vol_ML": _lookup(daily_volume, burst_date),
                "adv_preburst_vol_ML": _lookup(daily_volume, pre_date),
            })
        else:
            row["qualified"] = False
        rows.append(row)

    table = pd.DataFrame(rows).set_index("WaterYear")

    kept = table[table["qualified"]].copy()
    if len(kept):
        for col, stem in [("adv_burst_vol_ML", "adv_burst"),
                          ("adv_preburst_vol_ML", "adv_preburst")]:
            rank, pp = _cunnane(kept[col])
            table.loc[kept.index, f"rank_{stem}"] = rank
            table.loc[kept.index, f"cunnane_{stem}"] = pp.round(4)

    log.info("Antecedent storage: %d of %d water years yield a sample "
             "(%d dropped, no significant storm within %d days)",
             int(table["qualified"].sum()), len(table),
             int((~table["qualified"]).sum()), WINDOW_DAYS)
    return table


def _lookup(daily_volume, date):
    date = pd.Timestamp(date).normalize()
    if date in daily_volume.index:
        return round(float(daily_volume.loc[date]), 1)
    return np.nan


def antecedent_spill_check(result, table, fsl):
    """Did the lake reach full supply during any event's antecedent window?

    The gate operation only releases above full supply, so it can only affect an
    antecedent sample if the lake spilled between the pre-burst start and the
    moment the sample is read.  Returns the count of events that did and the
    highest level reached in any window, at the routing resolution of ``result``.

    The window closes at 09:00 on the day before ``burst_date``, which is where
    :func:`daily_volume_at_9am` reads the antecedent -- entering the rain-day the
    burst is labelled with.  Closing it at ``burst_date`` instead would run a day
    of the event's own rise past the reading, and for February 2015, whose burst
    is labelled the day after the peak it caused, it would sweep in the 217.29 m
    flood crest and report the sample as having spilled.
    """
    level = result["New_Level"]
    kept = table[table["qualified"]]
    maxima = []
    for _, row in kept.iterrows():
        window = level.loc[pd.Timestamp(row["preburst_start_date"]):
                           pd.Timestamp(row["burst_date"])
                           - pd.Timedelta(days=1) + pd.Timedelta(hours=9)]
        if len(window):
            maxima.append(float(window.max()))
    maxima = np.asarray(maxima, dtype=float)
    return {
        "n": len(maxima), "fsl": float(fsl),
        "n_spilled": int((maxima > fsl).sum()),
        "max_level": float(maxima.max()) if len(maxima) else float("nan"),
    }
