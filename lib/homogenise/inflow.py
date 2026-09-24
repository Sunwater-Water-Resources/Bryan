"""The inflow record: the net inflow stage 1 derives, as a series and as hydrographs.

Stage 1 of the homogenisation (``model.LakeModel.derive_inflow``) is a reverse
routing - inflow = change in storage + release + evaporation, with the release
from the rating in force at each step. This module turns it into what the inflow
flood frequency analysis and model calibration need:

* **the annual maximum series** - the peak inflow (averaged over a window, since
  a millimetre of gauge over a one-minute interval is hundreds of m3/s of noise)
  and the largest inflow volume over each of several durations, with each volume
  as a runoff depth over the catchment and the catchment rainfall over the same
  days beside it;
* **event hydrographs** for calibration and validation, native and smoothed.

Held against callide-fsl-reinstate's independent ``reverse_routing`` package on
Callide (23-24 September 2026): with the evaporation left out, as that package
leaves out every loss, the two agree to a median of 0.00% on the peak and on every
volume window over all 56 water years. What separates them otherwise: stage 1's
ratings are a shape-preserving curve through the published rows where the package
joins them with straight lines, which overstate a convex spillway rating between
rows (by up to 110 m3/s at Callide's 15-row tables) and so the flood volumes it
infers, by 1-4%.

**Recessions.** Above full supply a falling limb often derives *negative*
inflow: the gates were opened wider than the rating says, and the release the
rating misses turns up with an inflow's sign. The recession correction moves
it back into release (``Unaccounted_Release_ML``), which is right for the water
balance and harmless to peaks and burst volumes - they sit on rising limbs - but
leaves the recession of an inflow hydrograph at zero where it was really falling
gradually. Neither version is the true recession, because the actual release is
not recorded. So a hydrograph carries both, ``Inflow_m3s`` (corrected) and
``Inflow_uncorrected_m3s``, with ``Release_uncertain`` marking every interval the
correction touched: use the rising limb and the peak with confidence, and the
flagged part of the recession with care.

The averaging, the cumulative-volume windows and the 15-minute grid for the
burst volumes are the ``reverse_routing`` package's (``routing.py``,
``ams.py``), reimplemented here on stage 1's intervals.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import peaks

ML_TO_M3 = 1000.0
DURATIONS_H = (24, 36, 48, 72)
VOLUME_GRID = pd.Timedelta("15min")


def intervals(derived: pd.DataFrame, evaporation: bool = True,
              recession_correction: bool = True) -> pd.DataFrame:
    """Stage 1's steps as intervals: one row per interval, stamped at its midpoint.

    ``derived`` is ``LakeModel.derive_inflow()``'s frame, whose row at T holds the
    interval ending at T. ``Volume_ML`` is the inflow over the interval; it
    telescopes, so summed it is the storage change plus the release (plus the
    evaporation where kept) to the millilitre.
    """
    rows = derived.iloc[1:]
    start = derived.index[:-1]
    unaccounted = rows.get("Unaccounted_Release_ML", pd.Series(0.0, index=rows.index))
    corrected = rows["Inflow_ML"]
    uncorrected = corrected + unaccounted
    if not evaporation:
        corrected = corrected - rows["Evaporation_ML"]
        uncorrected = uncorrected - rows["Evaporation_ML"]
    volume = corrected if recession_correction else uncorrected
    dt = rows["dt_s"].to_numpy()
    frame = pd.DataFrame({
        "Interval_start": start,
        "dt_s": dt,
        "Volume_ML": volume.to_numpy(),
        "Inflow_native_m3s": volume.to_numpy() * ML_TO_M3 / dt,
        "Inflow_uncorrected_m3s": uncorrected.to_numpy() * ML_TO_M3 / dt,
        "Release_m3s": rows["Release_ML"].to_numpy() * ML_TO_M3 / dt,
        "Level_end": rows["Level"].to_numpy(),
        "Interpolated": rows["Interpolated"].to_numpy() if "Interpolated" in rows else False,
        "Release_uncertain": (unaccounted.to_numpy() > 0),
    }, index=start + (rows.index - start) / 2)
    frame.index.name = "Timestamp"
    frame["Inflow_m3s"] = frame["Inflow_native_m3s"]
    return frame


def _seconds(index, origin) -> np.ndarray:
    return (pd.DatetimeIndex(index) - origin).total_seconds().to_numpy()


def cumulative_volume(frame: pd.DataFrame) -> pd.Series:
    """Running total of inflow volume (ML), stamped at interval ends."""
    ends = pd.DatetimeIndex(frame["Interval_start"] + pd.to_timedelta(frame["dt_s"], unit="s"))
    running = pd.Series(frame["Volume_ML"].cumsum().to_numpy(), index=ends)
    first = pd.Series([0.0], index=pd.DatetimeIndex([frame["Interval_start"].iloc[0]]))
    series = pd.concat([first, running])
    return series[~series.index.duplicated(keep="last")].sort_index()


def smoothed(frame: pd.DataFrame, window) -> pd.Series:
    """Mean inflow (m3/s) over a window centred on each interval, off the cumulative
    volume - the exact time average, so a one-minute gauge step counts for one
    minute rather than as a sample equal to the six-hour interval beside it."""
    window = pd.Timedelta(window)
    cumulative = cumulative_volume(frame)
    origin = cumulative.index[0]
    stamps, volumes = _seconds(cumulative.index, origin), cumulative.to_numpy()
    centre = _seconds(frame.index, origin)
    half = window.total_seconds() / 2.0
    low = np.interp(centre - half, stamps, volumes)
    high = np.interp(centre + half, stamps, volumes)
    return pd.Series((high - low) * ML_TO_M3 / window.total_seconds(), index=frame.index)


def with_smoothing(frame: pd.DataFrame, window) -> pd.DataFrame:
    """``Inflow_m3s`` averaged over ``window``; native where the window is blank."""
    out = frame.copy()
    if window:
        out["Inflow_m3s"] = smoothed(frame, window)
    return out


def burst_volumes(frame: pd.DataFrame, durations_h=DURATIONS_H, start_month=10) -> pd.DataFrame:
    """The largest inflow volume (ML) over each duration, by water year."""
    cumulative = cumulative_volume(frame)
    grid = pd.date_range(cumulative.index[0].ceil(VOLUME_GRID),
                         cumulative.index[-1].floor(VOLUME_GRID), freq=VOLUME_GRID)
    origin = cumulative.index[0]
    running = np.interp(_seconds(grid, origin), _seconds(cumulative.index, origin),
                        cumulative.to_numpy())
    out = {}
    for hours in durations_h:
        steps = int(pd.Timedelta(hours=hours) / VOLUME_GRID)
        volume = pd.Series(running[steps:] - running[:-steps], index=grid[steps:])
        labels = peaks.water_year(volume.index, start_month)
        grouped = volume.groupby(labels)
        out[f"Volume_{hours}h_ML"] = grouped.max()
        out[f"Volume_{hours}h_end"] = grouped.idxmax()
    return pd.DataFrame(out)


def annual_maxima(frame: pd.DataFrame, durations_h=DURATIONS_H, start_month=10,
                  catchment_km2=None, rainfall=None, min_coverage=0.9) -> pd.DataFrame:
    """Peak inflow and burst volumes by water year, with depths and rainfall.

    A volume divided by the catchment area in km2 is a runoff depth in mm (1 mm
    over 1 km2 is 1 ML). With ``rainfall`` (daily ``rain_mm``, day D the 24 h to
    9 am on D), the rain over the days that volume's window spans - plus the day
    before it, since the storm precedes the runoff - is put beside it: a burst
    implying more runoff than rain fell is wrong however well the water balance
    closes, and it is the one check on the inflow that the inflow did not produce.
    """
    labels = peaks.water_year(frame.index, start_month)
    grouped = frame.groupby(labels)
    at = grouped["Inflow_m3s"].idxmax()
    table = pd.DataFrame({
        "Period": [peaks._period_label(int(year), start_month) for year in at.index],
        "Peak_inflow_m3s": grouped["Inflow_m3s"].max(),
        "Peak_time": at,
        "Release_at_peak_m3s": frame.loc[at, "Release_m3s"].to_numpy(),
        "Level_at_peak_m": frame.loc[at, "Level_end"].to_numpy(),
    }, index=at.index)
    table.index.name = "WaterYear"
    table = table.join(burst_volumes(frame, durations_h, start_month))

    span = grouped.apply(lambda block: (block.index.max() - block.index.min()).total_seconds())
    table["Coverage"] = (span / (365.25 * 86400)).clip(upper=1.0)
    table["Complete"] = table["Coverage"] >= min_coverage
    gaps = frame["dt_s"].groupby(labels)
    table["Median_gap_min"] = gaps.median() / 60.0
    table["Max_gap_d"] = gaps.max() / 86400.0

    for hours in durations_h:
        if catchment_km2:
            table[f"Depth_{hours}h_mm"] = table[f"Volume_{hours}h_ML"] / float(catchment_km2)
        if rainfall is not None and len(rainfall):
            ends = pd.to_datetime(table[f"Volume_{hours}h_end"])
            rain = []
            for end in ends:
                last = (end + pd.Timedelta(hours=15)).normalize()   # the 9 am day holding the end
                first = (end - pd.Timedelta(hours=hours)).normalize() - pd.Timedelta(days=1)
                window = rainfall.loc[first:last]
                rain.append(float(window.sum()) if len(window) else np.nan)
            table[f"Rain_{hours}h_mm"] = rain
    return table


# -- hydrographs --------------------------------------------------------------------

def hydrograph(frame: pd.DataFrame, start, end, window="1h") -> pd.DataFrame:
    """One event, for calibration: every interval between ``start`` and ``end``."""
    part = frame.loc[pd.Timestamp(start):pd.Timestamp(end)].copy()
    if part.empty:
        return part
    if window:
        part["Inflow_smoothed_m3s"] = smoothed(frame, window).loc[part.index]
    columns = ["Inflow_native_m3s", "Inflow_uncorrected_m3s", "Release_m3s", "Level_end",
               "Release_uncertain", "Interpolated", "dt_s"]
    if "Inflow_smoothed_m3s" in part:
        columns.insert(1, "Inflow_smoothed_m3s")
    return part[columns].rename(columns={"Level_end": "Level_m",
                                         "Inflow_native_m3s": "Inflow_m3s"})


def event_windows(ams: pd.DataFrame, top: int, before_days=3.0, after_days=7.0) -> list:
    """(name, start, end) around the ``top`` largest annual peaks."""
    chosen = ams.sort_values("Peak_inflow_m3s", ascending=False).head(int(top))
    out = []
    for year, row in chosen.iterrows():
        peak = pd.Timestamp(row["Peak_time"])
        out.append((f"{row['Period']}_peak_{peak:%Y%m%d}",
                    peak - pd.Timedelta(days=before_days), peak + pd.Timedelta(days=after_days)))
    return out
