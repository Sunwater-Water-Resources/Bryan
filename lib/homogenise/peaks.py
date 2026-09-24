"""Annual maxima and independent peaks from the homogenised record.

Two extractions, and they answer different questions.

The **annual maximum series** takes the largest level and storage in each water
year.  It is the input to a flood frequency analysis and it needs no
independence rule, because one value per year is taken by construction.

**Independent peaks** are for peak-over-threshold work and for counting events,
and there the rule matters.  The 2011 wet season alone puts seventeen separate
maxima within 0.3 m of each other between January and May; treating those as
seventeen floods would be wrong, and treating the whole season as one event
would be wrong too.  The rule used here is the standard one: decluster by runs
above the threshold, then merge neighbouring runs that are either closer than
``separation`` days or not separated by a fall of at least ``drop`` metres.
Both are reported alongside the peaks so the choice is visible rather than
buried.

Water years run October to September, labelled by the calendar year they end in,
so water year 2011 is October 2010 to September 2011.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("bryan.homogenise.peaks")

# 1 October. Chosen on the level series, because that is where the boundary does
# measurable damage to the range the curve is actually read over.
#
# In roughly 40% of years the annual maximum level is inherited -- the level the
# year opened at, carried over from the previous wet season, never rising during
# the year. That is real exposure and belongs in a hazard curve, but it means one
# flood supplies two years' maxima, which strains the independence a frequency
# analysis assumes. What matters is how many of those inherited maxima land high
# enough to matter: above 214 m AHD it is 4 under October against 7 under
# September, the extra three being 1979-80, 1981-82 and 2013-14 at 214.0-214.2 m,
# inside the 20-5% AEP band.
#
# The drawdown is why no month is clean. Median level falls monotonically through
# the dry season and past September -- Aug 208.37, Sep 208.23, Oct 208.07, Nov
# 207.82 -- bottoming in November, by which point the flood season has started.
# Irrigation drawdown runs through spring, so at Callide there is no month that is
# both post-drawdown and pre-flood.
#
# The inflow series would choose 1 September on its own evidence and this costs it
# one mid-ranking year: October splits the late-September 1996 event, taking its
# 225 m3/s shoulder as the 1996-97 maximum while the true 236 m3/s peak falls in
# 1995-96 and is masked. Section 2.6 of the reverse-routing write-up scores all
# twelve candidates on the inflow record and states that cost.
DEFAULT_WATER_YEAR_START = 10
DEFAULT_SEPARATION_DAYS = 5.0
DEFAULT_DROP_M = 0.5


def water_year(index, start_month=DEFAULT_WATER_YEAR_START):
    """Water year label: the calendar year the water year ends in."""
    index = pd.DatetimeIndex(index)
    if start_month == 1:
        return index.year.to_numpy()
    return (index.year + (index.month >= start_month)).to_numpy()


def _period_label(year, start_month):
    if start_month == 1:
        return str(year)
    return f"{year - 1}-{str(year)[-2:]}"


def annual_maxima(result, start_month=DEFAULT_WATER_YEAR_START,
                  min_coverage=0.9):
    """Annual maximum level and storage, recorded and homogenised.

    Storage is a monotone function of level, so the storage maximum falls at the
    same instant as the level maximum; both are reported because both are asked
    for, not because they are independent.

    Water years covered for less than ``min_coverage`` of their length are kept
    but flagged.  Dropping them silently would hide that the record starts on
    1970-10-08 and ends on 2026-02-09, and a partial year's maximum is still the
    true maximum of what was recorded -- it just cannot be assumed to be the
    year's.
    """
    labels = water_year(result.index, start_month)
    grouped = result.groupby(labels)

    rows = []
    for year, block in grouped:
        recorded_at = block["Level"].idxmax()
        new_at = block["New_Level"].idxmax()
        span_start = pd.Timestamp(year=year - (start_month > 1), month=start_month,
                                  day=1)
        span_end = span_start + pd.DateOffset(years=1)
        covered = (block.index.max() - block.index.min()).total_seconds()
        length = (span_end - span_start).total_seconds()
        rows.append({
            "WaterYear": year,
            "Period": _period_label(year, start_month),
            "Coverage": covered / length,
            "Observed_steps": int((~block["Interpolated"]).sum()),
            "FSL": float(block["FSL"].iloc[-1]),
            "Level_max": float(block["Level"].max()),
            "Level_max_at": recorded_at,
            "Volume_max_ML": float(block["Volume_ML"].max()),
            "New_FSL": float(block["New_FSL"].iloc[-1]),
            "New_Level_max": float(block["New_Level"].max()),
            "New_Level_max_at": new_at,
            "New_Volume_max_ML": float(block["New_Volume_ML"].max()),
            "New_Storage_above_FSL_max_ML": float(
                block["New_Storage_above_FSL_ML"].max()),
            "New_Release_peak_m3s": float(block["New_Release_flow_m3s"].max()),
            "Inflow_ML": float(block["Inflow_ML"].sum()),
        })

    ams = pd.DataFrame(rows).set_index("WaterYear")
    ams["Complete"] = ams["Coverage"] >= min_coverage
    incomplete = ams.index[~ams["Complete"]].tolist()
    if incomplete:
        log.info("Water years with less than %.0f%% coverage (kept, flagged): %s",
                 100 * min_coverage, incomplete)
    return ams


def _runs_above(values, threshold):
    """Start and stop positions of each run of ``values`` above ``threshold``."""
    above = values > threshold
    if not above.any():
        return np.empty((0, 2), dtype=int)
    edges = np.diff(above.astype(np.int8))
    starts = np.flatnonzero(edges == 1) + 1
    stops = np.flatnonzero(edges == -1) + 1
    if above[0]:
        starts = np.r_[0, starts]
    if above[-1]:
        stops = np.r_[stops, len(values)]
    return np.column_stack([starts, stops])


def independent_peaks(result, separation_days=DEFAULT_SEPARATION_DAYS,
                      drop=DEFAULT_DROP_M, threshold=None,
                      column="New_Level", start_month=DEFAULT_WATER_YEAR_START):
    """Declustered peaks of ``column`` above ``threshold``.

    Declustering is by runs above the threshold, then a merge pass over
    neighbouring peaks.  Two neighbours are treated as one event when they are
    closer than ``separation_days`` **or** when the level between them never
    falls ``drop`` metres below the lower of the two.  The merge repeats until
    nothing changes, keeping the higher peak of each merged pair.

    Working on runs rather than on every local maximum keeps this linear.  It
    also avoids a trap specific to this rating: ``CALLIDE_RFSL.sq`` releases
    142 m3/s at 10 mm of surcharge, so the homogenised level sits very close to full
    supply for long stretches, and every ripple there is a local maximum.
    """
    values = result[column].to_numpy()
    times = result.index
    if threshold is None:
        threshold = float(result["New_FSL"].iloc[-1])

    runs = _runs_above(values, threshold)
    if not len(runs):
        log.info("No peaks above %.3f m in %s", threshold, column)
        return pd.DataFrame(columns=["Peak", "Peak_at"])

    positions = [start + int(np.argmax(values[start:stop])) for start, stop in runs]
    log.info("%s: %d runs above %.3f m", column, len(positions), threshold)

    merged = True
    while merged and len(positions) > 1:
        merged = False
        kept = [positions[0]]
        for position in positions[1:]:
            previous = kept[-1]
            gap_days = (times[position] - times[previous]).total_seconds() / 86400.0
            trough = float(values[previous:position + 1].min())
            lower = min(values[position], values[previous])
            if gap_days < separation_days or trough > lower - drop:
                merged = True
                if values[position] > values[previous]:
                    kept[-1] = position
            else:
                kept.append(position)
        positions = kept

    index = np.asarray(positions)
    peaks = pd.DataFrame({
        "Peak_at": times[index],
        "Peak": values[index],
        "Level": result["Level"].to_numpy()[index],
        "FSL": result["FSL"].to_numpy()[index],
        "New_Level": result["New_Level"].to_numpy()[index],
        "New_Volume_ML": result["New_Volume_ML"].to_numpy()[index],
        "New_Storage_above_FSL_ML": result["New_Storage_above_FSL_ML"].to_numpy()[index],
        "New_Release_flow_m3s": result["New_Release_flow_m3s"].to_numpy()[index],
        "Interpolated": result["Interpolated"].to_numpy()[index],
    })
    peaks["WaterYear"] = water_year(peaks["Peak_at"], start_month)
    peaks = peaks.sort_values("Peak_at").reset_index(drop=True)
    log.info("%s: %d independent peaks above %.3f m (separation %.1f d, "
             "drop %.2f m)", column, len(peaks), threshold, separation_days, drop)
    return peaks


def format_ams(ams, limit=None):
    """Render the annual maximum series as a fixed-width table."""
    table = ams if limit is None else ams.tail(limit)
    lines = [
        "Annual maximum series (water year October-September)",
        "",
        f"{'Period':>9}  {'FSL':>7}  {'Recorded':>9}  {'FSL':>7}  "
        f"{'Homogenised':>11}  {'Volume ML':>10}  {'Above FSL':>10}  "
        f"{'Release':>8}  {'':>3}",
        f"{'':>9}  {'m AHD':>7}  {'m AHD':>9}  {'m AHD':>7}  {'m AHD':>11}  "
        f"{'':>10}  {'ML':>10}  {'m3/s':>8}  {'':>3}",
        "-" * 88,
    ]
    for _, row in table.iterrows():
        flag = "" if row["Complete"] else " *"
        lines.append(
            f"{row['Period']:>9}  {row['FSL']:>7.2f}  {row['Level_max']:>9.3f}  "
            f"{row['New_FSL']:>7.2f}  {row['New_Level_max']:>11.3f}  "
            f"{row['New_Volume_max_ML']:>10,.0f}  "
            f"{row['New_Storage_above_FSL_max_ML']:>10,.0f}  "
            f"{row['New_Release_peak_m3s']:>8,.0f}{flag}"
        )
    if (~table["Complete"]).any():
        lines.append("")
        lines.append("* water year not fully covered by the record")
    return "\n".join(lines)
