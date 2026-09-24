"""Open-water evaporation from SILO Data Drill pan evaporation.

Kept apart from the model because its timing convention is a property of the
data source, not of the routing, and it is the one input whose meaning is easy
to get wrong.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("bryan.homogenise.evaporation")

PAN_EVAPORATION_FACTORS = {
    1: 0.82, 2: 0.82, 3: 0.83, 4: 0.79, 5: 0.75, 6: 0.71,
    7: 0.76, 8: 0.81, 9: 0.79, 10: 0.81, 11: 0.82, 12: 0.83,
}

# SILO reads evaporation at 9am and shifts it to the previous day, so the value
# in row D is the evaporation over 09:00 D to 09:00 D+1.  See Evaporation.
SILO_READ_HOUR = 9

# Both the daily-maximum export and the GoldSim series are stamped a day late:
# the row dated D carries the maximum of day D-1.  Against the raw gauge record
# the daily maxima agree to 0.006 m rmse at this offset and 0.124 m without it.

# --------------------------------------------------------------------------
# Evaporation
# --------------------------------------------------------------------------

class Evaporation:
    """Open-water evaporation as a step function on 9am boundaries.

    The SILO Data Drill header states:

        As evaporation is read at 9am, it has been shifted to the day before
        ie The evaporation measured on 20 April is in row for 19 April

    So the value in row D is already the evaporation **over 09:00 D to 09:00
    D+1** -- an integral, not an instantaneous rate, and not aligned to midnight.
    Two things follow.

    The rate is a step function whose breaks fall at 9am, so ``E(t0, t1)`` is
    computed as the exact integral of that step function.  This is resolution
    independent: a day's steps sum to exactly the day's total no matter how the
    day is subdivided, and no interpolation is invented between readings.

    It also means the quantity must not be averaged across consecutive rows.
    Trapezoidally averaging an already-integrated daily total, as a rate would
    be averaged, smooths it a second time and shifts it half a day.
    """

    def __init__(self, daily_mm, read_hour=SILO_READ_HOUR):
        if daily_mm.isna().any():
            raise ValueError("evaporation series contains gaps")
        self.daily_mm = daily_mm
        self.read_hour = read_hour
        offset = pd.Timedelta(hours=read_hour)
        breaks = daily_mm.index + offset
        # One extra break closes the final day.
        self._breaks = breaks.append(
            pd.DatetimeIndex([breaks[-1] + pd.Timedelta(days=1)])
        )
        self._break_ns = self._breaks.astype("int64").to_numpy(dtype=float)
        self._cumulative = np.concatenate([[0.0], np.cumsum(daily_mm.to_numpy())])

    @property
    def start(self):
        return self._breaks[0]

    @property
    def end(self):
        return self._breaks[-1]

    def cumulative_mm(self, timestamps):
        """Evaporation accumulated since ``self.start``, in mm."""
        stamps = pd.DatetimeIndex(timestamps).astype("int64").to_numpy(dtype=float)
        return np.interp(stamps, self._break_ns, self._cumulative)

    def over(self, index):
        """Evaporation in each step of ``index``, in mm.

        The first entry is NaN: a step needs two endpoints.
        """
        cumulative = self.cumulative_mm(index)
        stepwise = np.empty(len(index))
        stepwise[0] = np.nan
        stepwise[1:] = np.diff(cumulative)
        return stepwise

    def covers(self, index):
        return index.min() >= self.start and index.max() <= self.end


def read_evaporation(path, pan_factors=None):
    """Read SILO Data Drill pan evaporation and apply the monthly pan factor.

    The header row is located by content rather than by a hard-coded
    ``skiprows``, so a change in SILO's preamble length is not silent.
    ``pan_factors`` maps month (1-12) to the pan-to-open-water factor; the
    default is the set Callide used.
    """
    pan_factors = {int(month): float(value)
                   for month, value in (pan_factors or PAN_EVAPORATION_FACTORS).items()}
    if sorted(pan_factors) != list(range(1, 13)):
        raise ValueError("pan factors need one value for each month, 1 to 12")
    header_row, names = None, None
    with open(path, "r", errors="replace") as handle:
        for i, line in enumerate(handle):
            fields = line.split()
            if fields[:2] == ["Date", "Day"] and "Evap" in fields:
                header_row, names = i, fields
                break
    if header_row is None:
        raise ValueError(f"{path}: could not find the SILO column header row")

    frame = pd.read_csv(
        path, sep=r"\s+", skiprows=header_row + 2, names=names, index_col=False,
    )
    pan = pd.to_numeric(frame["Evap"], errors="coerce")
    pan.index = pd.to_datetime(frame["Date"].astype(str), format="%Y%m%d")
    pan = pan.dropna().sort_index()

    if pan.index.to_series().diff().dt.days.dropna().ne(1).any():
        raise ValueError(f"{path}: the evaporation record is not continuous daily")

    factors = pan.index.month.map(pan_factors)
    open_water = pan * factors
    open_water.name = "Evaporation (mm/d)"
    log.info("Evaporation: %d days, %s to %s, mean pan %.2f mm/d, open water "
             "%.2f mm/d, applied over %02d:00 to %02d:00 the following day",
             len(pan), pan.index.min().date(), pan.index.max().date(),
             pan.mean(), open_water.mean(), SILO_READ_HOUR, SILO_READ_HOUR)
    return Evaporation(open_water)


