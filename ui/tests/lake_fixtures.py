"""Synthetic headwater level records, written the way Hydstra writes them.

Generated, never committed - Bryan's .gitignore excludes ``*.csv``. The record
is built to have the shape the lake level frequency analysis exists for: a
dry-season drawdown, wet-season floods, a lake held flat at full supply in the
wetter years with the occasional flood driven above it, and water years that
open at their highest because the flood came just before the boundary.

Shared by ``ui/tests`` and Bryan's own ``tests/``, which puts this folder on its
path.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

FSL = 215.0


def synthetic_level(years=40, start="1980-10-01", step="3h", seed=7, fsl=FSL) -> pd.Series:
    """A level record with a plateau at ``fsl`` and floods above it."""
    rng = np.random.default_rng(seed)
    index = pd.date_range(start, periods=int(years * 365.25 * 24 / 3), freq=step)
    hours = 3.0
    level = np.empty(len(index))
    current = fsl - 6.0
    surcharge = 0.0
    # Floods arrive mostly December to March, as a rise spread over a day.
    months = index.month.to_numpy()
    wet = np.isin(months, (12, 1, 2, 3))
    flood_start = rng.random(len(index)) < np.where(wet, 0.0025, 0.0002)
    rising = 0
    rise_per_step = 0.0
    for i in range(len(index)):
        if flood_start[i] and rising == 0:
            rising = 8
            rise_per_step = rng.exponential(2.6) / rising
        if rising:
            current += rise_per_step
            rising -= 1
        else:
            current -= 0.004 * hours / 24 * 3        # drawdown
        if current > fsl:
            surcharge = max(surcharge, (current - fsl - 0.5) * 0.9)
            current = fsl
        surcharge = max(0.0, surcharge - 0.03)
        current = max(current, fsl - 16.0)
        level[i] = current + surcharge
    return pd.Series(np.round(level, 3), index=index, name="level")


def stamp(when: pd.Timestamp) -> str:
    return f"{when.day}/{when.month:02d}/{when.year} {when.hour}:{when.minute:02d}"


def write_hydstra(path, series: pd.Series, site="999999A", kind="Point",
                  description="Test Ck at Test Dam HW") -> Path:
    meta = ["Sites:", f"{site} - {description}", "Variables:",
            "130.00 - Reservoir Water Level (Metres)", "Qualities:", "1 - Observed value"]
    lines = [f"Time,{site},", "and,130,", "Date,Storage Level (m),", f",{kind},"]
    for i, (when, value) in enumerate(series.items()):
        lines.append(f"{stamp(when)},{value:.3f},{meta[i] if i < len(meta) else ''}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def write_wmip(path, series: pd.Series, quality=None, site="999999A") -> Path:
    """The WMIP web export: seconds in the time, and a quality column."""
    quality = quality if quality is not None else [9] * len(series)
    lines = [f'"Time","{site}","",""', '"and","100.00","",""',
             '"Date","Level (Metres)","",""', '"Date and time","Point","Quality","Comments"']
    for i, ((when, value), code) in enumerate(zip(series.items(), quality)):
        comment = "Sites:" if i == 0 else " "
        lines.append(f"{when:%H:%M:%S %d/%m/%Y},     {value:.3f},      {code},{comment}")
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def write_mcdf(path, *, durations_shift=0.0, rows=4000, seed=3, fsl=FSL) -> Path:
    """A Monte Carlo database carrying just what the lake level page reads."""
    rng = np.random.default_rng(seed)
    aep = np.sort(10 ** rng.uniform(-6, 0, rows))
    z = -np.log10(aep)
    level = fsl - 8.0 + 8.0 * (1 - np.exp(-z * 1.2)) + 0.4 * z + durations_shift
    pd.DataFrame({"level": level, "level_aep": aep, "rain_aep": 1 / aep}).to_csv(
        path, index=False)
    return Path(path)
