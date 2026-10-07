"""The calibration table: a model's dam inflow against the reverse-routed inflow.

The **observed** hydrograph is one of the inflow record's event hydrographs
(the Lake record page's inflow step, ``lake_record/inflow/hydrographs``); the
**modelled** one is a CSV the analyst gives, of the model's dam inflow for the
same event. Per event the table gives both peaks and the peak ratio, both
volumes and the volume ratio, the Nash-Sutcliffe efficiency and, optionally,
the peak timing - each statistic shaded by the class its criteria give it
(``DEFAULT_CRITERIA``, the report's model performance table).

How each number is made - they change with each choice, so they are fixed here:

- **The window** is the modelled period. The observed side is taken from the
  whole inflow record (``inflow_intervals.csv.gz``) over it, not from the event
  hydrograph, whose window the inflow step chose and may be shorter than the
  model run. Where the model runs past the record, the window stops where the
  record does, and the table says so. Both volumes and the NSE are over it.
- **The observed series** is the recession-corrected inflow by default, or the
  uncorrected one (``observed_series``). Where the gates released more than the
  rating says, the corrected recession is zero; the share of the window where
  that happened is reported, because it lowers the observed volume and the NSE
  for something that is not a model error.
- **Observed volume** is exact: each interval's mean inflow times the part of the
  interval inside the window. **Modelled volume** is the trapezoidal integral of
  the modelled samples over the window.
- **Observed flow at a time** is the time average over the smoothing window
  (the inflow step's, 1 h by default) centred on it, off the cumulative volume -
  how the inflow record takes its peaks. With no smoothing it is the native
  interval's mean. The **observed peak** is the largest of those at the
  hydrograph's own stamps; the **NSE** compares them with the modelled samples
  at the modelled stamps.
- **PR** and **VR** are (modelled - observed) / observed, as a signed
  percentage; **timing** is the modelled peak's time less the observed one's,
  in hours. A class is judged on the absolute value of each.

pandas and numpy only: the UI imports this module.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import lakerecord
from .study import Study, resolve
from .wordtable import ReportTable

CORRECTED = "corrected"
UNCORRECTED = "uncorrected"
SERIES_COLUMNS = {CORRECTED: "Inflow_m3s", UNCORRECTED: "Inflow_uncorrected_m3s"}

ML_TO_M3 = 1000.0

# The report's model performance criteria (Callide, Table 15). Each statistic's
# classes are tried in order; a value past the last is the final class.
# ``limit`` bounds |PR|, |VR| (percent) and |timing| (hours) from above and NSE
# from below.
DEFAULT_CRITERIA = {
    "classes": [
        {"name": "Excellent", "score": 5, "fill": "#00B050"},
        {"name": "Good", "score": 4, "fill": "#92D050"},
        {"name": "Fair", "score": 3, "fill": "#FFFF00"},
        {"name": "Poor", "score": 2, "fill": "#FFC000"},
        {"name": "No Data/Exclude", "score": 0, "fill": "#A6A6A6"},
    ],
    "pr": [10, 15, 30, 50],
    "vr": [15, 25, 30, 50],
    "nse": [0.95, 0.90, 0.85, 0.50],
    "timing": [0.5, 1.0, 2.0, 3.0],
}
HIGHER_IS_BETTER = {"nse"}


# -- the observed side --------------------------------------------------------

def observed_events(study: Study) -> list:
    """The inflow step's event hydrographs: [{name, file, start, end, ...}]."""
    section = lakerecord.settings(study)
    summary = lakerecord.last_summary(study, section, "inflow")
    return [item for item in summary.get("hydrographs") or [] if item.get("name")]


def smoothing_of(study: Study) -> str:
    """The inflow step's smoothing window, as its last run recorded it."""
    section = lakerecord.settings(study)
    summary = lakerecord.last_summary(study, section, "inflow")
    return str((summary.get("settings") or {}).get("smoothing") or "")


def seconds(index) -> np.ndarray:
    """Seconds since 1970 - whatever unit pandas holds the times in."""
    return (pd.DatetimeIndex(index) - pd.Timestamp("1970-01-01")).total_seconds().to_numpy()


def _intervals(frame: pd.DataFrame, column: str):
    """(starts, ends, rates) of the observed intervals, in seconds since 1970."""
    mid = seconds(frame.index)
    dt = frame["dt_s"].to_numpy(dtype=float)
    return mid - dt / 2, mid + dt / 2, frame[column].to_numpy(dtype=float)


def _cumulative(starts, ends, rates):
    """Running volume (m3) at the interval boundaries."""
    times = np.concatenate([[starts[0]], ends])
    volumes = np.concatenate([[0.0], np.cumsum(np.nan_to_num(rates) * (ends - starts))])
    return times, volumes


def observed_at(frame: pd.DataFrame, column: str, times, window_s: float) -> np.ndarray:
    """Observed inflow (m3/s) at ``times`` (s): the centred time average over
    ``window_s``, or the native interval's mean when it is 0."""
    starts, ends, rates = _intervals(frame, column)
    times = np.asarray(times, dtype=float)
    if window_s > 0:
        stamps, volumes = _cumulative(starts, ends, rates)
        low = np.interp(times - window_s / 2, stamps, volumes)
        high = np.interp(times + window_s / 2, stamps, volumes)
        return (high - low) / window_s
    where = np.searchsorted(ends, times, side="left")
    where = np.clip(where, 0, len(rates) - 1)
    return rates[where]


def observed_volume(frame: pd.DataFrame, column: str, t0: float, t1: float) -> float:
    """ML of observed inflow between t0 and t1 (s): each interval's share inside."""
    starts, ends, rates = _intervals(frame, column)
    inside = np.clip(np.minimum(ends, t1) - np.maximum(starts, t0), 0, None)
    return float(np.nansum(rates * inside) / ML_TO_M3)


def uncertain_share(frame: pd.DataFrame, t0: float, t1: float) -> float:
    """The share of the window where the release above FSL is uncertain."""
    if "Release_uncertain" not in frame:
        return 0.0
    starts, ends, _ = _intervals(frame, "dt_s")
    inside = np.clip(np.minimum(ends, t1) - np.maximum(starts, t0), 0, None)
    flags = frame["Release_uncertain"].astype(str).str.lower().isin(["true", "1"]).to_numpy()
    total = inside.sum()
    return float((inside * flags).sum() / total) if total > 0 else 0.0


# -- the modelled side ----------------------------------------------------------

class CalibrationError(ValueError):
    pass


def read_modelled(path, column: str = "", start=None) -> pd.Series:
    """The modelled dam inflow (m3/s), indexed by time.

    The first column is the time: dates, or hours from ``start`` (the event's
    start, given with the event) when it is a plain number. The flow is
    ``column``, or the first numeric column after the time.
    """
    path = Path(path)
    if not path.is_file():
        raise CalibrationError(f"modelled inflow not found: {path}")
    frame = pd.read_csv(path)
    if frame.shape[1] < 2:
        raise CalibrationError(f"{path.name}: needs a time column and a flow column")
    time = frame.iloc[:, 0]
    if pd.api.types.is_numeric_dtype(time):
        if start in (None, ""):
            raise CalibrationError(f"{path.name}: the time column is in hours - give the "
                                   f"event's start time to place it")
        index = pd.Timestamp(start) + pd.to_timedelta(time.astype(float), unit="h")
    else:
        index = _dates(time, path.name)
    if column:
        if column not in frame.columns:
            raise CalibrationError(f"{path.name}: no column {column!r} (it has "
                                   f"{', '.join(map(str, frame.columns[1:]))})")
        flow = frame[column]
    else:
        numeric = [name for name in frame.columns[1:]
                   if pd.api.types.is_numeric_dtype(frame[name])]
        if not numeric:
            raise CalibrationError(f"{path.name}: no numeric flow column")
        flow = frame[numeric[0]]
    series = pd.Series(pd.to_numeric(flow, errors="coerce").to_numpy(),
                       index=pd.DatetimeIndex(index), name="modelled").dropna()
    return series[~series.index.duplicated()].sort_index()


def _dates(values: pd.Series, name: str) -> pd.DatetimeIndex:
    text = values.astype(str).str.strip()
    try:
        return pd.DatetimeIndex(pd.to_datetime(text, format="ISO8601"))
    except (ValueError, TypeError):
        pass
    try:
        return pd.DatetimeIndex(pd.to_datetime(text, dayfirst=True))
    except (ValueError, TypeError) as exc:
        raise CalibrationError(f"{name}: cannot read the times in the first column "
                               f"({text.iloc[0]!r})") from exc


# -- the statistics ---------------------------------------------------------------

def statistics(observed: pd.DataFrame, modelled: pd.Series, *,
               series: str = CORRECTED, smoothing: str = "1h") -> dict:
    """PR, VR, NSE and timing for one event, with what they were made from."""
    column = SERIES_COLUMNS[series]
    if column not in observed:
        raise CalibrationError(f"the observed hydrograph has no {column}")
    window_s = pd.Timedelta(smoothing).total_seconds() if smoothing else 0.0
    starts, ends, _ = _intervals(observed, column)
    m_times = seconds(modelled.index)
    t0, t1 = max(starts[0], m_times[0]), min(ends[-1], m_times[-1])
    if t1 <= t0:
        raise CalibrationError("the modelled and observed hydrographs do not overlap")

    # the observed peak, at the hydrograph's own stamps inside the window
    stamps = seconds(observed.index)
    keep = (stamps >= t0) & (stamps <= t1)
    o_rates = observed_at(observed, column, stamps[keep], window_s)
    o_peak_at = int(np.nanargmax(o_rates))
    o_peak, o_peak_time = float(o_rates[o_peak_at]), stamps[keep][o_peak_at]

    # the modelled samples inside the window, with the window's ends interpolated
    m_values = modelled.to_numpy(dtype=float)
    inside = (m_times >= t0) & (m_times <= t1)
    m_peak_at = int(np.nanargmax(np.where(inside, m_values, -np.inf)))
    m_peak, m_peak_time = float(m_values[m_peak_at]), m_times[m_peak_at]
    edge_t = np.concatenate([[t0], m_times[inside], [t1]])
    edge_q = np.interp(edge_t, m_times, m_values)
    m_volume = float(np.trapezoid(edge_q, edge_t) / ML_TO_M3) if hasattr(np, "trapezoid") \
        else float(np.trapz(edge_q, edge_t) / ML_TO_M3)
    o_volume = observed_volume(observed, column, t0, t1)

    o_at_m = observed_at(observed, column, m_times[inside], window_s)
    good = np.isfinite(o_at_m) & np.isfinite(m_values[inside])
    o_fit, m_fit = o_at_m[good], m_values[inside][good]
    spread = float(np.sum((o_fit - o_fit.mean()) ** 2)) if len(o_fit) else 0.0
    nse = 1.0 - float(np.sum((m_fit - o_fit) ** 2)) / spread if spread > 0 else math.nan

    return {
        "observed_peak": o_peak, "modelled_peak": m_peak,
        "pr": _ratio(m_peak, o_peak),
        "observed_volume": o_volume, "modelled_volume": m_volume,
        "vr": _ratio(m_volume, o_volume),
        "nse": nse,
        "timing": (m_peak_time - o_peak_time) / 3600.0,
        "window": (pd.Timestamp(t0, unit="s"), pd.Timestamp(t1, unit="s")),
        "uncertain": uncertain_share(observed, t0, t1),
        "points": int(good.sum()),
    }


def _ratio(modelled: float, observed: float) -> float:
    return (modelled - observed) / observed * 100.0 if observed else math.nan


def classify(statistic: str, value, criteria: dict | None = None) -> dict | None:
    """The class a statistic's value falls in, or None when there is no value."""
    criteria = criteria or DEFAULT_CRITERIA
    if value is None or not np.isfinite(value):
        return None
    classes = criteria["classes"]
    limits = criteria[statistic]
    for position, limit in enumerate(limits):
        if (value >= limit) if statistic in HIGHER_IS_BETTER else (abs(value) <= limit):
            return classes[position]
    return classes[min(len(limits), len(classes) - 1)]


# -- the table ----------------------------------------------------------------------

def observed_for(record, item, modelled: pd.Series, margin: pd.Timedelta) -> pd.DataFrame:
    """The observed intervals over the modelled period.

    From the whole inflow record where there is one, so the modelled period alone
    sets the window - an event hydrograph stops where the inflow step's event
    window does, which may be shorter than the model run. ``margin`` (half the
    smoothing window is enough; the whole is used) keeps the time average right
    at the window's ends. The event hydrograph is the fallback, for a record
    written before the intervals file was.
    """
    first, last = modelled.index[0] - margin, modelled.index[-1] + margin
    if record is not None:
        frame = lakerecord.read_record(record)
        ends = frame.index + pd.to_timedelta(frame["dt_s"] / 2, unit="s")
        starts = frame.index - pd.to_timedelta(frame["dt_s"] / 2, unit="s")
        part = frame[(ends > first) & (starts < last)]
        if part.empty:
            raise CalibrationError(f"the inflow record has nothing between "
                                   f"{modelled.index[0]:%Y-%m-%d %H:%M} and "
                                   f"{modelled.index[-1]:%Y-%m-%d %H:%M}")
        return part
    if item is None:
        raise CalibrationError("choose the observed hydrograph - the inflow record's "
                               "intervals file is not there to take it from")
    observed = lakerecord.hydrograph(item)
    if observed is None:
        raise CalibrationError(f"observed hydrograph file not found: {item.get('file')}")
    return observed


def _outside(modelled: pd.Series, window) -> str:
    """How much of the model run the window leaves out, as '6 h', or ''."""
    start, end = window
    hours = (max(start - modelled.index[0], pd.Timedelta(0))
             + max(modelled.index[-1] - end, pd.Timedelta(0))).total_seconds() / 3600
    return f"{hours:g} h" if hours > 0.01 else ""


def _fmt_flow(value) -> str:
    return f"{value:,.0f}" if np.isfinite(value) else "–"


def _fmt_percent(value) -> str:
    return f"{value:.0f}%" if np.isfinite(value) else "–"


def build(study: Study, spec: dict) -> ReportTable:
    mod, obs = spec.get("modelled_label") or "Mod", spec.get("observed_label") or "Rated"
    timing = bool(spec.get("timing"))
    header = ["Event", mod, obs, "PR", mod, obs, "VR", "Nash-Sutcliffe"]
    groups = [("", 1), ("Flow (m3/s)", 3), ("Volume (ML)", 3), ("", 1)]
    if timing:
        header.append("Timing (h)")
        groups.append(("", 1))
    table = ReportTable(header=header, align=["center"] * len(header))
    table.header_groups = groups
    criteria = spec.get("criteria") or DEFAULT_CRITERIA
    series = spec.get("observed_series") or CORRECTED
    smoothing = (smoothing_of(study) or "1h") if spec.get("smoothed", True) else ""
    summary = lakerecord.last_summary(study, lakerecord.settings(study), "inflow")
    events = {item["name"]: item for item in summary.get("hydrographs") or []
              if item.get("name")}
    record_path = lakerecord.intervals_path(summary)
    record = record_path if record_path is not None and record_path.is_file() else None
    if not events and record is None:
        table.problems.append("There is no inflow record yet - run the Lake record "
                              "page's inflow step first")
    margin = pd.Timedelta(smoothing) if smoothing else pd.Timedelta(0)
    for event in spec.get("events") or []:
        label = event.get("label") or event.get("observed") or "(unnamed)"
        item = events.get(event.get("observed") or "")
        try:
            path = resolve(study.folder, event.get("modelled") or "")
            if path is None:
                raise CalibrationError("no modelled inflow file given")
            modelled = read_modelled(path, event.get("column") or "",
                                     event.get("start") or (item or {}).get("start"))
            if modelled.empty:
                raise CalibrationError(f"{path.name} has no flows")
            observed = observed_for(record, item, modelled, margin)
            stats = statistics(observed, modelled, series=series, smoothing=smoothing)
        except (CalibrationError, ValueError, KeyError, OSError) as exc:
            table.problems.append(f"{label}: {exc}")
            table.add([label] + ["–"] * (len(header) - 1))
            continue
        short = _outside(modelled, stats["window"])
        if short:
            table.problems.append(f"{label}: the model runs {short} past the observed "
                                  f"record; the statistics stop where the record does")
        if stats["uncertain"] > 0 and series == CORRECTED:
            table.problems.append(
                f"{label}: the release is uncertain over {stats['uncertain']:.0%} of the "
                f"window, where the corrected inflow is held at zero - it lowers the "
                f"observed volume and the NSE")
        nse = stats["nse"]
        nse_text = ("–" if not np.isfinite(nse) else
                    f"{nse * 100:.0f}%" if spec.get("nse_percent", True) else f"{nse:.2f}")
        cells = [label, _fmt_flow(stats["modelled_peak"]), _fmt_flow(stats["observed_peak"]),
                 _fmt_percent(stats["pr"]), _fmt_flow(stats["modelled_volume"]),
                 _fmt_flow(stats["observed_volume"]), _fmt_percent(stats["vr"]), nse_text]
        shaded = {3: ("pr", stats["pr"]), 6: ("vr", stats["vr"]), 7: ("nse", nse)}
        if timing:
            cells.append(f"{stats['timing']:+.1f}")
            shaded[8] = ("timing", stats["timing"])
        row = table.add(cells)
        if spec.get("shade", True):
            row.fills = [""] * len(cells)
            for position, (name, value) in shaded.items():
                found = classify(name, value, criteria)
                if found:
                    row.fills[position] = found["fill"]
    return table


def criteria_table(spec: dict | None = None) -> ReportTable:
    """The criteria themselves, as the report's model performance table."""
    criteria = (spec or {}).get("criteria") or DEFAULT_CRITERIA
    table = ReportTable(header=["Class", "Score", "PR", "VR", "NS", "Timing (hrs)"],
                        align=["left"] + ["center"] * 5)
    count = len(criteria["pr"])
    for position, found in enumerate(criteria["classes"]):
        if position < count:
            cells = [found["name"], str(found["score"]),
                     f"+/- {criteria['pr'][position]:g}%", f"+/- {criteria['vr'][position]:g}%",
                     f">={criteria['nse'][position]:.2f}", f"<={criteria['timing'][position]:.2f}"]
        else:
            cells = [found["name"], str(found["score"]), f"> {criteria['pr'][-1]:g}%",
                     f"> {criteria['vr'][-1]:g}%", f"<{criteria['nse'][-1]:g}",
                     f">{criteria['timing'][-1]:g}"]
        row = table.add(cells)
        row.fills = [found["fill"]] + [""] * 5
    return table
