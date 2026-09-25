"""The Lake record page's settings, jobs and results - everything but the view.

Four steps, in the order they depend on each other, all kept in the study file
under ``"lake_record"`` with paths relative to it:

1. **Catchment rainfall** - ``core/awap.py``, in this process: a shapefile and
   the folder of AWAP/AWRA-L daily grids give ``date,rain_mm``. **The series is
   what the study keeps and ships**, in its own folder; the grids are tens of
   gigabytes not every user has, so their folder is each user's setting
   (``UiSettings.awap_folder``), never the study's, and steps 2-3 need only the
   series. On Callide the area-weighted series from ``C:\\AWAP`` reproduces the
   JPA project's to 0.008 mm rms over 42,126 days, dates and all.
2. **Homogenisation** - ``util/HomogeniseLakeLevels.py`` under Bryan's
   interpreter (scipy): the gauge exports, storage table, rating register,
   evaporation and target ratings give the lake re-routed through each target.
3. **Antecedent storage** - ``util/AntecedentStorage.py`` under Bryan's
   interpreter: the rainfall, the IFD and the homogenisation give the antecedent
   series and the ``lake_config.json`` files Bryan's Monte Carlo scheme samples.
4. **Inflow record** - ``util/InflowRecord.py`` under Bryan's interpreter: step
   2's derived inflow (stage 1, the rating in force at each step) as the annual
   maximum peak inflow and burst volumes, and event hydrographs for calibration.
   It shares step 2's inputs but not its target ratings, and by default leaves
   the evaporation out, so it is not clipped where the evaporation file ends.

This module writes the job files beside the outputs (so a run can be repeated
from a console exactly as the page ran it), runs the scripts, and reads back what
they wrote. No nicegui.
"""

from __future__ import annotations

import copy
import json
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import dam as dams
from .bryan import BRYAN_ROOT
from .paths import atomic_write_json, read_json
from .study import Study, portable, resolve
from .wordtable import ReportTable

KEY = "lake_record"
# Where the dam inputs are edited, which the step checks point to.
RECORD_CARD = "Dam inputs"
DAM_HINT = f"(These are the {RECORD_CARD.lower()}, set on the Study page.)"
HOMOGENISE_SCRIPT = BRYAN_ROOT / "util" / "HomogeniseLakeLevels.py"
ANTECEDENT_SCRIPT = BRYAN_ROOT / "util" / "AntecedentStorage.py"
INFLOW_SCRIPT = BRYAN_ROOT / "util" / "InflowRecord.py"
TIMEOUT_SECONDS = 3600

MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
CALLIDE_PAN_FACTORS = [0.82, 0.82, 0.83, 0.79, 0.75, 0.71, 0.76, 0.81, 0.79, 0.81, 0.82, 0.83]

DEFAULTS = {
    "rainfall": {"shapefile": "", "field": "", "value": "", "pattern": "rain_day_*.nc",
                 "variable": "", "weighting": "area", "output": "lake_record/rainfall.csv"},
    "homogenise": {"gauges": [], "overlay": None, "storage": "", "register": "",
                   "evaporation": "", "pan_factors": list(CALLIDE_PAN_FACTORS), "step": "1h",
                   "recession_correction": True, "water_year_start": 10,
                   "separation_days": 5, "drop_m": 0.5, "targets": [],
                   "out": "lake_record/homogenised"},
    "antecedent": {"rainfall": "", "ifd": "", "targets": [], "bases": ["burst", "storm"],
                   "settings": {"durations_d": [1, 2, 3, 4, 5],
                                "restriction": {"1": 1.15, "2": 1.11, "3": 1.07,
                                                "4": 1.05, "5": 1.04},
                                "window_days": 30, "threshold_fraction": 0.8,
                                "preburst_mm": 10, "cunnane_a": 0.4, "round_ml": 1000},
                   "out": "lake_record/antecedent"},
    "inflow": {"evaporation": False, "recession_correction": True, "smoothing": "1h",
               "durations_h": [24, 36, 48, 72], "catchment_km2": None, "rainfall": "",
               "top_events": 10, "before_days": 3, "after_days": 7, "events": [],
               "out": "lake_record/inflow",
               "window": {"start": "", "end": "", "step": ""}},
}


# -- settings ---------------------------------------------------------------------

def settings(study: Study) -> dict:
    """The study's section, completed with the defaults, one level deep, with the
    dam inputs (``core/dam.py``) laid over it where the jobs read them.

    The dam inputs are edited on the Study page and kept under ``"dam"``; the
    copies in this section are only what the jobs read, and what an older
    launcher finds.
    """
    stored = study.extra.get(KEY) or {}
    out = copy.deepcopy(DEFAULTS)
    for step, defaults in out.items():
        defaults.update(copy.deepcopy(stored.get(step) or {}))
    out["antecedent"]["settings"] = {**DEFAULTS["antecedent"]["settings"],
                                     **(out["antecedent"].get("settings") or {})}
    return dams.write_through(out, dams.settings(study))


def store(study: Study, section: dict) -> None:
    """Keep the section, with the dam inputs as they are now - a page opened before
    they were changed on the Study page must not write the old ones back."""
    study.extra[KEY] = dams.write_through(section, dams.settings(study))


def keep_path(study: Study, text) -> str:
    """A path as the study stores it: relative to the study file where it can be."""
    return portable(study.folder, text)


def path_of(study: Study, text) -> Path | None:
    return resolve(study.folder, text)


# -- the jobs -----------------------------------------------------------------------

def _absolute(study, value):
    path = path_of(study, value)
    return str(path) if path else ""


def homogenise_job(study: Study, section: dict) -> dict:
    """The homogenisation job, every path made absolute, ready to write."""
    h = section["homogenise"]
    job = {key: h[key] for key in ("step", "recession_correction", "water_year_start",
                                   "separation_days", "drop_m")}
    job["gauges"] = [_absolute(study, item) for item in h["gauges"] if str(item).strip()]
    for key in ("storage", "register", "evaporation"):
        job[key] = _absolute(study, h[key])
    factors = h.get("pan_factors") or CALLIDE_PAN_FACTORS
    job["pan_factors"] = {str(month + 1): float(value) for month, value in enumerate(factors)}
    overlay = h.get("overlay")
    job["overlay"] = ({"file": _absolute(study, overlay.get("file")),
                       "below": float(overlay["below"]),
                       "reconnect_margin": float(overlay.get("reconnect_margin", 0.10))}
                      if overlay and overlay.get("file") else None)
    job["targets"] = [{"name": str(t.get("name") or "").strip(),
                       "rating": _absolute(study, t.get("rating")),
                       "fsl": float(t["fsl"]) if t.get("fsl") not in (None, "") else None}
                      for t in h["targets"] if t.get("rating")]
    job["out"] = _absolute(study, h["out"])
    return job


def antecedent_job(study: Study, section: dict, homogenise_path: Path) -> dict:
    a = section["antecedent"]
    rainfall = a.get("rainfall") or section["rainfall"]["output"]
    return {"homogenise": str(homogenise_path), "rainfall": _absolute(study, rainfall),
            "ifd": _absolute(study, a["ifd"]), "targets": list(a.get("targets") or []),
            "bases": list(a.get("bases") or ["burst", "storm"]),
            "settings": copy.deepcopy(a["settings"]), "out": _absolute(study, a["out"])}


def inflow_job(study: Study, section: dict, homogenise_path: Path) -> dict:
    i = section["inflow"]
    rainfall = path_of(study, i.get("rainfall") or section["rainfall"]["output"])
    job = {key: copy.deepcopy(i[key]) for key in (
        "evaporation", "recession_correction", "smoothing", "durations_h", "catchment_km2",
        "top_events", "before_days", "after_days", "events")}
    job["homogenise"] = str(homogenise_path)
    job["rainfall"] = str(rainfall) if rainfall and rainfall.is_file() else ""
    job["out"] = _absolute(study, i["out"])
    return job


def parse_events(text) -> tuple:
    """``name, start, end`` a line -> (events, problems)."""
    events, problems = [], []
    for number, line in enumerate(str(text or "").splitlines(), start=1):
        if not line.strip():
            continue
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 3 or not parts[0]:
            problems.append(f"line {number}: give name, start, end")
            continue
        try:
            start, end = pd.Timestamp(parts[1]), pd.Timestamp(parts[2])
        except (ValueError, TypeError):
            problems.append(f"line {number}: the start or end is not a date")
            continue
        if end <= start:
            problems.append(f"line {number}: the end is not after the start")
            continue
        events.append({"name": parts[0], "start": parts[1], "end": parts[2]})
    return events, problems


def events_text(events) -> str:
    return "\n".join(f"{e['name']}, {e['start']}, {e['end']}" for e in events or [])


def problems_before_running(study: Study, section: dict, step: str, grids: str = "") -> list:
    """What would stop a step, said before anything is started.

    ``grids`` is the user's own folder of daily grids (the rainfall step only).
    """
    out = []

    def need(label, value, is_dir=False):
        path = path_of(study, value)
        if path is None:
            out.append(f"{label}: not given")
        elif not (path.is_dir() if is_dir else path.is_file()):
            out.append(f"{label}: not found ({path})")

    if step == "rainfall":
        r = section["rainfall"]
        need("Catchment shapefile", r["shapefile"])
        if out:
            out.append(DAM_HINT)
        if not str(grids or "").strip():
            out.append("Folder of daily grids on this computer: not given")
        elif not Path(grids).is_dir():
            out.append(f"Folder of daily grids on this computer: not found ({grids})")
    elif step in ("homogenise", "antecedent", "inflow"):
        h = section["homogenise"]
        if not [g for g in h["gauges"] if str(g).strip()]:
            out.append("Gauge exports: none given")
        for number, gauge in enumerate(g for g in h["gauges"] if str(g).strip()):
            need(f"Gauge export {number + 1}", gauge)
        need("Storage table (.els)", h["storage"])
        need("Rating register (.xlsx)", h["register"])
        if step != "inflow" or section["inflow"]["evaporation"]:
            need("Evaporation (SILO)", h["evaporation"])
        if h.get("overlay"):
            need("Overlay gauge", h["overlay"].get("file"))
        if out:
            out.append(DAM_HINT)
    if step in ("homogenise", "antecedent"):
        if not h["targets"]:
            out.append("Target ratings: none given")
        for target in h["targets"]:
            need(f"Target rating {target.get('name') or ''}".strip(), target.get("rating"))
        if len(set(str(t.get("name")) for t in h["targets"])) != len(h["targets"]):
            out.append("Target ratings: each needs its own name")
    if step == "antecedent":
        a = section["antecedent"]
        need("Rainfall series", a.get("rainfall") or section["rainfall"]["output"])
        need("IFD table", a["ifd"])
    return out


# -- running -------------------------------------------------------------------------

@dataclass
class RunResult:
    returncode: int
    output: str
    summary: dict

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and bool(self.summary)


def run_script(script: Path, job: dict, job_path: Path, python) -> RunResult:
    """Write the job, run the script on it, read its summary. Blocking."""
    atomic_write_json(job_path, job)
    environment = dict(os.environ)
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    try:
        finished = subprocess.run([str(python), "-u", str(script), str(job_path)],
                                  capture_output=True, text=True, timeout=TIMEOUT_SECONDS,
                                  env=environment, stdin=subprocess.DEVNULL,
                                  cwd=str(job_path.parent))
    except subprocess.TimeoutExpired:
        return RunResult(1, f"timed out after {TIMEOUT_SECONDS} s", {})
    except OSError as exc:
        return RunResult(1, f"could not start {script.name}: {exc}", {})
    output = (finished.stdout or "") + (finished.stderr or "")
    out_folder = Path(job.get("out") or job_path.parent)
    summary = read_json(out_folder / "summary.json", default={}) if finished.returncode == 0 else {}
    return RunResult(finished.returncode, output, summary or {})


def homogenise(study: Study, section: dict, python) -> RunResult:
    job = homogenise_job(study, section)
    return run_script(HOMOGENISE_SCRIPT, job, Path(job["out"]) / "homogenise_job.json", python)


def job_path_for_homogenisation(study: Study, section: dict) -> Path:
    return path_of(study, section["homogenise"]["out"]) / "homogenise_job.json"


def antecedent(study: Study, section: dict, python) -> RunResult:
    homogenise_path = job_path_for_homogenisation(study, section)
    atomic_write_json(homogenise_path, homogenise_job(study, section))
    job = antecedent_job(study, section, homogenise_path)
    return run_script(ANTECEDENT_SCRIPT, job, Path(job["out"]) / "antecedent_job.json", python)


def inflow(study: Study, section: dict, python) -> RunResult:
    homogenise_path = job_path_for_homogenisation(study, section)
    atomic_write_json(homogenise_path, homogenise_job(study, section))
    job = inflow_job(study, section, homogenise_path)
    return run_script(INFLOW_SCRIPT, job, Path(job["out"]) / "inflow_job.json", python)


def last_summary(study: Study, section: dict, step: str) -> dict:
    folder = path_of(study, section[step]["out"])
    return read_json(folder / "summary.json", default={}) if folder else {}


# -- reading the results back ---------------------------------------------------------

def annual_maxima(summary: dict) -> dict:
    """target -> AMS frame (recorded and homogenised maximum level by water year)."""
    out = {}
    for target in summary.get("targets") or []:
        path = Path(target.get("files", {}).get("ams", ""))
        if path.is_file():
            out[target["name"]] = pd.read_csv(path)
    return out


def inflow_maxima(summary: dict) -> pd.DataFrame | None:
    path = Path(summary.get("ams") or "")
    return pd.read_csv(path) if summary.get("ams") and path.is_file() else None


def hydrograph(item: dict) -> pd.DataFrame | None:
    """One event's hydrograph as the inflow step wrote it."""
    path = Path(item.get("file") or "")
    if not item.get("file") or not path.is_file():
        return None
    frame = pd.read_csv(path, index_col=0)
    frame.index = pd.DatetimeIndex(pd.to_datetime(frame.index, format="ISO8601"),
                                   name="Timestamp")
    return frame


def ams_table(ams: pd.DataFrame, durations) -> ReportTable:
    """The inflow AMS as a report table: peak, and each burst's volume, depth and rain."""
    header = ["Water year", "Peak inflow (m3/s)", "Peak time"]
    columns = []
    for hours in durations:
        for column, label in ((f"Volume_{hours}h_ML", f"{hours} h volume (ML)"),
                              (f"Depth_{hours}h_mm", f"{hours} h runoff (mm)"),
                              (f"Rain_{hours}h_mm", f"{hours} h rain (mm)")):
            if column in ams:
                columns.append((column, 0 if column.startswith("Volume") else 1))
                header.append(label)
    table = ReportTable(header=header, align=["left", "right", "left"])
    complete = ams["Complete"].astype(str).str.lower() == "true" if "Complete" in ams else None
    for index, row in ams.iterrows():
        cells = [row["Period"], f"{row['Peak_inflow_m3s']:,.0f}",
                 f"{pd.Timestamp(row['Peak_time']):%d %b %Y %H:%M}"]
        cells += ["" if pd.isna(row[column]) else f"{row[column]:,.{places}f}"
                  for column, places in columns]
        table.add(cells)
    if complete is not None and not complete.all():
        table.footnotes.append("Part years (less than 90% of the year recorded): "
                               + ", ".join(map(str, ams.loc[~complete, "Period"])) + ".")
    return table


# -- the inflow record between two dates ------------------------------------------------

_RECORD_CACHE: dict = {}


def intervals_path(summary: dict) -> Path | None:
    if summary.get("intervals"):
        return Path(summary["intervals"])
    return Path(summary["ams"]).parent / "inflow_intervals.csv.gz" if summary.get("ams") else None


def read_record(path: Path) -> pd.DataFrame:
    """Every interval of the inflow record, kept while the file is unchanged."""
    path = Path(path)
    stamp = path.stat().st_mtime_ns
    cached = _RECORD_CACHE.get(path)
    if cached is None or cached[0] != stamp:
        frame = pd.read_csv(path, index_col=0)
        # midpoints of odd intervals carry fractional seconds, so the formats are mixed
        frame.index = pd.DatetimeIndex(pd.to_datetime(frame.index, format="ISO8601"),
                                       name="Timestamp")
        frame["Interval_start"] = pd.to_datetime(frame["Interval_start"], format="ISO8601")
        _RECORD_CACHE.clear()
        _RECORD_CACHE[path] = (stamp, frame)
    return _RECORD_CACHE[path][1]


def parse_window(start, end, step="") -> tuple:
    """(start, end, step) as timestamps and a timedelta (None for native), or ValueError."""
    try:
        start, end = pd.Timestamp(str(start).strip()), pd.Timestamp(str(end).strip())
    except (ValueError, TypeError) as exc:
        raise ValueError("give the start and end as dates, e.g. 2013-01-20 or "
                         "2013-01-20 09:00") from exc
    if end <= start:
        raise ValueError("the end is not after the start")
    step = str(step or "").strip()
    if not step:
        return start, end, None
    try:
        delta = pd.Timedelta(step)
    except ValueError as exc:
        raise ValueError(f"'{step}' is not a time step - e.g. 15min, 1h, 1D") from exc
    if delta <= pd.Timedelta(0):
        raise ValueError("the time step must be positive")
    if (end - start) / delta > 2_000_000:
        raise ValueError("that time step gives more than two million rows")
    return start, end, delta


def record_window(frame: pd.DataFrame, start, end, step=None) -> pd.DataFrame:
    """The record between two times: its own intervals, or means over a regular step.

    Native, each row is one interval of the record stamped at its midpoint. On a
    regular step, each row is the step *ending* at its time stamp (as a model's
    hydrograph is), and every flow is the exact mean over the step - the change in
    the running volume across it - so the step's volume is kept whatever the gauge
    did inside it. The level is read at the stamp; ``Release_uncertain`` becomes the
    share of the step the recession correction touched.
    """
    if step is None:
        part = frame.loc[start:end]
        out = pd.DataFrame({
            "Inflow_m3s": part["Inflow_m3s"], "Inflow_native_m3s": part["Inflow_native_m3s"],
            "Inflow_uncorrected_m3s": part["Inflow_uncorrected_m3s"],
            "Release_m3s": part["Release_m3s"], "Level_m": part["Level_end"],
            "Volume_ML": part["Volume_ML"], "Release_uncertain": part["Release_uncertain"],
            "Interpolated": part["Interpolated"], "dt_s": part["dt_s"]})
        out.index.name = "Timestamp"
        return out

    starts = pd.DatetimeIndex(frame["Interval_start"])
    dt = frame["dt_s"].to_numpy(dtype=float)
    ends = starts + pd.to_timedelta(dt, unit="s")
    uncertain = frame["Release_uncertain"].astype(str).str.lower().eq("true").to_numpy()
    running = {
        "Inflow_m3s": frame["Volume_ML"].to_numpy(dtype=float) * 1000.0,
        "Inflow_uncorrected_m3s": frame["Inflow_uncorrected_m3s"].to_numpy(dtype=float) * dt,
        "Release_m3s": frame["Release_m3s"].to_numpy(dtype=float) * dt,
        "Release_uncertain": uncertain * dt,
    }
    origin = starts[0]
    knots = np.r_[0.0, (ends - origin).total_seconds().to_numpy()]
    first = max(start, starts[0])
    grid = pd.date_range(first.ceil(step), min(end, ends[-1]), freq=step)
    if len(grid) < 2:
        return pd.DataFrame()
    at = (grid - origin).total_seconds().to_numpy()
    seconds = np.diff(at)
    out = pd.DataFrame(index=pd.DatetimeIndex(grid[1:], name="Timestamp"))
    for column, amount in running.items():
        total = np.interp(at, knots, np.r_[0.0, np.cumsum(amount)])
        out[column] = np.diff(total) / seconds
    out["Volume_ML"] = out["Inflow_m3s"] * seconds / 1000.0
    out["Level_m"] = np.interp(at[1:], knots[1:], frame["Level_end"].to_numpy(dtype=float))
    out = out.rename(columns={"Release_uncertain": "Release_uncertain_share"})
    return out[["Inflow_m3s", "Inflow_uncorrected_m3s", "Release_m3s", "Level_m", "Volume_ML",
                "Release_uncertain_share"]]


def window_file(out_folder: Path, start, end, step=None) -> Path:
    tail = ""
    if step is not None:
        seconds = int(pd.Timedelta(step).total_seconds())
        tail = f"_{seconds // 3600}h" if seconds % 3600 == 0 else f"_{seconds // 60}min"
    return (Path(out_folder) / "extracts"
            / f"inflow_{pd.Timestamp(start):%Y%m%d%H%M}_{pd.Timestamp(end):%Y%m%d%H%M}{tail}.csv")


def write_window(table: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, float_format="%.4f", date_format="%Y-%m-%d %H:%M:%S")
    return path


def chart_points(table: pd.DataFrame, bins: int = 3000) -> pd.DataFrame:
    """At most ``bins`` rows for a chart: each flow the largest in its bin, so a long
    window still shows every peak; the level its mean. ``Uncertain`` (bool) is where
    the recession correction touched the release."""
    table = table.copy()
    if "Release_uncertain_share" in table:
        table["Uncertain"] = table["Release_uncertain_share"] > 0
    elif "Release_uncertain" in table:
        table["Uncertain"] = table["Release_uncertain"].astype(str).str.lower() == "true"
    if len(table) <= bins:
        return table
    labels = np.arange(len(table)) * bins // len(table)
    grouped = table.groupby(labels)
    flows = [c for c in ("Inflow_m3s", "Inflow_uncorrected_m3s", "Release_m3s", "Uncertain")
             if c in table]
    out = grouped[flows].max()
    out["Level_m"] = grouped["Level_m"].mean()
    out.index = pd.DatetimeIndex(pd.Series(table.index).groupby(labels).min().to_numpy())
    return out


def standard_normal(p: float) -> float:
    from statistics import NormalDist
    return NormalDist().inv_cdf(p) if 0 < p < 1 else math.nan


def scurve_points(target: dict, basis: str) -> tuple:
    """(samples as (z, ML), fitted curve as (z, ML)) for one target and basis."""
    table = pd.read_csv(target["table"])
    column = {"burst": "adv_burst_vol_ML", "storm": "adv_preburst_vol_ML"}[basis]
    position = {"burst": "cunnane_adv_burst", "storm": "cunnane_adv_preburst"}[basis]
    kept = table[table["qualified"].astype(str).str.lower() == "true"]
    samples = [(-standard_normal(p), v) for p, v in zip(kept[position], kept[column])
               if pd.notna(p) and pd.notna(v)]
    fit = target["fits"][basis]
    low, high = math.log10(fit["Vf"]), math.log10(fit["Vc"])
    curve = []
    for step in range(81):
        z = -3.0 + step * 0.075
        value = low + (high - low) / (fit.get("H", 1.0) + math.exp(-fit["k"] * (z - fit["z0"])))
        curve.append((z, 10 ** value))
    return samples, curve
