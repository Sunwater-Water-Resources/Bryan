"""The Lake record page's settings, jobs and results - everything but the view.

Three steps, in the order they depend on each other, all kept in the study file
under ``"lake_record"`` with paths relative to it:

1. **Catchment rainfall** - ``core/awap.py``, in this process: a shapefile and
   the folder of AWAP/AWRA-L daily grids give ``date,rain_mm``.
2. **Homogenisation** - ``util/HomogeniseLakeLevels.py`` under Bryan's
   interpreter (scipy): the gauge exports, storage table, rating register,
   evaporation and target ratings give the lake re-routed through each target.
3. **Antecedent storage** - ``util/AntecedentStorage.py`` under Bryan's
   interpreter: the rainfall, the IFD and the homogenisation give the antecedent
   series and the ``lake_config.json`` files Bryan's Monte Carlo scheme samples.

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

import pandas as pd

from .bryan import BRYAN_ROOT
from .paths import atomic_write_json, read_json
from .study import Study, portable, resolve

KEY = "lake_record"
HOMOGENISE_SCRIPT = BRYAN_ROOT / "util" / "HomogeniseLakeLevels.py"
ANTECEDENT_SCRIPT = BRYAN_ROOT / "util" / "AntecedentStorage.py"
TIMEOUT_SECONDS = 3600

MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
CALLIDE_PAN_FACTORS = [0.82, 0.82, 0.83, 0.79, 0.75, 0.71, 0.76, 0.81, 0.79, 0.81, 0.82, 0.83]

DEFAULTS = {
    "rainfall": {"shapefile": "", "field": "", "value": "", "grids": "", "pattern": "*.nc",
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
}


# -- settings ---------------------------------------------------------------------

def settings(study: Study) -> dict:
    """The study's section, completed with the defaults, one level deep."""
    stored = study.extra.get(KEY) or {}
    out = copy.deepcopy(DEFAULTS)
    for step, defaults in out.items():
        defaults.update(copy.deepcopy(stored.get(step) or {}))
    out["antecedent"]["settings"] = {**DEFAULTS["antecedent"]["settings"],
                                     **(out["antecedent"].get("settings") or {})}
    return out


def store(study: Study, section: dict) -> None:
    study.extra[KEY] = section


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


def problems_before_running(study: Study, section: dict, step: str) -> list:
    """What would stop a step, said before a subprocess is started."""
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
        need("Folder of rainfall grids", r["grids"], is_dir=True)
    elif step in ("homogenise", "antecedent"):
        h = section["homogenise"]
        if not [g for g in h["gauges"] if str(g).strip()]:
            out.append("Gauge exports: none given")
        for number, gauge in enumerate(g for g in h["gauges"] if str(g).strip()):
            need(f"Gauge export {number + 1}", gauge)
        need("Storage table (.els)", h["storage"])
        need("Rating register (.xlsx)", h["register"])
        need("Evaporation (SILO)", h["evaporation"])
        if h.get("overlay"):
            need("Overlay gauge", h["overlay"].get("file"))
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
