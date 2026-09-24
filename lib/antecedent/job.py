"""An antecedent storage job: inputs, method settings, the run and the lake configs.

    {
      "homogenise": "homogenise_job.json",       # or the job itself, inline
      "targets": ["rfsl", "fsl"],                 # which homogenised ratings; all if absent
      "rainfall": "rainfall/cld_rain.csv",        # date,rain_mm - day D is the 24 h to 9 am on D
      "ifd": "ifd/ifd_cld.csv",                   # duration_h x "1 in X" columns, depths in mm
      "settings": {"durations_d": [1, 2, 3, 4, 5],
                   "restriction": {"1": 1.15, "2": 1.11, "3": 1.07, "4": 1.05, "5": 1.04},
                   "window_days": 30, "threshold_fraction": 0.8, "preburst_mm": 10,
                   "cunnane_a": 0.4, "round_ml": 1000},
      "bases": ["burst", "storm"],
      "config_name": "lake_config_{basis}_{target}_01.json",
      "out": "antecedent"
    }

The antecedent volumes are the **homogenised** lake's - routed through the
target rating - so the homogenisation runs in this process, from its own job,
exactly as ``extract_antecedent`` ran callide's pipeline in process: reading the
volumes back from the trace file would round them to its four decimals.

``burst`` reads the volume on the day the main burst started, and is the
series for a design storm simulated without its pre-burst; ``storm`` reads it
where the pre-burst started, for a storm with the pre-burst prepended. Each is
fitted with the ceiling at the full supply volume of the target's FSL and
written as the sigmoid lake configuration Bryan's ``LakeConditions`` reads.
"""

from __future__ import annotations

import contextlib
import json
import logging
from pathlib import Path

import pandas as pd

from lib.homogenise import job as homogenise_jobs, peaks

from . import antecedent, scurve

log = logging.getLogger("bryan.antecedent.job")

SETTINGS = {"durations_d": [1, 2, 3, 4, 5],
            "restriction": {"1": 1.15, "2": 1.11, "3": 1.07, "4": 1.05, "5": 1.04},
            "window_days": 30, "threshold_fraction": 0.8, "preburst_mm": 10.0,
            "cunnane_a": 0.4, "round_ml": 1000.0}
BASES = {"burst": "adv_burst_vol_ML", "storm": "adv_preburst_vol_ML"}
CONFIG_NAME = "lake_config_{basis}_{target}_01.json"


class JobError(Exception):
    """An antecedent job that cannot be run, with what to change."""


# -- the settings ----------------------------------------------------------------

def settings_of(job: dict) -> dict:
    out = dict(SETTINGS)
    out.update(job.get("settings") or {})
    durations = [int(d) for d in out["durations_d"]]
    restriction = {int(k): float(v) for k, v in dict(out["restriction"]).items()}
    missing = [d for d in durations if d not in restriction]
    if missing:
        raise JobError(f"no restriction factor for {missing} day(s) - give one per duration "
                       f"(1.0 where no correction applies)")
    out.update(durations_d=durations, restriction=restriction)
    return out


@contextlib.contextmanager
def applied(settings: dict):
    """Bind the settings into the method modules for one run, then put them back.

    The methods read module constants (they were written as one study's code);
    this is the seam that makes them settings without rewriting them.
    """
    names = {(antecedent, "DURATIONS_D"): tuple(settings["durations_d"]),
             (antecedent, "RESTRICTION"): dict(settings["restriction"]),
             (antecedent, "WINDOW_DAYS"): int(settings["window_days"]),
             (antecedent, "THRESHOLD_FRACTION"): float(settings["threshold_fraction"]),
             (antecedent, "RAINFALL_THRESHOLD_MM"): float(settings["preburst_mm"]),
             (antecedent, "CUNNANE_A"): float(settings["cunnane_a"]),
             (scurve, "CUNNANE_A"): float(settings["cunnane_a"]),
             (scurve, "ROUND_ML"): float(settings["round_ml"])}
    saved = {key: getattr(*key) for key in names}
    try:
        for (module, name), value in names.items():
            setattr(module, name, value)
        yield
    finally:
        for (module, name), value in saved.items():
            setattr(module, name, value)


# -- the inputs --------------------------------------------------------------------

def _path(folder: Path, value, what: str) -> Path:
    if not value:
        raise JobError(f"the job names no {what}")
    path = Path(str(value))
    path = path if path.is_absolute() else folder / path
    if not path.is_file():
        raise JobError(f"{what} not found: {path}")
    return path


def read_rainfall(path) -> pd.Series:
    """``date,rain_mm``; '#' lines above it (the launcher's own series has them) ignored."""
    frame = pd.read_csv(path, comment="#")
    if not {"date", "rain_mm"} <= set(frame.columns):
        raise JobError(f"{Path(path).name} needs 'date' and 'rain_mm' columns "
                       f"(it has {', '.join(map(str, frame.columns))})")
    frame["date"] = pd.to_datetime(frame["date"])
    return frame.set_index("date")["rain_mm"].astype(float).sort_index()


def read_ifd(path, durations) -> pd.DataFrame:
    try:
        ifd = antecedent.load_ifd(path)
    except (KeyError, ValueError) as exc:
        raise JobError(f"{Path(path).name} is not an IFD table of duration_h against "
                       f"'1 in X' columns ({exc})") from exc
    missing = [d * 24 for d in durations if d * 24 not in ifd.index]
    if missing:
        raise JobError(f"the IFD has no rows for {missing} h - it needs one per duration")
    if 2.0 not in ifd.columns:
        raise JobError("the IFD needs a '1 in 2' column: the significance test is a "
                       "fraction of the 1 in 2 depth")
    return ifd


def homogenise_job(job: dict, folder: Path):
    source = job.get("homogenise")
    if isinstance(source, dict):
        return homogenise_jobs.from_dict(source, folder)
    return homogenise_jobs.load(_path(folder, source, "homogenisation job"))


# -- the run ------------------------------------------------------------------------

def coefficients(fit) -> dict:
    """The model's five sigmoid coefficients (reports/build_lake_config.py)."""
    return {"k": round(fit["k"], 4), "Vf": float(fit["floor_ML"]),
            "H": round(scurve.asym_parameters(fit)["H"], 4),
            "z0": round(fit["z0"], 4), "Vc": float(fit["ceiling_ML"])}


def lake_config(fit) -> dict:
    """One lake configuration, in the shape Bryan's LakeConditions reads."""
    return {"exceedance_layer_info": [{"lower_z": -99, "upper_z": 99, "type": "sigmoid",
                                       "coefficients": coefficients(fit)}],
            "volume_cap": "none"}


def run(job: dict, folder) -> dict:
    folder = Path(folder)
    settings = settings_of(job)
    rain = read_rainfall(_path(folder, job.get("rainfall"), "rainfall series"))
    ifd = read_ifd(_path(folder, job.get("ifd"), "IFD table"), settings["durations_d"])
    hjob = homogenise_job(job, folder)
    wanted = job.get("targets") or [t.get("name") for t in hjob["targets"]]
    targets = {str(t.get("name") or Path(str(t.get("rating"))).stem): t
               for t in hjob["targets"]}
    unknown = [name for name in wanted if name not in targets]
    if unknown:
        raise JobError(f"the homogenisation job has no target {unknown} (it has "
                       f"{', '.join(targets)})")
    out = Path(job.get("out") or "antecedent")
    out = out if out.is_absolute() else folder / out
    bases = job.get("bases") or list(BASES)
    name_pattern = job.get("config_name") or CONFIG_NAME

    inputs = homogenise_jobs.load_inputs(hjob)
    derived, summary = None, {"notes": list(inputs.notes), "targets": []}
    for name in wanted:
        result, rating, derived = homogenise_jobs.run_target(hjob, inputs, targets[name],
                                                            derived)
        fsv = float(inputs.volume_of_level(rating.fsl))
        ams = peaks.annual_maxima(result, start_month=int(hjob["water_year_start"]))
        daily_volume = antecedent.daily_volume_at_9am(result)
        with applied(settings):
            table = antecedent.antecedent_series(ams, rain, ifd, daily_volume)
            kept = table[table["qualified"]]
            fits, written = {}, {}
            for basis in bases:
                column = BASES[basis]
                fits[basis] = scurve.fit(kept[column].dropna().to_numpy(dtype=float),
                                         ceiling=fsv)
        target_out = out / name
        target_out.mkdir(parents=True, exist_ok=True)
        table.to_csv(target_out / "antecedent_storage.csv", float_format="%.4f")
        rows = []
        for basis, fit in fits.items():
            path = target_out / name_pattern.format(basis=basis, target=name)
            path.write_text(json.dumps(lake_config(fit), indent=2) + "\n", encoding="utf-8")
            written[basis] = str(path)
            rows.append({"basis": basis, "n": fit["n"], "floor_ML": fit["floor_ML"],
                         "ceiling_ML": fit["ceiling_ML"], "k": fit["k"], "z0": fit["z0"],
                         "mean_z0_ML": fit["mean_z0_ML"], "rmse_ML": fit["rmse_ML"]})
        pd.DataFrame(rows).to_csv(target_out / "scurve_params.csv", index=False,
                                  float_format="%.4f")
        no_rain = int(sum(1 for _, peak in ams.iterrows()
                          if pd.Timestamp(peak["New_Level_max_at"]).normalize()
                          < rain.index.min() + pd.Timedelta(days=settings["window_days"])
                          or pd.Timestamp(peak["New_Level_max_at"]) > rain.index.max()))
        summary["targets"].append({
            "name": name, "fsl": float(rating.fsl), "fsv_ML": fsv,
            "years": int(len(table)), "qualified": int(table["qualified"].sum()),
            "years_outside_rainfall": no_rain,
            "fits": {basis: {**coefficients(fit), "n": int(fit["n"]),
                             "mean_z0_ML": float(fit["mean_z0_ML"]),
                             "rmse_ML": float(fit["rmse_ML"])} for basis, fit in fits.items()},
            "configs": written,
            "table": str(target_out / "antecedent_storage.csv")})
    out.mkdir(parents=True, exist_ok=True)
    summary["settings"] = {**settings,
                           "restriction": {str(k): v for k, v in settings["restriction"].items()}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str),
                                      encoding="utf-8")
    return summary


def load_and_run(path) -> dict:
    path = Path(path)
    try:
        job = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise JobError(f"{path} could not be read ({exc})") from exc
    return run(job, path.parent)


