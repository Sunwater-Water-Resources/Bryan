"""The sub-burst neutrality check, and calibrated pattern weights, for a group of runs.

The run launcher's Results page (the Sub-bursts tab) shells out to this, because the
total probability theorem needs scipy and the launcher's environment has none. It can
equally be run by hand against a job file:

    python util/SubBurstCheck.py --job subburst_job.json --results subburst_results.json

For each storm duration of the tested group - and of a comparison group, if the job
names one - it applies the TPT to the sub-burst depths recorded in the Monte Carlo
database and divides by the same-AEP IFD depth: the neutrality margins of
Manual/SubDocs/sub_burst_check.md, computed exactly as ``Simulator.analyse_sub_bursts``
computes them. Nothing is re-run and nothing is written beside the runs.

With ``calibrate`` set it also calibrates temporal pattern weights to neutrality for each
duration of the tested group, with ``util/CalibrateTpWeights.py``'s own ``calibrate``,
and - where the database holds lake levels - returns the design level curve of the
tested group with those weights, beside the same group unweighted and the comparison
group: the exceedance probability at every level of a fine grid, so the page can read
off the AEP of any level (the dam crest, say) without running this again. The design
curve is the envelope over the durations. The weights for each duration are written to
``weights_folder`` in the form the simulation list's ``TP weights`` column reads.

The job is JSON; every path in it is absolute:

    {"scheme": {lower_aep, upper_aep, number_of_main_divisions,
                number_of_sub_divisions, number_of_temporal_patterns, aep_of_pmp},
     "tested":  {"label": "...", "databases": {"72h": "C:/.../x__mcdf.csv", ...}},
     "compare": {"label": "...", "databases": {...}} or null,
     "calibrate": false,
     "calibration": {"margin_target": 1.0, "weight_floor": 0.02,
                     "reduction_factor": 0.5, "max_iterations": 25, "aep_range": null},
     "weights_folder": "C:/.../_subburst/weights",
     "fingerprint": "..."}                  # echoed into the results
"""

import argparse
import contextlib
import importlib.util
import io
import json
import math
import os
import re
import sys

BRYAN_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BRYAN_ROOT not in sys.path:
    sys.path.insert(0, BRYAN_ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from lib.MCScheme import SampleScheme, TotalProbTheorem, neutrality_margin  # noqa: E402

CALIBRATION_DEFAULTS = {"margin_target": 1.0, "weight_floor": 0.02, "reduction_factor": 0.5,
                        "max_iterations": 25, "aep_range": None}
# The level grid for the design curves: round 5 mm steps, so that a level given to the
# centimetre (a crest, a spillway sill) is a grid point and its AEP is the TPT's own
# value rather than an interpolation along a curve that is stepped between realisations.
LEVEL_STEP = 0.005
LEVEL_GRID_MAX_POINTS = 20001


def _calibrate_module():
    """util/CalibrateTpWeights.py as a module: its functions, with its settings set per call."""
    path = os.path.join(BRYAN_ROOT, "util", "CalibrateTpWeights.py")
    spec = importlib.util.spec_from_file_location("CalibrateTpWeights", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scheme_of(job):
    s = job["scheme"]
    mc = SampleScheme(s["lower_aep"], s["upper_aep"], int(s["number_of_main_divisions"]),
                      int(s["number_of_sub_divisions"]), int(s.get("number_of_temporal_patterns", 10)))
    mc.aep_of_pmp = s.get("aep_of_pmp")
    return mc


def read_database(path):
    if str(path).lower().endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path, index_col=0)


def sub_columns(mcdf):
    return [c for c in mcdf.columns if c.startswith("subburst_") and f"ifd_{c[len('subburst_'):]}" in mcdf.columns]


def margins(mcdf, mc):
    """{sub-duration: {AEP: margin}} - Simulator.analyse_sub_bursts' arithmetic.

    On a copy: compute_std_quantiles adds a ``<column>_aep`` column per sub-duration to the
    frame it is given, and CalibrateTpWeights reads every ``subburst_`` column as one to
    calibrate against, so the database the calibration sees must be the one on disk.
    """
    mcdf = mcdf.copy()
    mc.df = mcdf
    out = {}
    for col in sub_columns(mcdf):
        sub = col[len("subburst_"):]
        with contextlib.redirect_stdout(io.StringIO()):
            mc.compute_std_quantiles(result_type=col)
        reference = mcdf[["rain_z", "storm_method", f"ifd_{sub}"]].dropna()
        series = neutrality_margin(mc.quantiles[col], reference, col, f"ifd_{sub}")
        out[sub] = {str(int(aep) if float(aep).is_integer() else aep): _num(value)
                    for aep, value in series.items()}
    return out


def exceedance(mcdf, mc, grid, weights=None):
    """P(level > x) on ``grid`` for one duration, by the TPT, optionally pattern-weighted."""
    data = mcdf[["m", "level"]].copy()
    data["tp_w"] = 1.0 if weights is None else weights
    tpt = TotalProbTheorem(mc.m, mc.n, mc.main_divisions, data, verbose=False)
    return tpt.assign_aep_all(grid, "level")


def level_grid(low, high):
    start, stop = math.floor(low / LEVEL_STEP) * LEVEL_STEP, math.ceil(high / LEVEL_STEP) * LEVEL_STEP
    if (stop - start) / LEVEL_STEP + 1 > LEVEL_GRID_MAX_POINTS:     # a very wide range: evenly spaced
        return np.linspace(start, stop, LEVEL_GRID_MAX_POINTS)
    return np.round(np.arange(start, stop + LEVEL_STEP / 2, LEVEL_STEP), 3)


def level_at(p_curve, grid, aep):
    order = np.argsort(p_curve, kind="mergesort")
    p = 1.0 / aep
    if p < p_curve[order][0] or p > p_curve[order][-1]:
        return None
    return float(np.interp(p, p_curve[order], grid[order]))


def _num(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) or math.isinf(value) else round(value, 4)


def _safe(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_") or "group"


def calibrate_duration(ctw, mcdf, mc, label, settings, folder):
    """CalibrateTpWeights.calibrate on one duration, its settings taken from the job."""
    for key, value in {**CALIBRATION_DEFAULTS, **(settings or {})}.items():
        setattr(ctw, key, value)
    os.makedirs(folder, exist_ok=True)
    ctw.output_folder = folder
    prepared = ctw.prepare_mcdf(mcdf.copy())
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        weights, _ = ctw.calibrate(prepared, mc, label)
    text = log.getvalue()
    iterations = [line for line in text.splitlines() if line.startswith("Iteration")]
    ctw.attach_weights(prepared, weights)
    floor = ctw.weight_floor
    return prepared["tp_w"].to_numpy(dtype=float), {
        "converged": "Converged" in text,
        "iterations": max(len(iterations) - 1, 0),
        "worst_before": _num(iterations[0].split("= ")[-1]) if iterations else None,
        "worst_after": _num(iterations[-1].split("= ")[-1]) if iterations else None,
        "weights": {group: [round(float(w), 4) for w in values] for group, values in weights.items()},
        # a weight within 1.5 x the floor is at it: renormalising lifts floored weights a little
        "at_floor": {group: int(sum(1 for w in values if w <= floor * 1.5)) for group, values in weights.items()},
        "weights_file": os.path.join(folder, f"{label}_tp_weights.csv"),
    }


def run(job, progress=print):
    mc = scheme_of(job)
    std_aeps = mc.get_standard_aeps()
    results = {"fingerprint": job.get("fingerprint"), "standard_aeps": [float(a) for a in std_aeps],
               "groups": {}, "problems": [], "calibration": None}
    frames = {}
    roles = [("tested", job["tested"])] + ([("compare", job["compare"])] if job.get("compare") else [])
    total = sum(len(spec["databases"]) for _, spec in roles)
    done = 0
    for role, spec in roles:
        entry = {"label": spec.get("label", role), "durations": {}}
        for duration, path in spec["databases"].items():
            done += 1
            try:
                mcdf = read_database(path)
            except (OSError, ValueError) as exc:
                results["problems"].append(f"{spec.get('label')} {duration}: could not read {path} ({exc})")
                continue
            if not sub_columns(mcdf):
                results["problems"].append(
                    f"{spec.get('label')} {duration}: no sub-burst columns - made by a Bryan version from "
                    "before sub-burst tracking; re-run the storms to record them")
                continue
            frames[(role, duration)] = mcdf
            entry["durations"][duration] = {"margins": margins(mcdf, mc), "database": path,
                                            "has_level": bool("level" in mcdf.columns and mcdf["level"].notna().any())}
            progress(f"[{done}/{total}] {spec.get('label')} {duration}: neutrality margins")
        results["groups"][role] = entry

    if job.get("calibrate"):
        ctw = _calibrate_module()
        label = _safe(job["tested"].get("label", "tested"))
        calibration = {"durations": {}}
        weights_by_duration = {}
        tested = [(d, f) for (role, d), f in frames.items() if role == "tested"]
        for n, (duration, mcdf) in enumerate(tested, start=1):
            progress(f"[{n}/{len(tested)}] calibrating the pattern weights for {duration}")
            try:
                weights, summary = calibrate_duration(ctw, mcdf, mc, f"{label}_{_safe(duration)}",
                                                      job.get("calibration"), job["weights_folder"])
            except Exception as exc:                                # noqa: BLE001 - report per duration
                results["problems"].append(f"calibration {duration}: {exc}")
                continue
            weights_by_duration[duration] = weights
            calibration["durations"][duration] = summary

        with_level = {key: f for key, f in frames.items() if "level" in f.columns and f["level"].notna().any()}
        if with_level and weights_by_duration:
            levels = np.concatenate([f["level"].dropna().to_numpy() for f in with_level.values()])
            grid = level_grid(float(levels.min()), float(levels.max()))
            cases = {"tested": [], "calibrated": [], "compare": []}
            for (role, duration), f in with_level.items():
                if role == "tested":
                    cases["tested"].append(exceedance(f, mc, grid))
                    if duration in weights_by_duration:
                        cases["calibrated"].append(exceedance(f, mc, grid, weights_by_duration[duration]))
                else:
                    cases["compare"].append(exceedance(f, mc, grid))
            curves = {}
            for case, parts in cases.items():
                if not parts:
                    continue
                envelope = np.max(np.vstack(parts), axis=0)
                curves[case] = {
                    "levels": {str(int(a) if float(a).is_integer() else a): _num(level_at(envelope, grid, a))
                               for a in std_aeps},
                    # the whole curve, for the AEP of any level: P(level > x), 6 significant figures
                    "p": [float(f"{value:.6g}") for value in envelope],
                }
            calibration["curves"] = curves
            calibration["grid"] = [round(float(x), 4) for x in grid]
            if len(weights_by_duration) < len(tested):
                results["problems"].append("the calibrated curve leaves out the durations whose calibration failed")
        results["calibration"] = calibration
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--job", required=True)
    parser.add_argument("--results", required=True)
    args = parser.parse_args(argv)
    with open(args.job, encoding="utf-8") as f:
        job = json.load(f)
    results = run(job, progress=lambda text: print(text, flush=True))
    temporary = args.results + ".tmp"
    with open(temporary, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    os.replace(temporary, args.results)
    print(f"wrote {args.results}", flush=True)
    for problem in results["problems"]:
        print("PROBLEM:", problem, flush=True)


if __name__ == "__main__":
    main()
