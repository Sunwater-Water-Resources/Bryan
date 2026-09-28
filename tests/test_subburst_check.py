"""util/SubBurstCheck.py: the script behind the launcher's Sub-bursts tab.

A synthetic Monte Carlo database stands in for a real one: ten temporal patterns, of
which three carry embedded bursts a third deeper than the IFD and the rest sit below
it. So the unfiltered ensemble must show margins above 1, calibration must take the
weight off those three patterns and nothing else, and the calibrated level curve must
come down - the offending patterns are also the ones that make the higher levels.

This is a **Bryan** test - it imports lib/MCScheme.py, which needs scipy. Run it with
Bryan's own interpreter:

    python -m pytest tests -q
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

from scipy.special import ndtri                                      # noqa: E402

SCHEME = {"lower_aep": 2, "upper_aep": 10000, "number_of_main_divisions": 6,
          "number_of_sub_divisions": 120, "number_of_temporal_patterns": 10}
OFFENDERS = (0, 1, 2)


def check_module():
    spec = importlib.util.spec_from_file_location("SubBurstCheck", BRYAN_ROOT / "util" / "SubBurstCheck.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def database(path, *, offending=True, seed=1):
    """One storm duration's mcdf: sub-burst depths per window, the IFD beside them, and a level."""
    rng = np.random.default_rng(seed)
    m, n = SCHEME["number_of_main_divisions"], SCHEME["number_of_sub_divisions"]
    edges = np.linspace(ndtri(1 - 1 / SCHEME["lower_aep"]), ndtri(1 - 1 / SCHEME["upper_aep"]), m + 1)
    rows = []
    for division in range(m):
        z = rng.uniform(edges[division], edges[division + 1], n)
        tp = rng.integers(0, 10, n)
        bursty = np.isin(tp, OFFENDERS) & offending
        factor = np.where(bursty, 1.35, 0.85)
        frame = pd.DataFrame({"m": division, "n": np.arange(n), "rain_z": z, "tp": tp,
                              "storm_method": "GSDM", "tp_frequency": ""})
        for window, scale in (("6h", 60.0), ("12h", 80.0)):
            ifd = scale * (1 + 0.35 * z)
            frame[f"ifd_{window}"] = ifd
            frame[f"subburst_{window}"] = ifd * factor
        frame["level"] = 200 + 0.8 * z + np.where(bursty, 0.6, 0.0)
        rows.append(frame)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(path)
    return path


def job(tmp_path, *, compare=True, calibrate=False):
    tested = {"24h": database(tmp_path / "a_24h__mcdf.csv", seed=1),
              "48h": database(tmp_path / "a_48h__mcdf.csv", seed=2)}
    other = {"24h": database(tmp_path / "b_24h__mcdf.csv", offending=False, seed=3)}
    return {"scheme": SCHEME,
            "tested": {"label": "unfiltered", "databases": {d: str(p) for d, p in tested.items()}},
            "compare": ({"label": "filtered", "databases": {d: str(p) for d, p in other.items()}}
                        if compare else None),
            "calibrate": calibrate, "calibration": {"max_iterations": 12},
            "weights_folder": str(tmp_path / "weights"), "fingerprint": "abc"}


def worst(margins, low=5, high=5000):
    return max(v for series in margins.values() for k, v in series.items()
               if v is not None and low <= float(k) <= high)


def test_embedded_bursts_that_are_too_deep_show_as_margins_above_one(tmp_path):
    results = check_module().run(job(tmp_path), progress=lambda _: None)

    tested = results["groups"]["tested"]["durations"]
    compared = results["groups"]["compare"]["durations"]
    assert set(tested) == {"24h", "48h"} and set(compared) == {"24h"}
    assert set(tested["24h"]["margins"]) == {"6h", "12h"}
    assert worst(tested["24h"]["margins"]) > 1.05
    assert worst(compared["24h"]["margins"]) < 1.0
    assert results["calibration"] is None and results["problems"] == []
    json.dumps(results, allow_nan=False)


def test_the_margins_leave_the_database_as_it_was(tmp_path):
    """compute_std_quantiles adds an '_aep' column per window to the frame it is given;
    the calibration would then read that as a window to calibrate."""
    module = check_module()
    frame = pd.read_csv(database(tmp_path / "x__mcdf.csv"), index_col=0)
    before = list(frame.columns)
    module.margins(frame, module.scheme_of({"scheme": SCHEME}))
    assert list(frame.columns) == before


def test_calibration_takes_the_weight_off_the_offending_patterns(tmp_path):
    results = check_module().run(job(tmp_path, calibrate=True), progress=lambda _: None)

    calibration = results["calibration"]
    entry = calibration["durations"]["24h"]
    weights = np.array(entry["weights"]["GSDM"])
    assert weights[list(OFFENDERS)].max() < weights[3:].min()
    assert entry["worst_before"] > entry["worst_after"]
    assert Path(entry["weights_file"]).is_file()

    curves = calibration["curves"]
    assert set(curves) == {"tested", "calibrated", "compare"}
    assert len(curves["tested"]["p"]) == len(calibration["grid"])
    # the offending patterns also make the higher levels, so down-weighting them lowers the curve
    assert curves["calibrated"]["levels"]["100"] < curves["tested"]["levels"]["100"]
    json.dumps(results, allow_nan=False)


def test_a_database_from_before_sub_burst_tracking_is_reported_not_fatal(tmp_path):
    spec = job(tmp_path, compare=False)
    old = tmp_path / "old__mcdf.csv"
    pd.read_csv(spec["tested"]["databases"]["48h"], index_col=0) \
        .drop(columns=["subburst_6h", "subburst_12h"]).to_csv(old)
    spec["tested"]["databases"]["48h"] = str(old)

    results = check_module().run(spec, progress=lambda _: None)
    assert list(results["groups"]["tested"]["durations"]) == ["24h"]
    assert any("before sub-burst tracking" in p for p in results["problems"])


def test_the_command_line_writes_the_results_file(tmp_path):
    module = check_module()
    job_path = tmp_path / "job.json"
    job_path.write_text(json.dumps(job(tmp_path, compare=False)), encoding="utf-8")
    results_path = tmp_path / "results.json"
    module.main(["--job", str(job_path), "--results", str(results_path)])
    assert json.loads(results_path.read_text(encoding="utf-8"))["fingerprint"] == "abc"
