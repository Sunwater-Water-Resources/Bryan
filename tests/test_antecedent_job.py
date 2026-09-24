"""Antecedent storage (lib/antecedent), generalised from callide-fsl-reinstate.

As with the homogenisation, the method is the callide modules copied and the
job layer is new - so the layer is tested here, and the whole is held to the
original: run on Callide's inputs, the job must write the four lake
configurations that were delivered to the design flood model, byte for byte.
"""

from __future__ import annotations

import filecmp
import json
from pathlib import Path

import pandas as pd
import pytest

from lib.antecedent import antecedent, job as jobs, scurve

BRYAN_ROOT = Path(__file__).resolve().parents[1]
CALLIDE = BRYAN_ROOT.parent / "callide-fsl-reinstate"


def test_every_duration_needs_a_restriction_factor():
    with pytest.raises(jobs.JobError, match=r"no restriction factor for \[6\]"):
        jobs.settings_of({"settings": {"durations_d": [1, 6]}})


def test_the_settings_are_bound_for_the_run_and_put_back():
    settings = jobs.settings_of({"settings": {"window_days": 45, "round_ml": 500,
                                              "threshold_fraction": 0.5}})
    with jobs.applied(settings):
        assert antecedent.WINDOW_DAYS == 45 and scurve.ROUND_ML == 500
        # the two defaults captured at definition time in the original now follow
        assert scurve.choose_bounds([1_234.0, 9_876.0]) == (1_000.0, 10_000.0)
    assert antecedent.WINDOW_DAYS == 30 and scurve.ROUND_ML == 1000.0
    assert scurve.choose_bounds([1_234.0, 9_876.0]) == (1_000.0, 10_000.0)


def test_the_significance_fraction_reaches_the_burst_search():
    rain = pd.Series([0.0, 60.0, 0.0], index=pd.date_range("2000-01-01", periods=3))
    ifd = pd.DataFrame({1.582: [50.0], 2.0: [100.0], 10.0: [200.0]},
                       index=pd.Index([24], name="duration_h"))
    with jobs.applied(jobs.settings_of({"settings": {"durations_d": [1],
                                                     "threshold_fraction": 0.5}})):
        assert antecedent.burst_frame(rain, ifd)["over"].any()          # 69 > 50
    with jobs.applied(jobs.settings_of({"settings": {"durations_d": [1]}})):
        assert not antecedent.burst_frame(rain, ifd)["over"].any()      # 69 < 80


def test_an_ifd_without_a_duration_or_the_1_in_2_column_is_refused(tmp_path):
    path = tmp_path / "ifd.csv"
    path.write_text("duration_h,1 in 5,1 in 10\n24,100,120\n", encoding="utf-8")
    with pytest.raises(jobs.JobError, match=r"no rows for \[48\] h"):
        jobs.read_ifd(path, [1, 2])
    with pytest.raises(jobs.JobError, match="'1 in 2' column"):
        jobs.read_ifd(path, [1])


def test_the_launcher_s_rainfall_file_is_read_past_its_header(tmp_path):
    path = tmp_path / "rain.csv"
    path.write_text("# catchment-average daily rainfall\n# from 1 file\n"
                    "date,rain_mm\n2000-01-01,1.5\n2000-01-02,0\n", encoding="utf-8")
    rain = jobs.read_rainfall(path)
    assert rain.index[0] == pd.Timestamp("2000-01-01") and rain.iloc[0] == 1.5
    path.write_text("day,mm\n2000-01-01,1\n", encoding="utf-8")
    with pytest.raises(jobs.JobError, match="needs 'date' and 'rain_mm'"):
        jobs.read_rainfall(path)


def test_a_lake_config_has_the_shape_lake_conditions_reads():
    config = jobs.lake_config({"k": 1.23456, "z0": -0.54321, "floor_ML": 12_000.0,
                               "ceiling_ML": 129_041.0})
    (layer,) = config["exceedance_layer_info"]
    assert layer["type"] == "sigmoid" and (layer["lower_z"], layer["upper_z"]) == (-99, 99)
    assert layer["coefficients"] == {"k": 1.2346, "Vf": 12_000.0, "H": 1.0, "z0": -0.5432,
                                     "Vc": 129_041.0}
    assert config["volume_cap"] == "none"


@pytest.mark.skipif(not (CALLIDE / "antecedent_storage" / "BryanLakeConfig").is_dir(),
                    reason="needs the callide-fsl-reinstate checkout beside Bryan")
def test_callide_s_delivered_lake_configs_come_back_byte_for_byte(tmp_path):
    data = CALLIDE / "data"
    homogenise = {
        "gauges": [str(data / "gauge" / "130314A.csv"), str(data / "gauge" / "130314C.csv")],
        "overlay": {"file": str(data / "gauge" / "130314B_hydstra.csv"), "below": 200.90,
                    "reconnect_margin": 0.10},
        "storage": str(data / "storage" / "CALLIDE_STORAGE.els"),
        "register": str(data / "ratings" / "RatingCurves.xlsx"),
        "evaporation": str(data / "climate" / "silo_-24.35_150.65.txt"),
        "targets": [{"name": "RFSL", "rating": str(data / "ratings" / "CALLIDE_RFSL.sq"),
                     "fsl": 215.5},
                    {"name": "FSL", "rating": str(data / "ratings" / "CALLIDE_FSL.sq"),
                     "fsl": 216.1}]}
    local = CALLIDE / "antecedent_storage" / "data"
    summary = jobs.run({"homogenise": homogenise,
                        "rainfall": str(local / "rainfall_cld_awap.csv"),
                        "ifd": str(local / "ifd_cld_sunwater.csv"),
                        "out": str(tmp_path / "antecedent")}, tmp_path)
    assert [t["qualified"] for t in summary["targets"]] == [24, 24]
    delivered = CALLIDE / "antecedent_storage" / "BryanLakeConfig"
    for condition in ("RFSL", "FSL"):
        for basis in ("burst", "storm"):
            name = f"lake_config_{basis}_{condition}_01.json"
            assert filecmp.cmp(tmp_path / "antecedent" / condition / name, delivered / name,
                               shallow=False), name
            json.loads((delivered / name).read_text(encoding="utf-8"))
