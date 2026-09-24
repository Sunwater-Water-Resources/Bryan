"""The Lake record page: its settings in the study, the jobs it writes, its results."""

from __future__ import annotations

import asyncio
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from core import lakerecord
from core import reporttables as rt
from core import study as studies
from report_fixtures import build_study


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def touch(path, text="x"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def configured(study):
    folder = study.folder
    for name in ("gauge/A.csv", "gauge/C.csv", "storage.els", "register.xlsx", "silo.txt",
                 "rfsl.sq", "ifd.csv", "lake_record/rainfall.csv"):
        touch(folder / name)
    section = lakerecord.settings(study)
    h = section["homogenise"]
    h.update(gauges=["gauge/A.csv", "gauge/C.csv"], storage="storage.els",
             register="register.xlsx", evaporation="silo.txt",
             targets=[{"name": "RFSL", "rating": "rfsl.sq", "fsl": 215.5}])
    section["antecedent"]["ifd"] = "ifd.csv"
    return section


def test_the_settings_complete_with_defaults_and_keep_what_was_stored(study):
    study.extra["lake_record"] = {"antecedent": {"settings": {"window_days": 45}}}
    section = lakerecord.settings(study)
    assert section["antecedent"]["settings"]["window_days"] == 45
    assert section["antecedent"]["settings"]["threshold_fraction"] == 0.8
    assert section["homogenise"]["pan_factors"] == lakerecord.CALLIDE_PAN_FACTORS
    assert section["rainfall"]["weighting"] == "area"


def test_the_homogenisation_job_has_absolute_paths_and_monthly_pan_factors(study):
    section = configured(study)
    section["homogenise"]["overlay"] = {"file": "gauge/C.csv", "below": 200.9}
    job = lakerecord.homogenise_job(study, section)
    assert job["gauges"] == [str((study.folder / "gauge/A.csv").resolve()),
                             str((study.folder / "gauge/C.csv").resolve())]
    assert job["pan_factors"]["1"] == 0.82 and len(job["pan_factors"]) == 12
    assert job["overlay"]["reconnect_margin"] == 0.10
    assert job["targets"] == [{"name": "RFSL", "fsl": 215.5,
                               "rating": str((study.folder / "rfsl.sq").resolve())}]


def test_the_antecedent_job_takes_step_one_s_rainfall_by_default(study, tmp_path):
    section = configured(study)
    job = lakerecord.antecedent_job(study, section, tmp_path / "h.json")
    assert job["rainfall"] == str((study.folder / "lake_record/rainfall.csv").resolve())
    assert job["homogenise"] == str(tmp_path / "h.json")
    assert job["settings"]["restriction"]["1"] == 1.15


def test_what_would_stop_a_step_is_said_before_anything_runs(study):
    section = lakerecord.settings(study)
    problems = lakerecord.problems_before_running(study, section, "antecedent")
    assert "Gauge exports: none given" in problems
    assert "Target ratings: none given" in problems
    assert "IFD table: not given" in problems
    assert not lakerecord.problems_before_running(study, configured(study), "antecedent")


def test_two_targets_of_one_name_are_refused(study):
    section = configured(study)
    section["homogenise"]["targets"].append({"name": "RFSL", "rating": "rfsl.sq", "fsl": 1})
    assert "Target ratings: each needs its own name" in \
        lakerecord.problems_before_running(study, section, "homogenise")


def test_a_run_writes_its_job_beside_the_outputs_and_reads_the_summary(tmp_path, monkeypatch):
    script = touch(tmp_path / "fake.py",
                   "import json, sys, pathlib\n"
                   "job = json.loads(pathlib.Path(sys.argv[1]).read_text())\n"
                   "out = pathlib.Path(job['out']); out.mkdir(parents=True, exist_ok=True)\n"
                   "(out / 'summary.json').write_text(json.dumps({'targets': [], 'notes': ['ok']}))\n")
    job = {"out": str(tmp_path / "out")}
    result = lakerecord.run_script(script, job, tmp_path / "out" / "job.json", sys.executable)
    assert result.ok and result.summary["notes"] == ["ok"]
    assert json.loads((tmp_path / "out" / "job.json").read_text()) == job


def test_a_failed_run_carries_its_output(tmp_path):
    script = touch(tmp_path / "fail.py", "print('ERROR: no gauges'); raise SystemExit(1)\n")
    result = lakerecord.run_script(script, {"out": str(tmp_path)}, tmp_path / "j.json",
                                   sys.executable)
    assert not result.ok and "ERROR: no gauges" in result.output


def test_the_s_curve_is_drawn_from_the_samples_and_the_fitted_coefficients(tmp_path):
    table = tmp_path / "antecedent_storage.csv"
    pd.DataFrame({"WaterYear": [2000, 2001, 2002], "qualified": [True, False, True],
                  "adv_burst_vol_ML": [50_000.0, None, 100_000.0],
                  "cunnane_adv_burst": [0.7, None, 0.3]}).to_csv(table, index=False)
    target = {"table": str(table), "fits": {"burst": {"k": 1.5, "z0": -0.5, "H": 1.0,
                                                       "Vf": 12_000.0, "Vc": 129_041.0}}}
    samples, curve = lakerecord.scurve_points(target, "burst")
    assert len(samples) == 2 and samples[1][0] > samples[0][0]     # rarer, higher z
    assert curve[0][1] == pytest.approx(12_000, rel=0.1)            # near the floor far left
    assert curve[-1][1] < 129_041 and curve[-1][1] > 100_000        # towards the ceiling


# -- the page ------------------------------------------------------------------------

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")


@pytest_asyncio.fixture
async def user():
    from nicegui.testing.user_simulation import user_simulation
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def opened(study, monkeypatch, tmp_path):
    import settings as settings_module
    from state import STATE
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")
    STATE.open_study(study.path)
    yield study
    STATE.study = None


@pytest.mark.asyncio
async def test_the_page_shows_the_three_steps_and_refuses_to_run_unready(user, opened):
    await user.open("/lake-record")
    for heading in ("1. Catchment rainfall", "2. Homogenisation", "3. Antecedent storage"):
        await user.should_see(heading)
    user.find(marker="run-homogenise").click()
    await user.should_see("Gauge exports: none given")


@pytest.mark.asyncio
async def test_a_path_typed_on_the_page_is_kept_relative_to_the_study(user, opened):
    await user.open("/lake-record")
    box = user.find(marker="ifd")
    box.clear().type(str(opened.folder / "ifd" / "cld.csv"))
    box.trigger("blur")
    await asyncio.sleep(0.2)
    stored = lakerecord.settings(studies.load_study(opened.path))
    assert stored["antecedent"]["ifd"] == "ifd/cld.csv"
