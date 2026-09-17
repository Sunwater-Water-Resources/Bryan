"""core/lakefreq.py: the Lake levels page's settings, job, chart and exports."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from conftest import MONTE_CARLO_COLUMNS, write_workbook
from core import events, lakefreq
from lake_fixtures import FSL, synthetic_level, write_hydstra, write_mcdf

DURATIONS = (12, 48)


@pytest.fixture
def project(tmp_path):
    """A project with a level record beside it and two Monte Carlo durations."""
    folder = tmp_path / "project"
    results = folder / "sims_mc" / "results"
    results.mkdir(parents=True)
    rows = []
    for hours in DURATIONS:
        name = f"DAM_mc_{hours}h_GWL0p3"
        write_mcdf(results / f"{name}__mcdf.csv", durations_shift=hours / 100)
        rows.append({"Include": "yes", "Method": "monte carlo", "Duration": hours,
                     "GWL": 0.3, "Output file": f"sims_mc\\results\\{name}",
                     "Run models": "no", "Analyse results": "yes"})
    write_workbook(folder / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    config = folder / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    write_hydstra(folder / "gauge" / "HW.csv", synthetic_level(years=25))

    from state import Project
    return Project.open(config)


def settings_for(project, **changes):
    settings = lakefreq.load_settings(project.config.config_path)
    settings["record"]["files"] = ["gauge/HW.csv"]
    settings["fsl"] = FSL
    for key, value in changes.items():
        settings[key] = value
    return settings


def test_settings_live_beside_the_config_and_keep_project_paths_relative(project, tmp_path):
    config = project.config.config_path
    outside = tmp_path / "elsewhere" / "HW.csv"
    settings = settings_for(project)
    settings["record"]["files"] = [str(config.parent / "gauge" / "HW.csv"), str(outside)]
    path = lakefreq.save_settings(config, settings)

    assert path == config.parent / "lake_frequency.json"
    stored = json.loads(path.read_text())
    assert stored["record"]["files"] == ["gauge/HW.csv", str(outside)]
    assert lakefreq.load_settings(config)["fsl"] == FSL


def test_the_job_resolves_paths_and_picks_the_durations_asked_for(project):
    config = project.config.config_path
    groups = events.sources_by_group(project)
    group = next(iter(groups))
    settings = settings_for(project, design={"include": True, "group": group,
                                             "durations": [48.0]})
    plan = lakefreq.build_job(config, settings, groups)

    assert plan.can_fit, plan.problems + plan.fit_problems
    assert plan.job["record"]["files"] == [str(config.parent / "gauge" / "HW.csv")]
    assert [source["duration"] for source in plan.job["design"]["sources"]] == [48.0]
    assert Path(plan.job["design"]["sources"][0]["path"]).name == "DAM_mc_48h_GWL0p3__mcdf.csv"

    everything = settings_for(project, design={"include": True, "group": group,
                                               "durations": None})
    assert len(lakefreq.build_job(config, everything, groups).job["design"]["sources"]) == 2


def test_what_stops_the_record_is_kept_apart_from_what_stops_only_the_curves(project):
    config = project.config.config_path
    plan = lakefreq.build_job(config, settings_for(project, fsl=None))
    assert plan.can_read and not plan.can_fit
    assert "full supply level" in plan.fit_problems[0]

    empty = lakefreq.load_settings(config)
    assert not lakefreq.build_job(config, empty).can_read

    missing = settings_for(project)
    missing["record"]["files"] = ["gauge/nope.csv"]
    assert "Not found" in lakefreq.build_job(config, missing).problems[0]


def test_the_view_is_read_once_per_file_and_rederived_per_setting(project, monkeypatch):
    config = project.config.config_path
    lakefreq.forget_records()
    job = lakefreq.build_job(config, settings_for(project)).job
    first = lakefreq.read_view(job)
    # 25 years of readings run a few hours into a 26th, which is kept and flagged.
    assert len(first.ams) == 26 and "no quality codes" in first.notes[0]
    assert "1 water year is not fully covered (2005-06)" in first.notes[-1]

    monkeypatch.setattr(lakefreq.RECORD, "read_record",
                        lambda *_: pytest.fail("read the export twice"))
    moved = lakefreq.read_view({**job, "water_year_start": 7})
    assert list(moved.ams["period"])[0].startswith("1980-")


def test_the_chart_shows_the_points_before_there_is_a_fit(project):
    config = project.config.config_path
    settings = settings_for(project, fsl_label="FSL",
                            reference_levels=[{"label": "crest", "level": 218.0}])
    view = lakefreq.read_view(lakefreq.build_job(config, settings).job)
    options = lakefreq.chart_options(view.positions, None, settings)

    names = [series["name"] for series in options["series"]]
    assert names == ["Storm-driven maxima", "Carried over the water year"]
    marks = options["series"][0]["markLine"]["data"]
    assert [mark["yAxis"] for mark in marks] == [218.0, FSL]
    labels = options["xAxis"]["axisLabel"][":formatter"]
    assert "1EY" in labels and "1 in 10" in labels


def results_for(positions):
    """The shape util/LakeLevelFrequency.py writes, without running it."""
    grid = [-0.3 + 0.1 * i for i in range(30)]
    block = {"form": "shouldered", "rmse": 0.2, "curve": [214.0] * 30,
             "band_lo": [213.0] * 30, "band_hi": [215.0] * 30, "draws_used": 380,
             "z_max": 1.0, "error": None}
    return {"grid": {"z": grid, "aep": [0.5] * 30},
            "fits": {"all": block, "storm": dict(block, rmse=0.3)},
            "design": {"durations": {"12": {"aep": [0.4, 0.01], "z": [0.25, 2.3],
                                            "level": [210.0, 216.0]}},
                       "envelope": [215.0] * 30}}


def test_the_chart_carries_the_fit_its_band_and_the_design_floods(project):
    config = project.config.config_path
    settings = settings_for(project)
    view = lakefreq.read_view(lakefreq.build_job(config, settings).job)
    results = results_for(view.positions)
    options = lakefreq.chart_options(view.positions, results, settings)
    names = [series["name"] for series in options["series"]]

    assert "90% band, all maxima (380 resamples)" in names
    assert "Design flood envelope" in names and "Monte Carlo durations" in names
    assert "Fit to all maxima (RMSE 0.20 m)" in names
    assert "Fit to storm-driven (RMSE 0.30 m)" in names
    assert "band" not in options["legend"]["data"]
    # The shouldered upper limb is not drawn past the rarest maximum it was fitted to.
    fit = next(s for s in options["series"] if s["name"].startswith("Fit to all"))
    assert max(point[0] for point in fit["data"]) <= 1.0

    hidden = lakefreq.chart_options(view.positions, results, settings, show_design=False)
    assert "Design flood envelope" not in [s["name"] for s in hidden["series"]]


def test_the_cached_results_are_only_the_ones_for_this_job(project):
    config = project.config.config_path
    job = lakefreq.build_job(config, settings_for(project)).job
    path = lakefreq.results_path(config, job)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"fingerprint": lakefreq.RECORD.fingerprint(job)}))

    assert lakefreq.cached_results(config, job) is not None
    assert lakefreq.cached_results(config, {**job, "fsl": FSL + 1}) is None


def test_the_command_runs_the_util_script_with_what_was_asked_for(tmp_path):
    argv = lakefreq.command("py", tmp_path / "job.json", results=tmp_path / "r.json",
                            png=tmp_path / "f.png", without_design=True)
    assert argv[:3] == ["py", "-u", str(lakefreq.SCRIPT)]
    assert "--without-design" in argv and str(tmp_path / "f.png") in argv
    assert "Interpreter not found" in lakefreq.interpreter_problem(str(tmp_path / "none"))


def test_export_names_say_which_figure_they_are(tmp_path):
    with_design = lakefreq.export_paths(tmp_path, "KRD", with_design=True)
    assert with_design.png.name == "KRD_validation.png"
    assert with_design.csv.name == "KRD_ams.csv"
    assert lakefreq.export_paths(tmp_path, "KRD", with_design=False).png.name == \
        "KRD_record.png"
    assert lakefreq.export_paths(tmp_path, "a/b", with_design=True).problems


def test_the_csv_export_is_in_water_year_order_and_reads_back(project, tmp_path):
    config = project.config.config_path
    job = lakefreq.build_job(config, settings_for(project)).job
    view = lakefreq.read_view(job)
    path = lakefreq.write_ams_csv(view, job, tmp_path / "out" / "ams.csv")

    text = path.read_text()
    assert "# sites: 999999A" in text
    back = lakefreq.RECORD.read_ams_csv(path)
    assert list(back["water_year"]) == sorted(back["water_year"])
    assert len(back) == len(view.positions)


def test_the_adopted_start_month_is_marked_in_the_scores(project):
    config = project.config.config_path
    view = lakefreq.read_view(lakefreq.build_job(config, settings_for(project)).job)
    scores = lakefreq.RECORD.water_year_scores(view.level, start_months=(9, 10))
    rows = lakefreq.score_rows(scores, 10)
    assert [row["adopted"] for row in rows] == [False, True]
    chart = lakefreq.seasonal_chart(lakefreq.RECORD.monthly_levels(view.level),
                                    lakefreq.RECORD.maxima_by_month(view.ams))
    assert len(chart["series"][-1]["data"]) == 12
