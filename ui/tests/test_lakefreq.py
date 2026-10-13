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


def test_the_quick_fit_is_the_same_job_without_resampling(project):
    config = project.config.config_path
    job = lakefreq.build_job(config, settings_for(project)).job
    quick = lakefreq.quick_job(job)
    assert quick["fit"]["draws"] == 0 and job["fit"]["draws"] > 0     # not mutated
    assert {**quick, "fit": None} == {**job, "fit": None}
    assert lakefreq.RECORD.fingerprint(quick) != lakefreq.RECORD.fingerprint(job)


def test_the_page_shows_the_bands_when_resampled_and_the_quick_fit_otherwise(project):
    config = project.config.config_path
    job = lakefreq.build_job(config, settings_for(project)).job
    assert lakefreq.shown_results(config, job) == (None, False)

    def save(which):
        path = lakefreq.results_path(config, which)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"fingerprint": lakefreq.RECORD.fingerprint(which),
                                    "kind": "banded" if which is job else "quick"}))

    save(lakefreq.quick_job(job))
    results, banded = lakefreq.shown_results(config, job)
    assert results["kind"] == "quick" and not banded
    save(job)
    results, banded = lakefreq.shown_results(config, job)
    assert results["kind"] == "banded" and banded


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


def test_the_curve_options_reach_the_job_and_old_settings_get_the_defaults(project):
    config = project.config.config_path
    old = settings_for(project)
    del old["fit"]["upper_degree"], old["fit"]["upper_join"]
    lakefreq.save_settings(config, old)
    loaded = lakefreq.load_settings(config)
    assert (loaded["fit"]["upper_degree"], loaded["fit"]["upper_join"]) == (1, "free")

    loaded["fit"].update(plateau_tolerance=0.0, upper_degree=2, upper_join="fsl")
    job = lakefreq.build_job(config, loaded).job
    assert job["fit"]["upper_join"] == "fsl" and job["fit"]["plateau_tolerance"] == 0.0
    assert set(lakefreq.UPPER_JOIN_LABELS) == {"free", "fsl"}


# -- the frequent end of the axis: where it starts, and where the fitted curves do ------------

AXIS_START, CURVE_START = "frequent_aep_1_in_x", "curve_frequent_aep_1_in_x"


def z_of(one_in_x):
    return lakefreq.RECORD.normal_variate_of_probability(1.0 / one_in_x)


def charted(project, results=True, whole_record=False, grid_from=None, **axes):
    """The chart for the synthetic record with these axis settings. ``grid_from``
    starts the results' grid there instead of at -0.3, as a start past 1EY does."""
    settings = settings_for(project)
    settings["axes"].update(axes)
    view = lakefreq.read_view(lakefreq.build_job(project.config.config_path, settings).job)
    fitted = results_for(view.positions) if results else None
    if grid_from is not None:
        fitted["grid"]["z"] = [grid_from + 0.1 * i for i in range(30)]
    return lakefreq.chart_options(view.positions, fitted, settings,
                                  whole_record=whole_record)


def first_z(options, name):
    """The most frequent point of the series whose name starts with ``name``."""
    return min(point["value"][0] if isinstance(point, dict) else point[0]
               for series in options["series"]
               if series["name"].startswith(name) for point in series["data"])


def test_both_starts_are_blank_by_default_and_a_blank_chart_is_the_one_from_1ey(project):
    axes = lakefreq.default_settings()["axes"]
    assert axes[AXIS_START] is None and axes[CURVE_START] is None

    blank = charted(project)
    ey1 = lakefreq.RECORD.normal_variate_of_probability(lakefreq.RECORD.EY1_AEP)
    assert blank["xAxis"]["min"] == ey1
    # The whole grid is drawn, from its first point: results_for's starts at -0.3.
    assert first_z(blank, "Fit to all") == first_z(blank, "Design flood envelope") == -0.3

    # Settings written before the keys existed draw the same chart.
    settings = settings_for(project)
    del settings["axes"][AXIS_START], settings["axes"][CURVE_START]
    view = lakefreq.read_view(lakefreq.build_job(project.config.config_path, settings).job)
    for whole_record in (False, True):
        old = lakefreq.chart_options(view.positions, results_for(view.positions), settings,
                                     whole_record=whole_record)
        assert old == charted(project, whole_record=whole_record)


def test_the_axis_start_moves_the_frame_and_the_design_floods_follow_it(project):
    blank, moved = charted(project), charted(project, **{AXIS_START: 5})
    assert moved["xAxis"]["min"] == pytest.approx(z_of(5))
    assert moved["xAxis"]["max"] == blank["xAxis"]["max"]
    # Nothing is cut: the curves still start at 1EY and the envelope with them,
    # and the frame clips both.
    for name in ("Fit to all", "90% band, all", "Design flood envelope"):
        assert first_z(moved, name) == -0.3, name
    # The level axis is set by what is inside the frame.
    assert moved["yAxis"]["min"] >= blank["yAxis"]["min"]
    # The ticks are the ones inside the new frame.
    ticks = list(json.loads(moved["xAxis"]["axisLabel"][":formatter"]
                            .split("const m = ")[1].split(";")[0]).values())
    assert ticks[:3] == ["1 in 5", "1 in 10", "1 in 20"]

    # More frequent than 1EY the axis reads in EY.
    wide = charted(project, **{AXIS_START: 1.1})
    assert wide["xAxis"]["min"] == pytest.approx(z_of(1.1))
    assert "2EY" in wide["xAxis"]["axisLabel"][":formatter"]
    assert "2EY" not in blank["xAxis"]["axisLabel"][":formatter"]


def test_a_blank_curve_start_is_1ey_wherever_the_axis_starts(project):
    """The grid runs out to an axis started past 1EY; the curves do not follow it
    there unless they are started there too. The design floods do."""
    ey1 = lakefreq.RECORD.normal_variate_of_probability(lakefreq.RECORD.EY1_AEP)
    out = z_of(1.1)
    wide = charted(project, grid_from=out, **{AXIS_START: 1.1})
    for name in ("Fit to all", "90% band, all", "band"):
        assert ey1 - 1e-9 <= first_z(wide, name) < ey1 + 0.1, name
    assert first_z(wide, "Design flood envelope") == pytest.approx(out)

    both = charted(project, grid_from=out, **{AXIS_START: 1.1, CURVE_START: 1.2})
    assert z_of(1.2) - 1e-9 <= first_z(both, "Fit to all") < z_of(1.2) + 0.1
    assert both["xAxis"] == wide["xAxis"]

    # Curves started past 1EY on a 1EY axis: there for the whole record view, and
    # the design floods are not carried out there with them.
    curves = charted(project, grid_from=out, whole_record=True, **{CURVE_START: 1.1})
    assert first_z(curves, "Fit to all") == pytest.approx(out)
    assert ey1 - 1e-9 <= first_z(curves, "Design flood envelope") < ey1 + 0.1


def test_the_curve_start_cuts_only_the_fitted_curves_and_their_bands(project):
    blank, cut = charted(project), charted(project, **{CURVE_START: 4})
    assert cut["xAxis"] == blank["xAxis"]
    for name in ("Fit to all", "90% band, all", "band", "Fit to storm-driven",
                 "90% band, storm-driven"):
        assert first_z(cut, name) >= z_of(4) - 1e-9, name
        assert first_z(cut, name) < z_of(4) + 0.11, name          # the next grid point
    assert first_z(cut, "Design flood envelope") == -0.3

    def untouched(options):
        return [series for series in options["series"]
                if not series["name"].lower().startswith(("fit to", "90% band", "band"))]
    assert untouched(cut) == untouched(blank)


def test_the_two_starts_are_independent_and_the_axis_clips_the_curves(project):
    both = charted(project, **{AXIS_START: 3, CURVE_START: 4})
    assert both["xAxis"]["min"] == pytest.approx(z_of(3))
    assert first_z(both, "Fit to all") >= z_of(4) - 1e-9

    # A curve start more frequent than the axis start is clipped by the frame,
    # which it does not move.
    ahead = charted(project, **{AXIS_START: 3, CURVE_START: 2})
    assert ahead["xAxis"] == both["xAxis"]
    assert z_of(2) - 1e-9 <= first_z(ahead, "Fit to all") < z_of(3)


def test_the_storm_driven_curve_still_starts_at_its_own_most_frequent_point(project):
    settings = settings_for(project)
    view = lakefreq.read_view(lakefreq.build_job(project.config.config_path, settings).job)
    floor = view.positions.loc[~view.positions["carried_over"], "storm_z"].min()
    # A grid reaching well past that point, as an axis started out there gives.
    results = results_for(view.positions)
    results["grid"]["z"] = [floor - 1.0 + 0.1 * i for i in range(30)]

    def storm_start(**axes):
        settings["axes"].update({AXIS_START: None, CURVE_START: None, **axes})
        return first_z(lakefreq.chart_options(view.positions, results, settings),
                       "Fit to storm-driven")

    assert z_of(1.01) < floor < storm_start()                     # blank: from 1EY
    assert storm_start(**{AXIS_START: 1.01}) == storm_start()     # the axis alone: still 1EY
    for axes in ({CURVE_START: 1.02}, {AXIS_START: 1.01, CURVE_START: 1.02}):
        assert floor <= storm_start(**axes) < floor + 0.1
    assert storm_start(**{CURVE_START: 4}) >= z_of(4) - 1e-9 > floor


def test_the_whole_record_still_shows_every_maximum_whatever_the_axis_start(project):
    blank = charted(project, whole_record=True)
    lowest = min(first_z(blank, "Storm-driven maxima"), first_z(blank, "Carried over"))
    assert blank["xAxis"]["min"] == pytest.approx(lowest - 0.1)
    moved = charted(project, whole_record=True, **{AXIS_START: 5})
    assert moved["xAxis"] == blank["xAxis"]
    # ...with the curves from 1EY, or from where they are started.
    assert first_z(moved, "Fit to all") == -0.3
    later = charted(project, whole_record=True, **{AXIS_START: 5, CURVE_START: 4})
    assert z_of(4) - 1e-9 <= first_z(later, "Fit to all") < z_of(4) + 0.1


def test_a_start_that_cannot_be_one_is_said_and_stops_only_the_curves(project):
    config = project.config.config_path

    def plan(**axes):
        settings = settings_for(project)
        settings["axes"].update(axes)
        return lakefreq.build_job(config, settings)

    assert plan().can_fit and plan(**{AXIS_START: 5, CURVE_START: 10}).can_fit
    for value in (1, 0.5, -3):
        bad = plan(**{AXIS_START: value})
        assert bad.can_read and not bad.can_fit
        assert "axis start" in bad.fit_problems[0] and "above 1" in bad.fit_problems[0]
    behind = plan(**{CURVE_START: 2000})
    assert "more frequent than where the axis reaches (1 in 2000)" in behind.fit_problems[0]
    assert plan(**{AXIS_START: 50, "rare_aep_1_in_x": 20}).fit_problems

    # The points are still drawn, in the frame a blank start gives.
    assert charted(project, results=False, **{AXIS_START: 0.5}) == \
        charted(project, results=False)
    assert charted(project, **{AXIS_START: 1e6, CURVE_START: 1e6}) == charted(project)
    assert lakefreq.RECORD.frequent_limits({AXIS_START: "x", CURVE_START: None}) == \
        (lakefreq.RECORD.EY1_AEP, lakefreq.RECORD.EY1_AEP)


def test_blank_starts_find_the_results_fitted_before_they_existed(project):
    config = project.config.config_path
    job = lakefreq.build_job(config, settings_for(project)).job
    before = json.loads(json.dumps(job))
    del before["axes"][AXIS_START], before["axes"][CURVE_START]
    fingerprint = lakefreq.RECORD.fingerprint
    assert fingerprint(job) == fingerprint(before)

    # Either start is part of the job, so curves fitted for another are not shown.
    path = lakefreq.results_path(config, job)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"fingerprint": fingerprint(job)}))
    assert lakefreq.shown_results(config, before)[0] is not None
    for key in (AXIS_START, CURVE_START):
        changed = {**job, "axes": {**job["axes"], key: 5}}
        assert fingerprint(changed) != fingerprint(job)
        assert lakefreq.shown_results(config, changed) == (None, False)


def test_the_starts_are_saved_and_an_old_settings_file_loads_blank(project):
    config = project.config.config_path
    path = lakefreq.settings_path(config)
    path.write_text(json.dumps({"record": {"files": ["gauge/HW.csv"]}, "fsl": FSL,
                                "axes": {"rare_aep_1_in_x": 500, "level_min": 210.0}}))
    loaded = lakefreq.load_settings(config)
    assert loaded["axes"] == {"rare_aep_1_in_x": 500, "level_min": 210.0, "level_max": None,
                              AXIS_START: None, CURVE_START: None}
    assert lakefreq.build_job(config, loaded).can_fit

    loaded["axes"].update({AXIS_START: 2, CURVE_START: 5})
    lakefreq.save_settings(config, loaded)
    assert json.loads(path.read_text())["axes"][CURVE_START] == 5
    job = lakefreq.build_job(config, lakefreq.load_settings(config)).job
    assert (job["axes"][AXIS_START], job["axes"][CURVE_START]) == (2, 5)


# -- the fitted curves as a table ---------------------------------------------------------

SHOULDERED = {"coef": [2.130803834, -2.958037438, -3.905986403, -1.373811495],
              "z1": 0.0, "z2": 0.5958401014, "upper_coef": [215.17714998, 0.4770768766],
              "upper_join": "free", "fsl": 215.0, "step": 0.17714998}
LOGISTIC = {"floor": 209.25, "k": 1.875, "z0": -0.8416212336, "ceiling": 216.125}


def fitted_results(form="shouldered", params=SHOULDERED, banded=True, storm=True):
    """Results as the script writes them for a fit that worked."""
    block = {"form": form, "rmse": 0.144836638, "params": params, "curve": [214.0],
             "band_lo": [213.0] if banded else None, "band_hi": [215.0] if banded else None,
             "draws_used": 338 if banded else 0, "draws": 400 if banded else 0,
             "z_min": -2.0, "z_max": 2.0, "error": None, "warning": None}
    fits = {"all": block}
    if storm:
        fits["storm"] = dict(block, rmse=0.318445864, draws_used=378 if banded else 0)
    return {"grid": {"z": [0.0], "aep": [0.5]}, "fits": fits, "design": None,
            "ams": [{"carried_over": index < 8} for index in range(31)]}


def cells(table, column="all"):
    return {row["quantity"]: row[column] for row in table["rows"]}


def test_the_shouldered_curve_is_tabled_one_column_a_fit_and_one_row_a_quantity():
    table = lakefreq.fit_table(fitted_results())
    assert table["columns"] == [("quantity", "Quantity"), ("all", "All maxima"),
                                ("storm", "Storm-driven maxima")]
    assert [row["quantity"] for row in table["rows"]] == [
        "Form", "RMSE (m)", "Maxima fitted", "Resamples fitted, of those drawn",
        "Full supply level (m AHD)", "Plateau starts, z1", "Plateau ends, z2",
        "Below FSL: coefficient of (z - z1)", "Below FSL: coefficient of (z - z1)^2",
        "Below FSL: coefficient of (z - z1)^3", "Below FSL: coefficient of (z - z1)^4",
        "Curve above FSL starts", "Above the plateau: level at z2 (m AHD)",
        "Above the plateau: coefficient of (z - z2)",
        "Step up at the end of the plateau (m)"]
    whole, storm = cells(table), cells(table, "storm")
    assert (whole["Form"], whole["RMSE (m)"], storm["RMSE (m)"]) == \
        ("Shouldered plateau", "0.1448", "0.3184")
    assert (whole["Maxima fitted"], storm["Maxima fitted"]) == ("31", "23")
    assert whole["Resamples fitted, of those drawn"] == "338 of 400"
    assert storm["Resamples fitted, of those drawn"] == "378 of 400"
    # z with the frequency it stands for, as the axis labels it.
    assert whole["Plateau starts, z1"] == "0 (1 in 2)"
    assert whole["Plateau ends, z2"].startswith("0.5958401 (1 in 4")
    # Enough figures to redraw the curve.
    assert whole["Below FSL: coefficient of (z - z1)^2"] == "-2.9580374"
    assert whole["Above the plateau: level at z2 (m AHD)"] == "215.17715"
    assert whole["Curve above FSL starts"] == "Free"
    assert "the level at z2 plus" in table["note"]


def test_an_upper_limb_starting_at_full_supply_has_no_level_of_its_own():
    params = dict(SHOULDERED, upper_join="fsl", upper_coef=[1.25, -0.125], step=0.0,
                  coef=[2.5], z1=-0.25, z2=-0.25)
    table = lakefreq.fit_table(fitted_results(params=params, banded=False, storm=False))
    assert table["columns"] == [("quantity", "Quantity"), ("all", "All maxima")]
    whole = cells(table)
    assert "Above the plateau: level at z2 (m AHD)" not in whole
    assert "Resamples fitted, of those drawn" not in whole          # no bands yet
    assert whole["Above the plateau: coefficient of (z - z2)"] == "1.25"
    assert whole["Above the plateau: coefficient of (z - z2)^2"] == "-0.125"
    assert "Below FSL: coefficient of (z - z1)^2" not in whole
    # More frequent than 50% AEP the frequency reads in EY.
    assert whole["Plateau starts, z1"].startswith("-0.25 (0.") and "EY)" in whole["Plateau starts, z1"]
    assert whole["Curve above FSL starts"] == "At full supply"
    assert "it is FSL plus the sum of the coefficients above the plateau" in table["note"]


def test_the_logistic_curve_is_tabled_with_its_four_parameters():
    table = lakefreq.fit_table(fitted_results("logistic", LOGISTIC))
    assert [row["quantity"] for row in table["rows"]][4:] == [
        "Floor (m AHD)", "Ceiling (m AHD)", "Steepness, k (per unit of z)", "Midpoint, z0"]
    whole = cells(table)
    assert (whole["Form"], whole["Floor (m AHD)"], whole["Ceiling (m AHD)"]) == \
        ("Logistic", "209.25", "216.125")
    assert whole["Steepness, k (per unit of z)"] == "1.875"
    assert whole["Midpoint, z0"] == "-0.84162123 (1.61EY)"
    assert "exp(-k (z - z0))" in table["note"]


def test_a_failed_fit_has_no_column_and_no_fit_has_no_table():
    failed = {"form": "shouldered", "rmse": None, "params": None, "curve": None,
              "error": "no annual maximum within 25 mm of full supply"}
    assert lakefreq.fit_table(None) is None
    assert lakefreq.fit_table({"fits": {}, "design": None}) is None      # form "none"
    assert lakefreq.fit_table({"fits": {"all": failed, "storm": failed}}) is None

    # One of the two failing leaves the other's column alone.
    results = fitted_results()
    results["fits"]["storm"] = failed
    assert lakefreq.fit_table(results)["columns"] == [("quantity", "Quantity"),
                                                      ("all", "All maxima")]

    # Results without the maxima (not what the script writes) still table the fit.
    del results["ams"]
    assert "Maxima fitted" not in cells(lakefreq.fit_table(results))


def test_the_table_is_copied_as_tab_separated_text_under_a_header_row():
    from core import wordtable
    table = lakefreq.fit_table(fitted_results())
    copied = lakefreq.fit_clipboard_table(table)
    lines = wordtable.to_text(copied).splitlines()

    assert lines[0].split("\t") == ["Quantity", "All maxima", "Storm-driven maxima"]
    assert lines[1].split("\t") == ["Form", "Shouldered plateau", "Shouldered plateau"]
    assert lines[2].split("\t") == ["RMSE (m)", "0.1448", "0.3184"]
    assert "Plateau ends, z2\t0.5958401 (1 in 4" in "\n".join(lines)
    # A line a row, every row the table's width, and the equation last.
    assert len(lines) == 1 + len(table["rows"]) + 1
    assert all(len(line.split("\t")) == 3 for line in lines[:-1])
    assert lines[-1] == table["note"]
    assert "<table" in wordtable.to_html(copied)           # and as a table for Word

    assert lakefreq.fit_clipboard_table(lakefreq.fit_table(None)) is None
    one = lakefreq.fit_clipboard_table(lakefreq.fit_table(fitted_results(storm=False)))
    assert wordtable.to_text(one).splitlines()[0] == "Quantity\tAll maxima"



# -- with a study open: the record, water year and curve settings are the study's ------------

@pytest.fixture
def with_study(project):
    """The project's own settings saved as the page left them, and an empty study."""
    from core import study as studies
    config = project.config.config_path
    settings = settings_for(project, water_year_start=7)
    settings["fit"]["plateau_tolerance"] = 0.35
    lakefreq.save_settings(config, settings)
    study = studies.new_study(config.parent.parent / "bryan_study.json", "Dam")
    return config, study


def test_without_a_study_the_page_keeps_its_own_settings(project):
    config = project.config.config_path
    lakefreq.save_settings(config, settings_for(project))
    assert lakefreq.load_settings(config, None)["record"]["files"] == ["gauge/HW.csv"]


def test_a_sims_list_that_names_a_record_keeps_it_as_its_own(with_study):
    """It may be a homogenised series (Callide's E011 reads _FFA/CLD_FSL_homogenised.csv),
    and swapping it for the study's raw gauges would change the analysis."""
    config, study = with_study
    settings = lakefreq.load_settings(config, study)
    assert settings["record_source"] == lakefreq.OWN_RECORD
    assert Path(settings["record"]["files"][0]) == (config.parent / "gauge" / "HW.csv").resolve()
    assert settings["fit"]["plateau_tolerance"] == 0.35 and settings["fsl"] == FSL
    # July was this list's; the study's default is October, so July is an override
    assert settings["water_year"] == 7 and settings["water_year_start"] == 7


def test_the_record_becomes_the_study_s_gauges_only_when_asked(with_study):
    from core import dam as dams
    config, study = with_study
    settings = lakefreq.load_settings(config, study)
    assert dams.gauges(dams.settings(study)) == []                 # nothing done unasked
    done = lakefreq.make_study_gauges(settings, study)
    assert "HW.csv" in done
    assert dams.settings(study)["gauges"] == ["project/gauge/HW.csv"]   # relative to the study
    assert settings["record_source"] == lakefreq.DAM_RECORD


def test_a_sims_list_without_a_record_reads_the_study_s_gauges(project):
    from core import dam as dams, study as studies
    config = project.config.config_path
    study = studies.new_study(config.parent.parent / "bryan_study.json", "Dam")
    dam = dams.settings(study)
    dam["gauges"] = ["project/gauge/HW.csv"]
    dams.store(study, dam)
    settings = lakefreq.load_settings(config, study)
    assert settings["record_source"] == lakefreq.DAM_RECORD
    assert Path(settings["record"]["files"][0]) == (config.parent / "gauge" / "HW.csv").resolve()


def test_saving_keeps_the_study_s_share_in_the_study_and_the_rest_here(with_study):
    from core import study as studies
    config, study = with_study
    settings = lakefreq.load_settings(config, study)
    settings["fit"]["degree"] = 3
    settings["water_year"], settings["water_year_start"] = None, 10     # the study's
    settings["export"]["name"] = "KRO"
    lakefreq.save_settings(config, settings, study)

    part = studies.load_study(study.path).extra[lakefreq.STUDY_KEY]
    assert part["fit"]["degree"] == 3 and part["water_year"] is None
    assert part["record_source"] == lakefreq.OWN_RECORD
    assert part["files"] == ["project/gauge/HW.csv"]                 # relative to the study
    assert "export" not in part and "design" not in part
    own = json.loads(lakefreq.settings_path(config).read_text(encoding="utf-8"))
    assert own["export"]["name"] == "KRO"
    assert own["fit"]["degree"] == 3                   # still whole, for older launchers
    assert own["record"]["files"] == ["gauge/HW.csv"]  # relative to the sims list


def test_the_study_s_settings_make_a_job_the_record_can_be_read_with(with_study, project):
    config, study = with_study
    plan = lakefreq.build_job(config, lakefreq.load_settings(config, study))
    assert plan.can_read and plan.job["water_year_start"] == 7
    assert lakefreq.read_view(plan.job).ams is not None
