"""The Sub-bursts tab: core/subburst.py, and the tab on the Results page.

The arithmetic runs in util/SubBurstCheck.py under Bryan's interpreter and has its own
test there (tests/test_subburst_check.py). What is pinned here is everything around it:
which groups can be tested, the job and its fingerprint, and reading a results file
back into the verdict, the worst-margin table, the chart and the level table - the
results file is written by hand, in the shape the script writes.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from conftest import MONTE_CARLO_COLUMNS, write_workbook
from test_results import level_curve, write_quantiles

AEPS = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
SCHEME = {"lower_aep": 2, "upper_aep": 10000, "number_of_main_divisions": 10,
          "number_of_sub_divisions": 100, "number_of_temporal_patterns": 10}


def mcdf(recorded=True):
    frame = pd.DataFrame({"m": [0, 1], "rain_z": [0.1, 2.0], "tp": [0, 1],
                          "storm_method": "GSDM", "level": [200.0, 201.0]})
    if recorded:
        frame["subburst_6h"] = [10.0, 30.0]
        frame["ifd_6h"] = [11.0, 28.0]
    return frame


def build(tmp_path, groups=("GWL1p3", "GWL1p3_noEBF"), unrecorded=()):
    """Two groups of two durations, analysed, each with its Monte Carlo database."""
    folder = tmp_path / "sims_mc" / "results"
    rows = []
    for group in groups:
        for duration in (24, 48):
            name = f"CLD_mc_{duration}h_{group}"
            write_quantiles(folder / f"{name}_level.csv", "level", level_curve(duration))
            mcdf(recorded=group not in unrecorded).to_csv(folder / f"{name}__mcdf.csv")
            rows.append({"Include": "yes", "Method": "monte carlo", "Duration": duration,
                         "Run models": "yes", "Analyse results": "yes",
                         "Output file": f"sims_mc\\results\\{name}", "Config file": "mc.json"})
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    (tmp_path / "mc.json").write_text(json.dumps({"scheme_config": SCHEME}), encoding="utf-8")
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return config


def open_project(config):
    from state import STATE
    STATE.open_project(config)
    return STATE.project


def margins(peak):
    """A margin series rising to ``peak`` at 1 in 1000 and falling away either side."""
    return {str(a): round(peak - 0.02 * abs(np.log10(a) - 3), 3) for a in AEPS}


def results_for(plan, *, tested_peak=1.2, compare_peak=0.97, calibrated=False):
    grid = list(np.linspace(199.0, 203.0, 41))
    out = {
        "fingerprint": plan.fingerprint, "standard_aeps": [float(a) for a in AEPS],
        "groups": {
            "tested": {"label": "tested", "durations": {
                "24h": {"margins": {"6h": margins(tested_peak - 0.2), "12h": margins(tested_peak - 0.05)},
                        "database": "a.csv", "has_level": True},
                "48h": {"margins": {"6h": margins(tested_peak - 0.08), "12h": margins(tested_peak)},
                        "database": "b.csv", "has_level": True}}},
        },
        "problems": [], "calibration": None,
    }
    if plan.job.get("compare"):
        out["groups"]["compare"] = {"label": "compared", "durations": {
            "24h": {"margins": {"6h": margins(compare_peak - 0.03), "12h": margins(compare_peak - 0.01)},
                    "database": "c.csv", "has_level": True},
            "48h": {"margins": {"6h": margins(compare_peak - 0.02), "12h": margins(compare_peak)},
                    "database": "d.csv", "has_level": True}}}
    if calibrated:
        def curve(shift):
            # P(level > x) falling from 1/2 at 199 m to 1/100000 at 203 m, shifted
            p = [float(10 ** (-0.3 - 4.7 * (x - 199.0 - shift) / 4.0)) for x in grid]
            p = [min(0.5, max(1e-6, v)) for v in p]
            levels = {str(a): round(float(np.interp(1 / a, p[::-1], grid[::-1])), 3) for a in AEPS}
            return {"levels": levels, "p": p}
        out["calibration"] = {
            "durations": {"24h": {"converged": True, "iterations": 4, "worst_before": 1.12,
                                  "worst_after": 0.99, "weights": {"GSDM": [0.02] * 3 + [0.137] * 7},
                                  "at_floor": {"GSDM": 3}, "weights_file": "w24.csv"},
                          "48h": {"converged": False, "iterations": 12, "worst_before": 1.2,
                                  "worst_after": 1.04, "weights": {"GSDM": [0.1] * 10},
                                  "at_floor": {"GSDM": 0}, "weights_file": "w48.csv"}},
            "curves": {"tested": curve(0.5), "calibrated": curve(0.0), "compare": curve(0.02)},
            "grid": [round(float(x), 4) for x in grid],
        }
    return out


def save(config, plan, results):
    from core import subburst
    path = subburst.results_path(config, plan.fingerprint)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results), encoding="utf-8")


# -- core ---------------------------------------------------------------------------

def test_candidates_are_the_groups_with_databases_and_say_which_recorded(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path, unrecorded=("GWL1p3_noEBF",)))
    found = subburst.candidates(project)
    recorded = {group: c.recorded for group, c in found.items()}
    assert len(recorded) == 2
    assert sorted(recorded.values()) == [False, True]
    for candidate in found.values():
        assert set(candidate.databases) == {"24h", "48h"}
        assert candidate.config_file.name == "mc.json"


def test_the_job_carries_the_scheme_and_both_groups(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested, compare = subburst.candidates(project).values()
    plan = subburst.build_job(project, tested, compare, calibrate=True)
    assert plan.problems == []
    assert plan.job["scheme"]["number_of_main_divisions"] == 10
    assert set(plan.job["tested"]["databases"]) == {"24h", "48h"}
    assert plan.job["compare"]["label"] == subburst.short_label(compare.group)
    assert plan.job["calibrate"] is True
    assert plan.fingerprint in plan.job["weights_folder"]


def test_the_fingerprint_changes_when_a_database_does(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested = next(iter(subburst.candidates(project).values()))
    before = subburst.build_job(project, tested).fingerprint
    assert subburst.build_job(project, tested).fingerprint == before
    assert subburst.build_job(project, tested, calibrate=True).fingerprint != before
    path = tested.databases["24h"]
    frame = pd.read_csv(path, index_col=0)
    pd.concat([frame, frame]).to_csv(path)
    assert subburst.build_job(project, tested).fingerprint != before


def test_a_group_without_sub_bursts_or_a_scheme_is_refused_with_a_reason(tmp_path):
    from core import subburst
    config = build(tmp_path, unrecorded=("GWL1p3_noEBF",))
    (tmp_path / "mc.json").write_text(json.dumps({"scheme_config": {"lower_aep": 2}}))
    project = open_project(config)
    unrecorded = next(c for c in subburst.candidates(project).values() if not c.recorded)
    problems = subburst.build_job(project, unrecorded).problems
    assert any("before sub-burst tracking" in p for p in problems)
    assert any("upper_aep" in p for p in problems)


def test_the_verdict_names_the_worst_window_and_reads_one_sided():
    from core import subburst

    class Plan:
        fingerprint, job = "x", {"compare": None}

    breach = results_for(Plan, tested_peak=1.2)
    aeps = subburst.window(breach)
    severity, message, _ = subburst.verdict(breach, aeps)
    assert severity == "warn"
    assert "12h windows of the 48h storms reach 1.20" in message and "1 in 1,000" in message

    severity, message, _ = subburst.verdict(results_for(Plan, tested_peak=1.03), aeps)
    assert (severity, message.split(":")[0]) == ("info", "Close to neutral")
    severity, message, _ = subburst.verdict(results_for(Plan, tested_peak=0.9), aeps)
    assert message.startswith("Neutral")


def test_the_window_drops_the_pinned_ends_unless_asked():
    from core import subburst
    results = {"standard_aeps": [float(a) for a in AEPS]}
    assert subburst.window(results) == [float(a) for a in AEPS[1:-1]]
    assert subburst.window(results, 100, 1000) == [100.0, 200.0, 500.0, 1000.0]


def test_the_matrix_shows_tested_then_compared_and_flags_breaches(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested, compare = subburst.candidates(project).values()
    plan = subburst.build_job(project, tested, compare)
    results = results_for(plan)
    subs, rows = subburst.matrix_rows(results, subburst.window(results))
    assert subs == ["6h", "12h"]
    assert [r["parent"] for r in rows] == ["24h", "48h"]
    assert rows[1]["12h"] == "1.20 → 0.97"
    assert rows[1]["12h_breach"] and not rows[0]["6h_breach"]


def test_the_chart_dashes_the_compared_group_and_marks_the_ifd(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested, compare = subburst.candidates(project).values()
    results = results_for(subburst.build_job(project, tested, compare))
    options = subburst.margin_chart(results, "48h")
    names = [s["name"] for s in options["series"]]
    assert names == ["6h", "6h (compared)", "12h", "12h (compared)"]
    assert options["series"][1]["lineStyle"]["type"] == "dashed"
    assert options["series"][0]["markLine"]["data"] == [{"yAxis": 1.0}]
    json.dumps(options, allow_nan=False)


def test_the_level_table_reads_the_aep_of_any_level_off_each_curve(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested, compare = subburst.candidates(project).values()
    results = results_for(subburst.build_job(project, tested, compare), calibrated=True)
    cases, rows = subburst.level_rows(results, 201.0)
    assert cases == ["compare", "tested", "calibrated"]
    assert rows[-1]["aep"] == "AEP of 201.00"
    # the tested curve is 0.5 m higher, so 201 m is more frequent on it
    tested_aep = subburst.curve_aep(results, "tested", 201.0)
    calibrated_aep = subburst.curve_aep(results, "calibrated", 201.0)
    assert tested_aep < calibrated_aep
    assert subburst.curve_aep(results, "tested", 250.0) is None
    assert rows[-1]["tested"].startswith("1 in ")


def test_the_calibration_rows_say_what_converged_and_what_hit_the_floor(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested = next(iter(subburst.candidates(project).values()))
    rows = subburst.calibration_rows(results_for(subburst.build_job(project, tested), calibrated=True))
    assert [(r["duration"], r["converged"], r["floored"]) for r in rows] == \
        [("24h", "yes", "GSDM: 3/10"), ("48h", "no", "none")]
    assert rows[0]["margin"] == "1.12 → 0.99"


def test_the_cached_result_is_found_and_the_calibrated_one_preferred(tmp_path):
    from core import subburst
    config = build(tmp_path)
    project = open_project(config)
    tested = next(iter(subburst.candidates(project).values()))
    assert subburst.shown(project, tested, None) is None
    plain = subburst.build_job(project, tested)
    save(config, plain, results_for(plain))
    assert subburst.shown(project, tested, None)["calibration"] is None
    calibrated = subburst.build_job(project, tested, calibrate=True)
    save(config, calibrated, results_for(calibrated, calibrated=True))
    assert subburst.shown(project, tested, None)["calibration"] is not None


# -- the tab ------------------------------------------------------------------------

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


async def open_tab(user, config, tab="Sub-bursts"):
    open_project(config)
    await user.open("/results")
    tabs = next(iter(user.find(marker="result-tabs").elements))
    tabs.value = tab
    return tabs


@pytest.mark.asyncio
async def test_the_tab_says_when_nothing_has_been_checked(user, tmp_path):
    await open_tab(user, build(tmp_path))
    await user.should_see("Check neutrality")
    await user.should_see(marker="subburst-none")


@pytest.mark.asyncio
async def test_the_tab_draws_a_saved_result(user, tmp_path):
    from core import subburst
    config = build(tmp_path)
    project = open_project(config)
    tested = subburst.candidates(project)[next(iter(subburst.candidates(project)))]
    plan = subburst.build_job(project, tested, calibrate=True)
    save(config, plan, results_for(plan, calibrated=True))

    await open_tab(user, config)
    await user.should_see("Embedded bursts occur more often than the IFD implies")
    await user.should_see(marker="subburst-matrix")
    chart = next(iter(user.find(marker="subburst-chart").elements))
    assert [s["name"] for s in chart.options["series"]] == ["6h", "12h"]
    await user.should_see(marker="subburst-weights")
    await user.should_see(marker="subburst-levels")


@pytest.mark.asyncio
async def test_a_project_whose_databases_predate_sub_bursts_says_how_to_get_them(user, tmp_path):
    await open_tab(user, build(tmp_path, groups=("GWL1p3",), unrecorded=("GWL1p3",)))
    await user.should_see("No group here has recorded sub-burst depths.")


@pytest.mark.asyncio
async def test_the_tab_opens_from_its_address(user, tmp_path):
    open_project(build(tmp_path))
    await user.open("/results?tab=Sub-bursts")
    await user.should_see("Check neutrality")
    await user.should_see(marker="subburst-none")


def test_groups_are_named_by_what_tells_them_apart_and_files_by_the_last_part(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    found = subburst.candidates(project)
    assert sorted(subburst.display_names(found).values()) == ["GWL1p3", "GWL1p3_noEBF"]
    tested = next(iter(found.values()))
    assert subburst.build_job(project, tested).job["tested"]["label"] == "CLD_mc_GWL1p3"
    assert subburst.short_label("sims_mc\results\CLD_mc_GWL1p3") == "CLD_mc_GWL1p3"
    assert subburst.short_label("sims_mc/results/x") == "x"


def test_the_margin_axis_always_takes_in_the_ifd(tmp_path):
    from core import subburst
    project = open_project(build(tmp_path))
    tested = next(iter(subburst.candidates(project).values()))
    options = subburst.margin_chart(results_for(subburst.build_job(project, tested)), "48h")
    assert "Math.min(v.min, 1)" in options["yAxis"][":min"]
    assert "Math.max(v.max, 1)" in options["yAxis"][":max"]
