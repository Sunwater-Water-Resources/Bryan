"""Render the results page against a project whose rows have all been analysed.

The core modules are tested without a browser stack; this exists because the
failure this feature is most exposed to - the page raising on load, or the
chart quietly rendering nothing - is invisible to those tests. It builds its
own NiceGUI simulation rather than enabling the plugin globally, so the rest of
the suite still runs in an environment without nicegui.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from conftest import MONTE_CARLO_COLUMNS, write_workbook          # noqa: E402
from test_results import level_curve, write_quantiles             # noqa: E402

DURATIONS = (12, 24, 48, 120)


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def project(tmp_path):
    """Four durations of one case, analysed for inflow, level and volumes."""
    results_folder = tmp_path / "sims_mc" / "results"
    rows = []
    for duration in DURATIONS:
        name = f"TFD_mc_{duration}h_GWL1p3"
        level = level_curve(duration)
        write_quantiles(results_folder / f"{name}_level.csv", "level", level)
        write_quantiles(results_folder / f"{name}_inflow.csv", "inflow",
                        {aep: value * 10 for aep, value in level.items()})
        write_quantiles(results_folder / f"{name}_inflowVol24h.csv", "Vol24h",
                        {aep: value * 100 for aep, value in level.items()})
        rows.append({
            "Include": "yes", "Method": "monte carlo", "Duration": duration,
            "Run models": "yes", "Analyse results": "yes", "GWL": 1.3,
            "Output file": f"sims_mc\\results\\{name}",     # a Windows path, as in the wild
            "Config file": "mc.json",
        })
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return config


def _open(config):
    from state import STATE
    STATE.open_project(config)
    return STATE


@pytest.mark.asyncio
async def test_the_page_renders_the_group_it_found(user, project):
    _open(project)
    await user.open("/results")
    await user.should_see("Durations")
    await user.should_see("Critical durations")
    await user.should_see("120h")


@pytest.mark.asyncio
async def test_the_pinned_end_warnings_reach_the_page(user, project):
    """The whole point of the tab: does the duration range need extending."""
    _open(project)
    await user.open("/results")
    await user.should_see("longest duration you ran")
    await user.should_see("shortest duration you ran")


@pytest.mark.asyncio
async def test_the_chart_carries_a_series_per_duration_and_an_envelope(user, project):
    _open(project)
    await user.open("/results")
    chart = next(iter(user.find(marker="duration-chart").elements))
    names = [series["name"] for series in chart.options["series"]]
    assert names == ["12h", "24h", "48h", "120h", "envelope"]
    # inflow leads TYPE_ORDER, so that is what opens - and flows are plotted
    # logarithmically, as UtilModule.plot_durations does it
    assert chart.options["yAxis"]["type"] == "log"
    assert chart.options["series"][-1]["markArea"]["data"]     # the duration bands
    json.dumps(chart.options, allow_nan=False)                # no NaN reaches the browser

    critical = next(iter(user.find(marker="critical-duration-chart").elements))
    durations = [point[1] for point in critical.options["series"][0]["data"]]
    assert durations[0] > durations[-1], (
        "the critical duration should shorten as the AEP rarens")


@pytest.mark.asyncio
async def test_a_project_with_no_analysed_rows_says_so(user, tmp_path):
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, [{
        "Include": "yes", "Method": "monte carlo", "Duration": 24,
        "Run models": "yes", "Analyse results": "no",
        "Output file": "results\\never_run", "Config file": "mc.json"}])
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    _open(config)
    await user.open("/results")
    await user.should_see("No analysed results found")


@pytest.mark.asyncio
async def test_the_files_expansion_survives_a_refresh(user, project):
    """Rebuilding an expansion closes it - the Run page console's lesson."""
    from nicegui import ui
    _open(project)
    await user.open("/results")
    expansion = next(iter(user.find(kind=ui.expansion).elements))
    expansion.value = True
    user.find(marker="curve-120h").click()
    assert expansion.value is True


@pytest.mark.asyncio
async def test_unticking_a_duration_takes_it_off_the_chart(user, project):
    _open(project)
    await user.open("/results")
    user.find(marker="curve-120h").click()

    chart = next(iter(user.find(marker="duration-chart").elements))
    names = [series["name"] for series in chart.options["series"]]
    assert names == ["12h", "24h", "48h", "envelope"]
    # and the warning follows the data rather than the sims list
    await user.should_see("longest duration you ran (48h)")


@pytest.mark.asyncio
async def test_the_export_dialog_previews_what_it_would_write(user, project):
    _open(project)
    await user.open("/results")
    user.find("Export critical durations").click()

    # the group's name without its duration, and the files it would produce
    await user.should_see("TFD_mc_GWL1p3")
    await user.should_see("TFD_mc_GWL1p3_inflow_critical.csv")
    await user.should_see("TFD_mc_GWL1p3_inflow_critical_durations.png")


@pytest.mark.asyncio
async def test_the_export_dialog_says_when_it_would_overwrite(user, project, tmp_path):
    existing = tmp_path / "sims_mc" / "results" / "TFD_mc_GWL1p3_inflow_critical.csv"
    existing.write_text("an earlier export\n")
    _open(project)
    await user.open("/results")
    user.find("Export critical durations").click()
    await user.should_see("exists, will be overwritten")


# -- the Groups tab ----------------------------------------------------------
#
# Two climate scenarios over the same four durations, so each group has an
# envelope and the difference between them is a known constant: +0.4 m on
# level, +8% on inflow. Anything else the overlay reports is a bug in the
# maths rather than an artefact of the fixture.

WARMINGS = {"GWL1p3": (0.0, 1.0), "GWL2p0": (0.4, 1.08)}


@pytest.fixture
def groups_project(tmp_path):
    results_folder = tmp_path / "sims_mc" / "results"
    rows = []
    for warming, (lift_m, lift_factor) in WARMINGS.items():
        for duration in DURATIONS:
            name = f"TFD_mc_{duration}h_{warming}"
            base = level_curve(duration)
            write_quantiles(results_folder / f"{name}_level.csv", "level",
                            {aep: value + lift_m for aep, value in base.items()})
            write_quantiles(results_folder / f"{name}_inflow.csv", "inflow",
                            {aep: value * 10 * lift_factor
                             for aep, value in base.items()})
            rows.append({
                "Include": "yes", "Method": "monte carlo", "Duration": duration,
                "Run models": "yes", "Analyse results": "yes",
                "Output file": f"sims_mc\\results\\{name}",
                "Config file": "mc.json",
            })
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return config


async def _open_groups_tab(user, config):
    """Switch tabs by value, not by clicking.

    A ``ui.tab`` has no python click handler - Quasar switches the panel in the
    browser - so a simulated click would leave the panel unbuilt.
    """
    _open(config)
    await user.open("/results")
    tabs = next(iter(user.find(marker="result-tabs").elements))
    tabs.value = "Groups"
    return tabs


def _chart(user, marker):
    return next(iter(user.find(marker=marker).elements))


@pytest.mark.asyncio
async def test_the_groups_tab_draws_one_envelope_per_group(user, groups_project):
    await _open_groups_tab(user, groups_project)
    chart = _chart(user, "overlay-chart")

    assert [series["name"] for series in chart.options["series"]] == \
        ["GWL1p3", "GWL2p0"]
    # no envelope of envelopes and no duration bands: groups are scenarios
    assert not any("markArea" in series for series in chart.options["series"])
    json.dumps(chart.options, allow_nan=False)


@pytest.mark.asyncio
async def test_the_group_key_is_trimmed_to_what_tells_them_apart(user,
                                                                 groups_project):
    """The full key is 'sims_mc\\results\\TFD_mc_GWL1p3'. The legend is not."""
    await _open_groups_tab(user, groups_project)
    await user.should_see("GWL2p0")


@pytest.mark.asyncio
async def test_flows_are_compared_as_percentages(user, groups_project):
    await _open_groups_tab(user, groups_project)
    delta = _chart(user, "delta-chart")

    assert [series["name"] for series in delta.options["series"]] == ["GWL2p0"]
    assert delta.options["yAxis"]["name"] == "Change from GWL1p3 (%)"
    values = [point[1] for point in delta.options["series"][0]["data"]]
    assert values == pytest.approx([8.0] * len(values))
    # a change is signed, so it is never drawn on a log axis
    assert delta.options["yAxis"]["type"] == "value"


@pytest.mark.asyncio
async def test_level_is_compared_in_metres(user, groups_project):
    """A percentage of a level on an arbitrary datum says nothing."""
    await _open_groups_tab(user, groups_project)
    toggle = _chart(user, "group-types")
    toggle.value = "level"

    delta = _chart(user, "delta-chart")
    assert delta.options["yAxis"]["name"] == "Change from GWL1p3 (m)"
    values = [point[1] for point in delta.options["series"][0]["data"]]
    assert values == pytest.approx([0.4] * len(values))


@pytest.mark.asyncio
async def test_a_group_keeps_its_colour_when_the_baseline_drops_out(
        user, groups_project):
    """The delta frame has no baseline column, so colouring off it would shift."""
    await _open_groups_tab(user, groups_project)
    overlay_colours = {series["name"]: series["itemStyle"]["color"]
                       for series in _chart(user, "overlay-chart").options["series"]}
    delta = _chart(user, "delta-chart").options["series"][0]

    assert delta["itemStyle"]["color"] == overlay_colours["GWL2p0"]


@pytest.mark.asyncio
async def test_unticking_a_group_takes_it_off_the_overlay(user, groups_project):
    await _open_groups_tab(user, groups_project)
    user.find(marker="group-GWL2p0").click()

    chart = _chart(user, "overlay-chart")
    assert [series["name"] for series in chart.options["series"]] == ["GWL1p3"]
    # one group left, so there is nothing to measure a change against: the
    # whole section is hidden rather than drawn blank, and user.find only
    # gathers what is visible
    with pytest.raises(AssertionError):
        _chart(user, "delta-chart")


@pytest.mark.asyncio
async def test_the_critical_duration_is_overlaid_too(user, groups_project):
    """Whether a warmer climate moves the critical duration at all."""
    await _open_groups_tab(user, groups_project)
    critical = _chart(user, "critical-overlay-chart")

    assert [series["name"] for series in critical.options["series"]] == \
        ["GWL1p3", "GWL2p0"]
    hours = [point[1] for point in critical.options["series"][0]["data"]]
    assert hours[0] > hours[-1], "level goes long, then short on the rare tail"


@pytest.mark.asyncio
async def test_a_pinned_envelope_is_reported_against_its_group(user,
                                                               groups_project):
    """An envelope pinned to the end of its range is a lower bound."""
    await _open_groups_tab(user, groups_project)
    await user.should_see("GWL1p3: The longest duration you ran")


@pytest.mark.asyncio
async def test_the_durations_tab_is_untouched_by_all_this(user, groups_project):
    _open(groups_project)
    await user.open("/results")
    chart = _chart(user, "duration-chart")
    names = [series["name"] for series in chart.options["series"]]
    assert names == ["12h", "24h", "48h", "120h", "envelope"]


@pytest.mark.asyncio
async def test_one_group_has_nothing_to_overlay(user, project):
    """The single-group sims list is the common case - say so, do not sulk."""
    await _open_groups_tab(user, project)
    await user.should_see("nothing to overlay it against")
