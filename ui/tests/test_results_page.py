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
