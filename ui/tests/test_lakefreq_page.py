"""Render the Lake levels page against a project with a level record beside it.

The core module is tested without a browser stack; this exists for what those
tests cannot see - the page raising on load, the record never reaching the
chart, or the fit finishing without the curves being drawn.
"""

from __future__ import annotations

import asyncio
import json

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from test_lakefreq import DURATIONS, project as lake_project      # noqa: E402,F401
from lake_fixtures import FSL                                     # noqa: E402
from test_real_bryan import BRYAN_PYTHON, needs_bryan             # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_settings(monkeypatch, tmp_path):
    """Opening a project records it in the user's settings; keep that in tmp."""
    import settings as settings_module
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")


@pytest.fixture
def opened(lake_project, private_settings):
    from core import lakefreq
    from state import STATE

    config = lake_project.config.config_path
    lakefreq.save_settings(config, {
        "record": {"files": ["gauge/HW.csv"]}, "fsl": FSL,
        "reference_levels": [{"label": "crest", "level": 218.0}],
        "fit": {"draws": 20},
        "design": {"include": True, "group": "", "durations": None}})
    STATE.open_project(config)
    return config


async def series_names(user, wanted, seconds=20.0):
    """The chart's series once ``wanted`` is among them - the record is read, and
    the fit run, off the event loop."""
    names = []
    for _ in range(int(seconds / 0.1)):
        chart = next(iter(user.find(marker="lake-chart").elements))
        names = [series["name"] for series in chart.options.get("series", [])]
        if any(name.startswith(wanted) for name in names):
            return names
        await asyncio.sleep(0.1)
    raise AssertionError(f"{wanted!r} never reached the chart; it has {names}")


@pytest.mark.asyncio
async def test_the_page_reads_the_record_and_plots_the_maxima(user, opened):
    await user.open("/lake-levels")
    await user.should_see("lake_frequency.json")
    names = await series_names(user, "Storm-driven maxima")
    assert "Fit to all maxima" not in " ".join(names)
    await user.should_see("press Fit curves")


@pytest.mark.asyncio
async def test_the_page_says_what_is_missing_instead_of_raising(user, lake_project,
                                                                  private_settings):
    from state import STATE
    STATE.open_project(lake_project.config.config_path)
    await user.open("/lake-levels")
    await user.should_see("Name a level record export")


@pytest.mark.asyncio
async def test_the_water_year_options_are_scored_with_the_adopted_month_marked(user, opened):
    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    user.find(marker="score-water-years").click()
    for _ in range(200):
        try:
            tables = user.find(marker="water-year-scores").elements
            break
        except AssertionError:          # find raises until the table exists
            await asyncio.sleep(0.1)
    rows = next(iter(tables)).rows
    assert len(rows) == 12
    assert [row["start_month"] for row in rows if row["adopted"]] == ["October"]


@needs_bryan
@pytest.mark.asyncio
async def test_fitting_draws_the_curves_and_the_design_floods(user, opened, monkeypatch):
    from state import STATE
    monkeypatch.setattr(STATE.settings, "bryan_python", BRYAN_PYTHON)
    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    user.find(marker="fit-curves").click()
    names = await series_names(user, "Fit to all maxima", seconds=90)
    assert "Design flood envelope" in names
    assert any(name.startswith("90% band, all maxima") for name in names)
    assert (opened.parent / "_lake_frequency").is_dir()


@needs_bryan
@pytest.mark.asyncio
async def test_the_figure_export_writes_the_report_png(user, opened, monkeypatch, tmp_path):
    from state import STATE
    monkeypatch.setattr(STATE.settings, "bryan_python", BRYAN_PYTHON)
    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    user.find("Export figure").click()
    await user.should_see("Export the report figure")
    user.find("Output folder").clear().type(str(tmp_path / "figures"))
    user.find("Base name").clear().type("DAM")
    user.find("Include the design floods").click()          # off: the record figure
    user.find(marker="export-run").click()
    for _ in range(900):
        if (tmp_path / "figures" / "DAM_record.png").is_file():
            break
        await asyncio.sleep(0.1)
    assert (tmp_path / "figures" / "DAM_record.png").stat().st_size > 10_000
    saved = json.loads((opened.parent / "lake_frequency.json").read_text())
    assert saved["export"]["name"] == "DAM"
