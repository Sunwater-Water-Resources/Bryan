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
    # Fit curves: the curves alone, quickly, and no band yet.
    user.find(marker="fit-curves").click()
    names = await series_names(user, "Fit to all maxima", seconds=90)
    assert "Design flood envelope" in names
    assert not any(name.startswith("90% band") for name in names)
    await user.should_see(marker="no-bands")
    assert (opened.parent / "_lake_frequency").is_dir()

    # Resample bands: the same curves, now with their bands.
    user.find(marker="resample-bands").click()
    await series_names(user, "90% band, all maxima", seconds=90)
    await user.should_not_see(marker="no-bands")


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
    figure = tmp_path / "figures" / "DAM_record.png"
    # Wait for the file to be written, not just created: the export runs in Bryan's
    # own interpreter, and a check the moment the name appeared read a size of 0
    # about one run in three.
    for _ in range(900):
        if figure.is_file() and figure.stat().st_size > 10_000:
            break
        await asyncio.sleep(0.1)
    assert figure.stat().st_size > 10_000
    saved = json.loads((opened.parent / "lake_frequency.json").read_text())
    assert saved["export"]["name"] == "DAM"


@pytest.mark.asyncio
async def test_a_band_from_too_few_resamples_is_warned_about(user, opened):
    from core import events, lakefreq
    from state import STATE
    from test_lakefreq import results_for

    project = STATE.project
    groups = events.sources_by_group(project)
    settings = lakefreq.load_settings(opened)
    settings["design"]["group"] = next(iter(groups))     # as the page picks it
    job = lakefreq.build_job(opened, settings, groups).job
    results = results_for(None)
    results["fingerprint"] = lakefreq.RECORD.fingerprint(job)
    results["fits"]["all"]["warning"] = "only 100 of 400 resamples could be fitted"
    path = lakefreq.results_path(opened, job)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results))

    await user.open("/lake-levels")
    await series_names(user, "Fit to all maxima")
    await user.should_see("only 100 of 400 resamples could be fitted")


@pytest.mark.asyncio
async def test_the_axis_and_curve_starts_are_blank_then_saved_and_move_the_chart(user, opened):
    from core import lakefreq

    def chart():
        return next(iter(user.find(marker="lake-chart").elements)).options

    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    boxes = {key: next(iter(user.find(marker=marker).elements))
             for key, marker in (("frequent_aep_1_in_x", "lake-axis-start"),
                                 ("curve_frequent_aep_1_in_x", "lake-curve-start"))}
    assert [box.value for box in boxes.values()] == [None, None]
    from_1ey = chart()["xAxis"]["min"]
    assert from_1ey == pytest.approx(lakefreq.RECORD.normal_variate_of_probability(
        lakefreq.RECORD.EY1_AEP))

    user.find(marker="lake-axis-start").type("5")
    user.find(marker="lake-curve-start").type("10")
    saved = json.loads((opened.parent / "lake_frequency.json").read_text())
    assert saved["axes"]["frequent_aep_1_in_x"] == 5
    assert saved["axes"]["curve_frequent_aep_1_in_x"] == 10
    assert chart()["xAxis"]["min"] == pytest.approx(
        lakefreq.RECORD.normal_variate_of_probability(0.2))

    # Cleared, the chart is back in the 1EY frame; a start that cannot be one is said.
    boxes["frequent_aep_1_in_x"].set_value(None)
    assert chart()["xAxis"]["min"] == from_1ey
    boxes["curve_frequent_aep_1_in_x"].set_value(0.5)
    await user.should_see("is not an AEP")
    assert chart()["xAxis"]["min"] == from_1ey


@needs_bryan
@pytest.mark.asyncio
async def test_fitted_curves_start_where_the_page_says(user, opened, monkeypatch):
    """Through the real script: the curves on the page start at the curve start,
    and the design flood envelope where it always did."""
    from core import lakefreq
    from state import STATE
    monkeypatch.setattr(STATE.settings, "bryan_python", BRYAN_PYTHON)
    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    user.find(marker="lake-curve-start").type("4")
    await user.should_see("press Fit curves")
    user.find(marker="fit-curves").click()
    await series_names(user, "Fit to all maxima", seconds=90)

    chart = next(iter(user.find(marker="lake-chart").elements)).options
    starts = {series["name"].split(" (")[0]: series["data"][0][0]
              for series in chart["series"] if series["type"] == "line"}
    z_of = lakefreq.RECORD.normal_variate_of_probability
    assert z_of(0.25) - 1e-9 <= starts["Fit to all maxima"] < z_of(0.25) + 0.05
    assert starts["Design flood envelope"] == pytest.approx(z_of(lakefreq.RECORD.EY1_AEP))


@pytest.mark.asyncio
async def test_with_a_study_the_list_s_record_is_kept_and_can_become_the_study_s(
        user, opened):
    """Kroombit's case: the gauge export was typed here before the study kept one.
    With the study open the page still reads it, as this analysis's own record,
    and one button makes it the study's gauge exports."""
    from core import dam as dams, study as studies
    from state import STATE
    study = studies.new_study(opened.parent.parent / "bryan_study.json", "Dam")
    STATE.study = study
    try:
        await user.open("/lake-levels")
        await user.should_see(marker="lake-record-source")
        await user.should_see(marker="lake-water-year")
        await series_names(user, "Storm-driven maxima")          # the record is read
        user.find(marker="lake-make-gauges").click()
        await user.should_see(marker="lake-edit-dam")
        assert dams.gauges(dams.settings(studies.load_study(study.path)))
    finally:
        STATE.study = None


@pytest.mark.asyncio
async def test_without_a_study_the_page_says_where_its_settings_are(user, opened):
    from state import STATE
    STATE.study = None
    await user.open("/lake-levels")
    await user.should_see(marker="lake-own-settings")
    await user.should_see("Open a study to share the record")
