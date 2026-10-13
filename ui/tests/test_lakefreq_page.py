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


def fit_tables(user):
    try:
        return list(user.find(marker="lake-fit-table").elements)
    except AssertionError:              # find raises while there is no table
        return []


def number(text):
    """'0.5958401 (1 in 4)' -> 0.5958401."""
    return float(text.split(" (")[0])


def curve_from_table(rows, column, z):
    """The level at ``z`` worked out from the table alone, by what its labels say."""
    cell = {row["quantity"]: row[column] for row in rows}
    fsl = number(cell["Full supply level (m AHD)"])
    z1, z2 = number(cell["Plateau starts, z1"]), number(cell["Plateau ends, z2"])

    def series(prefix, variable):
        total, order = 0.0, 1
        while f"{prefix} ({variable})" + ("" if order == 1 else f"^{order}") in cell:
            label = f"{prefix} ({variable})" + ("" if order == 1 else f"^{order}")
            total += number(cell[label]) * (z - (z1 if "z1" in variable else z2)) ** order
            order += 1
        return total

    if z < z1:
        return fsl + series("Below FSL: coefficient of", "z - z1")
    if z <= z2:
        return fsl
    start = number(cell.get("Above the plateau: level at z2 (m AHD)", str(fsl)))
    return start + series("Above the plateau: coefficient of", "z - z2")


@pytest.mark.asyncio
async def test_the_fitted_curves_are_tabled_from_the_results_the_chart_shows(
        user, opened, monkeypatch):
    """No table before a fit, the table of the cached fit - which can be copied -
    and none again once a setting has changed and those results are no longer
    the page's."""
    from nicegui import ui as nicegui_ui
    from core import events, lakefreq
    from state import STATE
    from test_lakefreq import fitted_results

    groups = events.sources_by_group(STATE.project)
    settings = lakefreq.load_settings(opened)
    settings["design"]["group"] = next(iter(groups))     # as the page picks it
    job = lakefreq.build_job(opened, settings, groups).job
    results = fitted_results()
    results["grid"] = {"z": [0.0, 0.5], "aep": [0.5, 0.3]}
    for block in results["fits"].values():
        block.update(curve=[214.0, 215.0], band_lo=[213.0, 214.0], band_hi=[215.0, 216.0])
    results["fingerprint"] = lakefreq.RECORD.fingerprint(job)
    path = lakefreq.results_path(opened, job)

    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    assert fit_tables(user) == []

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results))
    await user.open("/lake-levels")
    await series_names(user, "Fit to all maxima")
    await user.should_see("Fitted curves")
    table = fit_tables(user)[0]
    assert [column["label"] for column in table.columns] == \
        ["Quantity", "All maxima", "Storm-driven maxima"]
    rows = {row["quantity"]: row for row in table.rows}
    assert rows["RMSE (m)"]["all"] == "0.1448"
    assert rows["Resamples fitted, of those drawn"]["storm"] == "378 of 400"
    assert rows["Plateau starts, z1"]["all"] == "0 (1 in 2)"

    # Copy as text: the table as shown, tab-separated under its header row.
    copied = []
    monkeypatch.setattr(nicegui_ui.clipboard, "write", copied.append)
    await user.should_see(marker="copy-word-lake-fit")
    user.find(marker="copy-text-lake-fit").click()
    await user.should_see("Copied as text")
    lines = copied[0].splitlines()
    assert lines[0].split("\t") == ["Quantity", "All maxima", "Storm-driven maxima"]
    assert [line.split("\t")[0] for line in lines[1:-1]] == \
        [row["quantity"] for row in table.rows]
    assert "RMSE (m)\t0.1448\t0.3184" in lines

    # Another setting: these are not its results, and the table goes with the curves.
    user.find(marker="lake-curve-start").type("4")
    await user.should_see("press Fit curves")
    assert fit_tables(user) == []


@needs_bryan
@pytest.mark.asyncio
async def test_the_table_of_a_real_fit_redraws_the_curve_it_describes(user, opened,
                                                                    monkeypatch):
    """The labels say what each number multiplies; worked out from the table
    alone, the levels are the fitted curve's own."""
    from core import lakefreq
    from state import STATE
    monkeypatch.setattr(STATE.settings, "bryan_python", BRYAN_PYTHON)
    await user.open("/lake-levels")
    await series_names(user, "Storm-driven maxima")
    user.find(marker="fit-curves").click()
    await series_names(user, "Fit to all maxima", seconds=90)

    table = fit_tables(user)[0]
    rows = {row["quantity"]: row for row in table.rows}
    assert rows["Form"]["all"] == "Shouldered plateau"
    assert "Resamples fitted, of those drawn" not in rows           # Fit curves: no bands
    results = next(saved for saved in (json.loads(path.read_text()) for path in
                                       (opened.parent / "_lake_frequency").glob("*.json"))
                   if "fits" in saved)
    assert rows["Maxima fitted"]["all"] == str(len(results["ams"]))
    for name, block in results["fits"].items():
        for z, level in list(zip(results["grid"]["z"], block["curve"]))[::25]:
            assert curve_from_table(table.rows, name, z) == pytest.approx(level, abs=1e-5)


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
