"""Render the Downstream page against a project with a saved selection.

``core/downstream.py`` is tested against paths, so nothing there notices the
page handing it the open ``Project`` instead of the folder - which is a 500 on
load, and only on load. That is what this exists to catch.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from conftest import MONTE_CARLO_COLUMNS, write_workbook          # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


def _project(tmp_path):
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, [{
        "Include": "yes", "Method": "monte carlo", "Duration": 24,
        "Run models": "yes", "Analyse results": "yes", "GWL": 1.3,
        "Output file": "sims_mc\\results\\CLD_mc_24h_GWL1p3",
        "Config file": "mc.json"}])
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return config


def _selection(tmp_path, database):
    path = tmp_path / "GWL1p3_representative_events.json"
    path.write_text(json.dumps({"targets": [dict(
        kind="aep", value=1000.0, result_type="level", rain_aep=None, source="",
        output_file="CLD_mc_24h", database=str(database), count=10, picked=42,
        comment="")], "settings": {}}))
    return path


def _open(config):
    from state import STATE
    STATE.open_project(config)
    return STATE


@pytest.mark.asyncio
async def test_the_page_plans_the_selection_it_finds_under_the_project(user, tmp_path):
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    _open(config)
    await user.open("/downstream")
    await user.should_see("1 storm file would be written")
    from nicegui import ui as nicegui_ui
    table = next(iter(user.find(kind=nicegui_ui.table).elements))
    assert [row["file"] for row in table.rows] == ["000042_24h_GWL1p3.24"]


@pytest.mark.asyncio
async def test_a_project_with_no_saved_selection_says_so(user, tmp_path):
    _open(_project(tmp_path))
    await user.open("/downstream")
    await user.should_see("No saved event selections under this project.")


@pytest.mark.asyncio
async def test_the_suffix_starts_blank_and_a_clash_suggests_one(user, tmp_path):
    # the same Monte Carlo run re-routed at two full supply levels: one realisation, one
    # duration, one warming level, and so one storm file name between them. The page does
    # not fill a suffix in - it warns, and names the one that would keep them apart.
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    one = _selection(tmp_path, database)
    one.rename(tmp_path / "FSL_Design_GWL1p3_representative_events.json")
    two = _selection(tmp_path, database)
    two.rename(tmp_path / "RFSL_Design_GWL1p3_representative_events.json")
    _open(config)
    await user.open("/downstream")
    from nicegui import ui as nicegui_ui
    table = next(iter(user.find(kind=nicegui_ui.table).elements))
    file, = [row["file"] for row in table.rows]
    assert file == "000042_24h_GWL1p3.24"
    await user.should_see("one set replaces the other")
    await user.should_see("would keep them apart")
    # one click takes the suggestion, and the clash goes
    user.find(marker="use-suffix").click()
    await user.should_not_see("one set replaces the other")
    assert next(iter(user.find(marker="suffix").elements)).value in ("_FSL", "_RFSL")


@pytest.mark.asyncio
async def test_the_page_offers_urbs_and_leaves_the_dam_outflow_to_the_dam_failure_project(
        user, tmp_path):
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    _open(config)
    await user.open("/downstream")
    await user.should_see("Run URBS on each storm")
    await user.should_see("a dam failure project finds the runs by name", retries=3)
    await user.should_not_see("Write hydraulic model inflow files")


@pytest.mark.asyncio
async def test_a_run_shows_its_progress_and_how_it_ended(user, tmp_path, monkeypatch):
    import subprocess
    import sys
    from core import downstream
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model="urbs.json")
    script = ("import time; print('  wrote 000042_24h_GWL1p3.24', flush=True); time.sleep(1.5); "
              "print('1 storm files; what went into them: x.csv', flush=True)")

    def fake_launch(argv, cwd, log_path=None):
        with open(log_path, "w") as handle:
            return subprocess.Popen([sys.executable, "-c", script], stdout=handle,
                                    stderr=subprocess.STDOUT)

    monkeypatch.setattr(downstream, "launch", fake_launch)
    await user.open("/downstream")
    user.find("Check without writing").click()
    await user.should_see("wrote 000042_24h_GWL1p3.24", retries=40)
    await user.should_see("Finished. 1 storm files.", retries=60)


def _model(tmp_path, *, extra_dams=()):
    """A regional model config whose storms and results folders sit beside it."""
    (tmp_path / "regional" / "storms").mkdir(parents=True)
    (tmp_path / "regional" / "results").mkdir()
    config = tmp_path / "regional" / "urbs_config.json"
    config.write_text(json.dumps({
        "model_folder": "model", "storms_folder": "storms", "results_folder": "results",
        "vec_file": "x.vec", "additional_dams": [dict(location=d) for d in extra_dams]}))
    return config


@pytest.mark.asyncio
async def test_the_results_tab_says_which_runs_have_finished(user, tmp_path):
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    model = _model(tmp_path, extra_dams=["KROOMBIT"])
    regional = model.parent
    (regional / "storms" / "000042_24h_GWL1p3.24").write_text("")
    (regional / "results" / "000042_24h_GWL1p3.q").write_text("")
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model=str(model))
    await user.open("/downstream?tab=Results")
    await user.should_see("1 finished")
    await user.should_see("KROOMBIT from the event's lake_z")
    from nicegui import ui as nicegui_ui
    table = next(iter(user.find(marker="runs-table").elements))
    assert [(row["run"], row["status"]) for row in table.rows] == [("000042_24h_GWL1p3", "finished")]



@pytest.mark.asyncio
async def test_the_results_tab_holds_each_dam_against_its_crest(user, tmp_path):
    import shutil
    from pathlib import Path
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_72h_E013_GWL1p3_RFSL__mcdf.parquet"
    database.write_bytes(b"")
    path = tmp_path / "RFSL_Design_GWL1p3_representative_events.json"
    path.write_text(json.dumps({"targets": [dict(
        kind="aep", value=100.0, result_type="level", rain_aep=None, source="",
        output_file="x", database=str(database), count=10, picked=5044, name="1 in 100",
        comment="")], "settings": {}}))
    model = _model(tmp_path, extra_dams=["KROOMBIT"])
    data = json.loads(model.read_text())
    data["max_keys"] = {"level": "CALLIDE"}
    data["dam_levels"] = {"CALLIDE": {"DCL": 219.13, "PMF": {"1.3": 221.5}},
                          "KROOMBIT": {"DCL": 270.7}}
    model.write_text(json.dumps(data))
    real = Path(__file__).parent / "data" / "regional" / "005044_72h_GWL1p3_RFSL.h"
    shutil.copy(real, model.parent / "results" / "005044_72h_GWL1p3.h")
    (model.parent / "results" / "005044_72h_GWL1p3.q").write_text("")
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model=str(model))
    await user.open("/downstream?tab=Results")
    await user.should_see("0 of 1 events take every dam to its crest level")
    table = next(iter(user.find(marker="dams-table").elements))
    row, = table.rows
    assert row["KROOMBIT to DCL"] == "-2.67 m" and row["KROOMBIT to PMF"] == "not given"
    await user.should_see(marker="crest-chart")
    chart = next(iter(user.find(marker="crest-chart").elements))
    assert chart.options["yAxis"]["name"] == "Above the dam crest (m)"
    # The switch redraws the same event in m AHD, each dam against its own crest.
    user.find(marker="crest-absolute").click()
    await user.should_see(marker="crest-chart")
    chart = next(iter(user.find(marker="crest-chart").elements))
    assert chart.options["yAxis"]["name"] == "Level (m AHD)"
    names = [s["name"] for s in chart.options["series"]]
    assert "crest" not in names and "KROOMBIT" in names
    await user.should_see("not a probability")



@pytest.mark.asyncio
async def test_the_results_tab_draws_the_print_locations_and_checks_the_ends(user, tmp_path):
    import shutil
    from pathlib import Path
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_72h_E013_GWL1p3_RFSL__mcdf.parquet"
    database.write_bytes(b"")
    (tmp_path / "RFSL_Design_GWL1p3_representative_events.json").write_text(json.dumps(
        {"targets": [dict(kind="aep", value=100.0, result_type="level", rain_aep=None,
                          source="", output_file="x", database=str(database), count=10,
                          picked=5044, comment="")], "settings": {}}))
    model = _model(tmp_path)
    from test_urbsresults import FLOWS_CSV
    (model.parent / "results" / "005044_72h_GWL1p3_i.csv").write_text(FLOWS_CSV)
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model=str(model))
    await user.open("/downstream?tab=Results")
    await user.should_see("All 1 runs have receded")
    await user.should_see(marker="flows-GOOVIGEN")
    table = next(iter(user.find(marker="peaks-table").elements))
    row, = table.rows
    assert row["run"] == "005044_72h_GWL1p3" and row["GOOVIGEN"]
    await user.should_see(marker="copy-word-peaks")



@pytest.mark.asyncio
async def test_coming_back_finds_the_runs_under_the_suffix_they_were_written_with(user, tmp_path):
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    selection = _selection(tmp_path, database)
    (tmp_path / "GWL1p3_representative_events_downstream_storms.csv").write_text(
        "filename,run_name\n000042_24h_GWL1p3_RFSL.24,000042_24h_GWL1p3_RFSL\n")
    model = _model(tmp_path)
    (model.parent / "results" / "000042_24h_GWL1p3_RFSL.q").write_text("")
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model=str(model))
    await user.open("/downstream?tab=Results&selection=GWL1p3")
    await user.should_see("1 finished")
    assert next(iter(user.find(marker="suffix").elements)).value == "_RFSL"
    # the suffix box is above the tabs, so the Results tab shows it too, and the tab says
    # whose names it reads by
    await user.should_see("is named with _RFSL")
    await user.should_see("as the last run of this selection wrote them, with the suffix _RFSL")


@pytest.mark.asyncio
async def test_finished_runs_this_selection_never_wrote_are_called_out(user, tmp_path):
    # A selection never run, with the other selection's suffix in the box: its shared
    # realisations are found finished under that name, but they are the other's results.
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    model = _model(tmp_path)
    (model.parent / "results" / "000042_24h_GWL1p3_RFSL.q").write_text("")
    state = _open(config)
    state.settings.remember_downstream(config, config="storm.json", model=str(model))
    await user.open("/downstream?selection=GWL1p3")
    await user.should_see("No suffix: storm and run names")
    next(iter(user.find(marker="suffix").elements)).value = "_RFSL"
    await user.should_see("is named with _RFSL")
    from nicegui import ui as nicegui_ui
    next(iter(user.find(kind=nicegui_ui.tabs).elements)).value = "Results"
    await user.should_see("1 finished")
    await user.should_see("no run on record, so its runs are looked for under the suffix "
                          "in the box: _RFSL")
    await user.should_see("1 finished run this selection never wrote: 000042_24h_GWL1p3_RFSL")
