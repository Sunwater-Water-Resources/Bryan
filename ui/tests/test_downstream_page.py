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


@pytest.mark.asyncio
async def test_the_page_offers_urbs_and_the_dam_failure_outflow(user, tmp_path):
    config = _project(tmp_path)
    database = tmp_path / "CLD_mc_24h_E009_GWL1p3_RGN__mcdf.parquet"
    database.write_bytes(b"")
    _selection(tmp_path, database)
    _open(config)
    await user.open("/downstream")
    await user.should_see("Run URBS on each storm")
    await user.should_see("Dam outflow from a dam failure model")
    await user.should_see("Write hydraulic model inflow files")
