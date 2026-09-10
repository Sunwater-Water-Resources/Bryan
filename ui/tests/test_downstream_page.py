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
