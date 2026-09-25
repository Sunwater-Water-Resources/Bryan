"""Render the Ensemble page against a sims list with one analysed ensemble row."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from conftest import MONTE_CARLO_COLUMNS, monte_carlo_row, write_workbook  # noqa: E402
from test_enbresults import write_run                              # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


def _config(tmp_path, rows):
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    from state import STATE
    STATE.open_project(config)
    return config


@pytest.mark.asyncio
async def test_the_page_shows_the_ensemble_and_agrees_with_bryan(user, tmp_path):
    write_run(tmp_path, "CLD_enb_test")
    _config(tmp_path, [{"Include": "yes", "Method": "ensemble",
                        "Output file": "sims_enb\\results\\CLD_enb_test",
                        "Config file": "enb.json"}])
    await user.open("/ensemble")
    await user.should_see("Critical durations")
    await user.should_see("Agrees with Bryan's CLD_enb_test_critical.csv at all 2 AEPs")
    await user.should_see("'margin (m)'")                 # the level note, in metres
    await user.should_see(marker="ensemble-table")
    await user.should_see(marker="ensemble-box")
    await user.should_see("The patterns at one AEP")


@pytest.mark.asyncio
async def test_a_sims_list_without_an_ensemble_says_so(user, tmp_path):
    _config(tmp_path, [monte_carlo_row(**{"Output file": "sims_mc\\results\\mc"})])
    await user.open("/ensemble")
    await user.should_see("No ensemble results found")


@pytest.mark.asyncio
async def test_the_ensemble_table_is_copied_as_shown(user, tmp_path, monkeypatch):
    from nicegui import ui as nicegui_ui
    copied = []
    monkeypatch.setattr(nicegui_ui.clipboard, "write", copied.append)
    write_run(tmp_path, "CLD_enb_test")
    _config(tmp_path, [{"Include": "yes", "Method": "ensemble",
                        "Output file": "sims_enb\\results\\CLD_enb_test",
                        "Config file": "enb.json"}])
    await user.open("/ensemble")
    await user.should_see(marker="ensemble-table")
    user.find(marker="copy-text-ensemble").click()
    await user.should_see("Copied as text")
    assert copied and copied[0].splitlines()[0].startswith("AEP (1 in X)\t")
