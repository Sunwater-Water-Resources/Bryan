"""The settings panel behind the gear: saving, and Check setup."""

from __future__ import annotations

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_state(monkeypatch, tmp_path):
    import settings as settings_module
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    STATE.project = None
    yield STATE
    STATE.project = None


def _without(missing):
    def probe(python, modules):
        return {"python": "3.12.1", "executable": str(python),
                "packages": {name: {"error": "ModuleNotFoundError"} if name in missing
                             else {"version": "9.9"} for name in modules}}
    return probe


@pytest.mark.asyncio
async def test_the_gear_is_on_every_page(user, private_state):
    for route in ("/", "/study", "/select", "/report"):
        await user.open(route)
        await user.should_see(marker="open-settings")


@pytest.mark.asyncio
async def test_settings_are_saved_for_this_computer(user, private_state, tmp_path):
    import settings as settings_module
    grids = tmp_path / "awap"
    grids.mkdir()
    await user.open("/")
    user.find(marker="open-settings").click()
    await user.should_see(marker="settings-panel")
    user.find(marker="settings-awap").clear().type(str(grids))
    user.find(marker="save-settings").click()
    assert private_state.settings.awap_folder == str(grids)
    assert settings_module.UiSettings.load().awap_folder == str(grids)


@pytest.mark.asyncio
async def test_check_setup_says_what_is_missing_and_writes_it_down(user, private_state,
                                                                  monkeypatch, tmp_path):
    from core import setupcheck
    monkeypatch.setattr(setupcheck, "probe_interpreter", _without({"scipy"}))
    await user.open("/")
    user.find(marker="open-settings").click()
    user.find(marker="check-setup").click()
    await user.should_see(marker="setup-verdict")
    await user.should_see("Not ready")
    await user.should_see("Not installed in Bryan's interpreter")
    await user.should_see(marker="setup-report")
    report = tmp_path / setupcheck.REPORT_NAME
    assert "FAIL  scipy" in report.read_text(encoding="utf-8")


@pytest.mark.asyncio
async def test_the_simulations_page_no_longer_holds_the_settings(user, private_state):
    await user.open("/")
    await user.should_not_see("Python interpreter")
