"""Addresses that carry the view: the run, the group, the result type, the table."""

from __future__ import annotations

import asyncio
from urllib.parse import quote, urlencode

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import GROUP, build_study                    # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def opened(monkeypatch, tmp_path):
    import settings as settings_module
    from core import reporttables as rt, study as studies
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    for title in ("Table 26 near-term", "Table 27 long-term"):
        spec = rt.new_spec(rt.DESIGN_FLOODS)
        spec.update(title=title, source={"run": "E099 RFSL", "group": GROUP},
                    aeps=[5, 10], pmp_aep=None)
        study.put_table(spec)
    study.save()
    STATE.open_study(study.path)
    STATE.open_project(study.run_config_path("E099 RFSL"))
    yield study
    STATE.study, STATE.project = None, None


@pytest.fixture
def addresses(monkeypatch):
    """Every address the page writes, newest last."""
    import address
    written = []
    monkeypatch.setattr(address, "_replace", written.append)
    return written


def test_an_address_leaves_out_what_is_blank():
    import address
    assert address.url("/results", run="E013 RFSL", group="", type=None) \
        == "/results?run=E013+RFSL"
    assert address.url("/report") == "/report"


@pytest.mark.asyncio
async def test_the_run_in_the_address_is_opened(user, opened):
    from state import STATE
    await user.open(f"/select?run={quote('E099 PMF')}")
    assert STATE.project.config.config_path.resolve() \
        == opened.run_config_path("E099 PMF").resolve()
    header = next(iter(user.find(marker="header-run").elements))
    assert header.text == "E099 PMF"


@pytest.mark.asyncio
async def test_a_run_the_study_does_not_have_is_said(user, opened):
    await user.open("/select?run=E200")
    await user.should_see("names a run, E200, that the study does not have")


@pytest.mark.asyncio
async def test_every_page_puts_the_open_run_in_its_address(user, opened, addresses):
    await user.open("/select")
    assert any("run=E099+RFSL" in written for written in addresses)


@pytest.mark.asyncio
async def test_the_select_filters_come_from_the_address(user, opened, addresses):
    await user.open(f"/select?group={quote(GROUP)}&status={quote('up to date')}")
    from nicegui import ui
    selects = {element.props.get("label"): element.value
               for element in user.find(kind=ui.select).elements}
    assert selects["Group"] == GROUP
    assert selects["Status"] == "up to date"


@pytest.mark.asyncio
async def test_the_results_group_and_type_are_kept(user, opened, addresses):
    await user.open("/results?type=nonsense")      # a type it does not have: the first
    assert any(urlencode({"group": GROUP}) in written for written in addresses)


@pytest.mark.asyncio
async def test_a_report_table_named_in_the_address_is_opened(user, opened, addresses):
    await user.open("/report?table=table-27-long-term")
    for _ in range(30):
        chevron = next(iter(user.find(marker="chevron-table-27-long-term").elements))
        if chevron.name == "expand_more":
            break
        await asyncio.sleep(0.1)
    assert chevron.name == "expand_more"
    chevron = next(iter(user.find(marker="chevron-table-26-near-term").elements))
    assert chevron.name == "chevron_right"


@pytest.mark.asyncio
async def test_the_runs_panel_swaps_the_run_in_the_address(user, opened):
    from nicegui import ui
    gone_to = []
    await user.open(f"/results?run={quote('E099 RFSL')}&type=level")
    original = ui.navigate.to
    try:
        ui.navigate.to = gone_to.append
        user.find(marker="panel-run-E099 PMF").click()
    finally:
        ui.navigate.to = original
    assert gone_to and gone_to[0].startswith("/results?")
    assert "run=E099+PMF" in gone_to[0] and "type=level" in gone_to[0]
