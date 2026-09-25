"""Render the PMF page against a miniature study with one PMF group."""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import GROUP, PMF_GROUP, build_study        # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_settings(monkeypatch, tmp_path):
    import settings as settings_module
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")


@pytest.fixture
def opened(tmp_path, private_settings):
    from core import ensemble, reporttables as rt, study as studies
    from state import STATE

    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    section = ensemble.settings(study)
    section["groups"].append({"label": "Near-term RFSL",
                              "ensemble": {"run": "E099 PMF", "group": PMF_GROUP},
                              "mc": {"run": "E099 RFSL", "group": GROUP}})
    ensemble.store(study, section)
    study.save()
    STATE.open_study(study.path)
    yield study
    STATE.study = None


async def seen(user, marker, seconds=10.0):
    for _ in range(int(seconds / 0.1)):
        try:
            return next(iter(user.find(marker=marker).elements))
        except AssertionError:
            await asyncio.sleep(0.1)
    raise AssertionError(f"{marker} never appeared")


@pytest.mark.asyncio
async def test_without_a_study_the_page_sends_you_to_the_study_page(user, private_settings):
    from state import STATE
    STATE.study = None
    await user.open("/pmf")
    await user.should_see("No study is open")
    await user.should_see("Go to Study")


@pytest.mark.asyncio
async def test_the_pmf_and_its_notional_aep_are_shown(user, opened):
    await user.open("/pmf")
    highest = await seen(user, "pmf-highest")
    assert highest.text.startswith("221.38 m AHD, 9 h")
    answer = await seen(user, "pmf-aep")
    assert "is at 1 in" in answer.text
    await seen(user, "pmf-grid")
    await seen(user, "pmf-fit")


@pytest.mark.asyncio
async def test_a_changed_setting_is_kept_in_the_study(user, opened):
    from core import ensemble, study as studies
    await user.open("/pmf")
    await seen(user, "pmf-aep")
    user.find(marker="pmf-degree").click()      # the toggle; set the value directly
    toggle = next(iter(user.find(marker="pmf-degree").elements))
    toggle.set_value(2)
    await asyncio.sleep(0.3)
    stored = ensemble.settings(studies.load_study(opened.path))
    assert stored["groups"][0]["degree"] == 2
