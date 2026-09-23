"""Render the Figures page against a miniature study with one figure."""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import GROUP, build_study                   # noqa: E402


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
    from core import figures, reporttables as rt, study as studies
    from state import STATE

    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    figures.put(study, figures.new_spec(filename="levels", curves=[
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "GWL 1.3"}]))
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
async def test_without_a_study_the_page_sends_you_to_the_report_page(user, private_settings):
    from state import STATE
    STATE.study = None
    await user.open("/figures")
    await user.should_see("Go to Report")


@pytest.mark.asyncio
async def test_the_figure_is_previewed_with_its_label(user, opened):
    await user.open("/figures")
    chart = await seen(user, "preview-levels")
    assert [s["name"] for s in chart.options["series"]] == ["GWL 1.3"]


@pytest.mark.asyncio
async def test_a_label_edited_beside_the_preview_is_kept_for_this_figure(user, opened):
    from core import figures, study as studies
    await user.open("/figures")
    await seen(user, "preview-levels")
    box = next(iter(user.find(marker="label-levels-0").elements))
    box.set_value("Near-term")
    user.find(marker="label-levels-0").trigger("blur")
    for _ in range(50):
        stored = figures.figures(studies.load_study(opened.path))[0]
        if stored["curves"][0]["label"] == "Near-term":
            break
        await asyncio.sleep(0.1)
    assert stored["curves"][0]["label"] == "Near-term"
