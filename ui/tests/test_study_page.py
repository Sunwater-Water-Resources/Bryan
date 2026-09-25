"""The Study page, and reopening the last study when the launcher starts."""

from __future__ import annotations

import json
import time

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import build_study                           # noqa: E402


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
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    STATE.study, STATE.study_problem = None, ""
    yield STATE
    STATE.study, STATE.study_problem = None, ""


# -- the page ----------------------------------------------------------------------

@pytest.mark.asyncio
async def test_a_study_is_opened_and_what_it_holds_is_shown(user, private_settings,
                                                          tmp_path):
    study = build_study(tmp_path / "study")
    await user.open("/study")
    await user.should_see("Open a study")
    user.find(marker="study-path").clear().type(str(study.path))
    user.find(marker="open-study").click()
    await user.should_see(str(study.path))
    await user.should_see(marker="study-contents")
    await user.should_see(f"{len(study.runs)} runs")
    await user.should_see("Runs: " + ", ".join(run["name"] for run in study.runs))
    assert private_settings.study.path == study.path
    assert private_settings.settings.last_study == str(study.path)


@pytest.mark.asyncio
async def test_a_new_study_is_created_where_asked(user, private_settings, tmp_path):
    folder = tmp_path / "Kroombit"
    folder.mkdir()
    await user.open("/study")
    user.find(marker="study-path").clear().type(str(folder))
    user.find(marker="study-name").type("Kroombit Dam")
    user.find(marker="new-study").click()
    await user.should_see(marker="study-card")
    saved = json.loads((folder / "bryan_study.json").read_text(encoding="utf-8"))
    assert saved["name"] == "Kroombit Dam"


@pytest.mark.asyncio
async def test_closing_the_study_forgets_it_for_next_time(user, private_settings, tmp_path):
    private_settings.open_study(build_study(tmp_path / "study").path)
    await user.open("/study")
    user.find(marker="close-study").click()
    await user.should_see("Open a study")
    assert private_settings.study is None and private_settings.settings.last_study == ""


@pytest.mark.asyncio
async def test_why_the_last_study_did_not_reopen_is_said(user, private_settings):
    private_settings.study_problem = "The last study, F:\\x.json, did not answer"
    await user.open("/study")
    await user.should_see("did not answer")
    await user.open("/pmf")                  # and on the pages that need a study
    await user.should_see("did not answer")
    await user.should_see("Go to Study")


@pytest.mark.asyncio
async def test_the_study_is_first_in_the_menu(user, private_settings):
    from layout import NAV
    assert NAV[0] == ("Study", "/study")
    await user.open("/study")
    await user.should_see("Study")


# -- reopening at start-up ------------------------------------------------------------

def test_the_last_study_is_reopened_at_start_up(private_settings, tmp_path):
    study = build_study(tmp_path / "study")
    private_settings.settings.last_study = str(study.path)
    assert private_settings.reopen_last_study() is not None
    assert private_settings.study.path == study.path and not private_settings.study_problem


def test_a_study_that_is_not_there_leaves_none_open_and_says_why(private_settings,
                                                                tmp_path):
    private_settings.settings.last_study = str(tmp_path / "unplugged" / "bryan_study.json")
    assert private_settings.reopen_last_study() is None
    assert private_settings.study is None
    assert "could not be reopened" in private_settings.study_problem
    assert private_settings.settings.last_study      # kept, for when the disk is back


def test_a_study_that_does_not_answer_does_not_hold_up_the_launcher(private_settings,
                                                                   tmp_path, monkeypatch):
    from core import study as studies

    def slow(path):
        time.sleep(2.0)                                # a share that does not answer
        raise studies.StudyError("never reached")

    monkeypatch.setattr(studies, "load_study", slow)
    private_settings.settings.last_study = str(tmp_path / "bryan_study.json")
    started = time.perf_counter()
    assert private_settings.reopen_last_study(timeout=0.2) is None
    assert time.perf_counter() - started < 1.5
    assert "did not answer within 0.2 s" in private_settings.study_problem


def test_nothing_is_reopened_when_no_study_was_open(private_settings):
    private_settings.settings.last_study = ""
    assert private_settings.reopen_last_study() is None
    assert not private_settings.study_problem
