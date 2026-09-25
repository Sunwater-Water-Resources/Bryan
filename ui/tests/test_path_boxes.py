"""The path boxes on a page: what they say under themselves, and Browse."""

from __future__ import annotations

import json

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
def opened(monkeypatch, tmp_path):
    import settings as settings_module
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    study = build_study(tmp_path / "study")
    (study.folder / "survey").mkdir()
    (study.folder / "survey" / "callide.els").write_text("EL,A,V\n", encoding="utf-8")
    (study.folder / "survey" / "old.xlsx").write_text("", encoding="utf-8")
    STATE.open_study(study.path)
    yield study
    STATE.study, STATE.project = None, None


def _saved_dam(study):
    return json.loads(study.path.read_text(encoding="utf-8"))["dam"]


@pytest.mark.asyncio
async def test_a_box_says_whether_its_path_is_there(user, opened):
    await user.open("/study")
    box = user.find(marker="dam-storage")
    box.clear().type("survey/nothing.els")
    await user.should_see("not found")
    box.clear().type("survey/callide.els")
    await user.should_see("relative to the study folder")
    await user.should_not_see("not found")
    box.clear().type("survey")
    await user.should_see("is a folder, not a file")


@pytest.mark.asyncio
async def test_a_file_picked_with_browse_is_kept_relative(user, opened):
    await user.open("/study")
    user.find(marker="dam-storage-browse").click()
    await user.should_see(marker="browse-dialog")
    await user.should_see(marker="browse-folder-survey")
    user.find(marker="browse-folder-survey").click()
    await user.should_see(marker="browse-file-callide.els")
    await user.should_not_see(marker="browse-file-old.xlsx")      # not a storage table
    await user.should_see("1 other file not shown")
    user.find(marker="browse-file-callide.els").click()
    user.find(marker="browse-pick").click()
    assert _saved_dam(opened)["storage"] == "survey/callide.els"


@pytest.mark.asyncio
async def test_a_list_checks_each_line_and_browse_adds_one(user, opened):
    await user.open("/study")
    box = user.find(marker="dam-gauges")
    box.clear().type("survey/callide.els\nsurvey/gone.csv")
    await user.should_see("not found")
    box.trigger("blur")
    assert _saved_dam(opened)["gauges"] == ["survey/callide.els", "survey/gone.csv"]

    user.find(marker="dam-gauges-browse").click()        # opens where the last line is
    await user.should_see(marker="browse-file-callide.els")
    user.find(marker="browse-file-old.xlsx").click()
    user.find(marker="browse-pick").click()
    assert _saved_dam(opened)["gauges"][-1] == "survey/old.xlsx"


@pytest.mark.asyncio
async def test_browse_goes_up_and_to_the_study_folder(user, opened):
    await user.open("/study")
    user.find(marker="dam-storage-browse").click()
    user.find(marker="browse-folder-survey").click()
    user.find(marker="browse-up").click()
    await user.should_see(marker="browse-folder-survey")
    user.find(marker="browse-up").click()
    user.find(marker="browse-place-Study folder").click()
    await user.should_see(marker="browse-folder-survey")


@pytest.mark.asyncio
async def test_an_output_says_it_will_be_written(user, opened):
    await user.open("/lake-record")
    await user.should_see("Write the series to")
    await user.should_see("will be written")
