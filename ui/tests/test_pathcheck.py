"""Path boxes that check themselves, and the Browse dialog behind them."""

from __future__ import annotations

import pytest

from core import pathcheck as pc


@pytest.fixture
def folder(tmp_path):
    (tmp_path / "gauges").mkdir()
    (tmp_path / "gauges" / "130360A.csv").write_text("x", encoding="utf-8")
    (tmp_path / "gauges" / "notes.txt").write_text("x", encoding="utf-8")
    (tmp_path / "gauges" / "~$open.csv").write_text("x", encoding="utf-8")
    (tmp_path / "storage.els").write_text("x", encoding="utf-8")
    return tmp_path


def test_blank_says_nothing(folder):
    assert pc.check("  ", base=folder).status == pc.BLANK


def test_a_relative_path_says_what_it_is_relative_to(folder):
    found = pc.check("gauges/130360A.csv", base=folder, base_name="the study folder")
    assert found.status == pc.OK and found.resolved == folder / "gauges" / "130360A.csv"
    assert found.relative_note == "relative to the study folder"
    absolute = pc.check(str(folder / "storage.els"), base=folder)
    assert absolute.status == pc.OK and absolute.relative_note == ""


def test_a_copied_path_with_quotes_and_marks_is_found(folder):
    text = '‪"' + str(folder / "storage.els") + '"'
    assert pc.check(text).status == pc.OK


@pytest.mark.parametrize("text, expect, suffixes, status, words", [
    ("nowhere.els", pc.FILE, (), pc.MISSING, "not found"),
    ("gauges", pc.FILE, (), pc.WRONG_TYPE, "is a folder, not a file"),
    ("storage.els", pc.FILE, (".xlsx",), pc.WRONG_TYPE, "expected .xlsx, not .els"),
    ("storage.els", pc.FOLDER, (), pc.WRONG_TYPE, "is a file, not a folder"),
    ("gauges", pc.FOLDER, (), pc.OK, ""),
    ("gauges", pc.EITHER, (".json",), pc.OK, ""),
    ("storage.els", pc.OUTPUT, (), pc.WRITTEN, "there; it will be replaced"),
    ("new.csv", pc.OUTPUT, (), pc.WRITTEN, "will be written"),
    ("out/new.csv", pc.OUTPUT, (), pc.WRITTEN, "will be made"),
    ("gauges", pc.OUTPUT_FOLDER, (), pc.WRITTEN, "written into"),
    ("results", pc.OUTPUT_FOLDER, (), pc.WRITTEN, "will be made"),
    ("storage.els", pc.OUTPUT_FOLDER, (), pc.WRONG_TYPE, "is a file"),
])
def test_each_kind_of_box_says_what_it_finds(folder, text, expect, suffixes, status, words):
    found = pc.check(text, base=folder, expect=expect, suffixes=suffixes)
    assert found.status == status
    assert words in found.message
    assert str(found.resolved) in found.message          # the full path, always


def test_a_picked_file_under_the_base_is_kept_relative(folder, tmp_path_factory):
    assert pc.stored_text(folder / "gauges" / "130360A.csv", folder) == "gauges/130360A.csv"
    elsewhere = tmp_path_factory.mktemp("other") / "x.csv"
    assert pc.stored_text(elsewhere, folder) == str(elsewhere)
    assert pc.stored_text(elsewhere) == str(elsewhere)


def test_the_listing_keeps_folders_and_the_wanted_files(folder):
    found = pc.listing(folder / "gauges", suffixes=(".csv",))
    assert [path.name for path in found.files] == ["130360A.csv"]
    assert found.hidden_files == 1                       # notes.txt; ~$ lock files never
    everything = pc.listing(folder / "gauges", suffixes=(".csv",), show_all=True)
    assert [path.name for path in everything.files] == ["130360A.csv", "notes.txt"]
    top = pc.listing(folder)
    assert [path.name for path in top.folders] == ["gauges"]


def test_a_folder_that_cannot_be_listed_says_so(folder):
    assert pc.listing(folder / "gone").problem.startswith("cannot be listed")


def test_browse_opens_where_the_box_points(folder):
    assert pc.start_folder("gauges/130360A.csv", folder) == folder / "gauges"
    assert pc.start_folder("gauges", folder) == folder / "gauges"
    assert pc.start_folder("gone/x.csv", folder) == folder
    assert pc.start_folder("", None) == pc.Path.home()


def test_there_is_at_least_one_drive():
    assert pc.drives()
