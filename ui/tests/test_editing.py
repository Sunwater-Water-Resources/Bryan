"""Editing saves to a new file, and never damages the master."""

from __future__ import annotations

import pandas as pd
import pytest

from conftest import CALLIDE_COLUMNS, TINAROO_COLUMNS, reservoir_row, write_workbook
from core import editing
from core.simslist import audit_formulas, read_sims_list


def build(tmp_path, rows, columns=TINAROO_COLUMNS, **kwargs):
    return read_sims_list(write_workbook(tmp_path / "master.xlsx", columns, rows,
                                         **kwargs))


def test_saving_over_the_master_is_refused(tmp_path):
    """The whole point - saving it from Python would strip its formulas."""
    sims = build(tmp_path, [reservoir_row()])
    with pytest.raises(editing.EditError, match="never written by the UI"):
        editing.save_as(sims, sims.frame, sims.path)


def test_saving_to_a_new_file_preserves_shape(tmp_path):
    sims = build(tmp_path, [reservoir_row(Duration=d) for d in (18, 24)],
                 columns=CALLIDE_COLUMNS, sheet_name="SIMS",
                 extra_sheets=("Definitions",))
    destination = tmp_path / "copy.xlsx"

    edited = editing.set_cell(sims.frame.copy(), 0, "Comment", "changed")
    editing.save_as(sims, edited, destination)

    written = pd.read_excel(destination, sheet_name=0)
    assert list(written.columns) == list(sims.frame.columns)
    assert written.loc[0, "Comment"] == "changed"

    from openpyxl import load_workbook
    assert load_workbook(destination).worksheets[0].title == "SIMS"


def test_the_master_is_untouched_by_a_save(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    before = sims.path.read_bytes()

    edited = editing.set_cell(sims.frame.copy(), 0, "Comment", "changed")
    editing.save_as(sims, edited, tmp_path / "copy.xlsx")

    assert sims.path.read_bytes() == before


def test_a_saved_copy_has_no_uncached_formulas(tmp_path):
    """A copy must never be in the state TFD_SimsList_LongList_02.xlsx is in."""
    sims = build(tmp_path, [reservoir_row(Comment="=A1+1")])
    destination = tmp_path / "copy.xlsx"
    editing.save_as(sims, sims.frame, destination)

    assert not audit_formulas(destination).is_damaged
    assert pd.read_excel(destination, sheet_name=0).loc[0, "Comment"] == "=A1+1"


def test_overwriting_needs_confirmation(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    destination = tmp_path / "copy.xlsx"
    editing.save_as(sims, sims.frame, destination)

    with pytest.raises(editing.EditError, match="already exists"):
        editing.save_as(sims, sims.frame, destination)
    editing.save_as(sims, sims.frame, destination, overwrite=True)


def test_include_is_not_forced_when_editing(tmp_path):
    """Run copies force Include = yes; an edited copy must keep what it says."""
    sims = build(tmp_path, [reservoir_row(Include="no"),
                            reservoir_row(Include="yes", Duration=24)])
    destination = tmp_path / "copy.xlsx"
    editing.save_as(sims, sims.frame, destination)

    assert pd.read_excel(destination, sheet_name=0)["Include"].tolist() == ["no", "yes"]


def test_editing_a_wholly_blank_column_works(tmp_path):
    """A column that is blank throughout reads as float64.

    pandas 3 raises rather than upcasting when a string is assigned into it,
    so every edit to an empty Comment or Lake config would fail without the
    widening in set_cell.
    """
    sims = build(tmp_path, [reservoir_row(Comment=None)])
    assert sims.frame["Comment"].dtype == float

    edited = editing.set_cell(sims.frame.copy(), 0, "Comment", "a note")
    assert edited.loc[0, "Comment"] == "a note"

    destination = tmp_path / "copy.xlsx"
    editing.save_as(sims, edited, destination)
    assert pd.read_excel(destination, sheet_name=0).loc[0, "Comment"] == "a note"


def test_set_cell_keeps_numbers_numeric(tmp_path):
    """Widening must only happen when it is needed."""
    sims = build(tmp_path, [reservoir_row(Duration=18)])
    edited = editing.set_cell(sims.frame.copy(), 0, "Duration", 24)
    assert edited["Duration"].dtype == sims.frame["Duration"].dtype
    assert edited.loc[0, "Duration"] == 24


def test_set_cell_creates_a_missing_column(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    edited = editing.set_cell(sims.frame.copy(), 0, "Group", "my group")
    assert edited.loc[0, "Group"] == "my group"


def test_diff_reports_each_changed_cell(tmp_path):
    sims = build(tmp_path, [reservoir_row(), reservoir_row(Duration=24)])
    edited = editing.set_cell(sims.frame.copy(), 1, "Comment", "note")

    changes = editing.diff_frames(sims.frame, edited)
    assert any(c.column == "Comment" and c.after == "note" for c in changes)


def test_duplicate_rows_copies_rather_than_blanks(tmp_path):
    """A blank row would be missing every column pre-flight requires."""
    sims = build(tmp_path, [reservoir_row(**{"Output file": "orig"})])
    out = editing.duplicate_rows(sims.frame, [0], overrides={"GWL": 2.7})

    assert len(out) == 2
    assert out.loc[1, "Output file"] == "orig"
    assert out.loc[1, "GWL"] == 2.7
    assert out.loc[1, "SQ file"] == sims.frame.loc[0, "SQ file"]


def test_duplicate_with_no_selection_is_an_error(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    with pytest.raises(editing.EditError):
        editing.duplicate_rows(sims.frame, [])


def test_adding_a_group_column_writes_the_derived_keys(tmp_path):
    rows = [reservoir_row(Duration=d, **{"Output file": f"case_{d}h",
                                         "Output suffix": "A"})
            for d in (18, 24)]
    sims = build(tmp_path, rows)
    out = editing.add_group_column(sims.frame)

    assert "Group" in out.columns
    assert out["Group"].nunique() == 1, "both durations are one group"

    destination = tmp_path / "with_groups.xlsx"
    editing.save_as(sims, out, destination)
    written = pd.read_excel(destination, sheet_name=0)
    assert "Group" in written.columns


def test_added_rows_are_written(tmp_path):
    sims = build(tmp_path, [reservoir_row(**{"Output file": "orig"})])
    out = editing.duplicate_rows(sims.frame, [0],
                                 overrides={"Output file": "copy"})
    destination = tmp_path / "copy.xlsx"
    editing.save_as(sims, out, destination)

    written = pd.read_excel(destination, sheet_name=0)
    assert written["Output file"].tolist() == ["orig", "copy"]


def test_default_save_name_sits_beside_the_master(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    name = editing.default_save_name(sims.path)
    assert name.startswith("master_") and name.endswith("_01.xlsx")
