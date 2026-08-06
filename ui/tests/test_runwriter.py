"""The run copy must be exactly what the UI showed, and readable by Bryan."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from conftest import (CALLIDE_COLUMNS, TINAROO_COLUMNS, reservoir_row,
                      write_workbook)
from core import runplan, runwriter
from core.config import load_sims_config, run_log_path_for
from core.simslist import read_sims_list
from core.runwriter import WriteError, write_chunk_list


def build(tmp_path, rows, columns=TINAROO_COLUMNS, **kwargs):
    path = write_workbook(tmp_path / "sims.xlsx", columns, rows, **kwargs)
    return read_sims_list(path)


def test_round_trip_preserves_values_columns_and_order(tmp_path):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h"})
            for d in (18, 24, 36)]
    sims = build(tmp_path, rows)
    dest = tmp_path / "chunk.xlsx"

    write_chunk_list(sims.path, [0, 2], dest, source_row_numbers=False)
    written = pd.read_excel(dest, sheet_name=0)

    assert list(written.columns) == list(sims.frame.columns)
    assert written["Output file"].tolist() == ["out_18h", "out_36h"]
    assert written["FSL"].tolist() == [670.42, 670.42]


def test_include_is_forced_to_the_exact_string_bryan_compares(tmp_path):
    """Main.py:64 is `sim_df['Include'] == 'yes'` - not lowered."""
    sims = build(tmp_path, [reservoir_row(Include="no"),
                            reservoir_row(Include="Yes")])
    dest = tmp_path / "chunk.xlsx"
    write_chunk_list(sims.path, [0, 1], dest)
    assert pd.read_excel(dest, sheet_name=0)["Include"].tolist() == ["yes", "yes"]


def test_sheet_name_is_preserved_and_other_sheets_dropped(tmp_path):
    from openpyxl import load_workbook

    sims = build(tmp_path, [reservoir_row()], columns=CALLIDE_COLUMNS,
                 sheet_name="SIMS", extra_sheets=("Definitions",))
    dest = tmp_path / "chunk.xlsx"
    write_chunk_list(sims.path, [0], dest)

    workbook = load_workbook(dest)
    assert workbook.worksheets[0].title == "SIMS"
    assert len(workbook.worksheets) == 1


def test_a_value_that_looks_like_a_formula_stays_a_literal(tmp_path):
    """The exact mechanism that damaged TFD_SimsList_LongList_02.xlsx.

    openpyxl turns any string starting with '=' into a formula with no cached
    value, and then pandas - and Bryan - read it as blank.
    """
    sims = build(tmp_path, [reservoir_row(Comment="=SUM(A1:A2)")])
    dest = tmp_path / "chunk.xlsx"
    write_chunk_list(sims.path, [0], dest)

    written = pd.read_excel(dest, sheet_name=0)
    assert written.loc[0, "Comment"] == "=SUM(A1:A2)"

    from core.simslist import audit_formulas
    assert not audit_formulas(dest).is_damaged


def test_blank_cells_do_not_become_the_string_nan(tmp_path):
    sims = build(tmp_path, [reservoir_row(Comment=None, **{"Log file": None})])
    dest = tmp_path / "chunk.xlsx"
    write_chunk_list(sims.path, [0], dest)
    written = pd.read_excel(dest, sheet_name=0)
    assert pd.isna(written.loc[0, "Comment"])


def test_provenance_columns_are_added(tmp_path):
    sims = build(tmp_path, [reservoir_row(), reservoir_row()])
    dest = tmp_path / "chunk.xlsx"
    groups = pd.Series({0: "group A", 1: "group B"})
    write_chunk_list(sims.path, [1], dest, group_values=groups)

    written = pd.read_excel(dest, sheet_name=0)
    assert written["Group"].tolist() == ["group B"]
    assert written["Source row"].tolist() == [3]   # frame index 1 -> Excel row 3


def test_log_file_override_is_written(tmp_path):
    sims = build(tmp_path, [reservoir_row(), reservoir_row()])
    dest = tmp_path / "chunk.xlsx"
    write_chunk_list(sims.path, [0, 1], dest,
                     log_file_overrides={1: r"sims_mc\log\run_second.log"})
    written = pd.read_excel(dest, sheet_name=0)
    assert written["Log file"].tolist() == [r"sims_mc\log\run.log",
                                            r"sims_mc\log\run_second.log"]


def test_verification_catches_a_corrupted_copy(tmp_path, monkeypatch):
    """The read-back check is a production guarantee, not a test-only nicety."""
    sims = build(tmp_path, [reservoir_row()])

    real_set = runwriter._set

    def corrupt(cell, value):
        real_set(cell, "WRONG" if value == 670.42 else value)

    monkeypatch.setattr(runwriter, "_set", corrupt)
    with pytest.raises(WriteError, match="FSL"):
        write_chunk_list(sims.path, [0], tmp_path / "chunk.xlsx")


def test_stale_row_index_is_rejected(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    with pytest.raises(WriteError, match="changed since"):
        write_chunk_list(sims.path, [99], tmp_path / "chunk.xlsx")


# --- the whole run folder -------------------------------------------------

def test_run_log_lands_inside_the_run_folder(tmp_path, project):
    """The load-bearing one.

    lib/RunLog.py:29 builds the log path from the raw ``simulation_list``
    string, resolved against the CWD. Because the UI launches with
    cwd=project_folder and writes a relative simulation_list, each chunk's log
    lands in its own run folder - which is what makes per-chunk run logs work
    without any locking.
    """
    config_path = project(TINAROO_COLUMNS,
                          [reservoir_row(**{"Output file": f"out_{d}h"}, Duration=d)
                           for d in (18, 24)])
    config = load_sims_config(config_path)
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, [0, 1], n_chunks=2, keep_groups_together=False)
    folder = runwriter.write_run(sims, config, plan, run_id="20260805-101010")

    for chunk in folder.chunks:
        payload = json.loads(chunk.config.read_text())
        predicted = run_log_path_for(payload["simulation_list"],
                                     Path(payload["project_folder"]))
        assert predicted == chunk.run_log
        assert predicted.parent == folder.folder


def test_chunk_config_makes_mainpy_resolution_a_no_op(tmp_path, project):
    config_path = project(TINAROO_COLUMNS, [reservoir_row()])
    config = load_sims_config(config_path)
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, [0], n_chunks=1)
    folder = runwriter.write_run(sims, config, plan, run_id="20260805-101010")

    payload = json.loads(folder.chunks[0].config.read_text())
    assert Path(payload["project_folder"]) == config.project_folder
    for value in payload["filepaths"].values():
        assert Path(value).is_absolute()
    # Main.py:35 joins the two; it must find the chunk workbook.
    joined = Path(payload["project_folder"]) / payload["simulation_list"]
    assert joined.is_file()


def test_launch_scripts_are_written(tmp_path, project):
    config_path = project(TINAROO_COLUMNS, [reservoir_row()])
    config = load_sims_config(config_path)
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, [0], n_chunks=1)
    folder = runwriter.write_run(sims, config, plan, run_id="20260805-101010")

    assert (folder.folder / "launch.bat").is_file()
    text = (folder.folder / "launch.sh").read_text()
    assert "chunk_01.json" in text
    assert str(config.project_folder) in text


def test_manifest_records_provenance(tmp_path, project):
    config_path = project(CALLIDE_COLUMNS, [reservoir_row()],
                          sheet_name="SIMS", extra_sheets=("Definitions",))
    config = load_sims_config(config_path)
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, [0], n_chunks=1)
    folder = runwriter.write_run(sims, config, plan, run_id="20260805-101010")

    manifest = json.loads(folder.manifest.read_text())
    assert manifest["source"]["sha256"] == sims.sha256
    assert manifest["source"]["dropped_sheets"] == ["Definitions"]
    assert manifest["environment"]["pandas"] == pd.__version__


def test_run_id_must_not_break_the_run_log_path():
    """RunLog.log_filepath is a substring replace on '.xlsx'."""
    with pytest.raises(ValueError):
        runwriter.make_run_id(tag="bad.xlsx.name")
    with pytest.raises(ValueError):
        from core.paths import safe_run_id
        safe_run_id("has spaces")
