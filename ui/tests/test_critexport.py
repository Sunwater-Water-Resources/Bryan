"""Building the critical duration export - the command, and what it will write.

Nothing here runs the export; test_critical_export.py does that against the
real script. This is the part that has to be right before anything is launched:
the wrong ``--result-type`` reads the wrong column, and a name with a path
separator in it lands the files somewhere nobody will look.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

from core import critexport, results
from test_results import level_curve, mc_row, write_quantiles

PYTHON = sys.executable


def sources(tmp_path, durations, key="level", column=None):
    column = column or key
    rows = []
    for duration in durations:
        name = f"CLD_mc_{duration:g}h_E010_GWL0p3"
        write_quantiles(Path(f"{tmp_path / name}_{key}.csv"), column,
                        level_curve(duration))
        rows.append(mc_row(tmp_path, name, duration))
    frame = pd.DataFrame(rows)
    return results.sources_for_rows(frame, tmp_path)


def test_the_command_names_every_duration_and_its_file(tmp_path):
    found = sources(tmp_path, [24, 48, 72])
    built = critexport.plan(found, ["level"], folder=tmp_path,
                            base_name="CLD_mc_E010_GWL0p3", python=PYTHON)

    assert built.can_run
    command = list(built.jobs[0].command)
    assert command[0] == PYTHON
    assert command[2].endswith("CriticalDurationAnalysis.py")
    assert command.count("--sim") == 3
    assert "--result-type" in command and command[command.index("--result-type") + 1] == "level"
    # whole durations stay whole - they become the labels in the exported table
    assert "24" in command and "48" in command


def test_a_volume_export_names_the_column_not_the_file_tag(tmp_path):
    """The file is tagged inflowVol24h; the column inside it is Vol24h."""
    found = sources(tmp_path, [24, 48], key="inflowVol24h", column="Vol24h")
    built = critexport.plan(found, ["inflowVol24h"], folder=tmp_path,
                            base_name="CLD", python=PYTHON)
    command = list(built.jobs[0].command)

    assert command[command.index("--result-type") + 1] == "Vol24h"
    assert built.jobs[0].csv.name == "CLD_inflowVol24h_critical.csv"


def test_each_result_type_is_its_own_run(tmp_path):
    """The critical duration for level is not the critical duration for inflow."""
    found = sources(tmp_path, [24, 48])
    found.update(sources(tmp_path, [24, 48], key="inflow"))
    built = critexport.plan(found, ["inflow", "level"], folder=tmp_path,
                            base_name="CLD", python=PYTHON)

    assert [job.key for job in built.jobs] == ["inflow", "level"]
    assert [job.csv.name for job in built.jobs] == ["CLD_inflow_critical.csv",
                                                    "CLD_level_critical.csv"]


def test_the_outputs_are_the_csv_and_the_plot(tmp_path):
    found = sources(tmp_path, [24, 48])
    job = critexport.plan(found, ["level"], folder=tmp_path / "out",
                          base_name="CLD", python=PYTHON).jobs[0]
    assert job.csv == tmp_path / "out" / "CLD_level_critical.csv"
    assert job.png == tmp_path / "out" / "CLD_level_critical_durations.png"


def test_files_already_there_are_reported_before_anything_runs(tmp_path):
    found = sources(tmp_path, [24, 48])
    (tmp_path / "CLD_level_critical.csv").write_text("an earlier export\n")
    built = critexport.plan(found, ["level"], folder=tmp_path,
                            base_name="CLD", python=PYTHON)

    assert [path.name for path in built.existing] == ["CLD_level_critical.csv"]
    assert built.can_run          # reported, not blocked - the page confirms


def test_dropped_aeps_and_no_plot_reach_the_command(tmp_path):
    found = sources(tmp_path, [24, 48])
    command = list(critexport.plan(found, ["level"], folder=tmp_path,
                                   base_name="CLD", python=PYTHON,
                                   drop_aeps=(2, 2000000), plot=False
                                   ).jobs[0].command)
    assert command.count("--drop-aep") == 2
    assert "2000000" in command
    assert "--no-plot" in command


@pytest.mark.parametrize("name", ["", "sub\\folder\\CLD", "sub/CLD", "a:b"])
def test_a_name_that_is_not_a_filename_is_refused(tmp_path, name):
    found = sources(tmp_path, [24, 48])
    built = critexport.plan(found, ["level"], folder=tmp_path,
                            base_name=name, python=PYTHON)
    assert not built.can_run
    assert any("filename" in problem for problem in built.problems)


def test_a_missing_interpreter_is_refused_before_launching(tmp_path):
    found = sources(tmp_path, [24, 48])
    built = critexport.plan(found, ["level"], folder=tmp_path,
                            base_name="CLD", python="")
    assert not built.can_run
    assert any("interpreter" in problem for problem in built.problems)

    built = critexport.plan(found, ["level"], folder=tmp_path, base_name="CLD",
                            python=str(tmp_path / "not_a_python"))
    assert any("not found" in problem for problem in built.problems)


def test_a_curve_with_no_duration_cannot_be_analysed(tmp_path):
    """A hand-added file may have no duration; a critical duration needs one."""
    write_quantiles(Path(f"{tmp_path / 'odd'}_level.csv"), "level", {2: 1.0})
    frame = pd.DataFrame([mc_row(tmp_path, "odd", None)])
    found = results.sources_for_rows(frame, tmp_path)
    built = critexport.plan(found, ["level"], folder=tmp_path,
                            base_name="CLD", python=PYTHON)
    assert not built.can_run
    assert any("no duration" in problem for problem in built.problems)


def test_the_default_name_is_the_group_without_its_duration(tmp_path):
    found = sources(tmp_path, [24, 48, 72])
    base = critexport.default_base_name(found["level"])
    assert base == "CLD_mc_E010_GWL0p3"
    assert critexport.output_name(base, "level") == "CLD_mc_E010_GWL0p3_level_critical"


def test_the_default_folder_is_beside_the_quantile_files(tmp_path):
    found = sources(tmp_path, [24, 48])
    assert critexport.default_folder(found["level"]) == tmp_path
