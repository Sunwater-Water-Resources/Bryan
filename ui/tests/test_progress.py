"""Deriving progress, and pinning the Bryan output it is derived from."""

from __future__ import annotations

import pytest

from core import progress
from core.bryan import BRYAN_ROOT

MAIN_PY = BRYAN_ROOT / "Main.py"
RUNLOG_PY = BRYAN_ROOT / "lib" / "RunLog.py"


class FakeChunk:
    def __init__(self, index, rows, run_log, console_log):
        self.index = index
        self.rows = rows
        self.run_log = run_log
        self.console_log = console_log


def banner(position, total, name, method="ensemble", config="cfg.json"):
    return (f'\n{"=" * 90}\n'
            f'Simulation {position} of {total}: {name}\n'
            f'  method: {method}    method config: {config}\n'
            f'{"=" * 90}\n')


# --- the contract with Bryan ---------------------------------------------

def test_the_banner_format_still_matches_mainpy():
    """If Main.py's wording changes, fail here rather than silently show 0/N.

    Main.py:84 builds this string; core/progress.BANNER_RE parses it.
    """
    source = MAIN_PY.read_text(encoding="utf-8")
    assert "f'Simulation {position} of {total}: {sim[\"Output file\"]}'" in source
    rendered = "Simulation 3 of 12: TFD_mc_18h"
    match = progress.BANNER_RE.search(rendered)
    assert match and match.group("position") == "3" and match.group("total") == "12"


def test_the_summary_format_still_matches_runlog():
    """RunLog.print_summary:97,100 - what the UI reads as 'the chunk ended'."""
    source = RUNLOG_PY.read_text(encoding="utf-8")
    assert "All {total} simulations completed." in source
    assert "simulations did not complete:" in source
    assert progress.ALL_COMPLETED_RE.search("All 12 simulations completed.")
    assert progress.SOME_FAILED_RE.search("2 of 12 simulations did not complete:")


def test_run_log_columns_match_bryans():
    from core.bryan import run_log_columns

    columns = run_log_columns()
    for name in ("Simulation", "Start time", "End time", "Status", "Error"):
        assert name in columns


# --- parsing --------------------------------------------------------------

def test_parse_console_finds_the_running_simulation():
    text = banner(1, 3, "a") + banner(2, 3, "b")
    parsed = progress.parse_console(text)
    assert parsed["last"]["position"] == 2
    assert parsed["last"]["name"] == "b"
    assert parsed["total"] == 3
    assert not parsed["finished"]


def test_parse_console_detects_a_clean_finish():
    parsed = progress.parse_console(banner(1, 1, "a") + "\nAll 1 simulations completed.\n")
    assert parsed["all_completed"] and parsed["finished"]


def test_parse_console_detects_failures():
    text = banner(1, 2, "a") + "\n2 of 2 simulations did not complete:\n"
    parsed = progress.parse_console(text)
    assert parsed["finished"] and not parsed["all_completed"]
    assert "did not complete" in parsed["summary"]


# --- folding it together --------------------------------------------------

def test_run_log_rows_are_matched_positionally_not_by_name(tmp_path):
    """The Callide trap: many rows share one Output file.

    RunLog records Simulation = Output file (lib/RunLog.py:36) with no
    duration, so the Nth appended row is the Nth simulation and nothing else
    is reliable.
    """
    run_log = tmp_path / "chunk_01_log.csv"
    run_log.write_text(
        "Simulation,Start time,End time,Status,Error\n"
        "shared,2026-08-05 10:00,2026-08-05 10:10,completed,\n"
        "shared,2026-08-05 10:10,2026-08-05 10:20,FAILED,RuntimeError: bad\n",
        encoding="utf-8")
    console = tmp_path / "chunk_01_console.txt"
    console.write_text(banner(1, 3, "shared") + banner(2, 3, "shared")
                       + banner(3, 3, "shared"), encoding="utf-8")

    chunk = FakeChunk(1, (0, 1, 2), run_log, console)
    result = progress.read_chunk_progress(chunk, names=["shared"] * 3)

    assert [sim.status for sim in result.sims] == [
        "completed", "FAILED", progress.RUNNING]
    assert result.completed == 2
    assert result.current.position == 3
    assert result.failures[0].error.startswith("RuntimeError")


def test_a_first_long_simulation_shows_as_running_not_zero(tmp_path):
    """The run log only appends when a simulation finishes."""
    run_log = tmp_path / "log.csv"
    run_log.write_text("Simulation,Status\n", encoding="utf-8")
    console = tmp_path / "console.txt"
    console.write_text(banner(1, 5, "big"), encoding="utf-8")

    result = progress.read_chunk_progress(
        FakeChunk(1, tuple(range(5)), run_log, console))
    assert result.current is not None and result.current.position == 1
    assert result.completed == 0 and result.total == 5


def test_a_dead_process_marks_the_in_flight_simulation_cancelled(tmp_path):
    console = tmp_path / "console.txt"
    console.write_text(banner(1, 2, "a"), encoding="utf-8")
    result = progress.read_chunk_progress(
        FakeChunk(1, (0, 1), tmp_path / "missing.csv", console), alive=False)
    assert result.sims[0].status == progress.CANCELLED


def test_per_sim_log_size_is_reported(tmp_path):
    log = tmp_path / "sim.log"
    log.write_text("some output\n", encoding="utf-8")
    console = tmp_path / "console.txt"
    console.write_text(banner(1, 1, "a"), encoding="utf-8")

    result = progress.read_chunk_progress(
        FakeChunk(1, (0,), tmp_path / "missing.csv", console),
        log_paths=[log])
    assert result.sims[0].log_bytes > 0
    assert progress.stalled_for(result.sims[0], now=result.sims[0].log_mtime + 30) == 30


def test_missing_files_do_not_raise(tmp_path):
    result = progress.read_chunk_progress(
        FakeChunk(1, (0,), tmp_path / "nope.csv", tmp_path / "nope.txt"))
    assert result.total == 1 and result.sims[0].status == progress.PENDING


# --- explanations ---------------------------------------------------------

@pytest.mark.parametrize("code,expected", [
    (0, "completed"), (1, "at least one simulation failed"),
    (None, "still running"), (-9, "killed by signal"),
])
def test_returncode_explanations(code, expected):
    assert expected in progress.explain_returncode(code)


def test_eoferror_is_explained_in_terms_of_the_input_calls():
    """stdin=DEVNULL turns three hangs into three fast failures."""
    text = progress.explain_error("EOFError: EOF when reading a line")
    assert "open in Excel" in text and "URBS executable" in text


def test_missing_file_error_is_explained():
    assert "project folder" in progress.explain_error(
        "FileNotFoundError: no such file")


def test_unknown_errors_get_no_invented_explanation():
    assert progress.explain_error("ValueError: something") == ""


@pytest.mark.parametrize("seconds,expected", [
    (5, "5s"), (95, "1m 35s"), (3700, "1h 01m")])
def test_elapsed_text(seconds, expected):
    assert progress.elapsed_text(seconds) == expected
