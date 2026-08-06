"""Chunking, and the hazards that make a plan unsafe."""

from __future__ import annotations

import pandas as pd
import pytest

from conftest import (MONTE_CARLO_COLUMNS, TINAROO_COLUMNS, monte_carlo_row,
                      reservoir_row, write_workbook)
from core import runplan
from core.simslist import read_sims_list


def build(tmp_path, rows, columns=TINAROO_COLUMNS):
    return read_sims_list(write_workbook(tmp_path / "sims.xlsx", columns, rows))


def codes(plan):
    return {hazard.code for hazard in plan.hazards}


def test_rows_sharing_an_output_file_are_never_split(tmp_path):
    """The hard constraint - UrbsModel.__init__ rmtree's the working sub-folder.

    Reproduces the Callide shape, where the output name carries no duration.
    """
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (6, 12, 18, 24, 36, 48)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, list(range(6)), n_chunks=4,
                            keep_groups_together=False)

    assert len(plan.chunks) == 1, "one output name cannot be spread across chunks"
    assert "output-split" not in codes(plan)


def test_sharing_an_output_file_is_flagged_even_in_one_chunk(tmp_path):
    """Unsafe sequentially too: they all write the same results files."""
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (18, 24)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [0, 1], n_chunks=1)
    assert "output-collision" in codes(plan)
    assert plan.blocked


def test_distinct_outputs_split_across_chunks(tmp_path):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"opt{d}"})
            for d in (18, 24, 36, 48)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [0, 1, 2, 3], n_chunks=2,
                            keep_groups_together=False)

    assert len(plan.chunks) == 2
    assert sorted(plan.rows) == [0, 1, 2, 3]
    assert not plan.blocked


def test_groups_are_kept_together_by_default(tmp_path):
    rows = [reservoir_row(Duration=d, **{"Output file": f"caseA_{d}h",
                                         "Output suffix": "A"})
            for d in (18, 24)]
    rows += [reservoir_row(Duration=d, **{"Output file": f"caseB_{d}h",
                                          "Output suffix": "B"})
             for d in (18, 24)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [0, 1, 2, 3], n_chunks=2)

    assert len(plan.chunks) == 2
    for chunk in plan.chunks:
        suffixes = {sims.frame.loc[i, "Output suffix"] for i in chunk.rows}
        assert len(suffixes) == 1, "a group must not be split"


def test_rows_keep_sheet_order_within_a_chunk(tmp_path):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": "A"})
            for d in (18, 24, 36)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [2, 0, 1], n_chunks=1)
    assert plan.chunks[0].rows == (2, 0, 1)


def test_more_chunks_than_atoms_is_reported_not_failed(tmp_path):
    sims = build(tmp_path, [reservoir_row(**{"Output file": "only"})])
    plan = runplan.plan_run(sims, [0], n_chunks=5)
    assert len(plan.chunks) == 1
    assert "fewer-chunks" in codes(plan)


def test_empty_selection_blocks(tmp_path):
    sims = build(tmp_path, [reservoir_row()])
    plan = runplan.plan_run(sims, [], n_chunks=1)
    assert plan.blocked and "empty-selection" in codes(plan)


# --- cost and balance -----------------------------------------------------

def test_monte_carlo_is_estimated_far_above_reservoir_routing(tmp_path):
    mc = runplan.estimate_minutes(pd.Series(monte_carlo_row(Duration=24)))
    rr = runplan.estimate_minutes(pd.Series(reservoir_row(Duration=24)))
    assert mc > rr * 10


def test_analysis_only_rows_are_cheap(tmp_path):
    running = runplan.estimate_minutes(pd.Series(monte_carlo_row()))
    analysis = runplan.estimate_minutes(
        pd.Series(monte_carlo_row(**{"Run models": "no"})))
    assert analysis < running


def test_measured_history_overrides_the_heuristic():
    row = pd.Series(reservoir_row(**{"Output file": "known"}))
    assert runplan.estimate_minutes(row, {"known": 123.0}) == 123.0


def test_observed_minutes_reads_the_run_log(tmp_path, project):
    from core.config import load_sims_config

    config_path = project(TINAROO_COLUMNS, [reservoir_row()])
    config = load_sims_config(config_path)
    config.master_run_log.write_text(
        "Simulation,Start time,End time,Status\n"
        "out_18h,2026-08-05 10:00,2026-08-05 10:30,completed\n",
        encoding="utf-8")
    assert runplan.observed_minutes(config)["out_18h"] == pytest.approx(30.0)


def test_chunks_are_balanced_by_cost_not_row_count(tmp_path):
    """The expensive row goes on its own, the cheap ones share.

    Splitting five rows evenly by COUNT would put the 120 hour simulation with
    two short ones; by cost it goes alone. 120 is the floor here - a single
    simulation cannot be divided - so the test asserts the packing is optimal
    rather than merely even.
    """
    rows = [monte_carlo_row(Duration=120, **{"Output file": "big"})]
    rows += [monte_carlo_row(Duration=6, **{"Output file": f"small{i}"})
             for i in range(4)]
    sims = build(tmp_path, rows, columns=MONTE_CARLO_COLUMNS)
    plan = runplan.plan_run(sims, list(range(5)), n_chunks=2,
                            keep_groups_together=False)

    by_size = {chunk.size: chunk for chunk in plan.chunks}
    assert set(by_size) == {1, 4}, "the big row should be alone"
    assert by_size[1].rows == (0,)
    # Wall clock is the slowest chunk, and that is the lone big simulation.
    assert plan.estimated_minutes == by_size[1].estimated_minutes
    assert by_size[4].estimated_minutes < by_size[1].estimated_minutes


# --- log deduplication ----------------------------------------------------

def test_shared_log_paths_are_made_unique_across_chunks(tmp_path):
    """lib/LogFiles.resolve_duplicates only dedupes within one call.

    Two chunks each holding a row that names the same log would both open it
    with "w" and one run's log would be lost. Running Bryan's own dedupe once
    over the whole plan fixes it, and makes the UI's predicted paths exact.
    """
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}",
                                         "Log file": r"sims_mc\log\shared.log"})
            for d in (18, 24)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [0, 1], n_chunks=2,
                            keep_groups_together=False)

    resolved = list(plan.log_file_overrides.values())
    assert len(set(resolved)) == 2, f"logs still collide: {resolved}"
    assert "log-collision" not in codes(plan)
    assert plan.log_renames


def test_dedupe_matches_bryans_own_single_pass(tmp_path):
    """The renames must be the ones a single sequential run would produce."""
    from core.bryan import resolve_duplicate_logs

    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}",
                                         "Log file": r"sims_mc\log\shared.log"})
            for d in (18, 24, 36)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [0, 1, 2], n_chunks=3,
                            keep_groups_together=False)

    expected = resolve_duplicate_logs(sims.frame.loc[plan.rows])["Log file"]
    assert list(plan.log_file_overrides.values()) == list(expected)


# --- the duration dtype trap ---------------------------------------------

def test_dropping_the_fractional_duration_row_is_warned_about(tmp_path):
    """Live in CLD_RFSL_mc_sims_01.xlsx, which has a 4.5 hour duration.

    That one row makes the whole column float, so pandas gives 120.0 and
    URBSmodel.duration_string returns '120.0'. A chunk without it reads int and
    returns '120' - different storm filenames, and a simulation-period lookup
    that misses and silently doubles the duration (URBSmodel.py:287).
    """
    rows = [monte_carlo_row(Duration=d, **{"Output file": f"mc_{d}"})
            for d in (4.5, 24, 48)]
    sims = build(tmp_path, rows, columns=MONTE_CARLO_COLUMNS)
    assert sims.frame["Duration"].dtype == float

    plan = runplan.plan_run(sims, [1, 2], n_chunks=1, keep_groups_together=False)
    assert "duration-dtype" in codes(plan)
    assert not plan.blocked, "a warning, not a block"


def test_keeping_the_fractional_row_raises_no_warning(tmp_path):
    rows = [monte_carlo_row(Duration=d, **{"Output file": f"mc_{d}"})
            for d in (4.5, 24, 48)]
    sims = build(tmp_path, rows, columns=MONTE_CARLO_COLUMNS)
    plan = runplan.plan_run(sims, [0, 1, 2], n_chunks=1,
                            keep_groups_together=False)
    assert "duration-dtype" not in codes(plan)


def test_reservoir_routing_is_exempt_from_the_duration_warning(tmp_path):
    """Reservoir routing takes its durations from the stored results."""
    rows = [reservoir_row(Duration=d, **{"Output file": f"rr_{d}",
                                         "Output suffix": f"o{d}"})
            for d in (4.5, 24, 48)]
    sims = build(tmp_path, rows)
    plan = runplan.plan_run(sims, [1, 2], n_chunks=1,
                            keep_groups_together=False)
    assert "duration-dtype" not in codes(plan)


# --- parallel warnings ----------------------------------------------------

def test_mop_up_off_warns_when_running_in_parallel(tmp_path):
    rows = [monte_carlo_row(Duration=d, **{"Output file": f"mc_{d}",
                                           "Mop up files": "no"})
            for d in (18, 24)]
    sims = build(tmp_path, rows, columns=MONTE_CARLO_COLUMNS)
    plan = runplan.plan_run(sims, [0, 1], n_chunks=2,
                            keep_groups_together=False)
    assert "mop-up-off" in codes(plan)
    assert "memory" in codes(plan)


def test_no_parallel_warnings_for_a_single_chunk(tmp_path):
    rows = [monte_carlo_row(**{"Output file": "mc", "Mop up files": "no"})]
    sims = build(tmp_path, rows, columns=MONTE_CARLO_COLUMNS)
    plan = runplan.plan_run(sims, [0], n_chunks=1)
    assert not (codes(plan) & {"mop-up-off", "memory", "per-chunk-logs"})
