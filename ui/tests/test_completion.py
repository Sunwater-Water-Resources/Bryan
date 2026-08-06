"""The four completion states, and the Run models = no inversion."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from conftest import (MONTE_CARLO_COLUMNS, TINAROO_COLUMNS, monte_carlo_row,
                      reservoir_row)
from core import completion
from core.config import load_sims_config
from core.simslist import read_sims_list


def setup(project, rows, columns=TINAROO_COLUMNS):
    config = load_sims_config(project(columns, rows))
    return config, read_sims_list(config.sims_list_path)


def touch(path, when=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text("x", encoding="utf-8")
    if when is not None:
        os.utime(path, (when, when))


def test_no_results_means_not_run(project):
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.NOT_RUN
    assert state.primary is None


def test_results_newer_than_inputs_are_up_to_date(project):
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    results = config.project_folder / "sims_mc/results/out_18h__mcdf.csv"
    touch(results, time.time() + 60)

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.UP_TO_DATE
    assert state.primary == results
    assert state.needs_confirm_to_overwrite


def test_an_input_newer_than_the_results_is_stale(project):
    """The state that is invisible today - an edited .sq, say."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    results = config.project_folder / "sims_mc/results/out_18h__mcdf.csv"
    touch(results, time.time() - 3600)
    touch(config.project_folder / "reservoir/dam.sq", time.time())

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.STALE
    assert state.newest_input.name == "dam.sq"
    assert not state.needs_confirm_to_overwrite, "re-running a stale row is the point"


def test_a_changed_global_config_makes_results_stale(project):
    """A new climate config invalidates every result in the project."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf.csv",
          time.time() - 3600)
    touch(config.filepaths["climate_config"], time.time())

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.STALE
    assert state.newest_input.name == "climate_config.json"


def test_a_parquet_results_file_counts(project):
    """ReservoirRouting._read_indexed falls back between the extensions."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    results = config.project_folder / "sims_mc/results/out_18h__mcdf.parquet"
    touch(results, time.time() + 60)
    assert completion.assess(sims.frame.loc[0], config).state == completion.UP_TO_DATE


def test_a_truncated_mcdf_is_incomplete(project):
    """An mcdf short of main x sub divisions - a test_runs run, or a crash."""
    row = monte_carlo_row(**{"Output file": "results/mc_24h",
                             "Config file": "mc_config.json"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    (config.project_folder / "mc_config.json").write_text(json.dumps({
        "scheme_config": {"number_of_main_divisions": 10,
                          "number_of_sub_divisions": 10}}), encoding="utf-8")

    results = config.project_folder / "results/mc_24h__mcdf.csv"
    results.parent.mkdir(parents=True, exist_ok=True)
    results.write_text("index,inflow\n" + "".join(f"{i},1\n" for i in range(7)),
                       encoding="utf-8")
    os.utime(results, (time.time() + 60,) * 2)

    expected = completion.expected_realisations(sims.frame.loc[0], config)
    assert expected == 100

    state = completion.assess(sims.frame.loc[0], config,
                              check_truncation=True, mc_expected_rows=expected)
    assert state.state == completion.INCOMPLETE
    assert state.row_count == 7 and state.expected_rows == 100


def test_a_full_mcdf_is_up_to_date(project):
    row = monte_carlo_row(**{"Output file": "results/mc_24h",
                             "Config file": "mc_config.json"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    (config.project_folder / "mc_config.json").write_text(json.dumps({
        "scheme_config": {"number_of_main_divisions": 2,
                          "number_of_sub_divisions": 3}}), encoding="utf-8")

    results = config.project_folder / "results/mc_24h__mcdf.csv"
    results.parent.mkdir(parents=True, exist_ok=True)
    results.write_text("index,inflow\n" + "".join(f"{i},1\n" for i in range(6)),
                       encoding="utf-8")
    os.utime(results, (time.time() + 60,) * 2)

    state = completion.assess(sims.frame.loc[0], config, check_truncation=True,
                              mc_expected_rows=6)
    assert state.state == completion.UP_TO_DATE


def test_expected_realisations_falls_back_to_the_top_level(project):
    """A config written for a routing run holds them at the top level.

    lib/ReservoirRouting.py looks in scheme_config first, then the top level -
    reading only the top level was a real bug (change_log, 31 July 2026).
    """
    row = monte_carlo_row(**{"Config file": "mc_config.json"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    (config.project_folder / "mc_config.json").write_text(json.dumps({
        "number_of_main_divisions": 4, "number_of_sub_divisions": 5}),
        encoding="utf-8")
    assert completion.expected_realisations(sims.frame.loc[0], config) == 20


# --- the inversion --------------------------------------------------------

def test_run_models_no_without_results_is_an_error_not_not_run(project):
    """lib/ReservoirRouting.py:790-792 - an analysis-only row needs results."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Run models": "no"})])
    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.NEEDS_PRIOR
    assert "only re-analyses" in state.detail


def test_run_models_no_with_results_is_fine(project):
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Run models": "no"})])
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf.csv",
          time.time() + 60)
    assert completion.assess(sims.frame.loc[0], config).state == completion.UP_TO_DATE


def test_preflight_blocks_an_analysis_row_with_nothing_to_analyse(project):
    from core import preflight

    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Run models": "no"})])
    codes = {issue.code for issue in preflight.check(sims, config, [0])}
    assert "needs-prior-results" in codes


# --- the run log is provenance only --------------------------------------

def test_run_log_history_is_read_but_not_used_to_decide(project):
    """Simulation is Output file with no duration, so it cannot identify a row.

    In CLD_RFSL_mc_sims_01.xlsx twenty-four rows share one name.
    """
    config, sims = setup(project, [reservoir_row(**{"Output file": "shared"}),
                                   reservoir_row(**{"Output file": "shared",
                                                    "Duration": 24})])
    config.master_run_log.write_text(
        "Simulation,Start time,End time,Computer,Status,Error\n"
        "shared,2026-08-05 10:00,2026-08-05 10:30,PC1,completed,\n",
        encoding="utf-8")

    history = completion.last_run_from_logs(config)
    assert history["shared"][1] == "completed"
    # ...but both rows are still 'not run', because no results exist.
    states = completion.assess_frame(sims.frame, config)
    assert {state.state for state in states.values()} == {completion.NOT_RUN}


def test_rerun_worthy_covers_the_states_a_bulk_action_should_select():
    assert set(completion.RERUN_WORTHY) == {
        completion.NOT_RUN, completion.STALE, completion.INCOMPLETE}
