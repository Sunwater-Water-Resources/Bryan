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


def quantiles_of(config, output_file, suffix="FR-4C", kind="inflow"):
    """The path a routing row's own analysis writes - what says IT has run.

    Not the mcdf: lib/ReservoirRouting._output_base leaves the suffix off that,
    so every suffix over one Output file writes the same one.
    """
    return (config.project_folder / "sims_mc/results"
            / f"{output_file}__{kind}_quantiles_{suffix}.csv")


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
    results = quantiles_of(config, "out_18h")
    touch(results, time.time() + 60)

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.UP_TO_DATE
    assert state.primary == results
    assert state.needs_confirm_to_overwrite


def test_an_input_newer_than_the_results_is_stale(project):
    """The state that is invisible today - an edited .sq, say."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    touch(quantiles_of(config, "out_18h"), time.time() - 3600)
    touch(config.project_folder / "reservoir/dam.sq", time.time())

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.STALE
    assert state.newest_input.name == "dam.sq"
    assert not state.needs_confirm_to_overwrite, "re-running a stale row is the point"


def test_a_changed_global_config_makes_results_stale(project):
    """A new climate config invalidates every result in the project."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h"})])
    touch(quantiles_of(config, "out_18h"), time.time() - 3600)
    touch(config.filepaths["climate_config"], time.time())

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.STALE
    assert state.newest_input.name == "climate_config.json"


def test_a_parquet_results_file_counts(project):
    """ReservoirRouting._read_indexed falls back between the extensions.

    No Output suffix here, so the mcdf is this row's own - see
    test_a_suffixed_row_is_not_judged_by_the_shared_mcdf.
    """
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Output suffix": ""})])
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


def test_run_models_no_with_a_database_is_not_run_not_an_error(project):
    """A database to re-analyse is what lifts NEEDS_PRIOR - but re-analysing is
    exactly what has not happened yet, so the row reads as NOT_RUN until its own
    quantiles are there."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Run models": "no"})])
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf.csv",
          time.time() + 60)
    assert completion.assess(sims.frame.loc[0], config).state == completion.NOT_RUN

    touch(quantiles_of(config, "out_18h"), time.time() + 60)
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


def test_a_volume_row_lists_the_volume_table_it_writes(project):
    """The volume table is secondary - it says what the row writes, and never
    decides whether the row has run."""
    from core import outputs

    columns = TINAROO_COLUMNS + ["Analyse volumes"]
    row = reservoir_row(**{"Output file": "out", "Output suffix": "opt",
                           "Analyse volumes": "yes"})
    config, sims = setup(project, [row], columns=columns)

    written = outputs.outputs_for(sims.frame.loc[0], config.project_folder)
    names = [path.name for path in written.secondary]
    assert "out__inflow_volumes_opt.csv" in names
    assert not any("inflow_volumes" in path.name for path in written.primary)

    plain = reservoir_row(**{"Output file": "out", "Output suffix": "opt"})
    config, sims = setup(project, [plain], columns=columns)
    written = outputs.outputs_for(sims.frame.loc[0], config.project_folder)
    assert not any("inflow_volumes" in path.name for path in written.secondary)


# --- the suffix: what an Output suffix does and does not name --------------

def test_the_mcdf_carries_the_suffix(project):
    from core import outputs

    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Output suffix": "FR-4C"})])
    written = outputs.outputs_for(sims.frame.loc[0], config.project_folder)
    assert written.primary[0].name == "out_18h__mcdf_FR-4C.csv"
    assert written.databases[0].name == "out_18h__mcdf_FR-4C.csv"
    assert completion.assess(sims.frame.loc[0], config).state == completion.NOT_RUN

    touch(config.project_folder / "sims_mc/results/out_18h__mcdf_FR-4C.csv",
          time.time() + 60)
    assert completion.assess(sims.frame.loc[0], config).state == completion.UP_TO_DATE


def test_the_pre_suffix_mcdf_is_readable_but_proves_nothing(project):
    """The six databases thirty-six Tinaroo rows wrote between them.

    lib/ReservoirRouting._ensure_mcdf_loaded still falls back to that name, so
    it stays a database candidate - but it belongs to whichever suffix ran last,
    so it never says THIS row ran.
    """
    from core import outputs

    rows = [reservoir_row(**{"Output file": "out_18h", "Output suffix": name})
            for name in ("FR-4B", "FR-4C")]
    config, sims = setup(project, rows)
    legacy = config.project_folder / "sims_mc/results/out_18h__mcdf.csv"
    touch(legacy, time.time() + 60)

    written = outputs.outputs_for(sims.frame.loc[0], config.project_folder)
    assert [path.name for path in written.shared] == ["out_18h__mcdf.csv"]
    assert legacy not in written.primary
    assert legacy in written.databases, "an analysis-only row can still read it"

    states = completion.assess_frame(sims.frame, config)
    assert {state.state for state in states.values()} == {completion.NOT_RUN}

    # FR-4B has since been re-routed under the new naming; FR-4C has not.
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf_FR-4B.csv",
          time.time() + 60)
    states = completion.assess_frame(sims.frame, config)
    assert states[0].state == completion.UP_TO_DATE
    assert states[1].state == completion.NOT_RUN


def test_a_suffixed_row_dates_its_staleness_from_its_own_results(project):
    """A sibling's fresh database must not hide an input that changed after
    THIS row ran."""
    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Output suffix": "FR-4C"})])
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf_FR-4C.csv",
          time.time() - 3600)
    touch(config.project_folder / "reservoir/dam.sq", time.time() + 30)
    # the sibling that ran since, leaving a brand new database behind
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf_FR-4D.csv",
          time.time() + 120)
    touch(config.project_folder / "sims_mc/results/out_18h__mcdf.csv",
          time.time() + 120)

    state = completion.assess(sims.frame.loc[0], config)
    assert state.state == completion.STALE
    assert state.newest_input.name == "dam.sq"


def test_without_a_suffix_nothing_is_shared(project):
    from core import outputs

    config, sims = setup(project, [reservoir_row(**{"Output file": "out_18h",
                                                    "Output suffix": ""})])
    written = outputs.outputs_for(sims.frame.loc[0], config.project_folder)
    assert written.shared == ()
    assert written.primary[0].name == "out_18h__mcdf.csv"


def test_stored_hydrographs_identify_a_row_that_does_not_analyse(project):
    """Run models = yes, Analyse results = no, and an ensemble input, so no
    Monte Carlo database is written either: the hydrograph set is all there is."""
    config, sims = setup(project, [reservoir_row(**{
        "Output file": "out_18h", "Output suffix": "FR-4C",
        "Analyse results": "no", "Store hydrographs": "yes"})])

    assert completion.assess(sims.frame.loc[0], config).state == completion.NOT_RUN
    touch(config.project_folder / "sims_mc/hydrographs/out_18h_levels_FR-4C.csv",
          time.time() + 60)
    assert completion.assess(sims.frame.loc[0], config).state == completion.UP_TO_DATE


def test_store_hydrographs_reads_as_write_hydrographs_does():
    """ReservoirRouting._write_hydrographs guards with pd.notna, so an empty
    string is a value the user typed - it switches the hydrographs off - while
    an absent or NaN cell leaves the default 'yes' alone.

    A blank *cell* comes back from pandas as NaN, so this distinction only
    shows up where the value came from a formula. Tested on the Series
    directly for that reason.
    """
    import pandas as pd

    from core.columns import stores_hydrographs

    def row(value):
        return pd.Series({"Output file": "out_18h", "Store hydrographs": value})

    assert stores_hydrographs(row("yes"))
    assert stores_hydrographs(row(float("nan"))), "NaN leaves the default alone"
    assert stores_hydrographs(pd.Series({"Output file": "out_18h"})), "absent too"
    assert not stores_hydrographs(row("")), "an empty string is a value"
    assert not stores_hydrographs(row("no"))


def test_truncation_counts_the_database_not_a_quantile_table(project):
    """A quantile table holds one row per standard AEP, so counting the wrong
    file would call every routed row incomplete."""
    row = reservoir_row(**{"Output file": "out_18h", "Output suffix": "FR-4C",
                           "Config file": "mc_config.json"})
    config, sims = setup(project, [row])
    (config.project_folder / "mc_config.json").write_text(json.dumps({
        "scheme_config": {"number_of_main_divisions": 2,
                          "number_of_sub_divisions": 3}}), encoding="utf-8")

    mcdf = config.project_folder / "sims_mc/results/out_18h__mcdf_FR-4C.csv"
    mcdf.parent.mkdir(parents=True, exist_ok=True)
    mcdf.write_text("index,inflow\n" + "".join(f"{i},1\n" for i in range(6)),
                    encoding="utf-8")
    touch(quantiles_of(config, "out_18h", "FR-4C"), time.time() + 60)
    os.utime(mcdf, (time.time() + 60,) * 2)

    state = completion.assess(sims.frame.loc[0], config, check_truncation=True,
                              mc_expected_rows=6)
    assert state.state == completion.UP_TO_DATE

    mcdf.write_text("index,inflow\n0,1\n", encoding="utf-8")
    os.utime(mcdf, (time.time() + 60,) * 2)
    state = completion.assess(sims.frame.loc[0], config, check_truncation=True,
                              mc_expected_rows=6)
    assert state.state == completion.INCOMPLETE and state.row_count == 1
