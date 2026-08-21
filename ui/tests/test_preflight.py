"""One triggering and one clean case per pre-flight issue code."""

from __future__ import annotations

from conftest import (MONTE_CARLO_COLUMNS, TINAROO_COLUMNS, monte_carlo_row,
                      reservoir_row)
from core import preflight
from core.config import load_sims_config
from core.simslist import read_sims_list


def setup(project, rows, columns=TINAROO_COLUMNS, **kwargs):
    config = load_sims_config(project(columns, rows, **kwargs))
    return config, read_sims_list(config.sims_list_path)


def codes(sims, config, rows):
    return {issue.code for issue in preflight.check(sims, config, rows)}


def test_a_healthy_selection_raises_nothing(project):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}"})
            for d in (18, 24)]
    config, sims = setup(project, rows)
    assert preflight.check(sims, config, [0, 1]) == []


def test_empty_selection_blocks(project):
    config, sims = setup(project, [reservoir_row()])
    assert codes(sims, config, []) == {"empty-selection"}


def test_missing_input_file_blocks(project):
    config, sims = setup(project, [reservoir_row(**{"Output file": "out"})],
                         make_inputs=False)
    issues = preflight.check(sims, config, [0])
    assert "missing-input" in {issue.code for issue in issues}
    assert any("SQ file" in issue.message for issue in issues)


def test_run_only_inputs_are_skipped_for_an_analysis_only_row(project):
    """Run models = no returns early at Simulator.__init__'s do_runs guard.

    The rating curve and the inflows are never opened, so their absence is not
    a problem. The Config file still is - ReservoirRouting.py:817 reads it to
    get the TPT parameters even when only re-analysing.
    """
    row = reservoir_row(**{"Output file": "out", "Run models": "no"})
    config, sims = setup(project, [row], make_inputs=False)

    reported = " ".join(issue.message for issue in preflight.check(sims, config, [0])
                        if issue.code == "missing-input")
    for column in ("SQ file", "ELS file", "Inflow", "Input MCDF"):
        assert column not in reported, f"{column} is not read when Run models = no"


def test_run_only_inputs_are_checked_when_the_row_does_run(project):
    row = reservoir_row(**{"Output file": "out", "Run models": "yes"})
    config, sims = setup(project, [row], make_inputs=False)

    reported = " ".join(issue.message for issue in preflight.check(sims, config, [0])
                        if issue.code == "missing-input")
    assert "SQ file" in reported and "Inflow" in reported


def test_unrecognised_method_blocks(project):
    config, sims = setup(project, [reservoir_row(Method="Ensemble")])
    assert "bad-method" in codes(sims, config, [0])


def test_missing_required_column_blocks(project):
    columns = [c for c in TINAROO_COLUMNS if c != "Results folder"]
    config, sims = setup(project, [reservoir_row()], columns=columns)
    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "missing-column"]
    assert any("Results folder" in issue.message for issue in issues)


def test_blank_required_value_blocks(project):
    config, sims = setup(project, [reservoir_row(**{"Results folder": None})])
    assert "blank-required" in codes(sims, config, [0])


def test_a_blank_config_file_is_allowed_for_reservoir_routing(project):
    """sim_list.md says to leave it blank for ensemble input."""
    config, sims = setup(project, [reservoir_row(**{"Config file": None,
                                                    "Output file": "out"})])
    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "blank-required" and "Config file" in i.message]
    assert not issues


def test_output_collision_blocks(project):
    """The Callide case - twenty-four rows share one Output file there."""
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (18, 24, 36)]
    config, sims = setup(project, rows)
    issues = [i for i in preflight.check(sims, config, [0, 1, 2])
              if i.code == "output-collision"]
    assert issues and issues[0].rows == (0, 1, 2)
    assert "18" in issues[0].message and "36" in issues[0].message


def test_one_row_of_a_colliding_set_is_fine(project):
    """Which is how the Callide list is run today: one Include = yes."""
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (18, 24)]
    config, sims = setup(project, rows)
    assert "output-collision" not in codes(sims, config, [0])


def test_uncached_formulas_block_the_affected_rows(project):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}"})
            for d in (18, 24)]
    formulas = {(0, "Output file"): '="x"&C2'}
    config, sims = setup(project, rows, formulas=formulas, cached=False)

    issues = [i for i in preflight.check(sims, config, [0]) if
              i.code == "uncached-formulas"]
    assert issues and issues[0].rows == (0,)
    assert "Ctrl+Alt+F9" in issues[0].fix_hint


def test_uncached_formulas_elsewhere_only_warn(project):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}"})
            for d in (18, 24)]
    formulas = {(0, "Comment"): '="note"'}
    config, sims = setup(project, rows, formulas=formulas, cached=False)

    issues = [i for i in preflight.check(sims, config, [1])
              if i.code.startswith("uncached")]
    assert issues and issues[0].severity == preflight.WARN
    assert not preflight.blocking(issues)


def test_include_case_is_warned_not_blocked(project):
    """Main.py:64 does not lower - 'Yes' silently skips in a batch run."""
    config, sims = setup(project, [reservoir_row(Include="Yes")])
    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "include-case"]
    assert issues and issues[0].severity == preflight.WARN
    assert "does not affect runs launched from here" in issues[0].fix_hint


def test_the_documented_replicate_column_name_blocks(project):
    """Manual/SubDocs/sim_list.md:24 says 'Replication file'.

    Simulator.__init__ reads 'Replicate file'. A list built from the manual
    raises KeyError on every monte carlo row with Run models set.
    """
    columns = [("Replication file" if c == "Replicate file" else c)
               for c in MONTE_CARLO_COLUMNS]
    row = monte_carlo_row(**{"Output file": "mc"})
    row["Replication file"] = row.pop("Replicate file")
    config, sims = setup(project, [row], columns=columns)

    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "replicate-file-alias"]
    assert issues and "The manual is wrong" in issues[0].fix_hint


def test_blocked_rows_are_collected(project):
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (18, 24)]
    config, sims = setup(project, rows)
    issues = preflight.check(sims, config, [0, 1])
    assert preflight.blocked_rows(issues) == {0, 1}


def test_the_inflow_is_required_by_a_volume_row_that_does_not_run(project):
    """'Analyse volumes' reads the hydrographs whatever 'Run models' says.

    ReservoirRoutingSimulator._ensure_inflows_loaded measures the volumes off
    the inflow file, so the usual analysis-only exemption does not apply to it.
    """
    columns = TINAROO_COLUMNS + ["Analyse volumes"]
    row = reservoir_row(**{"Output file": "out", "Run models": "no",
                           "Analyse volumes": "yes"})
    config, sims = setup(project, [row], columns=columns, make_inputs=False)

    reported = " ".join(issue.message for issue in preflight.check(sims, config, [0])
                        if issue.code == "missing-input")
    assert "Inflow" in reported
    assert "SQ file" not in reported, "the rating curve is still not read"


def test_a_volume_row_that_does_not_run_needs_no_inflow_when_switched_off(project):
    columns = TINAROO_COLUMNS + ["Analyse volumes"]
    row = reservoir_row(**{"Output file": "out", "Run models": "no",
                           "Analyse volumes": "no"})
    config, sims = setup(project, [row], columns=columns, make_inputs=False)

    reported = " ".join(issue.message for issue in preflight.check(sims, config, [0])
                        if issue.code == "missing-input")
    assert "Inflow" not in reported
