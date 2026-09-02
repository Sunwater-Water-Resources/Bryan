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


def test_a_monte_carlo_row_that_replicates_nothing_is_clean(project):
    """Replication is opt-in: blank means sample everything afresh.

    set_replicates and set_exclusions both take the blank line as 'none found',
    and the replication file is only opened once a key is recognised, so all
    three columns are legitimately empty on an ordinary row.
    """
    row = monte_carlo_row(**{"Output file": "mc"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    blanks = [i for i in preflight.check(sims, config, [0])
              if i.code == "blank-required"]
    assert not blanks


def test_a_replicate_key_without_the_file_blocks(project):
    """pd.read_csv(nan) - the file is opened as soon as a key is recognised."""
    row = monte_carlo_row(**{"Output file": "mc", "Replicates": "rz, tp"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "replicates-without-file"]
    assert issues and issues[0].rows == (0,)


def test_a_replicate_key_with_a_file_is_clean(project):
    row = monte_carlo_row(**{"Output file": "mc", "Replicates": "rz,tp",
                             "Replicate file": r"sims_mc\mcdf.csv"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    assert "replicates-without-file" not in codes(sims, config, [0])


def test_an_unrecognised_key_warns(project):
    """set_replicates skips it in silence, so the run just does not replicate."""
    row = monte_carlo_row(**{"Output file": "mc", "Exclusions": "ebf,dd50"})
    config, sims = setup(project, [row], columns=MONTE_CARLO_COLUMNS)
    issues = [i for i in preflight.check(sims, config, [0])
              if i.code == "unknown-key"]
    assert len(issues) == 1
    assert "'dd50'" in issues[0].message
    assert issues[0].severity == preflight.WARN


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


def test_a_mis_capitalised_volume_header_is_still_recognised(project):
    """lib/Volumes.py accepts it, so pre-flight has to see it too."""
    columns = TINAROO_COLUMNS + ["Analyse Volumes"]
    row = reservoir_row(**{"Output file": "out", "Run models": "no",
                           "Analyse Volumes": "yes"})
    config, sims = setup(project, [row], columns=columns, make_inputs=False)

    reported = " ".join(issue.message for issue in preflight.check(sims, config, [0])
                        if issue.code == "missing-input")
    assert "Inflow" in reported


# --- triage: running what can be run ---------------------------------------

def good_rows(durations=(18, 24, 36)):
    return [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}"})
            for d in durations]


def test_triage_keeps_a_clean_selection_whole(project):
    config, sims = setup(project, good_rows())
    result = preflight.triage(sims, config, [0, 1, 2])
    assert result.can_run and result.is_whole_selection
    assert result.runnable == (0, 1, 2)


def test_triage_drops_the_row_with_a_missing_input(project):
    """The case this exists for: four bad paths in ninety-two rows."""
    rows = good_rows()
    rows[1]["SQ file"] = r"reservoir\gone.sq"
    config, sims = setup(project, rows)
    (config.project_folder / "reservoir/gone.sq").unlink()

    result = preflight.triage(sims, config, [0, 1, 2])
    assert result.runnable == (0, 2)
    assert list(result.skipped) == [1]
    assert result.skipped[1].code == "missing-input"
    assert result.can_run and not result.is_whole_selection
    assert result.codes() == ["missing-input"]
    assert result.reasons()[0][:2] == (1, "missing-input")


def test_triage_leaves_the_kept_rows_actually_runnable(project):
    """Whatever it hands back has to pass the same checks, or the run stops
    after the user has already been told it would not."""
    rows = good_rows()
    rows[0]["SQ file"] = r"reservoir\gone.sq"
    config, sims = setup(project, rows)
    (config.project_folder / "reservoir/gone.sq").unlink()

    result = preflight.triage(sims, config, [0, 1, 2])
    assert preflight.blocking(preflight.check(sims, config, result.runnable)) == []
    assert preflight.blocking(result.issues) == []


def test_triage_drops_every_member_of_a_collision(project):
    """Two rows writing the same output are both dropped - which of them the
    user wants is not the UI's to guess."""
    rows = good_rows((18, 24))
    rows.append(reservoir_row(Duration=36, **{"Output file": "out_18h",
                                              "Output suffix": "o18"}))
    config, sims = setup(project, rows)

    result = preflight.triage(sims, config, [0, 1, 2])
    assert result.runnable == (1,)
    assert set(result.skipped) == {0, 2}
    assert result.skipped[0].code == "output-collision"


def test_triage_reports_when_nothing_is_left(project):
    rows = good_rows((18,))
    rows[0]["SQ file"] = r"reservoir\gone.sq"
    config, sims = setup(project, rows)
    (config.project_folder / "reservoir/gone.sq").unlink()

    result = preflight.triage(sims, config, [0])
    assert not result.can_run and result.runnable == ()
    assert result.remaining[0].code == "nothing-left"
    assert list(result.skipped) == [0]


def test_triage_cannot_skip_past_a_problem_that_names_no_row(project):
    """A missing column or the wrong 'Replicate file' spelling is a property of
    the workbook, so no amount of deselecting fixes it."""
    columns = [("Replication file" if c == "Replicate file" else c)
               for c in MONTE_CARLO_COLUMNS]
    row = monte_carlo_row(**{"Output file": "mc"})
    row["Replication file"] = row.pop("Replicate file")
    config, sims = setup(project, [row], columns=columns)

    result = preflight.triage(sims, config, [0])
    assert not result.can_run and not result.skipped
    assert {issue.code for issue in result.remaining} == {"replicate-file-alias"}


def test_triage_folds_in_the_planners_own_hazards(project):
    """A hazard that only appears once the rows are split into chunks would
    otherwise stop a run the user had just been told was fine."""
    from core import runplan

    rows = good_rows((18, 24))
    config, sims = setup(project, rows)

    def plan_for(selection):
        return runplan.plan_run(sims, selection, n_chunks=1)

    clean = preflight.triage(sims, config, [0, 1], plan_for=plan_for)
    assert clean.can_run and clean.is_whole_selection

    invented = runplan.Hazard(runplan.BLOCK, "made-up", "row 1 is doomed",
                              rows=(1,))

    def plan_with_hazard(selection):
        plan = runplan.plan_run(sims, selection, n_chunks=1)
        if 1 in selection:
            return runplan.RunPlan(chunks=plan.chunks,
                                   hazards=plan.hazards + (invented,))
        return plan

    result = preflight.triage(sims, config, [0, 1], plan_for=plan_with_hazard)
    assert result.runnable == (0,)
    assert result.skipped[1].code == "made-up"


def test_triage_of_an_empty_selection_says_so(project):
    config, sims = setup(project, good_rows((18,)))
    result = preflight.triage(sims, config, [])
    assert not result.can_run
    assert {issue.code for issue in result.remaining} == {"empty-selection"}

