"""Finding the databases, and choosing the duration an event comes from.

The analysis itself is tested in test_representative_events.py against a
synthetic mcdf. What is pinned here is everything around it: which file belongs
to which row, that a level loading is read off the *envelope* rather than one
duration's curve, and that the suggested source is the duration that is
critical at the loading's AEP.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from conftest import MONTE_CARLO_COLUMNS, TINAROO_COLUMNS, write_workbook
from core import events
from test_representative_events import realisation
from test_results import level_curve, write_quantiles

DURATIONS = (12, 24, 48, 120)


def mcdf_frame():
    """The seven realisations of test_representative_events, as a file would hold."""
    rows = {
        0: realisation(0, 1000, 1000),
        1: realisation(1, 100, 1000, level=220.5),
        2: realisation(2, 1000, 100, level=214.0),
        3: realisation(3, 900, 1100, level=220.1,
                       embedded_bursts="Unfiltered embedded bursts: 2h burst exceeds by 20.0%",
                       subburst_2h=120.0),
        4: realisation(4, 900, 1100, level=220.1, preburst_p=0.88),
        5: realisation(5, 900, 1100, level=220.1, lake_z=2.0),
        6: realisation(6, 5_000_000, 1000, level=230.0),
    }
    return pd.DataFrame.from_dict(rows, orient="index")


def build(tmp_path, durations=DURATIONS, with_mcdf=True):
    """A monte carlo group: quantile tables per duration, and their databases."""
    folder = tmp_path / "sims_mc" / "results"
    rows = []
    for duration in durations:
        name = f"TFD_mc_{duration}h_GWL1p3"
        level = level_curve(duration)
        write_quantiles(folder / f"{name}_level.csv", "level", level)
        write_quantiles(folder / f"{name}_inflow.csv", "inflow",
                        {aep: value * 10 for aep, value in level.items()})
        if with_mcdf:
            mcdf_frame().to_csv(folder / f"{name}__mcdf.csv")
        rows.append({
            "Include": "yes", "Method": "monte carlo", "Duration": duration,
            "Run models": "yes", "Analyse results": "yes", "GWL": 1.3,
            "Output file": f"sims_mc\\results\\{name}",
            "Config file": "mc.json",
        })
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS, rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    return config


@pytest.fixture(autouse=True)
def clear_cache():
    events.forget_cached()
    yield
    events.forget_cached()


@pytest.fixture
def project(tmp_path):
    from state import Project
    return Project.open(build(tmp_path))


# -- finding the databases ---------------------------------------------------

def test_every_run_with_a_database_is_offered(project):
    sources = events.sources_for_rows(project)
    assert [source.label for source in sources] == ["12h", "24h", "48h", "120h"]
    assert all(source.path.name.endswith("__mcdf.csv") for source in sources)


def test_a_run_without_a_database_is_not_offered(tmp_path):
    from state import Project
    project = Project.open(build(tmp_path, with_mcdf=False))
    assert events.sources_for_rows(project) == []


def test_an_ensemble_row_is_not_offered(tmp_path):
    """It has no realisations to choose between - it ran every combination."""
    from state import Project
    config = build(tmp_path)
    project = Project.open(config)
    frame = project.frame
    frame.loc[0, "Method"] = "ensemble"
    assert "12h" not in [source.label for source in events.sources_for_rows(project)]


def test_a_reservoir_routing_row_resolves_to_its_suffixed_database(tmp_path):
    """<Output file>__mcdf<suffix>.csv - the routed database, not a sibling's."""
    from core.config import load_sims_config
    from core.simslist import read_sims_list
    from state import Project

    folder = tmp_path / "results"
    folder.mkdir(parents=True)
    mcdf_frame().to_csv(folder / "TFD_rr__mcdf_FR-4B.csv")
    rows = [{
        "Include": "yes", "Method": "reservoir routing", "Duration": 24,
        "Run models": "yes", "Analyse results": "yes",
        "Output file": "TFD_rr", "Output suffix": "FR-4B",
        "Results folder": "results", "Hydrographs folder": "results",
    }]
    write_workbook(tmp_path / "sims.xlsx", TINAROO_COLUMNS, rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))

    project = Project(config=load_sims_config(config),
                      sims=read_sims_list(tmp_path / "sims.xlsx"))
    sources = events.sources_for_rows(project)
    assert [source.path.name for source in sources] == ["TFD_rr__mcdf_FR-4B.csv"]


def test_two_rows_of_one_duration_get_distinct_labels(tmp_path):
    from state import Project
    config = build(tmp_path, durations=(24, 24, 48))
    project = Project.open(config)
    labels = [source.label for source in events.sources_for_rows(project)]
    assert len(set(labels)) == len(labels)


# -- the level frequency curve -----------------------------------------------

def test_a_level_target_is_read_off_the_envelope(project):
    """Not off one duration's curve - the envelope is the design quantile."""
    sources = events.sources_for_rows(project)
    curve = events.level_curve(project, sources)
    assert not curve.empty
    for duration in DURATIONS:
        single = pd.Series(level_curve(duration))
        assert (curve >= single.reindex(curve.index) - 1e-9).all()


# -- which duration the event comes from -------------------------------------

def test_the_suggested_source_is_the_critical_duration(project):
    """Level goes long at frequent AEPs and short on the rare tail.

    So the two ends must not suggest the same duration - if they did, the
    suggestion would be ignoring the physics the Results page exists to show.
    """
    sources = events.sources_for_rows(project)
    frequent, _ = events.critical_source(project, sources, "level", 10)
    rare, _ = events.critical_source(project, sources, "level", 100_000)
    assert frequent is not None and rare is not None
    assert frequent.duration > rare.duration


def test_the_suggestion_says_where_it_was_read(project):
    sources = events.sources_for_rows(project)
    _, note = events.critical_source(project, sources, "level", 1000)
    assert "critical for level at" in note


# -- evaluating a loading ----------------------------------------------------

def test_an_aep_loading_picks_the_neutral_event(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=3)
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert not outcome.problem
    assert outcome.picked_id == 0
    assert outcome.source.label == "24h"


def test_a_level_loading_becomes_an_aep_off_the_curve(project):
    sources = events.sources_for_rows(project)
    curve = events.level_curve(project, sources)
    level = float(curve.loc[1000])
    target = events.Target(kind="level", value=level, result_type="level",
                           source="24h")
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert outcome.aep == pytest.approx(1000, rel=0.05)


def test_a_level_above_the_curve_offers_the_highest_events(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="level", value=10_000, result_type="level",
                           source="24h", count=3)
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert not outcome.problem
    assert any("above the top of the curve" in note for note in outcome.notes)
    assert outcome.picked_id == 6            # level 230.0, the highest


def test_a_pick_overrides_the_ranking(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=7, picked=5)
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert outcome.picked_id == 5


def test_a_pick_that_is_no_longer_a_candidate_falls_back(project):
    """Tightening a filter must not leave the page showing a stale event."""
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=2, picked=6)
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert outcome.picked_id != 6
    assert outcome.picked_id in list(outcome.candidates.index)


def test_a_source_with_no_analysis_is_reported_not_raised(tmp_path):
    from state import Project
    config = build(tmp_path, durations=(24,))
    project = Project.open(config)
    source = events.sources_for_rows(project)[0]
    frame = mcdf_frame().drop(columns=["level_aep"])
    frame.to_csv(source.path)
    events.forget_cached()

    target = events.Target(kind="aep", value=1000, result_type="level", source="24h")
    outcome = events.evaluate(project, [source], target, events.Filters())
    assert "Analyse results" in outcome.problem
    assert outcome.candidates.empty


# -- the output --------------------------------------------------------------

def test_the_summary_names_the_hydrograph_column(project):
    """The whole point of the sim id: it is the column in the stored flows."""
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level", source="24h")
    outcome = events.evaluate(project, sources, target, events.Filters())
    row = events.summary_rows([outcome])[0]
    assert row["sim"] == 0
    assert row["hydrograph"] == "sim_00000"
    assert row["loading"] == "1 in 1,000"
    assert row["source"] == "24h"


def test_the_summary_carries_the_flags_of_a_flagged_pick(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=7, picked=4)
    row = events.summary_rows(
        [events.evaluate(project, sources, target, events.Filters())])[0]
    assert "pre-burst percentile" in row["flags"]


def test_the_table_rows_survive_json(project):
    """They are serialised to the browser - a numpy bool_ is a blank table."""
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=7)
    outcome = events.evaluate(project, sources, target, events.Filters())
    json.dumps(events.candidate_rows(outcome))
    json.dumps(events.summary_rows([outcome]))


def test_a_group_gets_its_own_selection_file():
    """Several groups commonly share one results folder - a GWL series does."""
    first = events.selection_path("/tmp/results", "TFD_mc_GWL1p3|exg")
    second = events.selection_path("/tmp/results", "TFD_mc_GWL1p7|exg")
    assert first != second
    assert "|" not in first.name


def test_targets_round_trip_through_the_selection_file(tmp_path):
    targets = [events.Target(kind="level", value=220.5, source="48h", picked=3)]
    path = events.selection_path(tmp_path, "group")
    events.save_targets(path, targets, {"result_type": "level"})

    loaded, settings = events.load_targets(path)
    assert loaded[0].value == 220.5
    assert loaded[0].picked == 3
    assert settings["result_type"] == "level"


# -- the cache ---------------------------------------------------------------

def test_the_database_is_read_once_and_re_read_when_it_changes(project):
    source = events.sources_for_rows(project)[0]
    first = events.prepared(source, "level")
    assert events.prepared(source, "level") is first

    frame = mcdf_frame()
    frame.loc[0, "level"] = 999.0
    frame.to_csv(source.path)
    assert events.prepared(source, "level") is not first


# -- what "closest" means ----------------------------------------------------

def test_ranking_on_the_result_reads_the_design_value_off_the_envelope(project):
    """A loading quoted as an AEP still has a level the event has to reach."""
    sources = events.sources_for_rows(project)
    curve = events.level_curve(project, sources)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=3)

    outcome = events.evaluate(project, sources, target, events.Filters(),
                              order=events.EVENTS.RESULT)
    assert outcome.target_value == pytest.approx(float(curve.loc[1000]), rel=1e-6)
    assert any("design value at 1 in" in note for note in outcome.notes)


def test_a_level_loading_ranks_on_the_level_it_names(project):
    """No curve read needed - the loading is already in metres."""
    sources = events.sources_for_rows(project)
    curve = events.level_curve(project, sources)
    level = float(curve.loc[1000])
    target = events.Target(kind="level", value=level, result_type="level",
                           source="24h", count=7)

    outcome = events.evaluate(project, sources, target, events.Filters(),
                              order=events.EVENTS.RESULT)
    assert outcome.target_value == pytest.approx(level)
    distances = (outcome.candidates["level"] - level).abs()
    assert distances.is_monotonic_increasing


def test_the_default_order_is_still_the_neutral_one(project):
    """The change is an option, not a new default."""
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=3)
    outcome = events.evaluate(project, sources, target, events.Filters())
    assert outcome.picked_id == 0
    assert outcome.target_value is None


def test_ranking_the_inflow_reads_the_inflow_curve(project):
    """The design value comes from the curve for the result being ranked."""
    sources = events.sources_for_rows(project)
    levels = events.level_curve(project, sources)
    inflows = events.envelope_curve(project, sources, "inflow")
    target = events.Target(kind="aep", value=1000, result_type="inflow",
                           source="24h", count=3)

    outcome = events.evaluate(project, sources, target, events.Filters(),
                              order=events.EVENTS.RESULT)
    assert outcome.target_value == pytest.approx(float(inflows.loc[1000]), rel=1e-6)
    assert outcome.target_value != pytest.approx(float(levels.loc[1000]), rel=1e-6)


def test_the_band_is_asked_for_in_millimetres_of_level(project):
    """The page asks in mm because that is what a level is argued about in."""
    assert events.band_units("level") == ("mm", 1000.0)
    assert events.band_in_result_units("level", 20) == pytest.approx(0.02)
    assert events.band_in_result_units("inflow", 20) == pytest.approx(20.0)
    assert events.band_in_result_units("level", None) == 0.0


def test_the_band_is_reported_against_the_loading(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=3)
    outcome = events.evaluate(project, sources, target, events.Filters(),
                              order=events.EVENTS.RESULT, band=0.02)
    assert any("within 20 mm" in note for note in outcome.notes)


def test_the_band_is_not_mentioned_when_ranking_on_neutrality(project):
    sources = events.sources_for_rows(project)
    target = events.Target(kind="aep", value=1000, result_type="level",
                           source="24h", count=3)
    outcome = events.evaluate(project, sources, target, events.Filters(),
                              band=0.02)
    assert not any("within" in note for note in outcome.notes)
