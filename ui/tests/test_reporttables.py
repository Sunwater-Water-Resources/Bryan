"""The Report page's core: the study file, the tables, and what goes on the clipboard.

What matters most, and is pinned here:

- a design flood table reads the inflow and outflow **at the level's critical
  duration**, not at their own - the report's footnote says so, and the old
  scripts behind Tables 33 and 34 got it the other way;
- an AEP where no duration spills is 0 outflow, not a missing row;
- a study file keeps its paths relative, so it survives being moved, and a run
  is named once so re-pointing the report at a re-run is one edit.
"""

from __future__ import annotations

import json
import math
import shutil

import pytest

from core import reporttables as rt
from core import study as studies
from core import wordtable as wt
from report_fixtures import GROUP, PMF_GROUP, build_study


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def design_spec(**overrides):
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Near-term", source={"run": "E099 RFSL", "group": GROUP},
                aeps=[5, 10, 100, 1000, 10000, 100000], pmp_aep=1900000)
    spec.update(overrides)
    return spec


def rows_of(table):
    return [row.cells for row in table.rows]


# -- the study file ------------------------------------------------------------

def test_a_study_round_trips_and_keeps_its_paths_relative(study):
    data = json.loads(study.path.read_text(encoding="utf-8"))
    assert data["bryan_study"] == 1
    assert data["runs"][0] == {"name": "E099 RFSL",
                               "sims_config": "runs/E099/CLD_RFSL_mc_sims_01.json"}
    again = studies.load_study(study.path)
    assert again.run_names() == ["E099 RFSL", "E099 PMF"]


def test_a_moved_study_still_opens_its_runs(study, tmp_path):
    study.put_table(design_spec())
    study.save()
    moved = tmp_path / "elsewhere"
    shutil.copytree(study.folder, moved)
    shutil.rmtree(study.folder)
    studies.forget_runs()
    rt.forget_cached()

    reopened = studies.load_study(moved / "bryan_study.json")
    table = rt.build(reopened, reopened.tables[0])
    assert not table.problems
    assert rows_of(table)[0][3] == "215.10"


def test_renaming_a_run_carries_every_table_with_it(study):
    study.put_table(design_spec())
    levels = rt.new_spec(rt.FLOOD_LEVELS)
    levels["sections"] = [{"heading": "", "rows": [
        {"label": "Near term", "run": "E099 RFSL", "group": GROUP}]}]
    study.put_table(levels)
    study.rename_run("E099 RFSL", "E100 RFSL")
    assert [s["run"] for t in study.tables for s in studies.table_sources(t)] \
        == ["E100 RFSL", "E100 RFSL"]


def test_removing_a_run_names_the_tables_left_without_it(study):
    study.put_table(design_spec(title="Table 26"))
    assert study.remove_run("E099 RFSL") == ["Table 26"]


def test_a_sims_config_is_not_mistaken_for_a_study(study):
    config = study.run_config_path("E099 RFSL")
    with pytest.raises(studies.StudyError, match="not a Bryan study file"):
        studies.load_study(config)


def test_a_newer_study_format_is_refused_rather_than_misread(tmp_path):
    path = tmp_path / "bryan_study.json"
    path.write_text(json.dumps({"bryan_study": 99, "runs": []}), encoding="utf-8")
    with pytest.raises(studies.StudyError, match="format 99"):
        studies.load_study(path)


def test_keys_this_version_does_not_know_are_kept(study):
    study.extra["figures"] = [{"id": "f1"}]
    study.save()
    assert studies.load_study(study.path).extra["figures"] == [{"id": "f1"}]


def test_table_ids_are_unique_and_duplicates_get_their_own(study):
    first = study.put_table(design_spec(title="Near-term"))
    second = study.put_table(design_spec(title="Near-term"))
    assert first["id"] == "near-term" and second["id"] == "near-term-2"


# -- design floods -------------------------------------------------------------

def test_flows_are_read_at_the_level_critical_duration(study):
    table = rt.build(study, design_spec())
    by_aep = {cells[0]: cells for cells in rows_of(table)}
    # 1 in 10: level critical at 72 h, so the 72 h inflow (950), not the 36 h
    # inflow that is the inflow envelope (1,100)
    assert by_aep["10"] == ["10", "950", "150", "216.00", "72"]
    # 1 in 1,000: level critical at 36 h
    assert by_aep["1,000"] == ["1,000", "5,000", "4,400", "217.30", "36"]


def test_an_aep_where_nothing_spills_is_zero_outflow(study):
    table = rt.build(study, design_spec())
    assert rows_of(table)[0] == ["5", "580", "0", "215.10", "72"]
    assert not table.problems


def test_the_pmp_row_carries_its_label(study):
    table = rt.build(study, design_spec())
    assert rows_of(table)[-1][0] == "1,900,000\n(PMPF)"
    assert rows_of(table)[-1][1:] == ["16,000", "13,100", "221.00", "36"]


def test_an_aep_the_run_did_not_produce_is_reported_not_invented(study):
    table = rt.build(study, design_spec(aeps=[5, 50]))
    assert rows_of(table)[1] == ["50", rt.NO_VALUE, rt.NO_VALUE, rt.NO_VALUE, rt.NO_VALUE]
    assert any("1 in 50" in problem for problem in table.problems)


def test_a_missing_group_says_so(study):
    table = rt.build(study, design_spec(source={"run": "E099 RFSL", "group": "nope"}))
    assert table.rows == []
    assert any("no group 'nope'" in problem for problem in table.problems)


def test_thousands_separators_can_be_left_off(study):
    table = rt.build(study, design_spec(thousands=False))
    assert rows_of(table)[3][1] == "5000"


def test_the_footnote_follows_the_table(study):
    assert rt.build(study, design_spec()).footnotes == [rt.INFLOW_FOOTNOTE]


# -- the multi-group kinds -----------------------------------------------------

def sections(*rows, heading="AEP for current RFSL of 215.5 m AHD"):
    return [{"heading": heading, "rows": list(rows)}]


def test_flood_levels_read_the_envelope_and_round_to_ten(study):
    spec = rt.new_spec(rt.FLOOD_LEVELS)
    spec["levels"] = [{"label": "DCF", "level": 219.13}, {"label": "Right", "level": 221.5}]
    spec["method"] = rt.ENVELOPE
    spec["sections"] = sections(
        {"label": "Baseline: Sunwater 2020*", "values": ["121,000", "", "24"]},
        {"label": "Near Term", "run": "E099 RFSL", "group": GROUP})
    table = rt.build(study, spec)
    assert table.rows[0].kind == wt.SECTION
    assert rows_of(table)[1] == ["Baseline: Sunwater 2020*", "121,000", "", "24"]

    near = rows_of(table)[2]
    envelope = rt.group_curves(study, "E099 RFSL", GROUP).envelope("level")
    expected = rt.EVENTS.aep_for_level(envelope, 219.13).aep
    assert near[1] == f"{round(expected / 10) * 10:,.0f}"
    assert near[2] == rt.NO_VALUE                     # 221.5 m is above the curve
    assert near[3] == "36"
    assert any("above the top of the level curve" in p for p in table.problems)


def test_peaks_at_one_aep_follow_the_level_unless_told_otherwise(study):
    spec = rt.new_spec(rt.PEAK_AT_AEP)
    spec["sections"] = sections({"label": "Near Term", "run": "E099 RFSL", "group": GROUP})
    assert rows_of(rt.build(study, spec))[1] == [
        "Near Term", "221.00", "16,000", "13,100", "36"]

    spec["aep"] = 100
    at_level = rows_of(rt.build(study, spec))[1]
    spec["flows_at"] = rt.AT_OWN
    own = rows_of(rt.build(study, spec))[1]
    assert at_level == ["Near Term", "216.80", "2,480", "2,250", "72"]
    assert own == ["Near Term", "216.80", "2,900", "2,250", "72"]


def test_the_pmf_is_the_highest_level_event_and_its_own_flows(study):
    spec = rt.new_spec(rt.ENSEMBLE_PEAK)
    spec["sections"] = sections({"label": "Near-Term", "run": "E099 PMF",
                                 "group": PMF_GROUP})
    assert rows_of(rt.build(study, spec))[1] == [
        "Near-Term", "221.38", "17,500", "16,310", "9"]
    spec["flows_at"] = rt.AT_OWN
    assert rows_of(rt.build(study, spec))[1][2] == "18,000"


def test_a_monte_carlo_group_offered_as_a_pmf_is_explained(study):
    spec = rt.new_spec(rt.ENSEMBLE_PEAK)
    spec["sections"] = sections({"label": "wrong", "run": "E099 RFSL", "group": GROUP})
    table = rt.build(study, spec)
    assert rows_of(table)[1][1] == rt.NO_VALUE
    assert table.problems


# -- formatting ----------------------------------------------------------------

@pytest.mark.parametrize("label, shown", [("72h", "72"), ("4.5h", "4.5"),
                                          (None, rt.NO_VALUE), ("custom", "custom")])
def test_durations_are_shown_in_hours_without_the_unit(label, shown):
    assert rt.fmt_duration(label) == shown


def test_aeps_and_levels_parse_from_the_page_text():
    assert rt.parse_aeps("5 10\n1,000 1_900_000 junk 1") == [5, 10, 1000, 1900000]
    assert rt.parse_levels("DCF = 219.13\nLeft overtopping=219.28\nno level") == [
        {"label": "DCF", "level": 219.13},
        {"label": "Left overtopping", "level": 219.28}]
    assert rt.parse_values("121,000 |  | 24") == ["121,000", "", "24"]


def test_a_group_without_results_is_reported_once_not_once_per_cell(study):
    spec = rt.new_spec(rt.FLOOD_LEVELS)
    spec["levels"] = [{"label": "a", "level": 216.0}, {"label": "b", "level": 217.0}]
    spec["sections"] = sections(
        {"label": "one", "run": "E099 RFSL", "group": "missing"},
        {"label": "two", "run": "E099 RFSL", "group": "missing"})
    problems = rt.build(study, spec).problems
    assert len(problems) == len(set(problems))


def test_a_long_title_makes_a_short_id():
    assert studies.slug("Table 29: Near-Term design hydrology results, reinstated FSL "
                        "216.1 m AHD (2021-2040)") == "table-29-near-term-design-hydrology"


# -- the level's AEP, and Table 1's extra rows ---------------------------------

def test_a_level_is_read_off_the_critical_duration_realisations(study):
    from core import ensemble
    from report_fixtures import MC_BASE, MC_SLOPE

    curves = rt.group_curves(study, "E099 RFSL", GROUP)
    aep, duration, problem = rt.aep_of_level(study, curves, 219.13)
    assert not problem and duration == "36h"
    # the fixture's realisations are level = 214 + 1.4 z, exactly
    assert aep == pytest.approx(ensemble.aep_of_variate((219.13 - MC_BASE) / MC_SLOPE),
                                rel=1e-3)
    envelope_aep, _, _ = rt.aep_of_level(study, curves, 219.13, rt.ENVELOPE)
    assert envelope_aep != pytest.approx(aep, rel=1e-3)     # the two methods differ


def test_without_a_database_the_level_falls_back_to_the_curve_and_says_so(study):
    for path in (study.folder / "runs" / "E099" / "sims_mc" / "results").glob("*__mcdf.csv"):
        path.unlink()
    curves = rt.group_curves(study, "E099 RFSL", GROUP)
    aep, _, problem = rt.aep_of_level(study, curves, 219.13)
    assert math.isfinite(aep) and "read off the curve" in problem


def test_table_1_carries_the_dcf_and_pmf_rows_in_aep_order(study):
    from core import ensemble

    section = ensemble.settings(study)
    section["adopted_aep"] = 8_000_000
    ensemble.store(study, section)
    spec = design_spec(dcf={"level": 219.13, "label": "DCF"},
                       pmf={"run": "E099 PMF", "group": PMF_GROUP, "label": "PMF"})
    table = rt.build(study, spec)
    labels = [row.cells[0] for row in table.rows]
    dcf = next(row for row in table.rows if row.cells[0].endswith("(DCF)"))
    assert dcf.bold and dcf.shaded and dcf.cells[1:4] == ["", "", "219.13"]
    order = [float(label.split()[0].replace(",", "")) for label in labels]
    assert order == sorted(order)                        # the DCF sits in AEP order
    assert labels.index(dcf.cells[0]) == labels.index("1,000") + 1
    assert table.rows[-1].cells == ["8,000,000 (PMF)", "17,500", "16,310", "221.38", "9"]


def test_a_pmf_row_without_an_adopted_aep_asks_for_one(study):
    spec = design_spec(pmf={"run": "E099 PMF", "group": PMF_GROUP})
    table = rt.build(study, spec)
    assert any("adopt one on the PMF page" in problem for problem in table.problems)


# -- representative events (Tables 35-36) --------------------------------------

def save_selection(study, targets):
    from core import events

    folder = study.folder / "runs" / "E099" / "sims_mc" / "results"
    events.save_targets(events.selection_path(folder, GROUP),
                        [events.Target.from_dict(target) for target in targets])


def representative_spec(**overrides):
    spec = rt.new_spec(rt.REPRESENTATIVE)
    spec["triggers"] = [{"label": "DCF", "level": 219.13}]
    spec["sections"] = [{"heading": "Near-Term", "run": "E099 RFSL", "group": GROUP,
                         "pmf": {"run": "E099 PMF", "group": PMF_GROUP}}]
    spec.update(overrides)
    return spec


def test_representative_events_come_from_the_saved_selection(study):
    save_selection(study, [
        {"kind": "level", "value": 219.13, "result_type": "level", "picked": 8096,
         "output_file": "CLD_mc_36h_GWL1p3_RFSL"},
        {"kind": "aep", "value": 100, "result_type": "level", "picked": 4980,
         "output_file": "CLD_mc_72h_GWL1p3_RFSL"},
        {"kind": "aep", "value": 1_900_000, "result_type": "level", "picked": 9728,
         "output_file": "CLD_mc_36h_GWL1p3_RFSL"},
        {"kind": "aep", "value": 1_000, "result_type": "level", "picked": None,
         "output_file": "CLD_mc_36h_GWL1p3_RFSL"}])
    table = rt.build(study, representative_spec())
    cells = [row.cells for row in table.rows]
    assert cells[0] == ["Near-Term"]
    assert cells[1] == ["100", "216.80", "", "4980", "72h"]         # AEP order
    assert cells[2][1:] == ["219.13", "DCF", "8096", "36h"]
    assert cells[3] == ["PMPF", "221.00", "", "9728", "36h"]
    assert cells[4] == ["PMF", "221.38", "", "2", "9h"]              # highest ensemble event
    assert any("no event picked" in problem for problem in table.problems)


def test_a_trigger_can_come_from_the_loading_s_own_comment(study):
    save_selection(study, [{"kind": "level", "value": 219.50, "result_type": "level",
                            "picked": 1, "comment": "Right embankment overtopping",
                            "output_file": "CLD_mc_36h_GWL1p3_RFSL"}])
    table = rt.build(study, representative_spec())
    assert table.rows[1].cells[2] == "Right embankment overtopping"


def test_a_group_without_a_saved_selection_says_where_to_make_one(study):
    table = rt.build(study, representative_spec())
    assert any("pick them on the Events page" in problem for problem in table.problems)


def test_renaming_a_run_carries_representative_sections_and_pmf_rows(study):
    study.put_table(representative_spec())
    study.put_table(design_spec(pmf={"run": "E099 PMF", "group": PMF_GROUP}))
    study.rename_run("E099 PMF", "E100 PMF")
    runs = [source["run"] for table in study.tables for source in studies.table_sources(table)]
    assert "E099 PMF" not in runs and runs.count("E100 PMF") == 2


# -- frequent levels (Table 37) --------------------------------------------------

def frequent_spec(**overrides):
    spec = rt.new_spec(rt.FREQUENT)
    spec["columns"] = ["GWL 1.3", "GWL 2.7"]
    spec["sections"] = [{"heading": "Current RFSL",
                         "groups": [{"run": "E099 RFSL", "group": GROUP}, {}]}]
    spec.update(overrides)
    return spec


def test_a_standard_aep_is_the_design_curve_and_1_ey_comes_off_the_realisations(study):
    from core import ensemble
    from report_fixtures import MC_BASE, MC_SLOPE
    table = rt.build(study, frequent_spec())
    one_in_two, one_ey = table.rows[1].cells, table.rows[2].cells
    assert one_in_two == ["1 in 2 AEP", "214.00", "", "72"]      # the envelope at 1 in 2
    # the fixture's realisations are level = 214 + 1.4 z, both durations alike
    z = ensemble.variate(1 / 1.582)
    assert one_ey[1] == f"{MC_BASE + MC_SLOPE * z:.2f}" and one_ey[2] == ""
    assert not table.problems


def test_a_duration_limit_leaves_the_other_runs_out(study):
    table = rt.build(study, frequent_spec(durations=[72]))
    assert table.rows[1].cells[-1] == "72" and table.rows[2].cells[-1] == "72"
    table = rt.build(study, frequent_spec(durations=[48]))           # nothing run at 48 h
    assert table.rows[1].cells[1] == rt.NO_VALUE and table.problems


def test_the_critical_duration_column_is_a_range_across_the_horizons():
    assert rt._duration_span(["24h", "72h", "24h"]) == f"24{rt.DASH}72"
    assert rt._duration_span(["24h", "24h"]) == "24"
    assert rt._duration_span([None]) == rt.NO_VALUE


def test_renaming_a_run_carries_a_grid_s_columns(study):
    study.put_table(frequent_spec())
    study.rename_run("E099 RFSL", "E100 RFSL")
    assert study.tables[0]["sections"][0]["groups"][0]["run"] == "E100 RFSL"


def test_a_save_waits_out_a_file_briefly_held_by_another_process(tmp_path, monkeypatch):
    """A virus scanner holding the study file for a moment must not fail the save."""
    import os as real_os
    from core import paths

    calls = []
    original = real_os.replace

    def flaky(source, target):
        calls.append(target)
        if len(calls) < 3:
            raise PermissionError(32, "being used by another process")
        return original(source, target)

    monkeypatch.setattr(paths.os, "replace", flaky)
    monkeypatch.setattr(paths, "REPLACE_PAUSE_SECONDS", 0.0)
    paths.atomic_write_json(tmp_path / "study.json", {"a": 1})
    assert len(calls) == 3 and json.loads((tmp_path / "study.json").read_text()) == {"a": 1}
