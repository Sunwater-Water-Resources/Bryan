"""Reading the regional runs back: where they are, and each one's status from what is on disk
and what the page knows - never from a clock."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from core import urbsresults as ur

RUN = "e00042_9h_GWL1p3_RFSL"
STORM = RUN + ".9"


@pytest.fixture
def found(tmp_path):
    (tmp_path / "storms").mkdir()
    (tmp_path / "results").mkdir()
    return ur.Folders(storms=tmp_path / "storms", results=tmp_path / "results")


def _touch(path: Path) -> None:
    path.write_text("")


def test_the_folders_are_the_ones_the_generator_ran_in(tmp_path):
    config = tmp_path / "Regional_E001" / "urbs_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"model_folder": "model", "storms_folder": "storms",
                                  "results_folder": "results"}))
    found = ur.folders(config)
    assert found.storms == config.parent / "storms" and found.results == config.parent / "results"
    # no results_folder: beside the model, as before the key was read
    config.write_text(json.dumps({"model_folder": "model", "storms_folder": "storms"}))
    assert ur.folders(config).results == config.parent / "model"


def test_a_config_that_cannot_be_read_says_why(tmp_path):
    with pytest.raises(ur.ModelConfigError, match="could not be read"):
        ur.folders(tmp_path / "missing.json")
    config = tmp_path / "c.json"
    config.write_text(json.dumps({"model_folder": "model"}))
    with pytest.raises(ur.ModelConfigError, match="storms_folder"):
        ur.folders(config)


def test_results_on_disk_is_finished_whatever_else_is_true(found):
    _touch(found.results / f"{RUN}.q")
    _touch(found.results / f"_run_{RUN}.bat")
    assert ur.run_status(found, RUN, STORM, generator_running=True, ran_urbs=True) == ur.FINISHED


def test_while_the_generator_runs_the_event_with_a_batch_file_is_running(found):
    assert ur.run_status(found, RUN, STORM, generator_running=True, ran_urbs=True) == ur.WAITING
    _touch(found.storms / STORM)
    _touch(found.results / f"_run_{RUN}.bat")
    assert ur.run_status(found, RUN, STORM, generator_running=True, ran_urbs=True) == ur.RUNNING
    # writing only: a batch file is not URBS running
    assert ur.run_status(found, RUN, STORM, generator_running=True, ran_urbs=False) == ur.WAITING


def test_after_it_stops_a_batch_file_left_is_failed_only_if_urbs_was_run(found):
    assert ur.run_status(found, RUN, STORM) == ur.NOT_WRITTEN
    _touch(found.storms / STORM)
    assert ur.run_status(found, RUN, STORM, ran_urbs=True) == ur.NOT_RUN
    _touch(found.results / f"_run_{RUN}.bat")
    assert ur.run_status(found, RUN, STORM, ran_urbs=True) == ur.FAILED
    assert ur.run_status(found, RUN, STORM, ran_urbs=False) == ur.NOT_RUN
    # a restarted launcher does not know which, and says so rather than guess
    assert ur.run_status(found, RUN, STORM, ran_urbs=None) == ur.UNKNOWN


def test_the_runs_table_carries_what_went_into_each_storm(found, tmp_path):
    selection = tmp_path / "PMF_representative_events.json"
    (tmp_path / "PMF_representative_events_downstream_storms.csv").write_text(
        "filename,initial_loss,continuing_loss,mean_depth_mm,temporal_pattern,embedded_bursts,"
        "duration_dip,run_name,antecedent\n"
        f"{STORM},0.0,1.8,366.0,as sampled,No embedded bursts,,{RUN},"
        "initial_lake_level 129041 ML from the event ADV; kroombit_lake_level full supply\n")
    _touch(found.results / f"{RUN}_i.csv")
    plan = SimpleNamespace(storms=[SimpleNamespace(filename=STORM, loading="PMF",
                                                   event_label="event 42")])
    row, = ur.runs_rows(plan, found, ur.storms_record(selection))
    assert (row["status"], row["CL (mm/h)"], row["mean depth (mm)"]) == ("finished", "1.8", "366.0")
    assert "kroombit_lake_level full supply" in row["starting storage"]
    assert row["flags"] == ""
    assert ur.counts([row]) == "1 finished"



# -- the dams, read off a real regional run (trimmed from Regional_E001 005044_72h_GWL1p3_RFSL.h)

DATA = Path(__file__).parent / "data" / "regional"
REAL = "005044_72h_GWL1p3_RFSL"
# its _i.csv, trimmed to seven of the 34 print locations and every 60th step; inline because
# the repository keeps no csv files
FLOWS_CSV = """Time_hrs,STAGCK_WEIR,MALAKOFF,RAINBOW_CK,96_KM,CLD_INFLOW,GERRARD,GOOVIGEN
0.000,0.00,0.00,0.00,0.00,0.00,0.00,0.00
14.750,25.85,33.74,65.83,318.24,378.83,290.20,712.31
29.750,77.77,99.88,204.72,713.64,757.48,576.03,1540.23
44.750,196.53,253.51,514.01,2280.92,2937.97,2574.65,11423.20
59.750,69.60,89.43,180.57,750.37,924.44,837.38,4870.54
74.750,18.76,24.26,46.56,268.88,383.62,342.77,2186.78
89.750,0.76,0.96,1.90,13.02,21.90,34.05,383.40
104.750,0.11,0.14,0.27,1.90,3.24,6.64,175.84
119.750,0.03,0.03,0.07,0.48,0.82,1.96,18.44
134.750,0.01,0.01,0.02,0.16,0.28,0.74,7.29
144.000,0.01,0.01,0.01,0.09,0.16,0.44,4.47
"""


@pytest.fixture
def flows_folder(tmp_path):
    (tmp_path / f"{REAL}_i.csv").write_text(FLOWS_CSV)
    return tmp_path
LEVELS = {"CALLIDE": {"DCL": 219.13, "PMF": {"1.3": 221.50, "2.7": 221.86}},
          "KROOMBIT": {"DCL": 270.70}}


@pytest.fixture(autouse=True)
def fresh():
    ur.forget_cached()
    yield
    ur.forget_cached()


def test_the_level_file_is_read_by_its_eight_character_names():
    frame = ur.read_levels(DATA, REAL)
    assert {"CALLIDE", "KROOMBIT", "CLD_OUTF", "GOOVIGEN"} <= set(frame.columns)
    assert list(frame.index) == [0.0, 15.0, 30.0, 45.0, 60.0, 75.0]
    assert frame.loc[45.0, "KROOMBIT"] == 268.03
    assert ur.read_levels(DATA, "no_such_run") is None


def test_a_name_that_appears_twice_keeps_both(tmp_path):
    text = (DATA / f"{REAL}.h").read_text().splitlines()
    text[8] = text[8].replace("LINKES  ", "CALLIDE ")
    (tmp_path / "twice.h").write_text("\n".join(text))
    frame = ur.read_levels(tmp_path, "twice")
    assert "CALLIDE" in frame and "CALLIDE_1" in frame


def test_each_dam_is_held_against_its_levels_at_the_runs_warming_level():
    frame = ur.read_levels(DATA, REAL)
    callide = ur.dam_peak(frame, "CALLIDE", LEVELS["CALLIDE"], 1.3)
    assert (callide.peak, callide.time, callide.start) == (216.52, 60.0, 200.97)
    assert callide.margins["DCL"] == pytest.approx(216.52 - 219.13)
    assert callide.margins["PMF"] == pytest.approx(216.52 - 221.50)
    # a GWL the per-GWL map does not give is "not given", not a guess
    assert ur.dam_peak(frame, "CALLIDE", LEVELS["CALLIDE"], 1.7).margins["PMF"] is None


def test_the_dams_table_says_not_given_rather_than_leave_a_blank():
    plan = SimpleNamespace(storms=[SimpleNamespace(filename=f"{REAL}.72", loading="1 in 100",
                                                   gwl=1.3)])
    row, = ur.dams_rows(plan, DATA, ["CALLIDE", "KROOMBIT"], LEVELS)
    assert row["CALLIDE to DCL"] == "-2.61 m" and row["CALLIDE to PMF"] == "-4.98 m"
    assert row["KROOMBIT peak"] == "268.03" and row["KROOMBIT to DCL"] == "-2.67 m"
    assert row["KROOMBIT to PMF"] == "not given"
    assert ur.reaching_crest([row], ["CALLIDE", "KROOMBIT"]) == 0


def test_dam_levels_and_the_dams_come_from_the_model_config(tmp_path):
    config = tmp_path / "urbs_config.json"
    config.write_text(json.dumps({"max_keys": {"level": "CALLIDE"},
                                  "additional_dams": [{"location": "KROOMBIT"}],
                                  "dam_levels": {"_note": "x", **LEVELS}}))
    dams, levels, problems = ur.dam_setup(config)
    assert dams == ["CALLIDE", "KROOMBIT"] and levels == LEVELS and problems == []


def test_the_scatter_takes_every_finished_run_of_every_selection(tmp_path):
    selection = tmp_path / "RFSL_Design_GWL1p3_representative_events.json"
    (tmp_path / "RFSL_Design_GWL1p3_representative_events_downstream_storms.csv").write_text(
        f"filename,gwl\n{REAL}.72,1.3\nnot_run_yet.72,1.3\n")
    points = ur.margin_points([selection], DATA, ["CALLIDE", "KROOMBIT"], LEVELS,
                              lambda path: "RFSL GWL1.3")
    (x, y, run), = points["RFSL GWL1.3"]
    assert run == REAL and x == pytest.approx(-2.61) and y == pytest.approx(-2.67)


def test_the_crest_chart_draws_height_above_each_crest():
    from core import damchart
    frame = ur.read_levels(DATA, REAL)
    series = {dam: ur.crest_series(frame, dam, LEVELS[dam], 1.3) for dam in LEVELS}
    options = damchart.crest_chart(series, LEVELS, 1.3)
    drawn = {s["name"]: s for s in options["series"]}
    assert drawn["KROOMBIT"]["data"][3] == [45.0, pytest.approx(268.03 - 270.70)]
    # Callide has a PMF level at GWL 1.3, drawn on the same basis; Kroombit has none
    assert [m["yAxis"] for m in drawn["CALLIDE"]["markLine"]["data"]] == \
        [pytest.approx(221.50 - 219.13)]
    assert "markLine" not in drawn["KROOMBIT"]
    assert "crest" in drawn and drawn["crest"]["markLine"]["data"] == [{"yAxis": 0}]
    assert options["yAxis"]["name"] == "Above the dam crest (m)"

    # In m AHD the same series carry their crests back, each dam gets its own crest
    # line and its PMF as given, and the shared zero line goes.
    absolute = damchart.crest_chart(series, LEVELS, 1.3, absolute=True)
    drawn = {s["name"]: s for s in absolute["series"]}
    assert drawn["KROOMBIT"]["data"][3] == [45.0, pytest.approx(268.03)]
    assert "crest" not in drawn
    marks = {m["label"]["formatter"]: m["yAxis"] for m in drawn["CALLIDE"]["markLine"]["data"]}
    assert marks == {"CALLIDE crest": 219.13, "CALLIDE PMF": 221.50}
    assert [m["label"]["formatter"] for m in drawn["KROOMBIT"]["markLine"]["data"]] == \
        ["KROOMBIT crest"]
    assert absolute["yAxis"]["name"] == "Level (m AHD)"



# -- the print locations, and whether a run stopped mid-flood ------------------------------

def test_the_flows_are_the_star_print_locations_by_name(flows_folder):
    flows = ur.read_flows(flows_folder, REAL)
    assert list(flows.columns) == ["STAGCK_WEIR", "MALAKOFF", "RAINBOW_CK", "96_KM",
                                   "CLD_INFLOW", "GERRARD", "GOOVIGEN"]
    assert flows.index[0] == 0.0 and flows["GOOVIGEN"].max() > 0
    # the real run receded everywhere before it stopped
    assert ur.end_of_run(flows) == []


def test_a_run_that_stops_before_its_flood_has_passed_is_cut_off():
    import pandas as pd
    times = [0.0, 10.0, 20.0, 30.0, 40.0]
    assert not ur.cut_off(pd.Series([0, 100, 50, 10, 5], index=times))      # receded
    assert ur.cut_off(pd.Series([0, 100, 80, 60, 40], index=times))        # still 40% of peak
    assert ur.cut_off(pd.Series([0, 1, 2, 3, 4], index=times))             # still rising
    assert not ur.cut_off(pd.Series([0, 0, 0, 0, 0], index=times))         # no flow at all


def test_the_first_location_shown_is_the_dams_outflow_or_the_last(tmp_path):
    config = tmp_path / "c.json"
    config.write_text(json.dumps({"max_keys": {"outflow": "CLD_OUTFLOW_RFSL"}}))
    assert ur.default_location(config, ["A", "CLD_OUTFLOW_RFSL", "B"]) == "CLD_OUTFLOW_RFSL"
    assert ur.default_location(config, ["A", "B"]) == "B"



def test_the_peak_table_has_events_down_and_locations_across(flows_folder):
    flows = {REAL: ur.read_flows(flows_folder, REAL)}
    row, = ur.peaks_rows(flows, {REAL: "1 in 100"})
    assert row["loading"] == "1 in 100" and row["run"] == REAL
    assert set(row) == {"loading", "run", "STAGCK_WEIR", "MALAKOFF", "RAINBOW_CK", "96_KM",
                        "CLD_INFLOW", "GERRARD", "GOOVIGEN"}
    peak = flows[REAL]["GOOVIGEN"].max()
    assert row["GOOVIGEN"] == f"{peak:,.0f}"
    timed, = ur.peaks_rows(flows, {}, ur.TIME)
    assert timed["GOOVIGEN"] == f"{flows[REAL]['GOOVIGEN'].idxmax():g}"



# -- coming back later: the runs as the last run named them ---------------------------------

def _planned(**overrides):
    base = dict(filename="e00042_9h_GWL1p3.9", realisation=42, duration=9.0, gwl=1.3,
                method="ensemble", loading="PMF", event_label="event 42")
    base.update(overrides)
    return SimpleNamespace(**base)


def test_a_selection_comes_back_with_the_suffix_its_last_run_used(tmp_path):
    selection = tmp_path / "PMF_representative_events.json"
    plan = SimpleNamespace(storms=[_planned()])
    assert ur.recorded_suffix(plan, selection) == ""            # never run: blank
    (tmp_path / "PMF_representative_events_downstream_storms.csv").write_text(
        "filename,run_name\ne00042_9h_GWL1p3_PMF_RFSL.9,e00042_9h_GWL1p3_PMF_RFSL\n")
    assert ur.recorded_suffix(plan, selection) == "_PMF_RFSL"


def test_results_are_read_under_the_names_the_run_wrote_whatever_the_box_says(tmp_path, found):
    selection = tmp_path / "PMF_representative_events.json"
    (tmp_path / "PMF_representative_events_downstream_storms.csv").write_text(
        "filename,run_name\ne00042_9h_GWL1p3_PMF_RFSL.9,e00042_9h_GWL1p3_PMF_RFSL\n")
    _touch(found.results / "e00042_9h_GWL1p3_PMF_RFSL.q")
    plan = SimpleNamespace(storms=[_planned()])               # planned with an empty suffix
    ran = ur.as_run(plan, selection)
    assert ran.storms[0].filename == "e00042_9h_GWL1p3_PMF_RFSL.9"
    assert plan.storms[0].filename == "e00042_9h_GWL1p3.9"     # the plan itself is untouched
    row, = ur.runs_rows(ran, found, ur.storms_record(selection))
    assert row["status"] == "finished"
    # a storm the record does not have keeps its planned name
    other = SimpleNamespace(storms=[_planned(filename="e00022_4p5h_GWL2p7.4p5", realisation=22,
                                             duration=4.5, gwl=2.7)])
    assert ur.as_run(other, selection).storms[0].filename == "e00022_4p5h_GWL2p7.4p5"


def test_a_record_written_before_run_names_were_recorded_still_gives_the_suffix(tmp_path):
    selection = tmp_path / "RFSL_Design_GWL1p3_representative_events.json"
    (tmp_path / "RFSL_Design_GWL1p3_representative_events_downstream_storms.csv").write_text(
        "filename\n005044_72h_GWL1p3_RFSL.72\n")
    plan = SimpleNamespace(storms=[_planned(filename="005044_72h_GWL1p3.72", realisation=5044,
                                            duration=72.0, method="monte carlo")])
    assert ur.recorded_suffix(plan, selection) == "_RFSL"
