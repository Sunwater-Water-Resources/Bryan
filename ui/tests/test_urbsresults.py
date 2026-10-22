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
    assert drawn["CALLIDE"]["markLine"]["data"] == [{"yAxis": pytest.approx(221.50 - 219.13)}]
    assert "markLine" not in drawn["KROOMBIT"]
