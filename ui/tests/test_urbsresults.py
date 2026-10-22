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
