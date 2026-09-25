"""Out-of-date results said beside what is drawn from them."""

from __future__ import annotations

import os
import time

from core import staleness
from report_fixtures import GROUP, build_study


def make_stale(study, run="E099 RFSL", name="dam.sq"):
    """An input edited after the run: the rating curve, written now, dated later."""
    folder = study.run_config_path(run).parent / "reservoir"
    folder.mkdir(exist_ok=True)
    path = folder / name
    path.write_text("edited", encoding="utf-8")
    later = time.time() + 3600
    os.utime(path, (later, later))
    return path


def test_a_group_whose_results_are_newer_than_its_inputs_is_current(tmp_path):
    study = build_study(tmp_path / "study")
    state = staleness.of_group(study, "E099 RFSL", GROUP)
    assert state.is_current and state.total == 2
    assert staleness.for_sources(study, [{"run": "E099 RFSL", "group": GROUP}]) == []


def test_an_input_changed_after_the_run_is_named(tmp_path):
    study = build_study(tmp_path / "study")
    make_stale(study)
    said = staleness.for_sources(study, [{"run": "E099 RFSL", "group": GROUP},
                                         {"run": "E099 RFSL", "group": GROUP}])
    assert said == [f"E099 RFSL {GROUP} is stale: dam.sq changed after the run."]


def test_a_duration_not_run_is_said(tmp_path):
    study = build_study(tmp_path / "study")
    for path in study.run_config_path("E099 RFSL").parent.rglob("*_36h_*"):
        path.unlink()                                    # every 36 h result gone
    said = staleness.for_sources(study, [{"run": "E099 RFSL", "group": GROUP}])
    assert said == [f"E099 RFSL {GROUP}: 1 of 2 durations not run yet, so what is "
                    f"drawn from it leaves them out."]


def test_a_source_that_cannot_be_read_is_left_to_the_table(tmp_path):
    study = build_study(tmp_path / "study")
    assert staleness.for_sources(study, [{"run": "no such run", "group": GROUP},
                                         {"run": "E099 RFSL", "group": "no such group"},
                                         {"run": "", "group": ""}]) == []
