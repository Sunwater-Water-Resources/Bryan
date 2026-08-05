"""The app state layer, driven exactly as the pages drive it.

No NiceGUI here - state.py is deliberately importable without it, so the whole
open -> select -> plan -> launch -> watch path can be tested headlessly.
"""

from __future__ import annotations

import sys
import time

import pytest

from conftest import TINAROO_COLUMNS, reservoir_row
from core import completion, runstate
from state import AppState


def wait_until(predicate, timeout=60.0, interval=0.05):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(interval)
    return False


@pytest.fixture
def app(project, fake_main, tmp_path, monkeypatch):
    """An AppState wired to the fake Bryan, with settings kept out of $HOME."""
    monkeypatch.setattr("settings.SETTINGS_PATH", tmp_path / "settings.json")
    state = AppState()
    state.settings.bryan_python = sys.executable
    state.settings.bryan_main = str(fake_main)
    state.settings.max_parallel = 1
    state.apply_settings()
    return state


def rows_for(durations=(18, 24)):
    return [reservoir_row(Duration=d,
                          **{"Output file": f"out_{d}h", "Output suffix": f"o{d}",
                             "Log file": rf"sims_mc\log\r{d}.log"})
            for d in durations]


def test_opening_seeds_the_selection_from_include(app, project):
    rows = rows_for()
    rows[1]["Include"] = "no"
    app.open_project(project(TINAROO_COLUMNS, rows))
    assert app.selected == {0}, "only Include = yes rows start selected"


def test_labels_and_groups_are_available(app, project):
    app.open_project(project(TINAROO_COLUMNS, rows_for()))
    assert app.project.label(0) == "out_18h"
    assert len(app.project.groups()) == 2


def test_completion_flows_through_to_the_bulk_actions(app, project):
    app.open_project(project(TINAROO_COLUMNS, rows_for()))
    app.refresh_completion()
    assert set(app.rows_needing_rerun()) == {0, 1}
    assert app.already_done([0, 1]) == []


def test_a_whole_run_end_to_end(app, project):
    """Open, select, plan, launch, watch it finish - the real path."""
    app.open_project(project(TINAROO_COLUMNS, rows_for()))
    app.set_selected([0, 1])

    assert app.preflight() == [], "the fixture project should be clean"

    plan = app.plan()
    assert not plan.blocked
    assert len(plan.chunks) == 1          # max_parallel is 1

    record, folder = app.launch(plan)
    assert wait_until(lambda: (app.manager.poll(), not record.is_live)[1], 90)
    assert record.status == runstate.COMPLETED

    # The results the fake Bryan wrote are found, so the rows now read as done.
    app.refresh_completion()
    states = {app.completion_of(index).state for index in (0, 1)}
    assert states == {completion.UP_TO_DATE}
    assert sorted(app.already_done([0, 1])) == [0, 1]


def test_rerunning_a_finished_row_asks_before_overwriting(app, project):
    app.open_project(project(TINAROO_COLUMNS, rows_for()))
    app.set_selected([0, 1])
    record, _ = app.launch(app.plan())
    assert wait_until(lambda: (app.manager.poll(), not record.is_live)[1], 90)

    app.refresh_completion()
    assert app.already_done(app.selected_in_order()), \
        "the confirm dialog must be triggered for up-to-date rows"


def test_editing_an_input_makes_the_results_stale(app, project):
    import os

    app.open_project(project(TINAROO_COLUMNS, rows_for((18,))))
    app.set_selected([0])
    record, _ = app.launch(app.plan())
    assert wait_until(lambda: (app.manager.poll(), not record.is_live)[1], 90)

    app.refresh_completion()
    assert app.completion_of(0).state == completion.UP_TO_DATE

    sq = app.project.config.project_folder / "reservoir/dam.sq"
    os.utime(sq, (time.time() + 120,) * 2)
    app.refresh_completion()

    state = app.completion_of(0)
    assert state.state == completion.STALE
    assert state.newest_input.name == "dam.sq"


def test_preflight_blocks_colliding_rows(app, project):
    rows = [reservoir_row(Duration=d, **{"Output file": "shared"})
            for d in (18, 24)]
    app.open_project(project(TINAROO_COLUMNS, rows))
    app.set_selected([0, 1])

    codes = {issue.code for issue in app.preflight()}
    assert "output-collision" in codes


def test_reload_keeps_the_selection_by_name(app, project):
    app.open_project(project(TINAROO_COLUMNS, rows_for((18, 24, 36))))
    app.set_selected([1])
    app.reload_project()
    assert app.selected == {1}, "the same simulation should still be selected"


def test_runs_are_rediscovered_for_the_history_page(app, project):
    app.open_project(project(TINAROO_COLUMNS, rows_for((18,))))
    app.set_selected([0])
    record, _ = app.launch(app.plan())
    assert wait_until(lambda: (app.manager.poll(), not record.is_live)[1], 90)

    fresh = AppState()
    fresh.settings.bryan_python = sys.executable
    fresh.settings.bryan_main = app.settings.bryan_main
    fresh.open_project(app.project.config.config_path)
    assert [r.run_id for r in fresh.runs()] == [record.run_id]
