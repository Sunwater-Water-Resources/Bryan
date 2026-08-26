"""Render the Select page and drive the skip-the-bad-rows path through it.

The core logic is tested in test_preflight.py and test_state.py; this exists
because the failure that path is most exposed to - a dialog that never opens,
or a button wired to nothing - is invisible to those. Same self-built NiceGUI
simulation as test_results_page.py, for the same reason.
"""

from __future__ import annotations

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from conftest import TINAROO_COLUMNS, reservoir_row              # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


def rows_for(durations=(18, 24, 36)):
    return [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}",
                                         "Log file": rf"sims_mc\log\r{d}.log"})
            for d in durations]


def _open(project, rows):
    from state import STATE
    config = project(TINAROO_COLUMNS, rows)
    STATE.open_project(config)
    STATE.set_selected(set(STATE.project.frame.index))
    return STATE


def _with_a_missing_rating_curve(project, bad=1):
    """One row of three pointing at a .sq that is not there."""
    rows = rows_for()
    rows[bad]["SQ file"] = r"reservoir\gone.sq"
    state = _open(project, rows)
    (state.project.config.project_folder / "reservoir/gone.sq").unlink()
    return state


@pytest.mark.asyncio
async def test_the_summary_counts_the_rows_that_would_stop_the_run(user, project):
    _with_a_missing_rating_curve(project)
    await user.open("/select")
    await user.should_see("1 would stop the run")
    await user.should_see("not found")


@pytest.mark.asyncio
async def test_deselecting_the_problem_rows_leaves_the_rest(user, project):
    state = _with_a_missing_rating_curve(project)
    await user.open("/select")
    user.find("Deselect problem rows").click()

    assert state.selected == {0, 2}
    await user.should_see("Deselected 1 row(s)")
    assert state.preflight() == [], "what is left has to actually run"


@pytest.mark.asyncio
async def test_deselecting_says_so_when_there_is_nothing_to_drop(user, project):
    state = _open(project, rows_for())
    await user.open("/select")
    user.find("Deselect problem rows").click()

    await user.should_see("No selected row is blocked")
    assert state.selected == {0, 1, 2}


@pytest.mark.asyncio
async def test_the_run_dialog_offers_to_skip_the_blocked_rows(user, project):
    _with_a_missing_rating_curve(project)
    await user.open("/select")
    user.find("Check and run").click()

    await user.should_see("must be fixed first")
    await user.should_see("2 of 3 rows are fine")
    await user.should_see("Skip 1 and run the rest")


@pytest.mark.asyncio
async def test_skipping_gets_through_to_the_run_dialog(user, project):
    state = _with_a_missing_rating_curve(project)
    await user.open("/select")
    user.find("Check and run").click()
    user.find("Skip 1 and run the rest").click()

    # the selection follows, so the page is not left claiming three rows
    assert state.selected == {0, 2}
    await user.should_see("1 row(s) are being skipped")
    await user.should_see("Estimated wall clock")


@pytest.mark.asyncio
async def test_a_selection_that_cannot_be_rescued_offers_no_skip(user, project):
    rows = rows_for((18,))
    rows[0]["SQ file"] = r"reservoir\gone.sq"
    state = _open(project, rows)
    (state.project.config.project_folder / "reservoir/gone.sq").unlink()

    await user.open("/select")
    user.find("Check and run").click()
    await user.should_see("nothing to run")
    with pytest.raises(AssertionError):
        user.find("Skip 1 and run the rest")
