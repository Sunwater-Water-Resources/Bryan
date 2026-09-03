"""Render the Events page against a project whose runs have databases.

The core modules are tested without a browser stack; this exists because the
failure this feature is most exposed to - the page raising on load, or a table
quietly rendering nothing - is invisible to those tests. Same arrangement as
test_results_page.py: its own NiceGUI simulation, so the rest of the suite still
runs in an environment without nicegui.
"""

from __future__ import annotations

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from test_events import build, mcdf_frame                        # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def project(tmp_path):
    from core import events
    events.forget_cached()
    return build(tmp_path)


@pytest.fixture
def bare(tmp_path):
    """Analysed runs, but no database - the empty state."""
    from core import events
    events.forget_cached()
    return build(tmp_path, with_mcdf=False)


def _open(config):
    from state import STATE
    STATE.open_project(config)
    return STATE


@pytest.mark.asyncio
async def test_the_page_offers_the_default_loadings(user, project):
    _open(project)
    await user.open("/events")
    await user.should_see("Loadings")
    await user.should_see("1 in 1,000")


@pytest.mark.asyncio
async def test_the_chosen_event_reaches_the_summary(user, project):
    _open(project)
    await user.open("/events")
    await user.should_see("Representative events")

    table = next(iter(user.find(marker="summary-table").elements))
    chosen = {row["loading"]: row["hydrograph"] for row in table.rows}
    assert set(chosen) == {"1 in 100", "1 in 1,000", "1 in 10,000"}
    # The AEP-neutral realisation of the fixture is the 1 in 1,000 one.
    assert chosen["1 in 1,000"] == "sim_00000"


@pytest.mark.asyncio
async def test_the_candidate_metrics_are_shown(user, project):
    """The flags are the reason to look at this page rather than sort a csv."""
    _open(project)
    await user.open("/events")

    table = next(iter(user.find(marker="candidates-0").elements))
    headings = [column["label"] for column in table.columns]
    assert {"Sub-burst", "Pre-burst p", "Lake z", "Flags"} <= set(headings)
    assert table.rows


@pytest.mark.asyncio
async def test_a_project_without_a_database_says_what_is_missing(user, bare):
    _open(bare)
    await user.open("/events")
    await user.should_see("No Monte Carlo database found")


@pytest.mark.asyncio
async def test_the_page_is_in_the_navigation(user, project):
    _open(project)
    await user.open("/results")
    await user.should_see("Events")


@pytest.mark.asyncio
async def test_typing_a_loading_value_does_not_rebuild_the_row(user, project):
    """The row must survive being typed into.

    Redrawing the loadings on every keystroke destroys the input the digits
    are going into, so the field took the first one and lost focus.
    """
    _open(project)
    await user.open("/events")

    field = next(iter(user.find(marker="loading-value-0").elements))
    identity = field.id
    for digits in (2, 20, 200, 2000):                 # typing "2000"
        field.set_value(digits)
        await user.should_see("Loadings")

    still = next(iter(user.find(marker="loading-value-0").elements))
    assert still.id == identity, "the input was replaced mid-edit"
    assert still.value == 2000


@pytest.mark.asyncio
async def test_adding_a_loading_does_redraw_the_rows(user, project):
    """The other half of it: a new row has to appear."""
    _open(project)
    await user.open("/events")
    assert len(user.find(marker="loading-value-2").elements) == 1

    user.find("Add loading").click()
    await user.should_see("Loadings")
    assert len(user.find(marker="loading-value-3").elements) == 1
