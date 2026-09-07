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


@pytest.mark.asyncio
async def test_picking_an_event_leaves_the_cards_where_they_were(user, project):
    """Every pick redraws the details, so the open cards have to survive it.

    They did not: the expansions were built with `value = position == 0`, so
    choosing an event in the third loading collapsed it and sprang the first
    one open underneath.
    """
    _open(project)
    await user.open("/events")

    def cards():
        return [next(iter(user.find(marker=f"target-{position}").elements))
                for position in range(3)]

    cards()[0].set_value(False)
    cards()[2].set_value(True)
    await user.should_see("Loadings")

    _pick_candidate(user, position=2, index=1)
    await user.should_see("Loadings")

    assert [card.value for card in cards()] == [False, False, True]


def _pick_candidate(user, position, index):
    """Choose a candidate the way the table's radio does."""
    table = next(iter(user.find(marker=f"candidates-{position}").elements))
    listener = next(listener for listener in table._event_listeners.values()
                    if listener.type == "pick")
    table._handle_event({"listener_id": listener.id, "args": table.rows[index]})


@pytest.mark.asyncio
async def test_the_rank_order_can_be_switched_to_the_result(user, project):
    """Reaching the loading is often what matters; neutrality is then a flag."""
    _open(project)
    await user.open("/events")

    toggle = next(iter(user.find(marker="rank-order").elements))
    assert toggle.value == "delta_z"

    toggle.set_value("result")
    await user.should_see("Loadings")
    await user.should_see("design value at 1 in")

    table = next(iter(user.find(marker="candidates-0").elements))
    assert "Δ target" in [column["label"] for column in table.columns]
    assert table.rows[0]["delta_value"] != "-"


@pytest.mark.asyncio
async def test_the_band_defaults_to_twenty_millimetres_of_level(user, project):
    _open(project)
    await user.open("/events")

    band = next(iter(user.find(marker="result-band").elements))
    assert band.value == 20
    await user.should_see("mm")


@pytest.mark.asyncio
async def test_the_band_follows_the_result_type(user, project):
    """20 mm of lake level is not 20 m3/s, so the two are kept apart."""
    _open(project)
    await user.open("/events")

    band = next(iter(user.find(marker="result-band").elements))
    band.set_value(50)
    await user.should_see("Loadings")

    next(iter(user.find(marker="result-type").elements)).set_value("inflow")
    await user.should_see("Loadings")
    band = next(iter(user.find(marker="result-band").elements))
    assert band.value == 10                       # the inflow default, in m3/s

    next(iter(user.find(marker="result-type").elements)).set_value("level")
    await user.should_see("Loadings")
    band = next(iter(user.find(marker="result-band").elements))
    assert band.value == 50                       # the level band, as it was left


@pytest.mark.asyncio
async def test_the_rounding_is_its_own_field(user, project):
    """The band and the grid answer different questions, so they are two boxes."""
    _open(project)
    await user.open("/events")

    band = next(iter(user.find(marker="result-band").elements))
    rounding = next(iter(user.find(marker="result-rounding").elements))
    assert (band.value, rounding.value) == (20, 10)

    rounding.set_value(25)
    next(iter(user.find(marker="result-type").elements)).set_value("inflow")
    await user.should_see("Loadings")
    assert next(iter(user.find(marker="result-rounding").elements)).value == 5

    next(iter(user.find(marker="result-type").elements)).set_value("level")
    await user.should_see("Loadings")
    assert next(iter(user.find(marker="result-rounding").elements)).value == 25
