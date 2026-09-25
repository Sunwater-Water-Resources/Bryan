"""The runs panel down the left of every page, and the run sums it shows."""

from __future__ import annotations

import asyncio
import os

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import build_design_run, build_study         # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_state(monkeypatch, tmp_path):
    import settings as settings_module
    from core import runstatus
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    STATE.study, STATE.study_problem, STATE.project = None, "", None
    runstatus.forget()
    yield STATE
    STATE.study, STATE.study_problem, STATE.project = None, "", None
    runstatus.forget()


# -- the sums ----------------------------------------------------------------------

def test_a_runs_rows_are_summed_by_state(tmp_path):
    from core import completion, runstatus
    study = build_study(tmp_path / "study")
    status = runstatus.judge(study, "E099 RFSL")
    assert status.total > 0 and not status.problem
    assert set(status.counts) <= set(runstatus.ORDER)
    assert status.describe().split(", ")[0].split(" ", 1)[1] in runstatus.ORDER
    worst_first = [state for state, _ in status.in_order()]
    assert worst_first == sorted(worst_first, key=runstatus.ORDER.index)
    assert completion.UP_TO_DATE not in worst_first[:-1]


def test_the_sums_are_kept_until_they_age_or_the_list_changes(tmp_path):
    from core import runstatus
    study = build_study(tmp_path / "study")
    clock = [1000.0]
    now = lambda: clock[0]                                     # noqa: E731
    runstatus.forget()
    assert runstatus.cached(study, "E099 RFSL", now=now) is None
    first = runstatus.status_of(study, "E099 RFSL", now=now)
    assert runstatus.cached(study, "E099 RFSL", now=now) == first

    clock[0] += runstatus.MAX_AGE + 1
    assert runstatus.cached(study, "E099 RFSL", now=now) is None

    key = runstatus._key(study.run_config_path("E099 RFSL"))
    runstatus.status_of(study, "E099 RFSL", now=now)          # judged afresh
    judged_at = clock[0]
    clock[0] += 1
    runstatus.status_of(study, "E099 RFSL", now=now)          # kept
    assert runstatus._CACHE[key][1] == judged_at

    sims = runstatus.load_sims_config(study.run_config_path("E099 RFSL")).sims_list_path
    stat = sims.stat()
    os.utime(sims, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5_000_000_000))
    clock[0] += 1
    again = runstatus.status_of(study, "E099 RFSL", now=now)   # the list changed
    assert runstatus._CACHE[key][1] == clock[0]
    assert again.counts == first.counts
    runstatus.forget()


def test_a_run_that_is_gone_says_so(tmp_path):
    from core import runstatus
    study = build_study(tmp_path / "study")
    study.run_config_path("E099 PMF").unlink()
    assert runstatus.judge(study, "E099 PMF").problem == "sims_config.json not found"


def test_the_launcher_keeps_what_it_judges_of_the_open_run(tmp_path, private_state):
    from core import runstatus
    study = build_study(tmp_path / "study")
    private_state.open_project(study.run_config_path("E099 RFSL"))
    private_state.refresh_completion()
    kept = runstatus.cached(study, "E099 RFSL")
    assert kept is not None and kept.total == len(private_state.project.frame)


# -- the panel ---------------------------------------------------------------------

async def _until(condition, seconds=5.0):
    for _ in range(int(seconds / 0.05)):
        if condition():
            return True
        await asyncio.sleep(0.05)
    return condition()


@pytest.mark.asyncio
async def test_the_panel_lists_the_studys_runs_and_sums_them(user, private_state, tmp_path):
    from core import runstatus
    study = build_study(tmp_path / "study")
    private_state.open_study(study.path)
    await user.open("/select")
    await user.should_see(marker="runs-panel")
    await user.should_see(marker="panel-study")
    await user.should_see("Dam inputs")
    await user.should_see(marker="panel-run-E099 RFSL")
    await user.should_see(marker="panel-run-E099 PMF")
    assert await _until(lambda: runstatus.cached(study, "E099 PMF") is not None)


@pytest.mark.asyncio
async def test_a_run_is_swapped_without_leaving_the_page(user, private_state, tmp_path):
    study = build_study(tmp_path / "study")
    private_state.open_study(study.path)
    private_state.open_project(study.run_config_path("E099 RFSL"))
    await user.open("/results")
    user.find(marker="panel-run-E099 PMF").click()
    assert private_state.project.config.config_path.resolve() \
        == study.run_config_path("E099 PMF").resolve()


@pytest.mark.asyncio
async def test_folding_is_remembered(user, private_state):
    import settings as settings_module
    await user.open("/")
    user.find(marker="fold-runs").click()
    assert private_state.settings.runs_panel_open is False
    assert settings_module.UiSettings.load().runs_panel_open is False
    drawer = user.find(marker="runs-panel").elements.pop()
    assert "mini" in drawer.props

    user.find(marker="header-run").trigger("click")         # the run name toggles it
    assert private_state.settings.runs_panel_open is True
    assert "mini" not in drawer.props


@pytest.mark.asyncio
async def test_without_a_study_the_panel_points_to_one(user, private_state):
    await user.open("/select")
    await user.should_see("No study is open.")
    await user.should_see(marker="panel-open-other")


@pytest.mark.asyncio
async def test_an_open_list_outside_the_study_can_be_added(user, private_state, tmp_path):
    study = build_study(tmp_path / "study")
    other = build_design_run(tmp_path / "study" / "runs" / "E100")
    private_state.open_study(study.path)
    private_state.open_project(other)
    await user.open("/select")
    await user.should_see(marker="panel-add-to-study")


@pytest.mark.asyncio
async def test_the_menu_bar_calls_the_open_run_by_its_study_name(user, private_state,
                                                                 tmp_path):
    study = build_study(tmp_path / "study")
    private_state.open_project(study.run_config_path("E099 RFSL"))
    await user.open("/select")
    await user.should_see("CLD_RFSL_mc_sims_01.json")         # no study: the file
    private_state.open_study(study.path)
    await user.open("/select")
    header = user.find(marker="header-run").elements.pop()
    assert header.text == "E099 RFSL"
