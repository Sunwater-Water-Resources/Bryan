"""The Simulations page: the study's runs opened with a click, others added to it."""

from __future__ import annotations

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import build_study                           # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_state(monkeypatch, tmp_path):
    import settings as settings_module
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    STATE.study, STATE.study_problem, STATE.project = None, "", None
    yield STATE
    STATE.study, STATE.study_problem, STATE.project = None, "", None


# -- guessing a run's name ---------------------------------------------------------

@pytest.mark.parametrize("path, name", [
    ("runs/E013/CLD_RFSL_mc_sims_01.json", "E013 RFSL"),
    ("runs\\E013\\CLD_FSL_mc_sims_01.json", "E013 FSL"),
    ("runs/E013/PMF_sims.json", "E013 PMF"),
    ("sims.json", "sims"),
    ("", ""),
])
def test_a_run_is_named_from_its_folder_and_state(path, name):
    from core.study import guess_run_name
    assert guess_run_name(path) == name


def test_the_study_knows_which_run_a_sims_list_is(tmp_path):
    study = build_study(tmp_path / "study")
    design = study.run_config_path("E099 RFSL")
    assert study.run_named_for(design) == "E099 RFSL"
    assert study.run_named_for(tmp_path / "elsewhere.json") is None


# -- the page ----------------------------------------------------------------------

@pytest.mark.asyncio
async def test_without_a_study_a_sims_config_is_opened_by_path(user, private_state,
                                                               tmp_path):
    study = build_study(tmp_path / "study")
    design = study.run_config_path("E099 RFSL")
    await user.open("/")
    await user.should_see("Open a sims_config.json")
    await user.should_not_see(marker="study-runs-card")
    user.find(marker="sims-config-path").clear().type(str(design))
    user.find(marker="open-sims").click()
    await user.should_see("Select")
    assert private_state.project.config.config_path.resolve() == design.resolve()


@pytest.mark.asyncio
async def test_the_studys_runs_are_listed_and_opened_with_a_click(user, private_state,
                                                                 tmp_path):
    study = build_study(tmp_path / "study")
    private_state.open_study(study.path)
    await user.open("/")
    await user.should_see(marker="study-runs-card")
    await user.should_see("E099 RFSL")
    await user.should_see("E099 PMF")
    await user.should_see("Open another sims_config.json")
    user.find(marker="open-run-E099 PMF").click()
    await user.should_see("Select")
    assert private_state.project.config.config_path.resolve() \
        == study.run_config_path("E099 PMF").resolve()

    await user.open("/")                     # now marked as the open one
    await user.should_see(marker="run-open-E099 PMF")
    await user.should_not_see(marker="add-to-study-card")


@pytest.mark.asyncio
async def test_a_run_whose_sims_config_is_gone_says_so(user, private_state, tmp_path):
    study = build_study(tmp_path / "study")
    study.run_config_path("E099 PMF").unlink()
    private_state.open_study(study.path)
    await user.open("/")
    await user.should_see("not found")
    await user.should_see(marker="open-run-E099 RFSL")
    await user.should_not_see(marker="open-run-E099 PMF")


@pytest.mark.asyncio
async def test_a_sims_list_outside_the_study_is_added_to_it(user, private_state,
                                                           tmp_path):
    from report_fixtures import build_design_run
    study = build_study(tmp_path / "study")
    other = build_design_run(tmp_path / "study" / "runs" / "E100")
    private_state.open_study(study.path)
    private_state.open_project(other)
    await user.open("/")
    await user.should_see(marker="add-to-study-card")
    assert user.find(marker="add-to-study-name").elements.pop().value == "E100 RFSL"
    user.find(marker="add-to-study").click()
    await user.should_see(marker="run-open-E100 RFSL")

    from core import study as studies
    saved = studies.load_study(study.path)
    assert saved.run_named_for(other) == "E100 RFSL"


@pytest.mark.asyncio
async def test_the_menu_and_empty_states_name_the_simulations_page(user, private_state):
    from layout import NAV
    assert ("Simulations", "/") in NAV
    assert not any(label == "Project" for label, _ in NAV)
    await user.open("/select")
    await user.should_see("Go to Simulations")
