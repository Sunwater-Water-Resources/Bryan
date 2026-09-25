"""Nothing is removed from a study without being asked first."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import GROUP, PMF_GROUP, build_study         # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def opened(monkeypatch, tmp_path):
    import settings as settings_module
    from core import ensemble, figures, reporttables as rt, study as studies
    from state import STATE
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Table 26 near-term", source={"run": "E099 RFSL", "group": GROUP},
                aeps=[5, 10], pmp_aep=None)
    study.put_table(spec)
    figures.put(study, figures.new_spec(filename="levels", curves=[
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "GWL 1.3"}]))
    section = ensemble.settings(study)
    section["groups"].append({"label": "Near-term RFSL",
                              "ensemble": {"run": "E099 PMF", "group": PMF_GROUP},
                              "mc": {"run": "E099 RFSL", "group": GROUP}})
    ensemble.store(study, section)
    study.save()
    STATE.open_study(study.path)
    yield study
    STATE.study, STATE.project = None, None


def _saved(study):
    return json.loads(study.path.read_text(encoding="utf-8"))


@pytest.mark.asyncio
async def test_a_run_is_removed_only_when_confirmed(user, opened):
    await user.open("/report")
    user.find(marker="remove-run-E099 RFSL").click()
    await user.should_see(marker="confirm-dialog")
    await user.should_see("Table 26 near-term")          # the tables left naming it
    user.find(marker="confirm-cancel").click()
    assert "E099 RFSL" in [run["name"] for run in _saved(opened)["runs"]]

    user.find(marker="remove-run-E099 RFSL").click()
    user.find(marker="confirm-yes").click()
    assert "E099 RFSL" not in [run["name"] for run in _saved(opened)["runs"]]


@pytest.mark.asyncio
async def test_a_table_is_removed_only_when_confirmed(user, opened):
    await user.open("/report")
    user.find(marker="remove-table-26-near-term").click()
    user.find(marker="confirm-cancel").click()
    assert len(_saved(opened)["tables"]) == 1
    user.find(marker="remove-table-26-near-term").click()
    user.find(marker="confirm-yes").click()
    assert _saved(opened)["tables"] == []


@pytest.mark.asyncio
async def test_a_figure_is_removed_only_when_confirmed(user, opened):
    await user.open("/figures")
    user.find(marker="remove-levels").click()
    user.find(marker="confirm-cancel").click()
    assert len(_saved(opened)["figures"]) == 1
    user.find(marker="remove-levels").click()
    user.find(marker="confirm-yes").click()
    assert _saved(opened)["figures"] == []


@pytest.mark.asyncio
async def test_a_pmf_group_is_removed_only_when_confirmed(user, opened):
    await user.open("/pmf")
    user.find(marker="remove-pmf-group").click()
    await user.should_see("Near-term RFSL")
    user.find(marker="confirm-cancel").click()
    assert len(_saved(opened)["pmf"]["groups"]) == 1
    user.find(marker="remove-pmf-group").click()
    user.find(marker="confirm-yes").click()
    assert _saved(opened)["pmf"]["groups"] == []
