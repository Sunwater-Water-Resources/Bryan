"""The dam inputs: kept once in the study, filled from where they used to be, and
still written there for launchers that predate them."""

from __future__ import annotations

import json

import pytest

from core import dam as dams
from core import lakerecord
from core import reporttables as rt
from core import study as studies
from report_fixtures import build_study


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def old_lake_record(study):
    """A study as the launcher before the dam section saved it."""
    study.extra["lake_record"] = {
        "rainfall": {"shapefile": "catchment/cld.shp", "field": "NAME", "value": "Callide",
                     "output": "lake_record/rainfall.csv"},
        "homogenise": {"gauges": ["gauge/A.csv", "gauge/C.csv"],
                       "overlay": {"file": "gauge/B.csv", "below": 200.9,
                                   "reconnect_margin": 0.1},
                       "storage": "storage.els", "register": "register.xlsx",
                       "evaporation": "silo.txt", "pan_factors": [0.8] * 12,
                       "water_year_start": 7, "step": "30min",
                       "targets": [{"name": "RFSL", "rating": "rfsl.sq", "fsl": 215.5}]},
        "inflow": {"catchment_km2": 519},
    }
    study.save()
    return study


def test_an_existing_study_fills_the_dam_inputs_from_its_lake_record(study):
    dam = dams.settings(old_lake_record(study))
    assert dam["gauges"] == ["gauge/A.csv", "gauge/C.csv"]
    assert dam["overlay"]["below"] == 200.9
    assert (dam["storage"], dam["register"], dam["evaporation"]) == \
        ("storage.els", "register.xlsx", "silo.txt")
    assert dam["pan_factors"] == [0.8] * 12 and dam["water_year_start"] == 7
    assert (dam["shapefile"], dam["field"], dam["value"]) == \
        ("catchment/cld.shp", "NAME", "Callide")
    assert dam["catchment_km2"] == 519


def test_the_jobs_are_unchanged_by_the_move(study):
    """Callide's study, opened in the new launcher, homogenises exactly as before."""
    old_lake_record(study)
    before = lakerecord.settings(study)            # dam section not yet written
    dams.store(study, dams.settings(study))        # the Study page saves it
    after = lakerecord.settings(study)
    assert lakerecord.homogenise_job(study, after) == lakerecord.homogenise_job(study, before)
    assert after["rainfall"]["shapefile"] == "catchment/cld.shp"
    assert after["inflow"]["catchment_km2"] == 519
    assert after["homogenise"]["step"] == "30min"            # the analysis's own, kept


def test_the_dam_section_wins_over_stale_lake_record_copies(study):
    old_lake_record(study)
    dam = dams.settings(study)
    dam["storage"] = "storage/new.els"
    dams.store(study, dam)
    study.extra["lake_record"]["homogenise"]["storage"] = "stale.els"   # an old launcher
    assert lakerecord.settings(study)["homogenise"]["storage"] == "storage/new.els"


def test_saving_the_lake_record_never_writes_old_dam_inputs_back(study):
    """The Lake record page opened before the storage table was changed."""
    old_lake_record(study)
    section = lakerecord.settings(study)                 # read now
    dam = dams.settings(study)
    dam["register"] = "ratings/new.xlsx"
    dams.store(study, dam)                               # changed on the Study page
    section["homogenise"]["recession_correction"] = False
    lakerecord.store(study, section)                     # the Lake record page saves
    assert study.extra["dam"]["register"] == "ratings/new.xlsx"
    assert study.extra["lake_record"]["homogenise"]["register"] == "ratings/new.xlsx"
    assert study.extra["lake_record"]["homogenise"]["recession_correction"] is False


def test_an_older_launcher_still_finds_the_inputs_where_it_looks(study):
    dam = dams.settings(study)
    dam.update(gauges=["gauge/A.csv"], storage="s.els", catchment_km2=520.0,
               shapefile="c.shp", water_year_start=9)
    dams.store(study, dam)
    study.save()
    saved = json.loads(study.path.read_text(encoding="utf-8"))
    old = saved["lake_record"]
    assert old["homogenise"]["gauges"] == ["gauge/A.csv"]
    assert old["homogenise"]["storage"] == "s.els"
    assert old["homogenise"]["water_year_start"] == 9
    assert old["rainfall"]["shapefile"] == "c.shp"
    assert old["inflow"]["catchment_km2"] == 520.0
    assert saved["dam"]["storage"] == "s.els"


def test_a_new_study_has_the_defaults(study):
    dam = dams.settings(study)
    assert dam == dams.DEFAULTS and dam is not dams.DEFAULTS


# -- the Study page ---------------------------------------------------------------------

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")


@pytest_asyncio.fixture
async def user():
    from nicegui.testing.user_simulation import user_simulation
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def opened(study, monkeypatch, tmp_path):
    import settings as settings_module
    from state import STATE
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")
    monkeypatch.setattr(STATE, "settings", settings_module.UiSettings.load())
    STATE.open_study(study.path)
    yield study
    STATE.study = None


async def _type(user, marker, text):
    import asyncio
    box = user.find(marker=marker)
    box.clear().type(text)
    box.trigger("blur")
    await asyncio.sleep(0.2)


@pytest.mark.asyncio
async def test_the_dam_inputs_are_entered_on_the_study_page(user, opened):
    await user.open("/study")
    await user.should_see(marker="dam-card")
    folder = opened.folder
    await _type(user, "dam-gauges", f"{folder / 'gauge' / 'A.csv'}\n{folder / 'gauge' / 'C.csv'}")
    await _type(user, "dam-register", str(folder / "ratings" / "RatingCurves.xlsx"))
    await _type(user, "dam-area", "519")
    stored = dams.settings(studies.load_study(opened.path))
    assert stored["gauges"] == ["gauge/A.csv", "gauge/C.csv"]      # relative to the study
    assert stored["register"] == "ratings/RatingCurves.xlsx"
    assert stored["catchment_km2"] == 519.0


@pytest.mark.asyncio
async def test_twelve_pan_factors_are_needed(user, opened):
    await user.open("/study")
    await _type(user, "dam-pan", "0.8 0.8 0.8")
    await user.should_see("Give twelve pan factors")
    assert dams.settings(studies.load_study(opened.path))["pan_factors"] == \
        dams.CALLIDE_PAN_FACTORS
