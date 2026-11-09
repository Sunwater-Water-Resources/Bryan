"""The dam's rating, as the design flood table's DCF row reads its outflow.

Pinned here: the rating is the one the group was routed through - the sims list's
for reservoir routing, the vec file's DAM ROUTE line for a run URBS routed - and
it is read the way lib/ReservoirRouting.py routes, with nothing extrapolated.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
import pytest

from core import damrating as dr
from report_fixtures import RATING_FSL, write_rating

VEC = """CALLIDE MODEL
{DAM ROUTE FSL=210.00 datafile=old.els il=210 location=CALLIDE file=old.sq}
DAM ROUTE FSL=215.50 datafile=dam.els il=initial_lake_level location=CALLIDE FILE=DAM.SQ
"""


@pytest.fixture(autouse=True)
def _fresh():
    dr.forget_cached()


def urbs_model(tmp_path, vec=VEC):
    """A model config, its vec and its ratings folder, laid out as the projects are."""
    (tmp_path / "urbs" / "model").mkdir(parents=True)
    (tmp_path / "urbs" / "model" / "callide.vec").write_text(vec, encoding="utf-8")
    write_rating(tmp_path / "urbs" / "ratings")
    config = tmp_path / "sim_options" / "urbs_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({
        "model_folder": "..\\urbs\\model", "vec_file": "callide.vec",
        "ratings_folder": "..\\urbs\\ratings", "max_keys": {"level": "CALLIDE"}}),
        encoding="utf-8")
    return config


def run_of(tmp_path, rows, model_config=None):
    frame = pd.DataFrame(rows)
    config = SimpleNamespace(filepaths={"model_config": model_config} if model_config else {},
                             config_path=tmp_path / "sims_config.json")
    return SimpleNamespace(frame=frame, config=config, project_folder=tmp_path)


def test_a_run_urbs_routed_takes_the_live_dam_route_line_not_a_commented_one(tmp_path):
    config = urbs_model(tmp_path)
    run = run_of(tmp_path, [{"Method": "monte carlo"}, {"Method": "Monte Carlo"}], config)
    rating = dr.rating_for(run, [0, 1])
    assert rating.fsl == RATING_FSL
    # FILE=DAM.SQ is dam.sq on disk: URBS runs on Windows, where case does not matter
    assert (rating.els.name, rating.sq.name.lower()) == ("dam.els", "dam.sq")
    assert rating.outflow_at(219.13) == pytest.approx(2_000 + 8_000 * (40_333.33 - 20_000)
                                                      / 40_000, abs=0.01)


def test_nothing_flows_at_or_below_full_supply(tmp_path):
    rating = dr.rating_for(run_of(tmp_path, [{"Method": "monte carlo"}],
                                  urbs_model(tmp_path)), [0])
    assert rating.outflow_at(RATING_FSL) == 0.0
    assert rating.outflow_at(212.0) == 0.0


def test_a_level_off_the_top_of_the_storage_curve_is_refused(tmp_path):
    rating = dr.rating_for(run_of(tmp_path, [{"Method": "monte carlo"}],
                                  urbs_model(tmp_path)), [0])
    with pytest.raises(dr.RatingError, match="above the top of the storage curve"):
        rating.outflow_at(226.0)


def test_a_dam_routed_by_volume_has_no_level_to_read_at(tmp_path):
    config = urbs_model(tmp_path, vec="DAM ROUTE VBF=ADV_below location=CALLIDE FILE=DAM.SQ\n")
    with pytest.raises(dr.RatingError, match="routed by volume"):
        dr.rating_for(run_of(tmp_path, [{"Method": "monte carlo"}], config), [0])


def test_a_vec_without_the_dam_says_which_location_it_looked_for(tmp_path):
    config = urbs_model(tmp_path, vec="DAM ROUTE FSL=215.5 datafile=dam.els location=OTHER "
                                      "file=dam.sq\n")
    with pytest.raises(dr.RatingError, match="no DAM ROUTE line for 'CALLIDE'"):
        dr.rating_for(run_of(tmp_path, [{"Method": "ensemble"}], config), [0])


def test_reservoir_routing_rows_that_name_different_ratings_are_refused(tmp_path):
    write_rating(tmp_path / "reservoir")
    write_rating(tmp_path / "reservoir", name="raised")
    rows = [{"Method": "reservoir routing", "ELS file": r"reservoir\dam.els",
             "SQ file": r"reservoir\dam.sq", "FSL": RATING_FSL},
            {"Method": "reservoir routing", "ELS file": r"reservoir\dam.els",
             "SQ file": r"reservoir\raised.sq", "FSL": RATING_FSL}]
    run = run_of(tmp_path, rows)
    assert dr.rating_for(run, [0]).sq.name == "dam.sq"
    with pytest.raises(dr.RatingError, match="different ratings: dam.sq .*raised.sq"):
        dr.rating_for(run, [0, 1])


def test_a_reservoir_routing_row_without_an_fsl_says_so(tmp_path):
    write_rating(tmp_path / "reservoir")
    run = run_of(tmp_path, [{"Method": "reservoir routing", "ELS file": r"reservoir\dam.els",
                             "SQ file": r"reservoir\dam.sq", "FSL": None}])
    with pytest.raises(dr.RatingError, match="no FSL"):
        dr.rating_for(run, [0])
