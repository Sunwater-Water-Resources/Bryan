"""The PMF page's core: which ensemble event, and the notional AEP of the PMF.

Pinned here: the two picking conventions (the highest event for the PMF, Bryan's
median pattern for anything else) with EnbAnalysis's exact median position; the
fit recovering a known AEP; and every way the fit can mislead - too few
realisations, a PMF above the sample, a window below the PMF - saying so.
"""

from __future__ import annotations

import math
from statistics import NormalDist

import numpy as np
import pandas as pd
import pytest

from core import ensemble, pmfchart
from core import reporttables as rt
from core import study as studies
from report_fixtures import (GROUP, MC_BASE, MC_SLOPE, MC_TOP_Z, PMF_GROUP,
                             build_study)


def events(rows):
    """An ensemble database from (duration, tp, level, inflow, outflow)."""
    return pd.DataFrame(rows, columns=["duration", "tp", "level", "inflow", "outflow"]) \
        .assign(storm_method="GSDM")


@pytest.fixture
def frame():
    # 4.5 h: ten patterns, levels 220.0 .. 220.9; 9 h: ten patterns 220.5 .. 221.4,
    # with the single highest event (221.4) at 9 h pattern 9.
    rows = [(4.5, tp, 220.0 + tp / 10, 14000 + tp, 12000 + tp) for tp in range(10)]
    rows += [(9.0, tp, 220.5 + tp / 10, 15000 + tp, 13000 + tp) for tp in range(10)]
    return events(rows)


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


# -- which event ---------------------------------------------------------------

@pytest.mark.parametrize("count, position", [(10, 5), (9, 4), (5, 2), (7, 4), (1, 0)])
def test_the_median_position_is_enbanalysis_own(count, position):
    # int(np.around(n / 2)) - numpy rounds half to even, so 5 patterns give 2
    assert ensemble.median_position(count) == position


def test_the_pmf_is_the_highest_event_with_its_own_flows(frame):
    picked = ensemble.pick(frame, ensemble.HIGHEST)
    assert (picked.duration, picked.level, picked.inflow, picked.outflow) == \
        (9.0, pytest.approx(221.4), 15009, 13009)
    assert picked.pattern == "GSDM: 9"


def test_the_median_pick_is_the_largest_per_duration_median(frame):
    picked = ensemble.pick(frame, ensemble.MEDIAN)
    # 9 h, sixth of ten ascending: pattern 5, 221.0 m
    assert (picked.duration, picked.level, picked.pattern) == \
        (9.0, pytest.approx(221.0), "GSDM: 5")
    assert picked.inflow == 15005


def test_the_spread_by_duration_carries_both_picks(frame):
    stats = ensemble.by_duration(frame)
    assert list(stats.index) == [4.5, 9.0]
    assert stats.loc[9.0, "median"] == pytest.approx(221.0)
    assert stats.loc[9.0, "max"] == pytest.approx(221.4)
    assert stats.loc[4.5, "median_pattern"] == "GSDM: 5"
    assert stats.loc[9.0, "n"] == 10


def test_an_empty_database_picks_nothing():
    assert not ensemble.pick(events([]), ensemble.HIGHEST).found


def test_the_pmf_table_can_take_the_median_instead(study):
    spec = rt.new_spec(rt.ENSEMBLE_PEAK)
    spec["sections"] = [{"heading": "", "rows": [
        {"label": "Near-Term", "run": "E099 PMF", "group": PMF_GROUP}]}]
    assert rt.build(study, spec).rows[0].cells[1] == "221.38"
    spec["pick"] = ensemble.MEDIAN
    # one pattern per duration, so each duration's median is its only event and
    # the largest of those is still the 9 h one
    assert rt.build(study, spec).rows[0].cells[1] == "221.38"


# -- the realisations ------------------------------------------------------------

def test_the_variate_keeps_its_digits_far_out_in_the_tail():
    p = 2.58e-8                                       # 1 in 38.8 million
    z = ensemble.variate(p)
    assert z == pytest.approx(-NormalDist().inv_cdf(p))
    assert ensemble.aep_of_variate(z) == pytest.approx(1 / p, rel=1e-9)
    assert math.isnan(ensemble.variate(0.0)) and math.isnan(ensemble.variate(1.0))


def test_realisations_are_read_as_probabilities_and_sorted(study):
    path = ensemble.mc_databases(study, "E099 RFSL", GROUP)["72h"]
    real = ensemble.read_realisations(path)
    assert np.all(np.diff(real.z) >= 0)
    assert real.top_value == pytest.approx(MC_BASE + MC_SLOPE * MC_TOP_Z)
    assert real.top_aep == pytest.approx(ensemble.aep_of_variate(MC_TOP_Z), rel=1e-6)


def test_an_mcdf_without_the_aep_column_is_explained(tmp_path):
    path = tmp_path / "x__mcdf.csv"
    pd.DataFrame({"level": [1.0]}).to_csv(path)
    with pytest.raises(studies.StudyError, match="level_aep"):
        ensemble.read_realisations(path)


def test_the_databases_are_found_by_duration_and_the_nearest_is_offered(study):
    databases = ensemble.mc_databases(study, "E099 RFSL", GROUP)
    assert list(databases) == ["36h", "72h"]
    assert ensemble.nearest_duration(databases, 9.0) == "36h"
    assert ensemble.nearest_duration(databases, 60) == "72h"


# -- the fit -----------------------------------------------------------------------

def realisations(top_z=MC_TOP_Z, count=5000):
    z = np.linspace(-1.0, top_z, count)
    return ensemble.Realisations(path=None, result="level", z=z,
                                 value=MC_BASE + MC_SLOPE * z)


def true_aep(level):
    return ensemble.aep_of_variate((level - MC_BASE) / MC_SLOPE)


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_the_fit_recovers_a_known_aep(degree):
    estimate = ensemble.fit(realisations(), 221.38, 500_000, None, degree)
    assert estimate.ok and not estimate.warnings
    assert estimate.aep == pytest.approx(true_aep(221.38), rel=0.02)


def test_a_pmf_above_every_realisation_is_called_an_extrapolation():
    estimate = ensemble.fit(realisations(top_z=5.0), 221.38, 500_000)
    assert estimate.ok
    assert any("above every realisation" in warning for warning in estimate.warnings)


def test_a_window_that_stops_below_the_pmf_says_so():
    # AEPofPMF.py's 1 in 4,000,000 upper bound did exactly this at Callide
    estimate = ensemble.fit(realisations(), 221.38, 500_000, 4_000_000)
    assert any("above the window" in warning for warning in estimate.warnings)


def test_too_few_realisations_refuses_rather_than_fits():
    estimate = ensemble.fit(realisations(count=200), 221.38, 50_000_000)
    assert not estimate.ok and "widen the window" in estimate.refused


def test_the_grid_is_three_windows_by_three_degrees():
    grid = ensemble.sensitivity(realisations(), 221.38, 500_000)
    assert list(grid.index) == [200_000.0, 500_000.0, 1_000_000.0]
    assert list(grid.columns) == [1, 2, 3]
    assert np.all(np.isfinite(grid.to_numpy()))


# -- the study keeps it --------------------------------------------------------------

def test_the_pmf_settings_survive_a_save_and_follow_a_renamed_run(study):
    section = ensemble.settings(study)
    section["adopted_aep"] = 8_000_000
    section["groups"].append({"label": "Near-term", "ensemble": {"run": "E099 PMF",
                              "group": PMF_GROUP}, "mc": {"run": "E099 RFSL", "group": GROUP}})
    ensemble.store(study, section)
    study.save()
    study = studies.load_study(study.path)
    study.rename_run("E099 RFSL", "E100 RFSL")
    reread = ensemble.settings(study)
    assert reread["adopted_aep"] == 8_000_000
    assert reread["groups"][0]["mc"]["run"] == "E100 RFSL"
    assert reread["groups"][0]["degree"] == 1 and reread["groups"][0]["duration"] == ""
    assert study.remove_run("E099 PMF") == ["PMF Near-term"]


# -- the charts --------------------------------------------------------------------

def test_the_box_chart_marks_the_highest_event(frame):
    highest = ensemble.pick(frame, ensemble.HIGHEST)
    options = pmfchart.box_chart(ensemble.by_duration(frame),
                                 ensemble.pattern_points(frame), highest=highest,
                                 reference_levels=[{"label": "DCL", "level": 219.13}])
    kinds = [series["type"] for series in options["series"]]
    assert kinds == ["boxplot", "scatter", "scatter"]
    assert options["series"][2]["data"][0][:2] == [1, pytest.approx(221.4)]
    assert options["series"][0]["markLine"]["data"][0]["yAxis"] == 219.13
    assert options["xAxis"]["data"] == ["4.5 h", "9 h"]


def test_the_fit_chart_draws_the_window_the_fit_and_the_answer():
    real = realisations()
    estimate = ensemble.fit(real, 221.38, 500_000)
    options = pmfchart.fit_chart(real, estimate, pmp_aep=1_900_000)
    names = [series["name"] for series in options["series"]]
    assert names == ["realisations", "in the fit window", "fit, degree 1"]
    marks = options["series"][0]["markLine"]["data"]
    assert len(marks) == 3 and marks[0]["yAxis"] == 221.38
    assert sum(len(s["data"]) for s in options["series"][:2]) <= pmfchart.MAX_POINTS
