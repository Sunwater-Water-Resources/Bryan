"""Rebuilding a realisation's rainfall from what the run wrote down.

The rebuild has to reproduce ``MonteCarloSimulator.run_models`` exactly - the
same filter buffers, the same order, the same pre-burst pattern - or the
hyetograph belongs to a storm nobody ran. The rainfall data itself lives outside
the repo, so what is pinned here is everything around it: the percent-to-mm
conversion, the pre-burst going on at negative times, the arguments the storm is
asked for, and the checks against the recorded depths actually failing when the
rebuild is wrong.

Bryan's own dependencies needed - lib/EventStorm.py drives StormGenerator.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

from lib import EventStorm                                          # noqa: E402

DURATION = 24.0
TIMESTEP = 1.0

# 24 hours of main burst, as a percentage of the burst depth - the units every
# temporal pattern in Bryan is in.
BURST = pd.Series([2.0] * 20 + [15.0, 20.0, 15.0, 10.0],
                  index=[float(hour) for hour in range(24)])


class StubStorm:
    """Stands in for StormBurst, and records how it was called."""

    def __init__(self, ave_rain=200.0, comment="No embedded bursts"):
        self.timesteps = TIMESTEP
        self.rainfall = self
        self.preburst_patterns = self
        self._ave_rain = ave_rain
        self._comment = comment
        self.calls = {}

    # -- the rainfall depths
    def get_depth_z(self, z, duration, storm_method, spatial_method=None):
        self.calls["get_depth_z"] = (z, duration, storm_method, spatial_method)
        return pd.Series([self._ave_rain], index=["sub_1"])

    def get_average_rain(self, depths, print_msg=True):
        return float(depths.iloc[0])

    # -- the pattern
    def get_temporal_pattern(self, storm_method, duration, tp_sample, rain_sample_z):
        self.calls["get_temporal_pattern"] = (storm_method, duration, tp_sample,
                                              rain_sample_z)
        return BURST.copy()

    def filter_temppat(self, pattern, z, duration, ave_rain, storm_method,
                       buffer=None, climate_adjustment=None):
        self.calls["filter_temppat"] = {"buffer": buffer,
                                        "climate_adjustment": climate_adjustment}
        return pattern, self._comment

    # -- the pre-burst
    def get_preburst_pattern(self, duration, proportion, timesteps, sample=None):
        self.calls["get_preburst_pattern"] = (duration, proportion, timesteps, sample)
        # Six hours of pre-burst, totalling proportion x 100 percent.
        hours = 6
        each = proportion * 100.0 / hours
        index = [float(-hours + step) for step in range(hours)]
        return pd.Series([each] * hours, index=index), sample

    def filter_embedded_bursts_in_preburst(self, preburst, pattern, duration, z,
                                           ave_rain, storm_method, buffer=None,
                                           climate_adjustment=None):
        self.calls["filter_preburst"] = {"buffer": buffer}
        return preburst, "No embedded bursts;"

    def uniform_preburst(self, duration, proportion, timesteps):
        self.calls["uniform_preburst"] = (duration, proportion, timesteps)
        index = [float(-4 + step) for step in range(4)]
        return pd.Series([proportion * 100.0 / 4] * 4, index=index)


def context(**overrides):
    settings = dict(storm_config="storm.json", focal_subcatchments="focal.csv",
                    duration=DURATION)
    settings.update(overrides)
    return EventStorm.StormContext(**settings)


def realisation(**overrides):
    row = {"rain_z": 3.0902, "tp": 4, "storm_method": "ARR point",
           "preburst_proportion": 0.3, "preburst_tp": 7,
           "mean_rain_mm": 200.0, "preburst_mm": 60.0,
           "embedded_bursts": "No embedded bursts"}
    row.update(overrides)
    return pd.Series(row)


def build(storm=None, ctx=None, sim=None):
    return EventStorm.hyetograph(storm or StubStorm(), None,
                                 ctx or context(), sim if sim is not None else realisation())


# -- the depths --------------------------------------------------------------

def test_a_percentage_pattern_becomes_millimetres():
    """Patterns are percentages of the burst depth - uniform_preburst says so."""
    hyeto = build()
    assert hyeto.burst_mm == pytest.approx(200.0)


def test_the_preburst_is_a_proportion_of_the_burst():
    hyeto = build()
    assert hyeto.preburst_mm == pytest.approx(60.0)          # 0.3 x 200 mm
    assert hyeto.total_mm == pytest.approx(260.0)


def test_the_preburst_runs_at_negative_times():
    """t = 0 is the start of the main burst, so events with different
    pre-burst durations can be plotted against each other."""
    hyeto = build()
    assert hyeto.depths.index.min() == -6.0
    assert hyeto.preburst_hours == 6.0
    assert hyeto.depths[hyeto.depths.index < 0].size == 6


def test_intensity_is_the_depths_over_the_timestep():
    hyeto = build()
    assert hyeto.intensity.max() == pytest.approx(hyeto.depths.max() / TIMESTEP)


def test_no_preburst_leaves_the_burst_alone():
    hyeto = build(sim=realisation(preburst_proportion=0.0, preburst_mm=0.0))
    assert hyeto.preburst_hours == 0.0
    assert hyeto.depths.index.min() >= 0
    assert hyeto.trustworthy


# -- replaying the run exactly ------------------------------------------------

def test_the_burst_filter_uses_the_run_s_buffer():
    """Simulator.run_models hard-codes 1.1. A different buffer is a different
    pattern with the same total, which no check on depths would catch."""
    storm = StubStorm()
    build(storm)
    assert storm.calls["filter_temppat"]["buffer"] == EventStorm.BURST_FILTER_BUFFER
    assert EventStorm.BURST_FILTER_BUFFER == 1.1


def test_the_preburst_filter_uses_a_different_buffer():
    """0.9, not 1.1 - and util/StormInstance.py omits it entirely."""
    storm = StubStorm()
    build(storm)
    assert storm.calls["filter_preburst"]["buffer"] == EventStorm.PREBURST_FILTER_BUFFER
    assert EventStorm.PREBURST_FILTER_BUFFER == 0.9


def test_the_recorded_pattern_number_is_the_one_rebuilt():
    storm = StubStorm()
    build(storm, sim=realisation(tp=9))
    assert storm.calls["get_temporal_pattern"][2] == 9


def test_the_recorded_preburst_pattern_is_the_one_rebuilt():
    storm = StubStorm()
    build(storm, sim=realisation(preburst_tp=3))
    assert storm.calls["get_preburst_pattern"][3] == 3


def test_excluding_the_filter_skips_it():
    storm = StubStorm()
    build(storm, ctx=context(exclusions="ebf,d50"))
    assert "filter_temppat" not in storm.calls


def test_excluding_the_preburst_leaves_the_burst_alone():
    storm = StubStorm()
    hyeto = build(storm, ctx=context(exclusions="pb"),
                  sim=realisation(preburst_mm=0.0))
    assert "get_preburst_pattern" not in storm.calls
    assert hyeto.preburst_hours == 0.0


def test_the_uniform_preburst_method_is_honoured():
    storm = StubStorm()
    hyeto = build(storm, ctx=context(preburst_method="uniform"))
    assert "uniform_preburst" in storm.calls
    assert "get_preburst_pattern" not in storm.calls
    assert hyeto.preburst_mm == pytest.approx(60.0)


# -- the checks against what the run recorded ---------------------------------

def test_a_matching_rebuild_is_trustworthy():
    hyeto = build()
    assert hyeto.trustworthy
    assert not hyeto.problems
    assert len(hyeto.checks) == 3          # depth, pre-burst, embedded bursts


def test_a_depth_that_does_not_match_is_reported():
    """The run recorded 250 mm; the rebuild gives 200. Something differs."""
    hyeto = build(sim=realisation(mean_rain_mm=250.0))
    assert not hyeto.trustworthy
    assert any("catchment average" in problem for problem in hyeto.problems)


def test_a_preburst_that_does_not_match_is_reported():
    hyeto = build(sim=realisation(preburst_mm=10.0))
    assert not hyeto.trustworthy
    assert any("pre-burst" in problem for problem in hyeto.problems)


def test_a_different_embedded_burst_outcome_is_reported():
    """The filter is where a rebuild most easily diverges, and the comment is
    the run's own record of what it did."""
    storm = StubStorm(comment="Embedded burst filtered")
    hyeto = build(storm)
    assert not hyeto.trustworthy
    assert any("embedded burst" in problem for problem in hyeto.problems)


def test_rounding_in_the_mcdf_does_not_count_as_a_mismatch():
    hyeto = build(sim=realisation(mean_rain_mm=200.4, preburst_mm=60.1))
    assert hyeto.trustworthy


def test_a_row_without_the_recorded_depths_makes_no_claim():
    sim = realisation()
    hyeto = build(sim=sim.drop(["mean_rain_mm", "preburst_mm", "embedded_bursts"]))
    assert hyeto.checks == []
    assert hyeto.trustworthy          # nothing to contradict, and nothing claimed


# -- the exclusion keys -------------------------------------------------------

def test_only_the_two_exclusions_that_shape_the_storm_are_read():
    parsed = EventStorm._parse_exclusions("ebf, pb, d50, ru")
    assert parsed == {"embedded_burst_filter": True, "preburst": True}


def test_a_blank_exclusions_cell_excludes_nothing():
    for blank in ("", None, float("nan")):
        parsed = EventStorm._parse_exclusions(blank)
        assert parsed == {"embedded_burst_filter": False, "preburst": False}
