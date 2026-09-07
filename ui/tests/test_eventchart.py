"""The AEP neutrality plot, and the two marks a level loading carries.

A level loading has two honest answers to 'where is it on this axis': the AEP
the design curve gives it, and the AEP this run's own realisations reached it
at. They differ - the curve is the envelope over the durations, read as a
straight line between the standard AEPs, and smoothed on the way - so the plot
draws both. Marking only the first is what makes the line look offset from the
level it was asked for.

The chart is plain dicts, so this builds an outcome by hand rather than running
a project through the page.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from core import eventchart
from core.events import Outcome
from lib.RepresentativeEvents import Target


def candidates():
    """Three realisations, with the columns the chart reads."""
    return pd.DataFrame(
        {"z_result": [3.0, 3.2, 2.8], "z_rain": [3.0, 2.6, 3.1],
         "rain_aep": [1000.0, 200.0, 1200.0],
         "result_aep": [1000.0, 1400.0, 800.0],
         "delta_z": [0.05, 0.6, 0.25],
         "level": [220.0, 220.4, 219.6],
         "flags": [(), ("high antecedent storage (z = +2.00)",), ()]},
        index=[0, 1, 2])


class _Ranking:
    def __init__(self, frame):
        self.candidates = frame
        self.excluded = {}
        self.notes = []


def outcome(**overrides):
    target = Target(kind="level", value=220.0, result_type="level")
    settings = dict(target=target, aep=1000.0, target_value=220.0,
                    ranking=_Ranking(candidates()))
    settings.update(overrides)
    return Outcome(**settings)


def series_named(chart, name):
    return next((series for series in chart["series"]
                 if series.get("name") == name), None)


def test_the_run_gets_its_own_mark_when_it_disagrees_with_the_curve():
    chart = eventchart.neutrality_chart(outcome(data_z=3.25), "level")
    drawn = series_named(chart, "220 m AHD in this run")
    assert drawn is not None
    assert [point[0] for point in drawn["data"]] == [3.25, 3.25]
    assert drawn["lineStyle"]["type"] == "dotted"


def test_the_design_aep_is_still_marked():
    """Both, not one instead of the other - the design AEP is the loading."""
    chart = eventchart.neutrality_chart(outcome(data_z=3.25), "level")
    marks = chart["series"][0]["markLine"]["data"]
    assert any(abs(mark.get("xAxis", 0) - 3.0902) < 1e-3 for mark in marks)


def test_the_run_mark_is_inside_the_axis():
    """Off the end of the range it would be drawn where nobody can see it."""
    chart = eventchart.neutrality_chart(outcome(data_z=4.6), "level")
    assert chart["xAxis"]["min"] <= 4.6 <= chart["xAxis"]["max"]


def test_an_aep_loading_carries_no_second_mark():
    """Nothing to disagree about: the loading is the AEP the line is drawn at."""
    plain = outcome(target=Target(kind="aep", value=1000, result_type="level"),
                    target_value=None, data_z=None)
    chart = eventchart.neutrality_chart(plain, "level")
    assert not [series for series in chart["series"]
                if str(series.get("name", "")).endswith("in this run")]
