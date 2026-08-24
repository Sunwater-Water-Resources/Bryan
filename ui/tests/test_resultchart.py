"""The ECharts options for the duration plot.

The one that really matters is JSON-serialisability: a NaN reaches the browser
as invalid JSON and the chart silently renders nothing, which looks exactly
like 'no results' to whoever is using it.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest

from core import resultchart, results
from test_results import STANDARD_AEPS, level_curve, mc_row, write_quantiles


def build(tmp_path, durations, key="level", curves=None):
    rows = []
    for duration in durations:
        name = f"TFD_mc_{duration:g}h"
        values = (curves[duration] if curves else level_curve(duration))
        write_quantiles(Path(f"{tmp_path / name}_{key}.csv"), key, values)
        rows.append(mc_row(tmp_path, name, duration))
    frame = pd.DataFrame(rows)
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)[key])
    return comparison, results.analyse(comparison)


def test_the_option_is_json_serialisable(tmp_path):
    """NaN is not JSON. A gap has to reach the browser as null."""
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    {2: 100.0, 100: 120.0})
    write_quantiles(Path(f"{tmp_path / 'b_48h'}_level.csv"), "level",
                    {2: 101.0, 100: 119.0, 10000: 150.0})
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    analysis = results.analyse(comparison)

    option = resultchart.duration_chart(comparison, analysis, "level")
    text = json.dumps(option, allow_nan=False)          # raises on NaN
    assert "NaN" not in text

    series = {entry["name"]: entry for entry in option["series"]}
    # 24h never reached 1 in 10,000, so its last point is a null, not a value
    assert series["24h"]["data"][-1][1] is None


def test_level_is_linear_and_inflow_is_logarithmic(tmp_path):
    comparison, analysis = build(tmp_path, [24, 48])
    assert resultchart.duration_chart(
        comparison, analysis, "level")["yAxis"]["type"] == "value"
    assert resultchart.duration_chart(
        comparison, analysis, "inflow")["yAxis"]["type"] == "log"


def test_the_axis_ticks_are_the_aeps_that_were_evaluated(tmp_path):
    comparison, analysis = build(tmp_path, [24, 48])
    option = resultchart.duration_chart(comparison, analysis, "level")
    expected = [results.normal_variate(aep) for aep in STANDARD_AEPS]
    assert option["xAxis"]["axisLabel"]["customValues"] == pytest.approx(expected)
    # and the formatter can label any z, for an ECharts that ignores those
    formatter = option["xAxis"]["axisLabel"][":formatter"]
    assert '"2.3263": "100"' in formatter
    assert "0.3989422804014327" in formatter        # the tail approximation


def test_the_envelope_and_its_markup_follow_their_switches(tmp_path):
    comparison, analysis = build(tmp_path, [12, 24, 48, 120])

    bare = resultchart.duration_chart(comparison, analysis, "level",
                                      show_envelope=False)
    assert not any(entry["name"] == "envelope" for entry in bare["series"])

    plain = resultchart.duration_chart(comparison, analysis, "level",
                                       show_markup=False)
    envelope = [e for e in plain["series"] if e["name"] == "envelope"][0]
    assert "markArea" not in envelope

    marked = resultchart.duration_chart(comparison, analysis, "level")
    envelope = [e for e in marked["series"] if e["name"] == "envelope"][0]
    bands = envelope["markArea"]["data"]
    assert [band[0]["name"] for band in bands] == [b.label for b in analysis.bands]
    # each band is shaded in its own duration's colour, not the envelope's
    colours = resultchart.colour_for(comparison.frame.columns)
    assert all(band[0]["itemStyle"]["color"] == colours[band[0]["name"]]
               for band in bands)
    assert "→" in envelope["markPoint"]["data"][0]["value"]


def test_a_noise_switch_gets_no_pin(tmp_path):
    """Pinning a hop that came out of sampling noise dresses it up as a finding."""
    comparison, analysis = build(tmp_path, [24, 48], curves={
        24: {2: 216.00, 100: 220.00},
        48: {2: 216.01, 100: 219.98},      # 10 mm apart - inside the floor
    })
    assert analysis.switches                       # idxmax did hop
    envelope = [e for e in resultchart.duration_chart(comparison, analysis, "level")
                ["series"] if e["name"] == "envelope"][0]
    assert "markPoint" not in envelope


def test_the_critical_duration_chart_falls_for_a_level_curve(tmp_path):
    comparison, analysis = build(tmp_path, [12, 24, 48, 120])
    option = resultchart.critical_duration_chart(comparison, analysis)
    durations = [point[1] for point in option["series"][0]["data"]]
    assert durations[0] == 120 and durations[-1] == 12
    assert durations == sorted(durations, reverse=True)
    json.dumps(option, allow_nan=False)


def test_an_empty_selection_makes_an_empty_chart():
    empty = results.Comparison()
    option = resultchart.duration_chart(empty, results.analyse(empty), "level")
    assert option == {"series": []}
    assert resultchart.critical_duration_chart(empty, results.analyse(empty)) == {}


def test_curves_without_a_duration_get_no_critical_duration_chart(tmp_path):
    """A file added by hand may have no duration - the second chart needs one."""
    write_quantiles(Path(f"{tmp_path / 'odd'}_level.csv"), "level", {2: 1.0, 10: 2.0})
    frame = pd.DataFrame([mc_row(tmp_path, "odd", None)])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    assert resultchart.critical_duration_chart(
        comparison, results.analyse(comparison)) == {}


def test_a_non_positive_value_does_not_break_a_log_axis(tmp_path):
    comparison, analysis = build(tmp_path, [24], key="inflow",
                                 curves={24: {2: 0.0, 100: 500.0}})
    option = resultchart.duration_chart(comparison, analysis, "inflow")
    json.dumps(option, allow_nan=False)
    assert option["series"][0]["data"][0][1] == 0.0
