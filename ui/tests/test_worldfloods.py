"""The world's maximum observed floods, and the PMF and PMPF laid against them.

Pinned: Table 1 of Li et al. (2013) as transcribed - its count, its surveyed and
modified flags, its extremes - and the WMO (2009) envelope, eq. (2), drawn as
written; the PMF is the largest ensemble inflow and the PMPF the inflow envelope
at the AEP of the PMP.
"""

from __future__ import annotations

import math

import pytest

from core import ensemble, pmfchart, reporttables as rt, study as studies
from core import worldfloods as wf
from report_fixtures import GROUP, PMF_GROUP, build_study


def test_table_1_is_all_there():
    floods = wf.WORLD_FLOODS
    assert len(floods) == 54
    assert [f.number for f in floods] == list(range(1, 55))
    assert sorted(f.number for f in floods if f.surveyed) == [1, 4, 6, 14, 17, 19, 29,
                                                              32, 35, 49, 50]
    amazon = floods[-1]
    assert amazon.modified and (amazon.area_km2, amazon.peak_m3s) == (4_640_300, 303_000)
    assert (floods[0].area_km2, floods[0].peak_m3s) == (0.9, 159)
    # the areas are in order, as the table lists them
    assert [f.area_km2 for f in floods] == sorted(f.area_km2 for f in floods)


def test_the_wmo_envelope_is_equation_2_as_written():
    assert wf.wmo_envelope(100) == pytest.approx(154 * 100 ** 0.738)
    assert wf.wmo_envelope(1_000) == pytest.approx(1830 * 1_000 ** 0.316)
    assert wf.wmo_envelope(299.9) < wf.wmo_envelope(300)          # it steps at 300 km²
    assert math.isnan(wf.wmo_envelope(3_000_000)) and math.isnan(wf.wmo_envelope(0))
    # a fitted envelope: Banqiao (1975) is under it, Sunjiazhai a little over
    assert wf.ratio_to_envelope(762, 13_100) == pytest.approx(0.88, abs=0.01)
    assert wf.ratio_to_envelope(0.9, 159) == pytest.approx(1.12, abs=0.01)


def test_the_herschy_envelope_is_equation_1():
    assert wf.herschy_envelope(10) == pytest.approx(100 * 10 ** 0.8)
    assert wf.herschy_envelope(1_000) == pytest.approx(500 * 1_000 ** 0.43)


def test_the_envelope_curve_keeps_its_step():
    curve = wf.envelope_curve(wf.wmo_envelope)
    areas = [a for a, _ in curve]
    assert areas == sorted(areas) and 300.0 in areas
    below = next(q for a, q in curve if a == 299.999)
    at = next(q for a, q in curve if a == 300.0)
    assert at > below


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    section = ensemble.settings(study)
    section["groups"].append({"label": "Near-term", "ensemble": {"run": "E099 PMF",
                                                                 "group": PMF_GROUP},
                              "mc": {"run": "E099 RFSL", "group": GROUP}})
    section["catchment_area_km2"] = 1_250
    ensemble.store(study, section)
    return study


def test_the_pmf_is_the_largest_inflow_and_the_pmpf_the_envelope_at_the_pmp(study):
    section = ensemble.settings(study)
    (flood,) = wf.dam_floods(study, section)
    assert not flood.problems
    assert flood.pmf == 18_000             # event 0, the largest inflow
    assert flood.pmf_level_event == 17_500  # event 2, which gives the PMF level
    assert flood.pmpf == 16_000            # the 36 h inflow at 1 in 1,900,000


def test_a_pmpf_at_an_aep_the_run_did_not_reach_says_so(study):
    section = ensemble.settings(study)
    section["pmp_aep"] = 3_000_000
    (flood,) = wf.dam_floods(study, section)
    assert math.isnan(flood.pmpf)
    assert any("no inflow quantile at 1 in 3,000,000" in p for p in flood.problems)


def test_the_settings_are_kept_in_the_study(study):
    study.save()
    stored = ensemble.settings(studies.load_study(study.path))
    assert stored["catchment_area_km2"] == 1_250 and stored["world_herschy"] is False


def test_the_chart_places_the_dam_floods_at_the_catchment_area(study):
    floods = wf.dam_floods(study, ensemble.settings(study))
    options = pmfchart.world_floods_chart(wf.WORLD_FLOODS, 1_250, floods)
    names = [s["name"] for s in options["series"]]
    assert names == ["World maximum floods", "World maximum floods (surveyed)",
                     wf.WMO_LABEL, "PMF", "PMPF"]
    assert options["xAxis"]["type"] == options["yAxis"]["type"] == "log"
    assert options["animation"] is False                 # a copy is the finished chart
    pmf = options["series"][3]["data"][0]["value"]
    assert pmf == [1_250, 18_000]
    assert len(options["series"][0]["data"]) + len(options["series"][1]["data"]) == 54
    assert options["series"][2]["markLine"]["data"] == [{"xAxis": 1_250}]
    # the legend stands beside the plot, never under the area axis, and names every entry
    assert options["legend"]["orient"] == "vertical" and options["legend"]["data"] == names
    assert options["grid"]["right"] == pmfchart.legend_width(names)


def test_a_legend_too_long_for_one_column_takes_the_room_of_two():
    one = pmfchart.legend_width(["x" * 10] * 19)
    two = pmfchart.legend_width(["x" * 10] * 20)
    assert two - pmfchart._LEGEND_PAD == 2 * (one - pmfchart._LEGEND_PAD)
    assert pmfchart.legend_width(["x" * 20]) > pmfchart.legend_width(["x" * 10])


def test_several_groups_are_told_apart_by_label_and_without_an_area_none_are_placed(study):
    one = wf.DamFlood("Near-term", pmf=18_000, pmpf=16_000)
    two = wf.DamFlood("Long-term", pmf=21_000, pmpf=19_000)
    names = [s["name"] for s in
             pmfchart.world_floods_chart(wf.WORLD_FLOODS, 1_250, [one, two])["series"]]
    assert names[3:] == ["Near-term PMF", "Near-term PMPF", "Long-term PMF", "Long-term PMPF"]
    bare = pmfchart.world_floods_chart(wf.WORLD_FLOODS, None, [one, two], herschy=True)
    assert [s["name"] for s in bare["series"]][3:] == [wf.HERSCHY_LABEL]


def test_the_copy_script_renders_the_chart_to_a_png():
    copy = pmfchart.copy_chart_script(42)
    assert "getElement(42)" in copy and "getDataURL" in copy
    assert "'image/png'" in copy and "ClipboardItem" in copy
    # entries switched off are cut from the legend for the render and put back after it
    assert "selected[nameOf(item)] !== false" in copy
    assert copy.index("chart.setOption(change)") < copy.index("getDataURL") \
        < copy.index("chart.setOption(restore)")
    save = pmfchart.copy_chart_script(42, "callide_world_floods.png")
    assert '"callide_world_floods.png"' in save and "link.download" in save
