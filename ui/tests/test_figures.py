"""Report figures: the data a figure is drawn from, and the PNG drawn from it.

Pinned: a curve's label belongs to the figure (the same group can be 'URBS' in
one and 'GWL 1.3 °C' in another); the three kinds of curve read what they claim
to; the preview and the export carry the same series; and the real script, where
Bryan's interpreter is available, draws the PNG and leaves the job beside it.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from core import figurechart, figures
from core import reporttables as rt
from core import study as studies
from report_fixtures import GROUP, LEVEL, build_study
from test_real_bryan import BRYAN_PYTHON, needs_bryan


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def write_ffa(path):
    """An RMC Bestfit export: x columns are exceedance probabilities."""
    probabilities = [0.5, 0.1, 0.01, 0.001]
    frame = pd.DataFrame({
        "Posterior Mode_x": probabilities, "Posterior Mode_y": [100, 500, 2000, 5000],
        "Posterior Mean_x": probabilities, "Posterior Mean_y": [110, 560, 2400, 6600],
        "90% Credible Intervals_x": probabilities, "90% Credible Intervals_y": [150, 800, 4000, 12000],
        "90% Credible Intervals_x2": probabilities, "90% Credible Intervals_y2": [80, 300, 1100, 2000],
        "Exact Data_x": [0.6, 0.2, 0.05, None], "Exact Data_y": [90, 400, 1500, None],
        "Interval Data_x": [0.004, None, None, None], "Interval Data_yLower": [3100, None, None, None],
    })
    frame.to_csv(path)
    return path


def level_figure(**overrides):
    spec = figures.new_spec(filename="levels", type="level", curves=[
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "GWL 1.3 °C"}],
        reference_levels=[{"label": "DCL", "level": 219.13, "colour": "red"}])
    spec.update(overrides)
    return spec


def test_a_group_curve_is_its_envelope_over_the_durations(study):
    data = figures.build(study, level_figure())
    assert not data.problems
    (series,) = data.series
    envelope = dict(series.points)
    assert envelope[5.0] == max(LEVEL[36][1], LEVEL[72][1])
    assert envelope[1_000.0] == max(LEVEL[36][4], LEVEL[72][4])


def test_the_label_belongs_to_the_figure_not_the_curve(study):
    figures.put(study, level_figure())
    figures.put(study, level_figure(filename="concordance", curves=[
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "URBS"}]))
    labels = [figures.build(study, spec).series[0].label for spec in figures.figures(study)]
    assert labels == ["GWL 1.3 °C", "URBS"]


def test_the_aep_range_trims_every_curve(study):
    data = figures.build(study, level_figure(min_aep=5, max_aep=10_000))
    aeps = [aep for aep, _ in data.series[0].points]
    assert min(aeps) == 5 and max(aeps) == 10_000


def test_a_file_curve_reads_the_first_column_as_the_aep(study, tmp_path):
    path = tmp_path / "adopted_2020.csv"
    pd.DataFrame({"": [5, 10, 100], "level": [216.38, 216.48, 217.0]}).to_csv(path, index=False)
    spec = level_figure(curves=[{"kind": figures.FILE, "file": str(path), "label": "Sunwater 2020"}])
    data = figures.build(study, spec)
    assert data.series[0].points == [(5.0, 216.38), (10.0, 216.48), (100.0, 217.0)]
    spec["curves"][0]["column"] = "outflow"
    assert "has no 'outflow' column" in figures.build(study, spec).problems[0]


def test_an_ffa_brings_its_posteriors_band_maxima_and_paleofloods(study, tmp_path):
    path = write_ffa(tmp_path / "ffa.csv")
    spec = figures.new_spec(type="inflow", max_aep=2000, curves=[
        {"kind": figures.FFA, "file": str(path), "label": "LPIII", "posterior": "Both"},
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "URBS"}])
    data = figures.build(study, spec)
    names = [(series.label, series.style) for series in data.series]
    assert names == [("AMS", "ams"), ("90% credible interval", "ffa-ci"),
                     ("90% credible interval", "ffa-ci"), ("LPIII: Posterior Mode", "ffa-mode"),
                     ("LPIII: Posterior Mean", "ffa-mean"),
                     ("Paleoflood (lower bound)", "paleo"), ("URBS", "line")]
    mode = dict(data.series[3].points)
    assert mode[10.0] == 500                                  # 1 in 10 from x = 0.1
    assert dict(data.series[5].points) == {250.0: 3100}       # 1 in 250 from x = 0.004
    # the 0.6 maximum is more frequent than 1 in 2 and trimmed with the axis
    assert [aep for aep, _ in data.series[0].points] == [5.0, 20.0]


def test_two_ffas_leave_the_band_off_because_two_cannot_be_told_apart(study, tmp_path):
    path = write_ffa(tmp_path / "ffa.csv")
    spec = figures.new_spec(type="inflow", curves=[
        {"kind": figures.FFA, "file": str(path), "label": "a"},
        {"kind": figures.FFA, "file": str(path), "label": "b"}])
    assert "ffa-ci" not in {series.style for series in figures.build(study, spec).series}


def test_a_file_that_is_not_an_ffa_says_so(study, tmp_path):
    path = tmp_path / "not.csv"
    pd.DataFrame({"a": [1]}).to_csv(path)
    spec = figures.new_spec(curves=[{"kind": figures.FFA, "file": str(path), "label": "x"}])
    assert "RMC Bestfit" in figures.build(study, spec).problems[0]


def test_the_preview_and_the_job_carry_the_same_series(study, tmp_path):
    data = figures.build(study, level_figure(aep_of_pmp=1_900_000))
    options = figurechart.preview(data)
    job = data.to_job(tmp_path / "x.png")
    assert [s["name"] for s in options["series"]] == [s["label"] for s in job["series"]]
    assert len(options["series"][0]["data"]) == len(job["series"][0]["points"])
    marks = options["series"][0]["markLine"]["data"]
    assert marks[0]["yAxis"] == 219.13 and "xAxis" in marks[1]
    assert options["yAxis"]["type"] == "value" and job["log_y"] is False
    json.dumps(options)                                       # JSON-safe for the browser


def test_figures_are_kept_in_the_study_and_named_uniquely(study):
    first = figures.put(study, level_figure())
    second = figures.put(study, level_figure())
    study.save()
    reread = studies.load_study(study.path)
    assert [spec["id"] for spec in figures.figures(reread)] == [first["id"], second["id"]]
    assert second["id"] == "levels-2"


def test_renaming_a_run_does_not_orphan_a_figure(study):
    figures.put(study, level_figure())
    study.rename_run("E099 RFSL", "E100 RFSL")
    assert figures.figures(study)[0]["curves"][0]["run"] == "E100 RFSL"


@pytest.mark.parametrize("text, parsed", [
    ("RFSL = 215.5 = yellowgreen", [{"label": "RFSL", "level": 215.5, "colour": "yellowgreen"}]),
    ("DCL = 219.13", [{"label": "DCL", "level": 219.13}]),
    ("no level here", [])])
def test_reference_levels_parse_with_or_without_a_colour(text, parsed):
    assert figures.parse_reference_levels(text) == parsed


def test_standard_aeps_step_as_v03_did():
    assert figures.standard_aeps(2, 2000) == [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]


@needs_bryan
def test_the_script_draws_the_png_and_leaves_the_job_beside_it(study, tmp_path):
    data = figures.build(study, level_figure(aep_of_pmp=1_900_000))
    result = figures.export(data, tmp_path / "out" / "levels.png", BRYAN_PYTHON)
    assert result.ok, result.output
    job = json.loads((tmp_path / "out" / "levels.json").read_text(encoding="utf-8"))
    assert job["series"][0]["label"] == "GWL 1.3 °C"
    assert (tmp_path / "out" / "levels.png").stat().st_size > 10_000
