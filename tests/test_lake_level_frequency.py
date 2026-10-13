"""The lake level frequency curves, and util/LakeLevelFrequency.py end to end.

The synthetic record comes from ``ui/tests/lake_fixtures.py``, shared with the
page's tests. The last test reproduces the Callide figure the analysis was
developed on, and runs only where the private Callide study checkout sits beside
this one.
"""

from __future__ import annotations

import glob
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
for path in (str(BRYAN_ROOT), str(BRYAN_ROOT / "ui" / "tests")):
    if path not in sys.path:
        sys.path.insert(0, path)

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

from lake_fixtures import FSL, synthetic_level, write_hydstra, write_mcdf  # noqa: E402
from lib import LakeLevelFrequency as frequency  # noqa: E402
from lib import LakeLevelRecord as record  # noqa: E402

CALLIDE = BRYAN_ROOT.parent / "callide-fsl-reinstate"


def cli():
    spec = importlib.util.spec_from_file_location(
        "LakeLevelFrequencyCli", BRYAN_ROOT / "util" / "LakeLevelFrequency.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def positions(years=40, seed=7):
    return record.with_positions(record.annual_maxima(synthetic_level(years, seed=seed)))


# -- the shouldered plateau ----------------------------------------------------------

def test_the_plateau_runs_from_the_first_maximum_on_it_to_the_first_real_jump():
    z = np.arange(8, dtype=float)
    level = np.array([210, 213, 215.0, 215.01, 215.06, 215.12, 215.6, 216.0])
    assert frequency.plateau_span(z, level, 215.0) == (2.0, 5.0)
    with pytest.raises(ValueError, match="no plateau"):
        frequency.plateau_span(z, level, 220.0)


def test_a_curve_built_from_the_form_is_recovered_exactly():
    params = {"coef": [2.0, -0.5, 0.0, 0.0], "z1": 0.0, "z2": 0.6,
              "upper_coef": [215.74, 0.9], "upper_join": "free", "fsl": 215.0}
    # The grid has to hold z1 and z2 themselves: the plateau's extent is read off
    # the maxima, never fitted.
    z = np.round(np.arange(-2.0, 2.55, 0.1), 10)
    level = frequency.shouldered_plateau(z, params)
    got, rmse = frequency.fit_shouldered_plateau(z, level, 215.0, degree=2)
    assert rmse < 1e-6
    assert got["upper_coef"] == pytest.approx([215.74, 0.9])
    assert got["step"] == pytest.approx(0.74)
    assert got["coef"][:2] == pytest.approx([2.0, -0.5], abs=1e-6)


def wet_dam(years=45, spilling_from=-0.9, fsl=215.0):
    """A dam that spills in most years: a steep shoulder below full supply, and
    above it a curve that flattens as the spillway takes over."""
    _, z = record.plotting_positions(np.arange(years, dtype=float))
    z = np.sort(z)
    d = z - spilling_from
    level = np.where(d < 0, fsl + 4.0 * d, fsl + 1.6 * d - 0.25 * d ** 2)
    return z, level


def test_no_plateau_meets_at_full_supply_between_the_maxima_either_side():
    z, level = wet_dam()
    params, rmse = frequency.fit_shouldered_plateau(
        z, level, 215.0, degree=1, tol=0, upper_degree=2, upper_join="fsl")
    assert params["z1"] == params["z2"]
    assert z[level <= 215.0].max() < params["z1"] < z[level > 215.0].min()
    assert rmse < 0.03
    with pytest.raises(ValueError, match="both sides"):
        frequency.fit_shouldered_plateau(z, level + 10, 215.0, tol=0)


def test_a_curved_upper_limb_is_only_allowed_where_enough_years_spill():
    z, level = wet_dam()
    linear = frequency.fit_shouldered_plateau(z, level, 215.0, degree=1, tol=0)[1]
    curved = frequency.fit_shouldered_plateau(z, level, 215.0, degree=1, tol=0,
                                              upper_degree=2)[1]
    assert curved < linear / 3

    # A dry-belt record: a handful above full supply.
    z, level = wet_dam(spilling_from=1.2)
    assert frequency.fit_shouldered_plateau(z, level, 215.0, degree=1, tol=0)[1] < 0.1
    with pytest.raises(ValueError, match="needs at least 12"):
        frequency.fit_shouldered_plateau(z, level, 215.0, degree=1, tol=0, upper_degree=2)


def test_the_upper_limb_can_start_at_full_supply_and_never_turns_back():
    z, level = wet_dam()
    params, _ = frequency.fit_shouldered_plateau(
        z, level, 215.0, degree=1, tol=0, upper_degree=3, upper_join="fsl")
    assert params["step"] == pytest.approx(0.0)
    grid = np.linspace(params["z2"], z.max(), 200)
    assert np.diff(frequency.shouldered_plateau(grid, params)).min() >= -1e-9


def test_too_few_maxima_above_the_plateau_is_a_reason_not_a_crash():
    z = np.linspace(-2, 2, 20)
    level = np.minimum(215.0, 210 + 3 * z)
    with pytest.raises(ValueError, match="above the plateau"):
        frequency.fit_shouldered_plateau(z, level, 215.0)


def test_both_forms_fit_the_synthetic_record():
    table = positions()
    _, shouldered = frequency.fit("shouldered", table["z"], table["level"], FSL)
    _, logistic = frequency.fit("logistic", table["z"], table["level"])
    assert 0 < shouldered < 0.5
    assert 0 < logistic < 0.6


def test_resampling_is_reproducible_and_says_how_many_draws_it_kept():
    table = positions()
    grid = np.linspace(-0.3, 1.8, 20)
    first = frequency.bootstrap_curves(table["level"], grid, "shouldered", FSL, draws=40)
    again = frequency.bootstrap_curves(table["level"], grid, "shouldered", FSL, draws=40)
    assert np.array_equal(first, again)
    assert 0 < len(first) <= 40

    storm = frequency.bootstrap_censored_curves(
        table["level"], table["carried_over"], grid, "logistic", draws=30)
    assert storm.shape[1] == len(grid)


def test_a_band_from_too_few_fitted_resamples_is_flagged():
    table = positions()
    out = frequency.analyse(table, fsl=FSL, form="shouldered", draws=40,
                            storm_driven=False, progress=lambda *_: None)
    block = out["fits"]["all"]
    assert block["draws"] == 40
    assert (block["warning"] is not None) == (block["draws_used"] < 28)


def test_no_resamples_fits_the_curves_without_bands():
    """The page's Fit curves: seconds, so settings can be tried on the curve alone."""
    table = positions()
    out = frequency.analyse(table, fsl=FSL, form="shouldered", draws=0,
                            progress=lambda *_: None)
    banded = frequency.analyse(table, fsl=FSL, form="shouldered", draws=20,
                               progress=lambda *_: None)
    for name in ("all", "storm"):
        block = out["fits"][name]
        assert block["curve"] == banded["fits"][name]["curve"]     # the same fit
        assert block["rmse"] == banded["fits"][name]["rmse"]
        assert block["band_lo"] is None and block["band_hi"] is None
        assert block["draws_used"] == 0 and block["warning"] is None


def test_the_analysis_records_a_form_it_could_not_fit_instead_of_failing():
    table = positions()
    out = frequency.analyse(table, fsl=FSL + 5, form="shouldered", draws=10,
                            progress=lambda *_: None)
    assert "no plateau" in out["fits"]["all"]["error"]
    assert out["fits"]["all"]["curve"] is None


def test_the_envelope_is_the_worst_duration_at_each_probability(tmp_path):
    low = frequency.read_design_curve(write_mcdf(tmp_path / "a.csv"))
    high = frequency.read_design_curve(write_mcdf(tmp_path / "b.csv", durations_shift=0.3))
    aep = np.array([0.5, 0.01, 1e-4])
    envelope = frequency.design_envelope({6: low, 24: high}, aep)
    assert envelope == pytest.approx(np.interp(np.log(aep), np.log(high[0]), high[1]))


# -- the script ---------------------------------------------------------------------

def job_file(tmp_path, **overrides):
    record_path = write_hydstra(tmp_path / "HW.csv", synthetic_level(years=30))
    mcdfs = [write_mcdf(tmp_path / f"run_{hours}h__mcdf.csv", durations_shift=shift)
             for hours, shift in ((12, 0.0), (48, 0.2))]
    job = {"record": {"files": [str(record_path)]}, "fsl": FSL, "fsl_label": "FSL",
           "reference_levels": [{"label": "crest", "level": 218.0}],
           "fit": {"draws": 30},
           "design": {"include": True,
                      "sources": [{"duration": 12, "path": str(mcdfs[0])},
                                  {"duration": 48, "path": str(mcdfs[1])}]}}
    job.update(overrides)
    path = tmp_path / "job.json"
    path.write_text(json.dumps(job))
    return path


def test_the_script_writes_results_a_figure_either_way_and_the_series(tmp_path, capsys):
    script = cli()
    job = job_file(tmp_path)
    results = tmp_path / "cache" / "results.json"
    assert script.main(["--job", str(job), "--results", str(results),
                        "--png", str(tmp_path / "with.png"),
                        "--ams-csv", str(tmp_path / "ams.csv")]) == 0

    saved = json.loads(results.read_text())
    assert saved["fits"]["all"]["rmse"] > 0 and saved["fits"]["all"]["draws_used"] > 0
    assert set(saved["design"]["durations"]) == {"12", "48"}
    assert len(saved["design"]["envelope"]) == len(saved["grid"]["z"])
    assert (tmp_path / "with.png").stat().st_size > 10_000
    settings = json.loads((tmp_path / "with.json").read_text())
    assert settings["job"]["fit"]["upper_join"] == "free"
    assert settings["design_floods"]["durations"] == [12.0, 48.0]
    assert settings["fits"]["all"]["rmse"] == saved["fits"]["all"]["rmse"]
    csv = (tmp_path / "ams.csv").read_text()
    assert csv.startswith("# source: ") and "water_year,period,level" in csv

    capsys.readouterr()
    assert script.main(["--job", str(job), "--results", str(results),
                        "--png", str(tmp_path / "without.png"), "--without-design"]) == 0
    assert "Reusing" in capsys.readouterr().out
    assert (tmp_path / "without.png").is_file()


def test_the_script_refuses_a_job_whose_record_is_missing(tmp_path, capsys):
    job = tmp_path / "job.json"
    job.write_text(json.dumps({"record": {"files": [str(tmp_path / "gone.csv")]}}))
    assert cli().main(["--job", str(job), "--results", str(tmp_path / "r.json")]) == 2
    assert "not found" in capsys.readouterr().out


def test_the_axis_reads_in_ey_then_one_in_x():
    positions_, labels = cli().frequency_ticks(frequency.EY1_AEP, 5e-4)
    assert labels == ["1EY", "1 in 2", "1 in 5", "1 in 10", "1 in 20", "1 in 50",
                      "1 in 100", "1 in 500", "1 in 2,000"]
    assert positions_ == sorted(positions_)


# -- the frequent end of the axis: where it starts, and where the fitted curves do ----

AXIS_START, CURVE_START = "frequent_aep_1_in_x", "curve_frequent_aep_1_in_x"


def test_the_axis_reads_correctly_from_wherever_it_starts():
    ticks = cli().frequency_ticks
    assert ticks(1 / 5, 5e-4)[1] == ["1 in 5", "1 in 10", "1 in 20", "1 in 50", "1 in 100",
                                     "1 in 500", "1 in 2,000"]
    assert ticks(1 / 3, 5e-4)[1][0] == "1 in 5"          # no tick is made up at the limit
    assert ticks(1 / 1.1, 0.1)[1] == ["2EY", "1EY", "1 in 2", "1 in 5", "1 in 10"]


def test_the_starts_are_read_as_blank_unless_they_can_be_drawn():
    limits, ey1 = record.frequent_limits, frequency.EY1_AEP
    assert record.EY1_AEP == ey1
    assert limits({}) == limits({AXIS_START: None, CURVE_START: ""}) == (ey1, ey1)
    assert limits({AXIS_START: 5}) == (0.2, ey1)                    # each is its own
    assert limits({CURVE_START: 10}) == (ey1, 0.1)
    assert limits({AXIS_START: 1.25, CURVE_START: 10}) == (0.8, 0.1)
    assert limits({AXIS_START: 5, CURVE_START: 2}) == (0.2, 0.5)     # the frame clips it
    for bad in (1, 0.5, 0, -4, "x", float("nan"), float("inf")):
        assert limits({AXIS_START: bad, CURVE_START: bad}) == (ey1, ey1)
    # Not more frequent than the end the axis reaches: never a frame drawn backwards.
    assert limits({AXIS_START: 50, CURVE_START: 100}, rare_aep=0.05) == (ey1, ey1)
    assert limits({AXIS_START: 5, CURVE_START: 100}, rare_aep=0.05) == (0.2, ey1)

    assert record.axis_problems({}) == record.axis_problems({AXIS_START: 5}) == []
    assert "above 1" in record.axis_problems({CURVE_START: 0.5})[0]
    assert "1 in 500" in record.axis_problems({AXIS_START: 500, "rare_aep_1_in_x": 500})[0]


def analysed(tmp_path, name, **axes):
    """The script's results for the synthetic job with these axis settings."""
    folder = tmp_path / name
    folder.mkdir()
    job = job_file(folder, **({"axes": axes} if axes else {}))
    results = folder / "results.json"
    assert cli().main(["--job", str(job), "--results", str(results)]) == 0
    return json.loads(results.read_text())


def test_neither_start_changes_the_fit_and_only_a_wider_axis_changes_the_grid(tmp_path):
    blank = analysed(tmp_path, "blank")
    assert blank["job"]["axes"][AXIS_START] is None
    assert blank["grid"]["aep"][0] == pytest.approx(frequency.EY1_AEP)

    def fitted(results):
        return {name: {key: block[key] for key in ("rmse", "params", "draws_used", "z_min",
                                                   "z_max", "error", "warning")}
                for name, block in results["fits"].items()}

    # Starts inside the 1EY frame are a matter of what is drawn: the results are
    # the blank ones to the last digit.
    for name, axes in (("axis", {AXIS_START: 5}), ("curves", {CURVE_START: 10}),
                       ("both", {AXIS_START: 3, CURVE_START: 10})):
        got = analysed(tmp_path, name, **axes)
        for key in ("grid", "fits", "design", "ams"):
            assert got[key] == blank[key], (name, key)

    # An axis, or curves, started further out are evaluated further out, on the
    # same fit.
    for name, axes in (("wide", {AXIS_START: 1.1, CURVE_START: 3}),
                       ("ahead", {CURVE_START: 1.1})):
        wide = analysed(tmp_path, name, **axes)
        assert wide["grid"]["aep"][0] == pytest.approx(1 / 1.1)
        assert wide["grid"]["aep"][-1] == pytest.approx(blank["grid"]["aep"][-1])
        assert fitted(wide) == fitted(blank)
        assert len(wide["design"]["envelope"]) == len(wide["grid"]["z"])


def drawn(monkeypatch, tmp_path, results, with_design=True):
    """The axes of the report figure, kept open to be read."""
    import matplotlib.pyplot as plt
    kept = []
    with monkeypatch.context() as patched:
        patched.setattr(plt, "close", kept.append)
        cli().figure(results, str(tmp_path / "figure.png"), with_design=with_design)
    axis = kept[0].axes[0]
    lines = {line.get_label(): line.get_xdata() for line in axis.get_lines()}
    lines.update({patch.get_label(): patch.get_paths()[0].vertices[:, 0]
                  for patch in axis.collections
                  if patch.get_label().startswith("90% band") and patch.get_paths()})
    plt.close(kept[0])

    def start(prefix):
        return min(float(np.min(x)) for label, x in lines.items()
                   if label.startswith(prefix) and len(x))
    return axis.get_xlim(), [tick.get_text() for tick in axis.get_xticklabels()], start


def test_the_figure_starts_the_axis_and_the_curves_where_the_job_says(tmp_path, monkeypatch):
    script = cli()
    blank = analysed(tmp_path, "blank")
    ey1 = script._z(frequency.EY1_AEP)
    xlim, labels, start = drawn(monkeypatch, tmp_path, blank)
    assert xlim[0] == ey1 and labels[0] == "1EY"
    assert start("Fit to all") == start("Design flood envelope") == pytest.approx(ey1)
    assert start("90% band, all") == pytest.approx(ey1)

    # Results written before the keys existed draw the same frame.
    old = json.loads(json.dumps(blank))
    del old["job"]["axes"][AXIS_START], old["job"]["axes"][CURVE_START]
    assert drawn(monkeypatch, tmp_path, old)[:2] == (xlim, labels)

    def redrawn(**axes):
        results = json.loads(json.dumps(blank))
        results["job"]["axes"].update(axes)
        return drawn(monkeypatch, tmp_path, results)

    # The axis start alone: the frame moves and clips the curves and the design
    # floods, which are not cut.
    moved, labels, start = redrawn(**{AXIS_START: 5})
    assert moved == pytest.approx((script._z(0.2), xlim[1])) and labels[0] == "1 in 5"
    assert start("Fit to all") >= script._z(0.2) - 1e-9
    assert start("Design flood envelope") == pytest.approx(ey1)

    # The curve start alone: the frame stays, the curves and their bands start later.
    kept, labels, start = redrawn(**{CURVE_START: 10})
    assert kept == xlim and labels[0] == "1EY"
    for name in ("Fit to all", "90% band, all"):
        assert script._z(0.1) - 1e-9 <= start(name) < script._z(0.1) + 0.05
    assert start("Design flood envelope") == pytest.approx(ey1)

    # Both, each doing its own; and a curve start out past the axis is clipped by it.
    both, _, start = redrawn(**{AXIS_START: 3, CURVE_START: 10})
    assert both[0] == pytest.approx(script._z(1 / 3))
    assert start("Fit to all") >= script._z(0.1) - 1e-9
    _, _, start = redrawn(**{AXIS_START: 3, CURVE_START: 1.2})
    assert script._z(1 / 3) - 1e-9 <= start("Fit to all") < script._z(1 / 3) + 0.05

    # An axis started past 1EY: the design floods go with it, and the curves stay
    # at 1EY until they are started out there too.
    wide = analysed(tmp_path, "wide", **{AXIS_START: 1.1})
    out, _, start = drawn(monkeypatch, tmp_path, wide)
    assert out[0] == pytest.approx(script._z(1 / 1.1))
    assert start("Design flood envelope") == pytest.approx(out[0])
    for name in ("Fit to all", "90% band, all"):
        assert ey1 - 1e-9 <= start(name) < ey1 + 0.05
    wide["job"]["axes"][CURVE_START] = 1.1
    assert drawn(monkeypatch, tmp_path, wide)[2]("Fit to all") == pytest.approx(out[0])

    # The record alone reaches only as far as its rarest maximum, and a start
    # past that end is read as blank instead of turning the axis round.
    short, labels, _ = redrawn(**{AXIS_START: 5000})
    assert short == xlim
    record_only, _, _ = drawn(monkeypatch, tmp_path, blank, with_design=False)
    results = json.loads(json.dumps(blank))
    results["job"]["axes"][AXIS_START] = 1 / (script._rare_limit(blank, False) * 0.9)
    assert drawn(monkeypatch, tmp_path, results, with_design=False)[0] == record_only


def test_blank_starts_reuse_results_and_a_start_that_cannot_be_one_is_refused(tmp_path,
                                                                             capsys):
    script = cli()
    job = job_file(tmp_path)
    results = tmp_path / "results.json"
    assert script.main(["--job", str(job), "--results", str(results)]) == 0

    # A job written before the keys existed, and one carrying them blank, are one job.
    saved = json.loads(job.read_text())
    assert "axes" not in saved
    saved["axes"] = {AXIS_START: None, CURVE_START: None}
    job.write_text(json.dumps(saved))
    capsys.readouterr()
    assert script.main(["--job", str(job), "--results", str(results)]) == 0
    assert "Reusing" in capsys.readouterr().out

    # A start is part of the job: the curves are not reused from another one's.
    saved["axes"][CURVE_START] = 10
    job.write_text(json.dumps(saved))
    assert script.main(["--job", str(job), "--results", str(results)]) == 0
    assert "Reusing" not in capsys.readouterr().out

    for axes in ({AXIS_START: 0.5}, {CURVE_START: 5000}):
        saved["axes"] = axes
        job.write_text(json.dumps(saved))
        assert script.main(["--job", str(job), "--results", str(results)]) == 2
        assert "ERROR" in capsys.readouterr().out


# -- the figure it was developed on ------------------------------------------------

@pytest.mark.skipif(not (CALLIDE / "out" / "rfsl-215-5-rapid" / "ams.csv").is_file(),
                    reason="needs Sunwater's private Callide study checkout (callide-fsl-reinstate) and its outputs")
def test_the_callide_validation_figure_is_reproduced():
    """curve_review_validation_shouldered_a4.png: RMSE 0.40 and 0.56 m, 338 resamples."""
    ams = record.read_ams_csv(CALLIDE / "out" / "rfsl-215-5-rapid" / "ams.csv")
    table = record.with_positions(ams)
    assert int(table["carried_over"].sum()) == 19

    sources = []
    for path in sorted(glob.glob(str(CALLIDE / "data" / "design_flood_modelling"
                                     / "*_GWL0p3_RFSL__mcdf.csv"))):
        sources.append((float(re.search(r"_mc_([\d.]+)h_", path).group(1)), path))
    out = frequency.analyse(table, fsl=215.5, design_sources=sources,
                            progress=lambda *_: None)

    assert out["fits"]["all"]["rmse"] == pytest.approx(0.4002290, abs=1e-6)
    assert out["fits"]["storm"]["rmse"] == pytest.approx(0.5588990, abs=1e-6)
    assert out["fits"]["all"]["draws_used"] == 338
    assert out["fits"]["storm"]["draws_used"] == 378
