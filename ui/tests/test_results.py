"""The results viewer's reading and critical-duration maths.

The physics being pinned here, because it decides what counts as a finding:
for lake level the critical duration is long at frequent AEPs (it takes volume
to fill and charge the storage) and short on the rare tail (the dam behaving as
a conveyance). Peak inflow does not do that. So the pinned-end warnings must
fire in whichever direction the data actually goes - never a hard-coded one.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from core import results
from core.outputs import quantile_files

STANDARD_AEPS = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                 10000, 20000, 50000, 100000]


def write_quantiles(path: Path, column: str, values: dict) -> Path:
    """A quantile file in the shape lib/MCScheme.py:353 writes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame({
        "aep (1 in x)": list(values),
        "probability": [1.0 / aep for aep in values],
        column: list(values.values()),
    })
    frame.to_csv(path, index=False)
    return path


def level_curve(duration: float) -> dict:
    """A level curve whose critical duration goes long-then-short.

    The long durations win while the dam is filling, the short ones once it is
    spilling and behaving as a conveyance.
    """
    curve = {}
    for aep in STANDARD_AEPS:
        z = results.normal_variate(aep)
        # filling: rewards volume, and volume comes with duration
        filling = 8.0 * math.log10(duration) * math.exp(-z)
        # conveyance: rewards peak intensity, which the short storms have
        conveyance = 2.5 * z ** 2 / math.log10(duration)
        curve[aep] = round(100.0 + filling + conveyance, 3)
    return curve


def mc_row(folder: Path, name: str, duration) -> pd.Series:
    return pd.Series({
        "Method": "monte carlo",
        "Output file": str(folder / name),
        "Duration": duration,
        "Output suffix": "",
        "Results folder": "",
    })


# -- finding the files -------------------------------------------------------

def test_monte_carlo_quantile_files_are_found(tmp_path):
    base = tmp_path / "results" / "TFD_mc_24h"
    write_quantiles(Path(f"{base}_inflow.csv"), "inflow", {2: 100.0})
    write_quantiles(Path(f"{base}_level.csv"), "level", {2: 670.0})
    write_quantiles(Path(f"{base}_inflowVol24h.csv"), "Vol24h", {2: 50000.0})
    (base.parent / "TFD_mc_24h_volume.csv").write_text("not a quantile file\n")

    found = quantile_files(mc_row(tmp_path / "results", "TFD_mc_24h", 24), tmp_path)

    assert set(found) == {"inflow", "level", "inflowVol24h"}
    assert found["level"].column == "level"
    # the file is tagged inflowVol24h, the column inside it is Vol24h
    assert found["inflowVol24h"].column == "Vol24h"
    assert found["inflowVol24h"].is_volume


def test_outflow_that_was_never_analysed_is_simply_absent(tmp_path):
    base = tmp_path / "TFD_mc_24h"
    write_quantiles(Path(f"{base}_inflow.csv"), "inflow", {2: 100.0})
    found = quantile_files(mc_row(tmp_path, "TFD_mc_24h", 24), tmp_path)
    assert set(found) == {"inflow"}


def test_reservoir_routing_naming_is_understood(tmp_path):
    results_folder = tmp_path / "rr"
    row = pd.Series({
        "Method": "reservoir routing",
        "Output file": "TFD_rr_24h",
        "Output suffix": "optB",
        "Results folder": str(results_folder),
        "Analyse volumes": "yes",
    })
    write_quantiles(results_folder / "TFD_rr_24h__level_quantiles_optB.csv",
                    "level", {2: 670.0})
    write_quantiles(results_folder / "TFD_rr_24h__inflowVol72h_quantiles_optB.csv",
                    "Vol72h", {2: 90000.0})

    found = quantile_files(row, tmp_path)

    assert set(found) == {"level", "inflowVol72h"}
    assert found["inflowVol72h"].column == "Vol72h"


def test_the_ensemble_method_contributes_nothing(tmp_path):
    """lib/EnbAnalysis.py writes its own critical duration analysis."""
    row = pd.Series({"Method": "ensemble", "Output file": str(tmp_path / "x")})
    assert quantile_files(row, tmp_path) == {}


# -- labelling ---------------------------------------------------------------

def _frame(tmp_path, durations, column="level", key="level"):
    rows = []
    for duration in durations:
        name = f"TFD_mc_{duration:g}h"
        write_quantiles(Path(f"{tmp_path / name}_{key}.csv"), column,
                        level_curve(duration))
        rows.append(mc_row(tmp_path, name, duration))
    return pd.DataFrame(rows)


def test_curves_are_labelled_and_ordered_by_duration(tmp_path):
    frame = _frame(tmp_path, [72, 24, 120])
    sources = results.sources_for_rows(frame, tmp_path)["level"]
    assert [source.label for source in sources] == ["24h", "72h", "120h"]


def test_two_rows_of_one_duration_keep_separate_labels(tmp_path):
    """Otherwise the second curve silently replaces the first in the frame."""
    for name in ("TFD_mc_24h_a", "TFD_mc_24h_b"):
        write_quantiles(Path(f"{tmp_path / name}_level.csv"), "level",
                        level_curve(24))
    frame = pd.DataFrame([mc_row(tmp_path, "TFD_mc_24h_a", 24),
                          mc_row(tmp_path, "TFD_mc_24h_b", 24)])
    sources = results.sources_for_rows(frame, tmp_path)["level"]
    labels = [source.label for source in sources]
    assert len(set(labels)) == 2
    assert all(label.startswith("24h") for label in labels)

    comparison = results.compare(sources)
    assert len(comparison.frame.columns) == 2


def test_a_windows_path_in_output_file_does_not_become_the_label(tmp_path):
    """Sims lists hold Windows paths; 'a\\b\\c'.name is the lot on POSIX.

    It reaches the exported filename, where it is refused as a name, so the
    basename has to be taken after normalising the separators.
    """
    folder = tmp_path / "sims_mc" / "results"
    write_quantiles(folder / "TFD_mc_18h_level.csv", "level", level_curve(18))
    row = pd.Series({
        "Method": "monte carlo", "Duration": None, "Output suffix": "",
        "Results folder": "",
        "Output file": "sims_mc\\results\\TFD_mc_18h",
    })
    source = results.sources_for_rows(pd.DataFrame([row]), tmp_path)["level"][0]
    assert source.output_name == "TFD_mc_18h"
    assert "\\" not in source.label


def test_duration_comes_off_the_name_when_the_column_is_blank(tmp_path):
    write_quantiles(Path(f"{tmp_path / 'TFD_mc_18h'}_level.csv"), "level",
                    level_curve(18))
    row = mc_row(tmp_path, "TFD_mc_18h", None)
    frame = pd.DataFrame([row])
    sources = results.sources_for_rows(frame, tmp_path)["level"]
    assert sources[0].duration == 18


# -- reading -----------------------------------------------------------------

def test_mismatched_aep_sets_outer_join_without_interpolation(tmp_path):
    """get_standard_aeps extends to the AEP of the PMP, so the sets can differ."""
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    {2: 100.0, 10: 110.0, 100: 120.0})
    write_quantiles(Path(f"{tmp_path / 'b_48h'}_level.csv"), "level",
                    {2: 101.0, 10: 109.0, 100: 119.0, 1000: 130.0})
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])

    assert list(comparison.frame.index) == [2, 10, 100, 1000]
    assert pd.isna(comparison.frame.loc[1000, "24h"])
    analysis = results.analyse(comparison)
    assert analysis.critical.loc[1000] == "48h"      # the only one that reaches


def test_a_damaged_file_is_reported_not_raised(tmp_path):
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    level_curve(24))
    # the right name, the wrong contents
    (tmp_path / "b_48h_level.csv").write_text("something,else\n1,2\n")
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])

    assert comparison.labels == ["24h"]
    assert len(comparison.problems) == 1
    assert "aep (1 in x)" in comparison.problems[0][1]


# -- the analysis ------------------------------------------------------------

def test_envelope_is_the_per_aep_maximum(tmp_path):
    frame = _frame(tmp_path, [24, 48, 120])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    analysis = results.analyse(comparison)
    for aep in comparison.frame.index:
        assert analysis.envelope[aep] == pytest.approx(comparison.frame.loc[aep].max())
        assert (comparison.frame.loc[aep, analysis.critical[aep]]
                == pytest.approx(analysis.envelope[aep]))


def test_level_goes_long_at_the_frequent_end_and_short_on_the_rare_tail(tmp_path):
    frame = _frame(tmp_path, [12, 24, 48, 120])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    analysis = results.analyse(comparison)

    assert analysis.critical.loc[2] == "120h"
    assert analysis.critical.loc[100000] == "12h"
    durations = [comparison.durations[owner] for owner in analysis.critical]
    assert durations == sorted(durations, reverse=True), (
        "the critical duration should shorten monotonically as the AEP rarens")
    assert [band.label for band in analysis.bands][0] == "120h"
    assert any(switch.before == "120h" for switch in analysis.switches)


def test_both_ends_warn_when_the_critical_duration_is_pinned(tmp_path):
    frame = _frame(tmp_path, [12, 24, 48, 120])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    warnings = results.analyse(comparison).warnings

    assert any("longest duration you ran (120h)" in text for text in warnings)
    assert any("shortest duration you ran (12h)" in text for text in warnings)


def test_no_pinned_warning_when_the_range_brackets_the_critical_duration(tmp_path):
    """The interior durations own every AEP, so nothing needs extending."""
    write_quantiles(Path(f"{tmp_path / 'a_12h'}_level.csv"), "level",
                    {2: 100.0, 100: 100.0, 10000: 100.0})
    write_quantiles(Path(f"{tmp_path / 'b_24h'}_level.csv"), "level",
                    {2: 105.0, 100: 130.0, 10000: 160.0})
    write_quantiles(Path(f"{tmp_path / 'c_48h'}_level.csv"), "level",
                    {2: 101.0, 100: 101.0, 10000: 101.0})
    frame = pd.DataFrame([mc_row(tmp_path, "a_12h", 12),
                          mc_row(tmp_path, "b_24h", 24),
                          mc_row(tmp_path, "c_48h", 48)])
    analysis = results.analyse(
        results.compare(results.sources_for_rows(frame, tmp_path)["level"]))

    assert set(analysis.critical) == {"24h"}
    assert not any("critical at" in text for text in analysis.warnings)
    assert analysis.never_critical == ["12h", "48h"]


def test_a_switch_inside_the_noise_floor_is_called_out(tmp_path):
    """Nearly coincident curves make idxmax hop; that is not a crossover."""
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    {2: 216.00, 100: 220.00})
    write_quantiles(Path(f"{tmp_path / 'b_48h'}_level.csv"), "level",
                    {2: 216.01, 100: 219.98})
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    analysis = results.analyse(
        results.compare(results.sources_for_rows(frame, tmp_path)["level"]))

    assert analysis.margin.loc[2] == pytest.approx(0.01)     # metres, not percent
    assert analysis.switches[0].strength < analysis.noise_floor
    assert any("inside sampling noise" in text for text in analysis.warnings)
    # the 24h band's peak margin, in metres - 220.00 against 219.98
    assert any("0.020 m clear" in text for text in analysis.warnings)


def test_level_margins_are_metres_and_flows_are_percentages(tmp_path):
    """Level is an interval scale on an arbitrary datum.

    A percentage of 217 m AHD says nothing: at Callide the durations separate
    by 0.01-0.10 m, which is 0.005-0.05% - so on a percentage floor every real
    level crossover is dismissed as noise, on the very result type where the
    critical duration is most interesting.
    """
    assert results.margin_scale("level") == (results.ABSOLUTE,
                                             results.MARGIN_NOISE_METRES,
                                             "margin (m)")
    for key in ("inflow", "outflow", "inflowVol24h"):
        kind, floor, label = results.margin_scale(key)
        assert (kind, floor, label) == (results.PERCENT,
                                        results.MARGIN_NOISE_PERCENT, "margin %")

    # the same numbers, judged both ways
    for kind, values in (("level", {2: 216.0, 100: 220.0}),
                         ("inflow", {2: 216.0, 100: 220.0})):
        write_quantiles(Path(f"{tmp_path / kind}_a_24h_{kind}.csv"), kind, values)
        write_quantiles(Path(f"{tmp_path / kind}_b_48h_{kind}.csv"), kind,
                        {2: 216.1, 100: 219.0})
        frame = pd.DataFrame([mc_row(tmp_path, f"{kind}_a_24h", 24),
                              mc_row(tmp_path, f"{kind}_b_48h", 48)])
        analysis = results.analyse(
            results.compare(results.sources_for_rows(frame, tmp_path)[kind]))
        if kind == "level":
            # 0.1 m clear - worth acting on, and above the metre floor
            assert analysis.margin.loc[2] == pytest.approx(0.1)
            assert not any("sampling noise" in t for t in analysis.warnings)
        else:
            # the same 0.1 in 216 is 0.05% - inside the noise floor for a flow
            assert analysis.margin.loc[2] == pytest.approx(0.0463, abs=1e-3)
            assert any("sampling noise" in t for t in analysis.warnings)


def test_the_noise_floor_can_be_overridden(tmp_path):
    """Different dams and different result types justify different floors."""
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    {2: 216.00, 100: 220.00})
    write_quantiles(Path(f"{tmp_path / 'b_48h'}_level.csv"), "level",
                    {2: 216.02, 100: 219.98})
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])

    assert any("sampling noise" in t
               for t in results.analyse(comparison).warnings)
    loose = results.analyse(comparison, noise_floor=0.005)
    assert loose.noise_floor == 0.005
    assert not any("sampling noise" in t for t in loose.warnings)


def test_a_crossover_is_judged_over_its_range_not_at_the_crossing_point(tmp_path):
    """The margin AT a crossover is ~0 by definition - the curves are equal.

    So strength has to come from how convincingly the new owner wins over the
    range it then holds, or every real crossover would be dismissed as noise.
    """
    frame = _frame(tmp_path, [12, 24, 48, 120])
    comparison = results.compare(results.sources_for_rows(frame, tmp_path)["level"])
    analysis = results.analyse(comparison)

    switch = analysis.switches[0]
    assert analysis.margin.loc[switch.aep] < switch.strength / 10, (
        "the curves are nearly equal at the crossing point, whatever the "
        "crossover turns out to mean")
    assert switch.strength > 10 * analysis.noise_floor
    assert not any("sampling noise" in text for text in analysis.warnings)


def test_a_real_crossover_is_not_called_noise(tmp_path):
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    {2: 100.0, 100: 160.0})
    write_quantiles(Path(f"{tmp_path / 'b_48h'}_level.csv"), "level",
                    {2: 130.0, 100: 120.0})
    frame = pd.DataFrame([mc_row(tmp_path, "a_24h", 24),
                          mc_row(tmp_path, "b_48h", 48)])
    analysis = results.analyse(
        results.compare(results.sources_for_rows(frame, tmp_path)["level"]))

    assert not any("sampling noise" in text for text in analysis.warnings)
    assert [(s.before, s.after) for s in analysis.switches] == [("48h", "24h")]


def test_bands_are_contiguous_and_cover_the_axis(tmp_path):
    frame = _frame(tmp_path, [12, 24, 48, 120])
    analysis = results.analyse(
        results.compare(results.sources_for_rows(frame, tmp_path)["level"]))
    bands = analysis.bands
    assert len(bands) > 1
    for earlier, later in zip(bands, bands[1:]):
        assert earlier.z_to == pytest.approx(later.z_from)


# -- the axis ----------------------------------------------------------------

@pytest.mark.parametrize("aep, expected", [
    (2, 0.0), (10, 1.2815515655), (100, 2.3263478740), (100000, 4.2648907939)])
def test_normal_variate_matches_ndtri(aep, expected):
    """The same values scipy.special.ndtri gives UtilModule.plot_durations."""
    assert results.normal_variate(aep) == pytest.approx(expected, abs=1e-9)


def test_an_impossible_aep_does_not_raise():
    assert math.isnan(results.normal_variate(1))
    assert math.isnan(results.normal_variate("not a number"))


def test_level_is_linear_and_flows_are_logarithmic():
    assert results.y_axis("level") == ("Peak lake level (m AHD)", False)
    assert results.y_axis("inflow")[1] is True
    label, log = results.y_axis("inflowVol24h")
    assert "24 h" in label and log is True
