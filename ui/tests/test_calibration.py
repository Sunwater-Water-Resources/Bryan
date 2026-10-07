"""The calibration table: modelled dam inflow against the reverse-routed inflow."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from core import calibration, lakerecord, reporttables, wordtable
from core import reporttables as rt
from core import study as studies
from report_fixtures import build_study

START = pd.Timestamp("2013-01-26 00:00")


def observed_frame(rates, step_h=1.0, start=START, uncertain=None):
    """An event hydrograph as the inflow step writes it: one row per interval,
    stamped at its midpoint, with the interval's length in dt_s."""
    rates = np.asarray(rates, dtype=float)
    dt = step_h * 3600.0
    mids = start + pd.to_timedelta((np.arange(len(rates)) + 0.5) * step_h, unit="h")
    frame = pd.DataFrame({"Inflow_m3s": rates, "Inflow_uncorrected_m3s": rates,
                          "Release_uncertain": uncertain if uncertain is not None
                          else [False] * len(rates),
                          "dt_s": dt}, index=pd.DatetimeIndex(mids, name="Timestamp"))
    return frame


def triangle(n=25, peak=1000.0):
    half = n // 2
    return np.concatenate([np.linspace(0, peak, half + 1), np.linspace(peak, 0, n - half)[1:]])


def as_modelled(frame, factor=1.0, shift_h=0.0):
    return pd.Series(frame["Inflow_m3s"].to_numpy() * factor,
                     index=frame.index + pd.Timedelta(hours=shift_h))


def test_a_model_that_matches_the_record_scores_perfectly():
    observed = observed_frame(triangle())
    stats = calibration.statistics(observed, as_modelled(observed), smoothing="")
    assert stats["pr"] == pytest.approx(0)
    assert stats["nse"] == pytest.approx(1)
    assert stats["timing"] == pytest.approx(0)
    # the modelled samples are the interval means at the midpoints, so their
    # trapezoid misses only the two half-intervals at the ends, which are zero here
    assert stats["vr"] == pytest.approx(0, abs=0.5)


def test_a_model_ten_percent_low_has_ratios_of_minus_ten_and_the_nse_that_follows():
    observed = observed_frame(triangle())
    stats = calibration.statistics(observed, as_modelled(observed, 0.9), smoothing="")
    assert stats["pr"] == pytest.approx(-10)
    o = observed["Inflow_m3s"].to_numpy()
    expected = 1 - np.sum((0.1 * o) ** 2) / np.sum((o - o.mean()) ** 2)
    assert stats["nse"] == pytest.approx(expected)


def test_a_late_peak_is_timed_in_hours_and_the_window_is_the_overlap():
    observed = observed_frame(triangle())
    stats = calibration.statistics(observed, as_modelled(observed, shift_h=2), smoothing="")
    assert stats["timing"] == pytest.approx(2)
    start, end = stats["window"]
    assert start == observed.index[0] + pd.Timedelta(hours=2)    # the model starts late
    assert end == observed.index[-1] + pd.Timedelta(hours=0.5)   # the record's last interval ends


def test_the_observed_volume_counts_only_the_part_of_each_interval_inside_the_window():
    observed = observed_frame([100.0, 100.0, 100.0], step_h=1.0)
    t0 = calibration.seconds([START + pd.Timedelta(minutes=30)])[0]
    t1 = calibration.seconds([START + pd.Timedelta(hours=2)])[0]
    # 1.5 h at 100 m3/s = 540,000 m3 = 540 ML
    assert calibration.observed_volume(observed, "Inflow_m3s", t0, t1) == pytest.approx(540)


def test_smoothing_takes_the_peak_off_the_time_average_as_the_inflow_record_does():
    rates = np.zeros(48)
    rates[24] = 1200.0                                   # one spike of 15 minutes
    observed = observed_frame(rates, step_h=0.25)
    modelled = as_modelled(observed)
    raw = calibration.statistics(observed, modelled, smoothing="")
    smooth = calibration.statistics(observed, modelled, smoothing="1h")
    assert raw["observed_peak"] == pytest.approx(1200)
    assert smooth["observed_peak"] == pytest.approx(300)    # a quarter of the hour


def test_the_uncertain_release_share_is_measured_over_the_window():
    flags = [False] * 20 + [True] * 5
    observed = observed_frame(triangle(), uncertain=flags)
    stats = calibration.statistics(observed, as_modelled(observed), smoothing="")
    assert stats["uncertain"] == pytest.approx(5 / 25, abs=0.03)


@pytest.mark.parametrize("statistic, value, expected", [
    # Callide's Table 16 against its Table 15, as the report grades them
    ("pr", -14, "Good"), ("pr", -21, "Fair"), ("pr", -12, "Good"),
    ("vr", -20, "Good"), ("vr", 5, "Excellent"), ("vr", -11, "Excellent"),
    ("nse", 0.70, "Poor"), ("nse", 0.87, "Fair"), ("nse", 0.93, "Good"),
    ("pr", 10, "Excellent"), ("pr", -55, "No Data/Exclude"), ("nse", 0.4, "No Data/Exclude"),
    ("timing", -0.5, "Excellent"), ("timing", 2.5, "Poor"),
])
def test_the_classes_are_the_report_s(statistic, value, expected):
    assert calibration.classify(statistic, value)["name"] == expected


def test_no_value_has_no_class():
    assert calibration.classify("nse", float("nan")) is None


def test_a_modelled_csv_may_give_dates_or_hours_from_the_event_start(tmp_path):
    dated = tmp_path / "dated.csv"
    dated.write_text("Time,Q\n26/01/2013 00:00,1\n26/01/2013 01:00,2\n", encoding="utf-8")
    series = calibration.read_modelled(dated)
    assert series.index[1] == pd.Timestamp("2013-01-26 01:00")

    hours = tmp_path / "hours.csv"
    hours.write_text("hour,Inflow,Outflow\n0,1,0\n1.5,2,0\n", encoding="utf-8")
    series = calibration.read_modelled(hours, "Inflow", "2013-01-26 00:00")
    assert series.index[1] == pd.Timestamp("2013-01-26 01:30")
    assert list(series) == [1, 2]

    with pytest.raises(calibration.CalibrationError, match="hours"):
        calibration.read_modelled(hours)
    with pytest.raises(calibration.CalibrationError, match="no column 'Q'"):
        calibration.read_modelled(hours, "Q", "2013-01-26")


# -- the table, built from a study --------------------------------------------------

@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def with_event(study, rates, *, uncertain=None):
    out = lakerecord.path_of(study, lakerecord.settings(study)["inflow"]["out"])
    (out / "hydrographs").mkdir(parents=True, exist_ok=True)
    frame = observed_frame(rates, uncertain=uncertain)
    path = out / "hydrographs" / "2013_peak.csv"
    frame.to_csv(path)
    (out / "summary.json").write_text(json.dumps({
        "settings": {"smoothing": ""},
        "hydrographs": [{"name": "2013_peak", "file": str(path),
                         "start": f"{START:%Y-%m-%d %H:%M}"}]}), encoding="utf-8")
    return frame


def modelled_csv(study, frame, factor):
    path = study.folder / "urbs" / "jan13.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"Time": frame.index.strftime("%Y-%m-%d %H:%M"),
                  "Dam inflow": frame["Inflow_m3s"].to_numpy() * factor}) \
        .to_csv(path, index=False)
    return "urbs/jan13.csv"


def spec(**events):
    out = reporttables.new_spec(reporttables.CALIBRATION)
    out["smoothed"] = False
    out["events"] = [dict(label="Jan-13", observed="2013_peak", column="", start="", **events)]
    return out


def test_the_table_reads_the_event_and_the_csv_and_shades_each_statistic(study):
    frame = with_event(study, triangle())
    table = reporttables.build(study, spec(modelled=modelled_csv(study, frame, 0.86)))
    assert table.problems == []
    assert table.header == ["Event", "Mod", "Rated", "PR", "Mod", "Rated", "VR",
                            "Nash-Sutcliffe"]
    assert table.header_groups == [("", 1), ("Flow (m3/s)", 3), ("Volume (ML)", 3), ("", 1)]
    row = table.rows[0]
    assert row.cells[:4] == ["Jan-13", "860", "1,000", "-14%"]
    assert row.fills[3] == "#92D050"                    # -14% is Good
    assert row.fills[0] == "" and row.fills[1] == ""


def test_what_is_missing_is_said_per_event_and_the_row_is_kept(study):
    frame = with_event(study, triangle())
    table = reporttables.build(study, spec(modelled="urbs/none.csv"))
    assert table.rows[0].cells[0] == "Jan-13" and table.rows[0].cells[1] == "–"
    assert any("not found" in problem for problem in table.problems)

    # no intervals file to take the observed side from, and no event chosen
    table = reporttables.build(study, {**spec(modelled=""), "events": [
        {"label": "Feb-15", "observed": "1999_peak",
         "modelled": modelled_csv(study, frame, 1.0)}]})
    assert any("choose the observed hydrograph" in p for p in table.problems)


def test_an_uncertain_release_in_the_window_is_reported(study):
    frame = with_event(study, triangle(), uncertain=[False] * 20 + [True] * 5)
    table = reporttables.build(study, spec(modelled=modelled_csv(study, frame, 1.0)))
    assert any("uncertain" in problem for problem in table.problems)


def test_the_timing_column_is_optional(study):
    frame = with_event(study, triangle())
    table = reporttables.build(study, {**spec(modelled=modelled_csv(study, frame, 1.0)),
                                       "timing": True})
    assert table.header[-1] == "Timing (h)"
    assert table.rows[0].cells[-1] == "+0.0"


def test_the_criteria_table_is_the_report_s(study):
    table = reporttables.build(study, reporttables.new_spec(reporttables.CRITERIA))
    assert [row.cells[0] for row in table.rows] == ["Excellent", "Good", "Fair", "Poor",
                                                     "No Data/Exclude"]
    assert table.rows[0].cells[2:] == ["+/- 10%", "+/- 15%", ">=0.95", "<=0.50"]
    assert table.rows[-1].cells[2:] == ["> 50%", "> 50%", "<0.5", ">3"]


# -- in Word --------------------------------------------------------------------------

def test_a_grouped_header_pastes_as_two_rows_and_a_fill_as_a_cell_background():
    table = wordtable.ReportTable(header=["Event", "Mod", "Rated", "NSE"])
    table.header_groups = [("", 1), ("Flow (m3/s)", 2), ("", 1)]
    row = table.add(["Jan-13", "903", "1,056", "70%"])
    row.fills = ["", "", "", "#FFC000"]
    html = wordtable.to_html(table)
    assert html.count("<tr>") == 3
    assert 'rowspan="2"' in html and 'colspan="2"' in html
    assert html.count(">Event<") == 1                   # not repeated in the second row
    assert "background:#FFC000" in html
    text = wordtable.to_text(table)
    assert text.splitlines()[:2] == ["\tFlow (m3/s)\t\t", "Event\tMod\tRated\tNSE"]


def test_a_table_without_groups_or_fills_pastes_as_before():
    table = wordtable.ReportTable(header=["A", "B"])
    table.add(["1", "2"])
    html = wordtable.to_html(table)
    assert "rowspan" not in html and html.count("<tr>") == 2


def with_record(study, rates, start=START):
    """The whole inflow record, as the inflow step writes inflow_intervals.csv.gz."""
    out = lakerecord.path_of(study, lakerecord.settings(study)["inflow"]["out"])
    frame = observed_frame(rates, start=start)
    frame.insert(0, "Interval_start", frame.index - pd.Timedelta(minutes=30))
    path = out / "inflow_intervals.csv.gz"
    frame.to_csv(path)
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    summary["intervals"] = str(path)
    (out / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return frame


def test_the_window_is_the_modelled_period_even_past_the_event_hydrograph(study):
    # the record runs 2 days; the event hydrograph only its first 10 hours
    record = np.concatenate([triangle(), np.zeros(23)])
    with_event(study, record[:10])
    full = with_record(study, record)
    csv = modelled_csv(study, full.iloc[:40], 1.0)

    table = reporttables.build(study, spec(modelled=csv))
    assert table.problems == []
    row = table.rows[0].cells
    assert row[1] == row[2] == "1,000"                  # the peak at hour 12 is in it
    assert row[3] == "0%"

    from core import calibration as c
    observed = c.observed_for(c.lakerecord.intervals_path(c.lakerecord.last_summary(
        study, lakerecord.settings(study), "inflow")), None,
        c.read_modelled(study.folder / csv), pd.Timedelta(0))
    assert len(observed) == 40                          # clipped to the model run


def test_a_model_run_past_the_record_is_reported(study):
    with_event(study, triangle())
    full = with_record(study, triangle())
    longer = observed_frame(np.concatenate([triangle(), np.zeros(6)]))
    table = reporttables.build(study, spec(modelled=modelled_csv(study, longer, 1.0)))
    assert any("past the observed record" in problem for problem in table.problems)
    assert table.rows[0].cells[3] == "0%"
