"""lib/LakeLevelRecord.py: reading level exports and deriving the annual maxima.

Here rather than in Bryan's tests/ because the module is on the UI's allow-list
and must stay importable without scipy - which this environment proves.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from core.bryan import lake_level_record
from lake_fixtures import FSL, synthetic_level, write_hydstra, write_wmip

R = lake_level_record()


def series(values, start="2000-09-28 00:00", freq="6h"):
    return pd.Series(values, index=pd.date_range(start, periods=len(values), freq=freq),
                     name="level")


# -- reading -------------------------------------------------------------------------

def test_a_hydstra_export_reads_with_its_site_and_says_it_was_not_screened(tmp_path):
    record = series([200.001, 200.002, 200.5, 201.0, 200.9])
    path = write_hydstra(tmp_path / "x.csv", record, site="130314C",
                         description="Callide Ck at Callide Dam HW (Dam)")
    got = R.read_export(path)

    assert got.sites == [("130314C", "130314C - Callide Ck at Callide Dam HW (Dam)")]
    assert got.value_kinds == ["Point"]
    assert not got.screened
    assert any("no quality codes" in note for note in got.notes)
    pd.testing.assert_series_equal(got.level, record, check_names=False,
                                   check_index_type=False, check_freq=False)


def test_hydstra_writes_the_day_and_hour_unpadded_and_it_still_parses(tmp_path):
    record = series([1.0, 2.0], start="1989-10-01 08:00", freq="1h")
    path = write_hydstra(tmp_path / "x.csv", record)
    assert "1/10/1989 8:00" in path.read_text()
    assert list(R.read_export(path).level.index) == list(record.index)


def test_duplicate_timestamps_keep_the_last_and_say_so(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("Time,A,\nand,130,\nDate,Storage Level (m),\n,Point,\n"
                    "1/10/2000 0:00,1.000,Sites:\n1/10/2000 0:00,2.000,\n"
                    "1/10/2000 1:00,3.000,\n")
    got = R.read_export(path)
    assert list(got.level) == [2.0, 3.0]
    assert any("1 duplicate" in note for note in got.notes)


def test_a_wmip_export_drops_the_unusable_quality_codes(tmp_path):
    record = series([200.0, 250.0, 201.0])
    got = R.read_export(write_wmip(tmp_path / "x.csv", record, quality=[9, 255, 59]))
    assert got.screened
    assert list(got.level) == [200.0, 201.0]


def test_averaged_values_are_flagged_because_they_understate_peaks(tmp_path):
    got = R.read_export(write_hydstra(tmp_path / "x.csv", series([1.0, 2.0]), kind="Mean"))
    assert any("'Mean'" in note for note in got.notes)


def test_a_file_with_a_changed_layout_is_refused_rather_than_half_read(tmp_path):
    path = tmp_path / "x.csv"
    rows = "\n".join(f"2000-10-{d:02d} 00:00,1.0," for d in range(1, 20))
    path.write_text("Time,A,\nand,130,\nDate,Storage Level (m),\n,Point,\n"
                    "1/10/2000 0:00,1.0,\n" + rows + "\n")
    with pytest.raises(ValueError, match="layout may have changed"):
        R.read_export(path)


def test_gauges_hand_over_at_the_later_gauges_first_reading(tmp_path):
    early = series([1.0, 1.1, 1.2, 1.3], start="2000-01-01", freq="1D")
    late = series([5.0, 5.1, 5.2], start="2000-01-03", freq="1D")
    a = write_hydstra(tmp_path / "a.csv", early, site="A")
    b = write_hydstra(tmp_path / "b.csv", late, site="B")

    got = R.read_record([a, b])
    assert list(got.level) == [1.0, 1.1, 5.0, 5.1, 5.2]
    assert [site for site, _ in got.sites] == ["A", "B"]
    assert any("jump" in note for note in got.notes)
    with pytest.raises(ValueError, match="earliest first"):
        R.read_record([b, a])


# -- water years and maxima ------------------------------------------------------------

def test_water_years_are_labelled_by_the_year_they_end_in():
    index = pd.DatetimeIndex(["2010-09-30 23:00", "2010-10-01 00:00", "2011-01-01"])
    assert list(R.water_year(index, 10)) == [2010, 2011, 2011]
    assert R.period_label(2011, 10) == "2010-11"
    assert list(R.water_year(index, 1)) == [2010, 2010, 2011]
    assert R.period_label(2011, 1) == "2011"
    assert R.water_year_start(2011, 10) == pd.Timestamp("2010-10-01")


def test_a_maximum_the_year_opened_at_is_carried_over():
    # Falling from the start of water year 2001, then a smaller storm in February.
    index = pd.date_range("2000-10-01", "2001-09-30 18:00", freq="6h")
    level = pd.Series(210.0 - np.arange(len(index)) * 0.001, index=index)
    level.loc["2001-02-01":"2001-02-02"] = 209.9
    ams = R.annual_maxima(level, 10)

    row = ams.iloc[0]
    assert row["carried_over"] and row["days_into_year"] == 0.0
    assert row["complete"]


def test_a_lake_back_on_the_plateau_later_is_not_carried_over():
    """Held at full supply on 1 October and put back there by a flood in March:
    read to the millimetre the two are equal, and the flood is the maximum."""
    index = pd.date_range("2000-10-01", "2001-09-30 18:00", freq="6h")
    level = pd.Series(FSL - 0.5, index=index)
    level.iloc[:4] = FSL
    level.loc["2001-03-10"] = FSL
    row = R.annual_maxima(level, 10).iloc[0]

    assert not row["carried_over"]
    assert row["level_at"] == pd.Timestamp("2001-03-10")
    # ...but not when it only comes back to within a few millimetres,
    level.loc["2001-03-10"] = FSL - 0.004
    assert R.annual_maxima(level, 10).iloc[0]["carried_over"]
    # ...nor when it was simply held on the plateau past the window and let go.
    level.loc["2001-03-10"] = FSL - 0.5
    level.iloc[:60] = FSL
    row = R.annual_maxima(level, 10).iloc[0]
    assert row["carried_over"] and row["days_into_year"] == 0.0


def test_a_partly_covered_year_is_kept_and_flagged_with_its_reading_interval():
    index = pd.date_range("2000-10-01", "2001-01-01", freq="1D")
    ams = R.annual_maxima(pd.Series(np.linspace(200, 201, len(index)), index=index), 10)
    assert not ams.iloc[0]["complete"]
    assert ams.iloc[0]["coverage"] == pytest.approx(92 / 365, abs=0.01)
    assert ams.iloc[0]["reading_interval_min"] == 1440.0


def test_cunnane_positions_and_censored_positions_share_a_denominator():
    p, z = R.plotting_positions([3.0, 1.0, 2.0])
    assert p == pytest.approx([(1 - 0.4) / 3.2, (3 - 0.4) / 3.2, (2 - 0.4) / 3.2])
    assert z[0] == pytest.approx(R.normal_variate_of_probability(0.6 / 3.2))

    zc, values = R.censored_positions([3.0, 2.0], n_years=10)
    assert list(values) == [2.0, 3.0]
    assert zc[1] == pytest.approx(R.normal_variate_of_probability(0.6 / 10.2))


def test_the_curve_table_carries_storm_positions_only_for_storm_driven_years():
    ams = R.annual_maxima(synthetic_level(years=20), 10)
    table = R.with_positions(ams)
    assert list(table["z"]) == sorted(table["z"])
    assert table.loc[table["carried_over"], "storm_z"].isna().all()
    assert table.loc[~table["carried_over"], "storm_z"].notna().all()
    fewer = R.with_positions(ams.assign(complete=[True] * 19 + [False]),
                             include_incomplete=False)
    assert len(fewer) == 19


# -- the CSV --------------------------------------------------------------------------

def test_a_homogenised_series_is_read_from_its_own_column_names(tmp_path):
    path = tmp_path / "ams.csv"
    pd.DataFrame({
        "WaterYear": [1971, 1972], "Level_max": [209.7, 208.4],
        "New_Level_max": [215.5179, 215.3728], "Coverage": [0.97, 0.99],
        "New_Level_max_at": ["1971-02-09 19:00:00.000", "1971-10-01 00:00:00.000"],
    }).to_csv(path, index=False)
    got = R.read_ams_csv(path)

    assert list(got["level"]) == [215.5179, 215.3728]
    assert list(got["carried_over"]) == [False, True]
    assert list(got["period"]) == ["1970-71", "1971-72"]
    assert list(R.read_ams_csv(path, level_column="Level_max")["level"]) == [209.7, 208.4]


def test_the_written_series_reads_back_the_same(tmp_path):
    job = {"record": {"files": ["x.csv"]}}
    ams = R.annual_maxima(synthetic_level(years=12), 10)
    text = R.ams_csv_text(R.with_positions(ams), R.ams_metadata(job))
    assert text.startswith("# source: x.csv\n# water_year_start: October\n")
    path = tmp_path / "ams.csv"
    path.write_text(text)

    back = R.read_ams_csv(path).sort_values("water_year").reset_index(drop=True)
    ams = ams.sort_values("water_year").reset_index(drop=True)
    assert list(back["level"]) == pytest.approx(list(ams["level"]), abs=5e-4)
    assert list(back["carried_over"]) == list(ams["carried_over"])


# -- the job ------------------------------------------------------------------------

def test_a_job_is_completed_without_touching_what_was_given():
    given = {"fit": {"form": "logistic"}, "fsl": 215.5}
    job = R.complete_job(given)
    assert job["fit"]["form"] == "logistic" and job["fit"]["draws"] == 400
    assert given == {"fit": {"form": "logistic"}, "fsl": 215.5}


def test_the_fingerprint_moves_with_the_settings_and_with_the_input_file(tmp_path):
    path = write_hydstra(tmp_path / "x.csv", series([1.0, 2.0]))
    job = {"record": {"files": [str(path)]}}
    first = R.fingerprint(job)
    assert R.fingerprint(job) == first
    assert R.fingerprint({**job, "water_year_start": 9}) != first
    write_hydstra(path, series([1.0, 2.0, 3.0]))
    assert R.fingerprint(job) != first


# -- choosing the water year ------------------------------------------------------------

def test_every_start_month_is_scored_and_a_cut_through_a_flood_counts():
    # Three years of a lake drawing down, with a flood rising across 1 February
    # each year and peaking a day after it.
    index = pd.date_range("2000-10-01", "2003-09-30 21:00", freq="3h")
    level = pd.Series(205.0, index=index)
    for year in (2001, 2002, 2003):
        rise = (index >= f"{year}-01-31") & (index <= f"{year}-02-02")
        level[rise] = 205.0 + np.linspace(0, 3, rise.sum())
    scores = R.water_year_scores(level, high_level=204.0, full_supply=207.0)
    by_month = scores.set_index("start_month")

    assert list(scores["start_month"]) == list(R.MONTH_NAMES)
    assert (scores["carried_over_high"].fillna(0) <= scores["carried_over"]).all()
    assert by_month.loc["February", "cuts_rising"] == 3
    assert by_month.loc["February", "cuts_above_fsl"] == 0      # still below it at midnight
    # A flood peaking a day into the year reads as carried over - the reason the
    # rising cuts are counted beside it.
    assert by_month.loc["February", "carried_over"] == 3
    assert by_month.loc["February", "closest_storm_max_d"] <= 1.0
    assert by_month.loc["October", ["cuts_rising", "carried_over"]].tolist() == [0, 0]
    assert by_month.loc["October", "closest_storm_max_d"] > 100


def test_the_monthly_levels_do_not_depend_on_the_water_year():
    level = synthetic_level(years=6)
    monthly = R.monthly_levels(level)
    assert list(monthly["month"]) == [name[:3] for name in R.MONTH_NAMES]
    assert (monthly["p10"] <= monthly["median"]).all()
    assert (monthly["median"] <= monthly["p90"]).all()
