"""The inflow record: stage 1's derived inflow as annual maxima and hydrographs."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lib.homogenise import inflow

BRYAN_ROOT = Path(__file__).resolve().parents[1]
CALLIDE = BRYAN_ROOT.parent / "callide-fsl-reinstate"


def derived(rates, step="1h", start="2000-01-01", evaporation_ml=0.0, unaccounted=None,
            level=215.0):
    """A ``derive_inflow`` frame: row T holds the interval ending at T."""
    stamps = pd.date_range(start, periods=len(rates) + 1, freq=step)
    seconds = pd.Timedelta(step).total_seconds()
    inflow_ml = np.r_[0.0, np.asarray(rates, float) * seconds / 1000.0]
    frame = pd.DataFrame({
        "dt_s": np.r_[0.0, np.full(len(rates), seconds)],
        "Inflow_ML": inflow_ml + np.r_[0.0, np.full(len(rates), evaporation_ml)],
        "Evaporation_ML": np.r_[0.0, np.full(len(rates), evaporation_ml)],
        "Release_ML": 0.0, "Level": level, "Interpolated": False,
        "Unaccounted_Release_ML": 0.0}, index=stamps)
    if unaccounted is not None:
        frame["Unaccounted_Release_ML"] = np.r_[0.0, unaccounted]
    return frame


def test_intervals_are_stamped_at_their_midpoints_and_keep_the_volume():
    frame = inflow.intervals(derived([10.0, 20.0, 30.0]))
    assert list(frame.index) == list(pd.date_range("2000-01-01 00:30", periods=3, freq="1h"))
    assert frame["Inflow_m3s"].tolist() == pytest.approx([10.0, 20.0, 30.0])
    assert frame["Volume_ML"].sum() == pytest.approx(60.0 * 3600 / 1000)


def test_the_evaporation_is_taken_back_out_where_the_inflow_leaves_it_out():
    source = derived([10.0, 10.0], evaporation_ml=3.6)                 # 1 m3/s of evaporation
    assert inflow.intervals(source)["Inflow_m3s"].tolist() == pytest.approx([11.0, 11.0])
    assert inflow.intervals(source, evaporation=False)["Inflow_m3s"].tolist() == \
        pytest.approx([10.0, 10.0])


def test_the_recession_keeps_both_inflows_and_flags_where_the_release_is_uncertain():
    source = derived([50.0, 0.0, 0.0], unaccounted=[0.0, 72.0, 0.0])   # 20 m3/s booked as release
    frame = inflow.intervals(source)
    assert frame["Inflow_m3s"].tolist() == pytest.approx([50.0, 0.0, 0.0])
    assert frame["Inflow_uncorrected_m3s"].tolist() == pytest.approx([50.0, 20.0, 0.0])
    assert frame["Release_uncertain"].tolist() == [False, True, False]
    raw = inflow.intervals(source, recession_correction=False)
    assert raw["Inflow_m3s"].tolist() == pytest.approx([50.0, 20.0, 0.0])


def test_the_peak_is_a_time_average_so_a_one_minute_step_counts_for_a_minute():
    rates = [10.0] * 18
    source = derived(rates, step="10min")
    source.loc[source.index[9], "Inflow_ML"] += 60.0                   # 100 m3/s for ten minutes
    frame = inflow.with_smoothing(inflow.intervals(source), "1h")
    assert frame["Inflow_native_m3s"].max() == pytest.approx(110.0)
    assert frame["Inflow_m3s"].max() == pytest.approx(10.0 + 100.0 / 6, rel=1e-6)


def test_burst_volumes_are_the_largest_over_each_duration_in_each_water_year():
    rates = [0.0] * 48 + [100.0] * 24 + [0.0] * 48                     # a day at 100 m3/s
    frame = inflow.intervals(derived(rates, start="2000-01-01"))
    volumes = inflow.burst_volumes(frame, durations_h=[24, 48], start_month=10)
    assert volumes.loc[2000, "Volume_24h_ML"] == pytest.approx(8640.0)
    assert volumes.loc[2000, "Volume_48h_ML"] == pytest.approx(8640.0)


def test_annual_maxima_give_depths_and_the_rain_over_the_burst_and_the_day_before():
    rates = [0.0] * 48 + [100.0] * 24 + [0.0] * 48
    frame = inflow.intervals(derived(rates, start="2000-01-01"))
    rain = pd.Series(10.0, index=pd.date_range("1999-12-25", "2000-01-15", freq="D"))
    ams = inflow.annual_maxima(frame, [24], 10, catchment_km2=86.4, rainfall=rain)
    row = ams.loc[2000]
    assert row["Peak_inflow_m3s"] == pytest.approx(100.0)
    assert row["Depth_24h_mm"] == pytest.approx(100.0)                  # 8,640 ML over 86.4 km2
    assert row["Rain_24h_mm"] == pytest.approx(30.0)                    # two days and the one before
    assert not row["Complete"]                                          # five days of a year


def test_a_hydrograph_carries_the_flags_and_its_events_are_named_by_year_and_day():
    source = derived([50.0, 0.0, 0.0], unaccounted=[0.0, 72.0, 0.0])
    frame = inflow.with_smoothing(inflow.intervals(source), "1h")
    part = inflow.hydrograph(frame, "2000-01-01", "2000-01-02", "1h")
    assert {"Inflow_m3s", "Inflow_smoothed_m3s", "Inflow_uncorrected_m3s", "Release_m3s",
            "Level_m", "Release_uncertain"} <= set(part.columns)
    ams = inflow.annual_maxima(frame, [1], 10)
    [(name, start, end)] = inflow.event_windows(ams, 1, before_days=1, after_days=2)
    assert name.endswith("_peak_20000101")
    assert end - start == pd.Timedelta(days=3)


# -- held to the independent reverse routing ---------------------------------------------

@pytest.mark.skipif(not (CALLIDE / "reverse_routing" / "run.py").is_file(),
                    reason="needs the callide-fsl-reinstate checkout beside Bryan")
def test_callide_s_inflow_agrees_with_the_reverse_routing_package(tmp_path):
    data = CALLIDE / "data"
    homogenise = {
        "gauges": [str(data / "gauge" / "130314A.csv"), str(data / "gauge" / "130314C.csv")],
        "overlay": {"file": str(data / "gauge" / "130314B_hydstra.csv"), "below": 200.90},
        "storage": str(data / "storage" / "CALLIDE_STORAGE.els"),
        "register": str(data / "ratings" / "RatingCurves.xlsx"),
        "evaporation": str(data / "climate" / "silo_-24.35_150.65.txt")}
    job = tmp_path / "inflow.json"
    job.write_text(json.dumps({"homogenise": homogenise, "evaporation": False,
                               "catchment_km2": 519, "top_events": 2,
                               "out": str(tmp_path / "inflow")}), encoding="utf-8")
    subprocess.run([sys.executable, str(BRYAN_ROOT / "util" / "InflowRecord.py"), str(job)],
                   check=True, capture_output=True, timeout=900)
    summary = json.loads((tmp_path / "inflow" / "summary.json").read_text(encoding="utf-8"))
    assert summary["record"]["end"] > "2026-02-10"          # past the SILO file's last day
    assert len(summary["hydrographs"]) == 2

    python = CALLIDE / ".venv" / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    theirs = tmp_path / "rr.csv"
    subprocess.run([str(python if python.exists() else sys.executable), "-c",
                    "import sys; from reverse_routing import sources, run; "
                    "run.build(sources.load_level_record())['ams_primary']"
                    f".to_csv(sys.argv[1])", str(theirs)],
                   cwd=CALLIDE, check=True, capture_output=True, timeout=900)
    rr = pd.read_csv(theirs, index_col=0)
    ours = pd.read_csv(tmp_path / "inflow" / "inflow_ams.csv").set_index("Period").reindex(rr.index)
    for column in ("Peak_inflow_m3s", "Volume_24h_ML", "Volume_72h_ML"):
        ratio = ours[column] / rr[column] - 1
        assert abs(ratio.median()) < 1e-4, column
    big = rr["Peak_inflow_m3s"] > 500                        # the floods: within the rating shape
    assert (ours.loc[big, "Volume_72h_ML"] / rr.loc[big, "Volume_72h_ML"] - 1).abs().max() < 0.05
