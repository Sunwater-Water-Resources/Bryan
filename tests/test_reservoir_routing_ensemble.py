"""Reservoir routing with ensemble input.

The method was built around a Monte Carlo database, and since 4 August 2026 it
takes an ensemble one too: it tells the two apart by their columns (_detect_scheme),
takes one antecedent volume from the sims list for every event, and re-analyses the
routed peaks with lib/EnbAnalysis.py - the median pattern of each duration, then the
critical duration at each AEP. Checked on Callide against URBS (the E007 ensemble to
5 mm, and the E013 gate failure PMFs), but nothing here pinned it until now.

The ensemble below is synthetic. Every pattern of a duration has the same triangular
hydrograph scaled by its own factor, so the patterns' order - and so the median pattern
- is known without routing anything, and the 48 h storm carries the most volume, so it
must be critical. A linear storage-level curve and a linear rating keep the routing
simple enough to reason about.

This is a **Bryan** test - it imports lib/ReservoirRouting.py, which needs scipy and
matplotlib. Run it with Bryan's own interpreter:

    python -m pytest tests -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from lib.ReservoirRouting import ReservoirRoutingSimulator                    # noqa: E402

AEPS = (100, 1000, 10000)
DURATIONS = (12, 24, 48)
PATTERNS = 10
# each pattern's scale factor: shuffled, so the median pattern is not simply the middle id
SCALE = dict(zip(range(PATTERNS), (1.00, 1.30, 0.80, 1.10, 0.90, 1.25, 0.95, 1.05, 1.20, 0.85)))
FSL, SOURCE_ADV = 110.0, 90_000.0          # the source run started below full supply
STEP_HOURS = 1.0


def fsv():
    return 10_000.0 * (FSL - 100.0)        # the ELS below: 10,000 ML per metre from 100 m


def hydrograph(aep, duration, pattern, hours):
    """A triangle rising over the storm and falling over twice its length. Rarer and longer
    storms are bigger; the pattern scales the whole thing."""
    peak = 200.0 * np.log10(aep) * (duration / 12.0) ** 0.5 * SCALE[pattern]
    t = np.arange(0.0, hours + STEP_HOURS, STEP_HOURS)
    rise = np.clip(t / duration, 0, 1)
    fall = np.clip(1 - (t - duration) / (2 * duration), 0, 1)
    return np.where(t <= duration, rise, fall) * peak


def build(tmp_path, *, adv="fsv", suffix="test", run=True, analyse=True):
    """An ensemble database and its inflow hydrographs, the ELS and SQ, and the sims row."""
    events, columns = [], {}
    longest = 3 * max(DURATIONS) + 24
    for aep in AEPS:
        for duration in DURATIONS:
            for pattern in range(PATTERNS):
                name = f"sim_{len(events)}"
                flow = hydrograph(aep, duration, pattern, 3 * duration + 12)
                series = np.full(int(longest / STEP_HOURS) + 1, np.nan)
                series[:len(flow)] = flow          # shorter storms end early, as URBS writes them
                columns[name] = series
                events.append({"rain_aep": aep, "duration": duration, "tp": pattern,
                               "storm_method": "GSDM", "ADV": SOURCE_ADV,
                               # the source run's own peaks, which the re-route must replace
                               "inflow": -1.0, "level": -1.0, "outflow": -1.0})
    inflows = pd.DataFrame(columns, index=np.arange(0.0, longest + STEP_HOURS, STEP_HOURS))
    inflows.index.name = "time"
    database = pd.DataFrame(events, index=list(columns))

    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    inflows.to_csv(inputs / "enb_inflows.csv")
    database.to_csv(inputs / "enb_database.csv")
    pd.DataFrame({"EL": np.arange(100.0, 131.0), "V": 10_000.0 * np.arange(0.0, 31.0)}) \
        .to_csv(inputs / "dam.els", index=False)
    pairs = [(s, 0.02 * s) for s in range(0, 200_001, 10_000)]     # 200 m3/s per 10,000 ML
    (inputs / "dam.sq").write_text("TEST DAM\n*\n*\n*\n" + f"{len(pairs)} PAIRS:\n"
                                   + "".join(f"{s} {q}\n" for s, q in pairs))

    row = pd.Series({
        "Include": "yes", "Method": "reservoir routing", "Output file": "CLD_enb_rr",
        "Output suffix": suffix, "Run models": "yes" if run else "no",
        "Analyse results": "yes" if analyse else "no", "Store hydrographs": "yes",
        "Input MCDF": str(inputs / "enb_database.csv"), "Inflow": str(inputs / "enb_inflows.csv"),
        "ELS file": str(inputs / "dam.els"), "SQ file": str(inputs / "dam.sq"), "FSL": FSL,
        "ADV": adv, "Hydrographs folder": str(tmp_path / "hydrographs"),
        "Results folder": str(tmp_path / "results"), "Config file": np.nan, "Log file": np.nan,
    })
    return row


def route(row):
    """Run the simulator, then put stdout back and close its log, as Main.py does."""
    try:
        return ReservoirRoutingSimulator(row, {})
    finally:
        logger, sys.stdout = sys.stdout, sys.__stdout__
        if hasattr(logger, "close"):
            logger.close()


def routed(tmp_path, suffix="test"):
    return pd.read_csv(tmp_path / "results" / f"CLD_enb_rr_{suffix}.csv", index_col=0)


def critical(tmp_path, suffix="test"):
    return pd.read_csv(tmp_path / "results" / "csv" / f"CLD_enb_rr_{suffix}_critical.csv", index_col=0)


# -- which scheme, and where it is written --------------------------------------------

def test_an_ensemble_database_is_routed_as_an_ensemble(tmp_path):
    simulator = route(build(tmp_path))
    assert simulator.scheme == "ensemble"
    results = tmp_path / "results"
    # the ensemble's own naming, beside plots/ and csv/ - not the Monte Carlo '__mcdf'
    assert (results / "CLD_enb_rr_test.csv").is_file()
    assert not list(results.glob("*__mcdf*"))
    assert (results / "csv" / "CLD_enb_rr_test_critical.csv").is_file()
    assert (results / "plots" / "CLD_enb_rr_test_1000_bp.png").is_file()
    for series in ("outflows", "levels", "volumes"):
        assert (tmp_path / "hydrographs" / f"CLD_enb_rr_{series}_test.csv").is_file()


def test_the_source_runs_peaks_are_replaced_and_its_other_columns_kept(tmp_path):
    route(build(tmp_path))
    frame = routed(tmp_path)
    assert (frame[["inflow", "level", "outflow"]] > 0).all().all()
    assert list(frame["tp"]) == [e % PATTERNS for e in range(len(frame))]
    # the inflow peak is the hydrograph's own, trailing NaNs of the short storms and all
    first = frame.iloc[0]
    expected = hydrograph(first["rain_aep"], first["duration"], first["tp"], 3 * first["duration"] + 12).max()
    assert first["inflow"] == pytest.approx(expected)


def test_routing_attenuates_and_the_lake_never_falls_below_where_it_started(tmp_path):
    route(build(tmp_path))
    frame = routed(tmp_path)
    assert (frame["outflow"] < frame["inflow"]).all()
    assert (frame["level"] > FSL).all()


# -- the antecedent volume ------------------------------------------------------------

def test_fsv_starts_every_event_at_the_routed_curves_full_supply(tmp_path):
    route(build(tmp_path, adv="fsv"))
    frame = routed(tmp_path)
    assert (frame["ADV"] == fsv()).all()
    assert (frame["ADV_input"] == SOURCE_ADV).all()      # the source run's, kept for checking


def test_database_starts_every_event_where_the_source_run_did(tmp_path):
    route(build(tmp_path, adv="database"))
    frame = routed(tmp_path)
    assert (frame["ADV"] == SOURCE_ADV).all()
    # starting 1 m below full supply, the same storms peak lower
    route(build(tmp_path, adv="fsv", suffix="fsv"))
    assert (frame["level"].to_numpy() < routed(tmp_path, "fsv")["level"].to_numpy()).all()


def test_a_volume_in_ml_is_used_as_given(tmp_path):
    route(build(tmp_path, adv=95_000))
    assert (routed(tmp_path)["ADV"] == 95_000).all()


def test_a_varying_adv_is_refused_for_an_ensemble(tmp_path):
    with pytest.raises(ValueError, match="cannot be used with ensemble input"):
        route(build(tmp_path, adv="varying"))


# -- the analysis ---------------------------------------------------------------------

def median_pattern():
    """The pattern EnbAnalysis takes as the median: position round(n / 2) in ascending order.
    Every pattern of a duration is the same hydrograph scaled, so the routed peaks sort as
    the scale factors do."""
    ranked = sorted(SCALE, key=SCALE.get)
    return ranked[int(np.around(PATTERNS / 2, 0))]


def test_the_critical_event_is_the_median_pattern_of_the_critical_duration(tmp_path):
    route(build(tmp_path))
    table = critical(tmp_path)
    assert sorted(table.index) == sorted(AEPS)
    for result in ("inflow", "level", "outflow"):
        # the 48 h storm has the largest peak and the most volume at every AEP
        assert (table[f"{result}_duration"] == max(DURATIONS)).all(), result
        assert (table[f"{result}_tp"] == f"GSDM: {median_pattern()}").all(), result


def test_the_critical_values_are_the_routed_peaks_of_those_events(tmp_path):
    route(build(tmp_path))
    frame, table = routed(tmp_path), critical(tmp_path)
    for aep in AEPS:
        event = frame[(frame["rain_aep"] == aep) & (frame["duration"] == max(DURATIONS))
                      & (frame["tp"] == median_pattern())].iloc[0]
        for result in ("inflow", "level", "outflow"):
            assert table.loc[aep, result] == pytest.approx(event[result])


def test_an_analysis_only_row_finds_the_routed_ensemble_again(tmp_path):
    route(build(tmp_path))
    before = critical(tmp_path)
    (tmp_path / "results" / "csv" / "CLD_enb_rr_test_critical.csv").unlink()
    simulator = route(build(tmp_path, run=False))
    assert simulator.scheme == "ensemble"
    pd.testing.assert_frame_equal(critical(tmp_path), before)
