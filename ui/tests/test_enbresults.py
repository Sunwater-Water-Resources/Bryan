"""The Ensemble page's core: median patterns by AEP and duration, and the checks.

Pinned here: EnbAnalysis's median position and critical-duration rule (a tie
goes to the first duration listed), each result type's own critical duration,
the margin in metres for lake level, which rows count as ensemble results, and
the check against Bryan's ``_critical.csv``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import pytest

from conftest import MONTE_CARLO_COLUMNS, write_workbook
from core import enbresults, results

DURATIONS = (12.0, 24.0, 48.0)
PATTERNS = range(10)

# Each AEP and duration's base value; the patterns add tp / 100 (level) or tp * 10
# (flows), so the median - the sixth of ten ascending - is base + 0.05 or base + 50.
LEVEL = {100: {12: 215.0, 24: 215.3, 48: 215.2},       # 24h critical, 0.10 m clear
         1000: {12: 216.4, 24: 216.2, 48: 216.0}}      # 12h critical, 0.20 m clear
INFLOW = {100: {12: 800, 24: 700, 48: 600},             # inflow's own: 12h throughout
          1000: {12: 1500, 24: 1300, 48: 1100}}
OUTFLOW = {100: {12: 0, 24: 0, 48: 0},                  # a tie: goes to 12h, listed first
           1000: {12: 900, 24: 950, 48: 700}}

# What lib/EnbAnalysis.py writes for that database.
BRYAN_CRITICAL = pd.DataFrame(
    {"inflow": [850.0, 1550.0], "inflow_duration": [12.0, 12.0],
     "inflow_tp": ["ARR areal: 5"] * 2,
     "level": [215.35, 216.45], "level_duration": [24.0, 12.0],
     "level_tp": ["ARR areal: 5"] * 2,
     "outflow": [0.0, 1000.0], "outflow_duration": [12.0, 24.0],
     "outflow_tp": ["ARR areal: 5"] * 2},
    index=[100, 1000])


def database(level=LEVEL):
    rows = []
    for aep in (100, 1000):
        for duration in DURATIONS:
            for tp in PATTERNS:
                flows = 0 if OUTFLOW[aep][duration] == 0 else tp * 10
                rows.append({"rain_aep": aep, "duration": duration, "tp": tp,
                             "storm_method": "ARR areal",
                             "inflow": INFLOW[aep][duration] + tp * 10,
                             "level": level[aep][duration] + tp / 100,
                             "outflow": OUTFLOW[aep][duration] + flows})
    return pd.DataFrame(rows)


def write_run(folder: Path, name="CLD_enb_test", *, critical=True) -> Path:
    path = folder / "sims_enb" / "results" / f"{name}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    database().to_csv(path)
    if critical:
        (path.parent / "csv").mkdir(exist_ok=True)
        BRYAN_CRITICAL.to_csv(path.parent / "csv" / f"{name}_critical.csv")
    return path


def source(path) -> enbresults.Source:
    return enbresults.Source(label=path.stem, path=path)


# -- the medians -----------------------------------------------------------------

def test_the_median_is_the_sixth_of_ten_ascending():
    found = enbresults.medians(database(), "level")
    frame = found.comparison.frame
    assert list(frame.columns) == ["12h", "24h", "48h"]
    assert frame.loc[100.0, "24h"] == pytest.approx(215.35)
    assert frame.loc[1000.0, "12h"] == pytest.approx(216.45)
    assert found.patterns.loc[100.0, "24h"] == "ARR areal: 5"
    assert found.counts.loc[100.0, "24h"] == 10
    assert found.comparison.key == "level"
    assert found.comparison.durations == {"12h": 12.0, "24h": 24.0, "48h": 48.0}


def test_the_critical_duration_and_the_level_margin_in_metres():
    found = enbresults.medians(database(), "level")
    analysis = results.analyse(found.comparison)
    assert analysis.critical.to_dict() == {100.0: "24h", 1000.0: "12h"}
    assert analysis.margin_label == "margin (m)"
    assert analysis.margin.loc[100.0] == pytest.approx(0.10)
    assert analysis.margin.loc[1000.0] == pytest.approx(0.20)


def test_each_result_type_has_its_own_critical_duration():
    """As EnbAnalysis: the inflow's critical duration is not the level's."""
    found = enbresults.medians(database(), "inflow")
    analysis = results.analyse(found.comparison)
    assert analysis.critical.to_dict() == {100.0: "12h", 1000.0: "12h"}
    assert analysis.margin_label == "margin %"


def test_a_tie_goes_to_the_first_duration_listed():
    """np.argmax over the medians in the order the database lists the durations."""
    found = enbresults.medians(database(), "outflow")
    analysis = results.analyse(found.comparison)
    assert analysis.critical.loc[100.0] == "12h"
    assert analysis.critical.loc[1000.0] == "24h"


def test_the_highest_event_at_each_aep():
    top = enbresults.highest(database(), "level")
    assert top.loc[100.0, "highest"] == pytest.approx(215.39)
    assert top.loc[100.0, "highest duration"] == "24h"
    assert top.loc[1000.0, "highest pattern"] == "ARR areal: 9"


def test_the_table_carries_the_median_pattern_and_the_highest_event():
    found = enbresults.medians(database(), "level")
    analysis = results.analyse(found.comparison)
    table = enbresults.critical_table(found, analysis, enbresults.highest(database()))
    assert table.loc[100.0, "critical duration"] == "24h"
    assert table.loc[100.0, "median pattern"] == "ARR areal: 5"
    assert table.loc[1000.0, "highest"] == pytest.approx(216.49)


# -- reading -------------------------------------------------------------------------

def test_failed_events_and_overlapping_databases_are_reported(tmp_path):
    first = write_run(tmp_path, "a", critical=False)
    second = write_run(tmp_path, "b", critical=False)
    frame = pd.read_csv(first, index_col=0)
    frame.loc[0, "level"] = float("nan")
    frame.to_csv(first)
    loaded = enbresults.load([source(first), source(second)])
    notes = " ".join(loaded.notes)
    assert "1 event(s) have no result" in notes
    assert "more than one of the group's databases" in notes
    assert loaded.aeps == [100.0, 1000.0]


def test_a_monte_carlo_database_is_refused(tmp_path):
    path = tmp_path / "mc.csv"
    pd.DataFrame({"m": [0], "n": [0], "level": [215.0]}).to_csv(path)
    with pytest.raises(enbresults.EnsembleError, match="is this an ensemble run"):
        enbresults.load([source(path)])


# -- which rows ----------------------------------------------------------------------

def project(tmp_path, rows):
    write_workbook(tmp_path / "sims.xlsx", MONTE_CARLO_COLUMNS + ["Results folder"], rows)
    config = tmp_path / "sims_config.json"
    config.write_text(json.dumps({"simulation_list": "sims.xlsx", "filepaths": {}}))
    from core.config import load_sims_config
    from core.simslist import read_sims_list

    class Project:                                      # what sources_by_group reads
        def __init__(self):
            self.config = load_sims_config(config)
            self.frame = read_sims_list(self.config.sims_list_path).frame

        def group_keys(self):
            from core import grouping
            return grouping.add_group_keys(self.frame)

    return Project()


def test_ensemble_rows_and_re_routed_ensembles_are_found(tmp_path):
    write_run(tmp_path, "CLD_enb_test")
    routed = tmp_path / "routed" / "CLD_enb_rerouted.csv"
    routed.parent.mkdir()
    database().to_csv(routed)
    mc = tmp_path / "routed" / "CLD_mc_rerouted__mcdf.csv"
    pd.DataFrame({"m": [0], "n": [0], "level": [215.0]}).to_csv(mc)
    rows = [
        {"Include": "yes", "Method": "ensemble", "Output file":
         "sims_enb\\results\\CLD_enb_test", "Config file": "enb.json"},
        {"Include": "yes", "Method": "reservoir routing", "Output file":
         "CLD_enb_rerouted", "Results folder": "routed", "Config file": ""},
        {"Include": "yes", "Method": "reservoir routing", "Output file":
         "CLD_mc_rerouted", "Results folder": "routed", "Config file": "mc.json"},
    ]
    found = enbresults.sources_by_group(project(tmp_path, rows))
    names = sorted(s.path.name for sources in found.values() for s in sources)
    assert names == ["CLD_enb_rerouted.csv", "CLD_enb_test.csv"]


# -- the check against Bryan ---------------------------------------------------------

def _analysis(path, result="level"):
    loaded = enbresults.load([source(path)])
    return results.analyse(enbresults.medians(loaded.frame, result).comparison)


@pytest.mark.parametrize("result", ["level", "inflow", "outflow"])
def test_the_page_agrees_with_bryans_critical_csv(tmp_path, result):
    path = write_run(tmp_path)
    found = enbresults.check([source(path)], _analysis(path, result), result)
    assert found.agrees and found.compared == 2


def test_a_difference_from_bryan_is_listed_and_a_stale_csv_said_so(tmp_path):
    path = write_run(tmp_path)
    csv = enbresults.critical_csv(source(path))
    changed = BRYAN_CRITICAL.copy()
    changed.loc[100, "level_duration"] = 48.0
    changed.to_csv(csv)
    os.utime(csv, (path.stat().st_mtime - 60, path.stat().st_mtime - 60))
    found = enbresults.check([source(path)], _analysis(path), "level")
    assert not found.agrees
    assert found.differences == ["1 in 100: Bryan 215.35 at 48h, here 215.35 at 24h"]
    assert found.older


def test_no_check_without_the_csv_or_with_two_databases(tmp_path):
    path = write_run(tmp_path, critical=False)
    assert enbresults.check([source(path)], _analysis(path), "level") is None
    other = write_run(tmp_path, "other")
    assert enbresults.check([source(path), source(other)], _analysis(path), "level") is None
