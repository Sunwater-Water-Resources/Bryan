"""Representative events from an ensemble: the PMF is the highest event, a design AEP the
median pattern at the critical duration, and the event is named by its row in its own file."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from core import enbevents

PMP = 1_900_000.0


def ensemble_frame(aeps=(PMP,), durations=(4.5, 9.0, 12.0), patterns=10, start=0) -> pd.DataFrame:
    """Every pattern of every duration at each AEP, as an ensemble database holds them.

    Levels rise with the pattern number, and the 9 h duration is the highest at every
    pattern, so the median and the highest are both 9 h - except that one 4.5 h pattern
    is the single highest event of all, which is what makes the two rules differ."""
    rows = []
    for aep in aeps:
        for duration in durations:
            for tp in range(patterns):
                level = 220.0 + tp * 0.01 + (0.5 if duration == 9.0 else 0.0)
                if duration == 4.5 and tp == 7:
                    level = 222.0
                rows.append({"rain_aep": aep, "duration": duration, "tp": tp,
                             "storm_method": "GSDM", "level": level, "inflow": level * 10,
                             "outflow": level * 9, "ADV": 129041.0})
    frame = pd.DataFrame(rows)
    frame.index = range(start, start + len(frame))
    return frame


def sources(tmp_path, *frames) -> list:
    out = []
    for number, frame in enumerate(frames):
        path = tmp_path / f"CLD_enb_{number}.csv"
        frame.to_csv(path)
        out.append(enbevents.Source(output_name=f"sims_enb\\results\\CLD_enb_{number}", path=path))
    return out


@pytest.fixture(autouse=True)
def clear():
    enbevents.forget_cached()
    yield
    enbevents.forget_cached()


def test_the_pmf_is_the_highest_event_of_any_duration(tmp_path):
    found = sources(tmp_path, ensemble_frame())
    frame = enbevents.load(found)
    choice = enbevents.choose(frame, enbevents.EVENTS.Target(kind="pmf"), found)
    assert choice.found and choice.rule == "highest"
    assert (choice.duration, choice.event, choice.values["level"]) == (4.5, 7, 222.0)


def test_a_design_aep_is_the_median_pattern_at_the_critical_duration(tmp_path):
    found = sources(tmp_path, ensemble_frame())
    frame = enbevents.load(found)
    choice = enbevents.choose(frame, enbevents.EVENTS.Target(kind="aep", value=PMP), found)
    # 9 h is critical; of ten patterns the median is position round(10 / 2) = 5, the sixth
    assert choice.rule == "median" and choice.duration == 9.0
    assert choice.event == 10 + 5 and choice.pattern == "GSDM: 5"
    rows = enbevents.duration_rows(choice)
    assert [row["duration"] for row in rows if row["chosen"]] == [9.0]
    # numbers, so the page's table sorts 4.5 h before 12 h rather than after it as text
    assert all(isinstance(row[key], float) for row in rows
               for key in ("duration", "median", "highest"))
    assert all(column.get("sortable") for column in enbevents.DURATION_COLUMNS)


def test_a_duration_overrides_the_critical_one(tmp_path):
    found = sources(tmp_path, ensemble_frame())
    frame = enbevents.load(found)
    target = enbevents.EVENTS.Target(kind="aep", value=PMP, source="12")
    choice = enbevents.choose(frame, target, found)
    assert choice.duration == 12.0 and choice.event == 20 + 5


def test_an_aep_the_ensemble_did_not_run_is_reported(tmp_path):
    found = sources(tmp_path, ensemble_frame())
    choice = enbevents.choose(enbevents.load(found),
                              enbevents.EVENTS.Target(kind="aep", value=1000), found)
    assert not choice.found and "did not run 1 in 1,000" in choice.problem


def test_the_event_is_its_row_in_its_own_file_when_a_group_has_several(tmp_path):
    first = ensemble_frame(aeps=(1000.0,))
    second = ensemble_frame(aeps=(PMP,), start=0)
    found = sources(tmp_path, first, second)
    frame = enbevents.load(found)
    choice = enbevents.choose(frame, enbevents.EVENTS.Target(kind="pmf"), found)
    # both files have the 222.0 event at row 7; the first wins the tie, and is named as such
    assert choice.event == 7 and choice.source == found[0]
    pmpf = enbevents.choose(frame, enbevents.EVENTS.Target(kind="aep", value=PMP), found)
    assert pmpf.event == 15 and pmpf.source == found[1]


def test_a_new_pmf_group_starts_with_the_pmf(tmp_path):
    # the PMPF is the Monte Carlo design flood at the AEP of the PMP, not an ensemble pick
    defaults = enbevents.default_targets(ensemble_frame())
    assert [(t.kind, t.name) for t in defaults] == [("pmf", "PMF")]
    several = enbevents.default_targets(ensemble_frame(aeps=(100.0, 1000.0, 1e4, 1e5)))
    assert [t.value for t in several] == [1000.0, 1e4, 1e5]


def test_what_is_saved_names_the_event_its_duration_and_its_rule(tmp_path):
    found = sources(tmp_path, ensemble_frame())
    choice = enbevents.choose(enbevents.load(found),
                              enbevents.EVENTS.Target(kind="pmf", notional_aep=8e6), found)
    stored = enbevents.stored_target(choice, lambda path: Path(path).name)
    assert (stored.picked, stored.duration, stored.rule) == (7, 4.5, "highest")
    assert stored.database == "CLD_enb_0.csv"
    assert stored.output_file == "sims_enb\\results\\CLD_enb_0"
    row, = enbevents.summary_rows([choice])
    assert row["hydrograph"] == "sim_00007" and row["notional AEP"] == "1 in 8,000,000"
