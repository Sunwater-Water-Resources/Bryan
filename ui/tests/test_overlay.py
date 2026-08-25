"""Overlaying one group's envelope against another's.

The Durations tab compares storm durations inside one group; this compares
groups - climate scenarios, dam options, antecedent storage assumptions - and
each contributes a single line, its envelope over the durations. What the
tests here mostly pin is the honesty of that picture: an envelope hides how
many durations went into it and where it stopped, so those have to be said.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from core import overlay, results
from test_results import level_curve, mc_row, write_quantiles


def build(tmp_path, groups, *, key="level", column=None, curve=None,
          methods=None):
    """``({group key: {result key: [CurveSource]}}, sims frame)``.

    The group keys are given rather than derived: ``grouping`` is tested on its
    own, and building them from output names here would tie every assertion
    about labelling to the naming convention of a fixture.
    """
    column = column or key
    rows, membership = [], {}
    for group, durations in groups.items():
        for duration in durations:
            name = f"{group}_{duration:g}h"
            values = (curve(group, duration) if curve
                      else level_curve(duration))
            write_quantiles(Path(f"{tmp_path / name}_{key}.csv"), column, values)
            row = mc_row(tmp_path, name, duration)
            if methods is not None:
                row["Method"] = methods(group, duration)
            membership.setdefault(group, []).append(len(rows))
            rows.append(row)

    frame = pd.DataFrame(rows)
    available = {}
    for group, indices in membership.items():
        found = results.sources_for_rows(frame, tmp_path, indices)
        if found:
            available[group] = found
    return available, frame


# -- labelling ---------------------------------------------------------------

def test_the_label_is_the_part_that_tells_the_groups_apart():
    labels = overlay.distinguishing_labels(
        ["TFD_mc_C030_ebf_L20-4_GWL1p3|exg", "TFD_mc_C030_ebf_L20-4_GWL1p7|exg"])
    assert list(labels.values()) == ["GWL1p3", "GWL1p7"]


def test_a_shared_prefix_is_never_cut_mid_token():
    """'GWL1p3' and 'GWL1p7' share 'GWL1p'. They must not become '3' and '7'."""
    labels = overlay.distinguishing_labels(["a_GWL1p3", "a_GWL1p7"])
    assert list(labels.values()) == ["GWL1p3", "GWL1p7"]


def test_one_group_keeps_its_whole_key():
    key = "TFD_mc_C030_GWL1p3|exg"
    assert overlay.distinguishing_labels([key]) == {key: key}


def test_keys_with_nothing_in_common_are_left_alone():
    keys = ["callide_E010", "tinaroo_H012"]
    assert overlay.distinguishing_labels(keys) == {key: key for key in keys}


def test_a_longer_key_keeps_the_part_the_shorter_one_lacks():
    labels = overlay.distinguishing_labels(["TFD_mc_exg", "TFD_mc_exg_raised"])
    assert list(labels.values()) == ["exg", "exg_raised"]


def test_trimming_that_would_leave_nothing_falls_back_to_the_full_keys():
    """Head and tail overlap here, so one key would trim away to nothing."""
    keys = ["run_exg", "run_raised_exg"]
    assert overlay.distinguishing_labels(keys) == {key: key for key in keys}


def test_the_trailing_part_comes_off_too():
    labels = overlay.distinguishing_labels(["run_GWL1p3_v2", "run_GWL2p0_v2"])
    assert list(labels.values()) == ["GWL1p3", "GWL2p0"]


# -- building ----------------------------------------------------------------

def test_each_group_contributes_one_envelope(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48, 120],
                                        "a_GWL1p7": [24, 48, 120]})
    built = overlay.build(available, list(available), "level", frame)

    assert built.labels == ["GWL1p3", "GWL1p7"]
    # and the envelope really is the maximum over that group's durations
    one = results.compare(available["a_GWL1p3"]["level"])
    expected = results.analyse(one).envelope
    pd.testing.assert_series_equal(built.frame["GWL1p3"], expected,
                                   check_names=False)


def test_a_group_without_that_result_type_is_reported_not_dropped(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48]})
    available["a_GWL1p7"] = {"inflow": []}

    built = overlay.build(available, list(available), "level", frame)

    # the label still comes off the whole selected set, so the surviving
    # curve is named the way it would be if both had drawn
    assert built.labels == ["GWL1p3"]
    assert built.problems == [("GWL1p7", "has no level results")]


def test_an_unreadable_file_is_reported_against_its_group(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48],
                                        "a_GWL1p7": [24, 48]})
    # the file exists, so it is found - but it does not hold a level column
    write_quantiles(Path(f"{tmp_path / 'a_GWL1p7_48h'}_level.csv"),
                    "something else", {2: 1.0})

    built = overlay.build(available, list(available), "level", frame)

    assert any(who.startswith("GWL1p7") for who, _why in built.problems)
    assert built.labels == ["GWL1p3", "GWL1p7"]      # 24h still carries it


def test_groups_that_stop_at_different_aeps_outer_join(tmp_path):
    def curve(group, duration):
        aeps = [2, 100, 10000] if group.endswith("1p3") else [2, 100]
        return {aep: 100.0 + aep for aep in aeps}

    available, frame = build(tmp_path, {"a_GWL1p3": [24], "a_GWL1p7": [24]},
                             curve=curve)
    built = overlay.build(available, list(available), "level", frame)

    assert list(built.frame.index) == [2, 100, 10000]
    assert pd.isna(built.frame.loc[10000, "GWL1p7"])
    assert any("do not reach the same AEP" in note for note in built.notes)


# -- the notes that keep it honest -------------------------------------------

def test_uneven_duration_coverage_is_noted(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48, 120],
                                        "a_GWL1p7": [24]})
    built = overlay.build(available, list(available), "level", frame)

    note = next(note for note in built.notes
                if "same number of durations" in note)
    assert "GWL1p3: 3" in note and "GWL1p7: 1" in note


def test_even_coverage_says_nothing_about_it(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48],
                                        "a_GWL1p7": [24, 48]})
    built = overlay.build(available, list(available), "level", frame)
    assert not any("same number of durations" in note for note in built.notes)


def test_an_envelope_pinned_to_the_end_of_its_range_is_noted(tmp_path):
    """A lower bound, so a difference measured against it is biased."""
    available, frame = build(tmp_path, {"a_GWL1p3": [12, 24, 48, 120],
                                        "a_GWL1p7": [12, 24, 48, 120]})
    built = overlay.build(available, list(available), "level", frame)

    assert any("GWL1p3: The longest duration you ran" in note
               for note in built.notes)
    assert any("GWL1p3: The shortest duration you ran" in note
               for note in built.notes)


def test_a_group_mixing_two_methods_is_noted(tmp_path):
    """Grouping usually separates them, but nothing guarantees it."""
    write_quantiles(Path(f"{tmp_path / 'a_24h'}_level.csv"), "level",
                    level_curve(24))
    routed = tmp_path / "rr"
    write_quantiles(routed / "a_48h__level_quantiles.csv", "level",
                    level_curve(48))
    frame = pd.DataFrame([
        mc_row(tmp_path, "a_24h", 24),
        pd.Series({"Method": "reservoir routing", "Output file": "a_48h",
                   "Output suffix": "", "Results folder": str(routed),
                   "Duration": 48}),
    ])
    available = {"a": results.sources_for_rows(frame, tmp_path)}

    built = overlay.build(available, ["a"], "level", frame)

    assert built.curves[0].duration_count == 2
    assert any("more than one method" in note for note in built.notes)


def test_without_the_sims_frame_no_method_is_claimed(tmp_path):
    available, _frame = build(tmp_path, {"a_GWL1p3": [24, 48]})
    built = overlay.build(available, list(available), "level")
    assert built.curves[0].methods == ()
    assert not any("more than one method" in note for note in built.notes)


# -- against a baseline ------------------------------------------------------

def _two_groups(tmp_path, key, column, base_value, other_value):
    def curve(group, duration):
        value = base_value if group.endswith("base") else other_value
        return {2: value, 100: value * 1.5 if value else 0.0}

    return build(tmp_path, {"a_base": [24], "a_warm": [24]},
                 key=key, column=column, curve=curve)


def test_level_changes_are_metres(tmp_path):
    available, frame = _two_groups(tmp_path, "level", "level", 100.0, 100.4)
    built = overlay.build(available, list(available), "level", frame)
    change = overlay.deltas(built, "base")

    assert change.kind == results.ABSOLUTE
    assert change.label == "change (m)"
    assert change.frame.loc[2, "warm"] == pytest.approx(0.4)


def test_flow_changes_are_percentages(tmp_path):
    available, frame = _two_groups(tmp_path, "inflow", "inflow", 1000.0, 1200.0)
    built = overlay.build(available, list(available), "inflow", frame)
    change = overlay.deltas(built, "base")

    assert change.kind == results.PERCENT
    assert change.frame.loc[2, "warm"] == pytest.approx(20.0)


def test_a_zero_baseline_is_reported_rather_than_infinite(tmp_path):
    """A dam that does not spill has no frequent outflow quantile."""
    def curve(group, duration):
        return ({2: 0.0, 100: 500.0} if group.endswith("base")
                else {2: 0.0, 100: 600.0})

    available, frame = build(tmp_path, {"a_base": [24], "a_warm": [24]},
                             key="outflow", curve=curve)
    built = overlay.build(available, list(available), "outflow", frame)
    change = overlay.deltas(built, "base")

    assert pd.isna(change.frame.loc[2, "warm"])
    assert change.undefined == (2,)
    assert change.frame.loc[100, "warm"] == pytest.approx(20.0)


def test_the_baseline_is_not_compared_with_itself(tmp_path):
    available, frame = _two_groups(tmp_path, "level", "level", 100.0, 100.4)
    built = overlay.build(available, list(available), "level", frame)
    assert list(overlay.deltas(built, "base").frame.columns) == ["warm"]


def test_one_group_alone_has_nothing_to_compare_against(tmp_path):
    available, frame = build(tmp_path, {"a_base": [24]})
    built = overlay.build(available, list(available), "level", frame)
    assert overlay.deltas(built, "a_base").is_empty


# -- the critical duration, per group ----------------------------------------

def test_the_critical_duration_frame_carries_one_column_per_group(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [12, 24, 48, 120],
                                        "a_GWL1p7": [12, 24, 48, 120]})
    built = overlay.build(available, list(available), "level", frame)
    critical, skipped = overlay.critical_frame(built)

    assert list(critical.columns) == ["GWL1p3", "GWL1p7"]
    assert not skipped
    # level goes long while the storage fills and short on the rare tail
    assert critical["GWL1p3"].iloc[0] > critical["GWL1p3"].iloc[-1]


def test_a_group_whose_durations_are_unknown_is_left_out(tmp_path):
    available, frame = build(tmp_path, {"a_GWL1p3": [24, 48],
                                        "a_GWL1p7": [24, 48]})
    curves = overlay.build(available, list(available), "level", frame).curves
    built = overlay.Overlay(
        frame=pd.DataFrame({"GWL1p3": curves[0].envelope}),
        curves=(curves[0],
                overlay.GroupCurve(key="x", label="GWL1p7",
                                   envelope=curves[1].envelope,
                                   critical=curves[1].critical,
                                   durations={"24h": None, "48h": None})),
        key="level")
    critical, skipped = overlay.critical_frame(built)

    assert list(critical.columns) == ["GWL1p3"]
    assert skipped == ["GWL1p7"]


# -- the table ---------------------------------------------------------------

def test_the_table_carries_the_envelopes_then_the_changes(tmp_path):
    available, frame = _two_groups(tmp_path, "level", "level", 100.0, 100.4)
    built = overlay.build(available, list(available), "level", frame)
    table = overlay.table(built, overlay.deltas(built, "base"))

    assert list(table.columns) == ["base", "warm", "warm change (m)"]
