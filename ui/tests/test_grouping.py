"""Deriving groups, and the guards that keep the derivation honest.

The three worked examples are the real conventions from the Tinaroo and
Callide projects - see the module docstring in core/grouping.py.
"""

from __future__ import annotations

import pandas as pd
import pytest

from core.grouping import (FROM_COLUMN, FROM_FALLBACK, FROM_NAME,
                           FROM_NAME_STRIPPED, add_group_keys,
                           derived_group_report, group_key,
                           strip_duration_token)


def row(**values):
    return pd.Series(values)


# --- the duration token ---------------------------------------------------

def test_token_is_stripped_when_it_is_the_rows_duration():
    name, stripped = strip_duration_token("TFD_mc_C030_ebf_18h_L20-4_GWL1p3", 18)
    assert (name, stripped) == ("TFD_mc_C030_ebf_L20-4_GWL1p3", True)


def test_token_survives_when_it_is_not_the_rows_duration():
    """A fixed 24h in a model name must not be mistaken for the varying part."""
    assert strip_duration_token("model_24h_fixed_GWL1p7", 18) == (
        "model_24h_fixed_GWL1p7", False)


def test_a_name_with_no_token_is_untouched():
    """The Callide convention: the duration is simply not in the name."""
    name = "E010_no-pbp_L105-1.7_GWL0p3_RGN"
    assert strip_duration_token(name, 120) == (name, False)


def test_gwl_and_layout_numbers_are_not_durations():
    for name in ("GWL1p7", "L20-4", "C030", "L105-1.7"):
        assert strip_duration_token(name, 20) == (name, False)


@pytest.mark.parametrize("text,duration", [
    ("run_1.5h_x", 1.5), ("run_1_5h_x", 1.5), ("run_24 hr_x", 24),
    ("run_24hrs_x", 24), ("run_24hour_x", 24),
])
def test_token_spellings(text, duration):
    _, stripped = strip_duration_token(text, duration)
    assert stripped


def test_float_and_int_durations_match():
    assert strip_duration_token("a_120h_b", 120.0)[1]
    assert strip_duration_token("a_120h_b", 120)[1]


# --- whole rows -----------------------------------------------------------

def test_group_column_wins_verbatim():
    key = group_key(row(Group="my group", **{"Output file": "x_18h"}, Duration=18))
    assert key.key == "my group"
    assert key.provenance == FROM_COLUMN


def test_tinaroo_reservoir_routing_convention():
    keys = {
        group_key(row(**{"Output file": f"TFD_mc_C030_ebf_{d}h_L20-4_GWL1p3",
                         "Output suffix": "exg", "Duration": d,
                         "Results folder": r"sims_mc\results"})).key
        for d in (18, 24, 36, 48, 72, 96)
    }
    assert len(keys) == 1, "six durations must collapse to one group"


def test_callide_convention_without_a_duration_in_the_name():
    keys = {
        group_key(row(Basename="E010_no-pbp_L105-1.7_GWL0p3_RGN",
                      **{"Output file": "E010_no-pbp_L105-1.7_GWL0p3_RGN",
                         "Output suffix": "RFSL", "Duration": d,
                         "Results folder": "sims_mc/results"})).key
        for d in (6, 9, 12, 18, 24, 36, 48, 72, 96, 120)
    }
    assert len(keys) == 1


def test_blank_output_file_falls_back_to_the_option_name():
    """Every row of TFD_SimsList_LongList_02 lands here - uncached formulas."""
    key = group_key(row(**{"Output file": None, "Output suffix": "FR-4C-no_raise",
                           "FSL": 670.42, "Method": "reservoir routing",
                           "Duration": 18}))
    assert key.key == "FR-4C-no_raise|670.42|reservoir routing"
    assert key.provenance == FROM_FALLBACK


def test_two_options_over_one_base_name_do_not_collapse():
    """The Tinaroo case: every option re-routes the same source run."""
    def key_for(suffix):
        return group_key(row(**{"Output file": "TFD_mc_18h", "Duration": 18,
                                "Output suffix": suffix,
                                "Results folder": r"sims_mc\results"})).key

    assert key_for("FR-4C") != key_for("FR-4D")


def test_provenance_distinguishes_stripped_from_untouched():
    stripped = group_key(row(**{"Output file": "a_18h_b", "Duration": 18}))
    plain = group_key(row(**{"Output file": "a_b", "Duration": 18}))
    assert stripped.provenance == FROM_NAME_STRIPPED
    assert plain.provenance == FROM_NAME


# --- the report -----------------------------------------------------------

def test_report_flags_one_group_per_row():
    frame = pd.DataFrame([
        {"Output file": f"unrelated_{i}", "Duration": 18} for i in range(4)
    ])
    report = derived_group_report(frame)
    assert not report.is_useful
    assert "its own group" in report.degenerate


def test_report_flags_everything_in_one_group():
    frame = pd.DataFrame([{"Output file": "same", "Duration": 18}] * 4)
    report = derived_group_report(frame)
    assert not report.is_useful
    assert "one group" in report.degenerate


def test_report_is_happy_with_a_real_grouping():
    frame = pd.DataFrame([
        {"Output file": f"case{c}_{d}h", "Duration": d}
        for c in (1, 2) for d in (18, 24, 36)
    ])
    report = derived_group_report(frame)
    assert report.is_useful and report.group_count == 2


def test_keys_align_with_the_frame_index():
    frame = pd.DataFrame([{"Output file": "a_18h", "Duration": 18},
                          {"Output file": "b_24h", "Duration": 24}])
    keys = add_group_keys(frame)
    assert list(keys.index) == list(frame.index)
