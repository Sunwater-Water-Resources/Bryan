"""Which storm durations a table's or figure's critical duration is taken over.

The study keeps a default; a table or figure follows it, takes every duration, or
chooses its own (core/study.py). The fixture's level is critical at 72 h up to
1 in 100 and at 36 h from 1 in 1,000, so leaving either out changes both the
critical duration and the flows read at it - which is what is checked.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from core import figures, results
from core import reporttables as rt
from core import study as studies
from report_fixtures import GROUP, INFLOW, LEVEL, build_study


@pytest.fixture
def study(tmp_path):
    studies.forget_runs()
    rt.forget_cached()
    return build_study(tmp_path / "study")


def design_spec(**overrides):
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Near-term", source={"run": "E099 RFSL", "group": GROUP},
                aeps=[10, 1000], pmp_aep=None)
    spec.update(overrides)
    return spec


def by_aep(table):
    return {row.cells[0]: row.cells for row in table.rows}


# -- the filter ------------------------------------------------------------------

def _source(label, duration):
    return results.CurveSource(label=label, duration=duration, path=Path(f"{label}.csv"),
                               column="level", key="level", output_name=label, row_index=0)


def test_no_limit_keeps_every_duration():
    sources = [_source("36h", 36.0), _source("72h", 72.0), _source("odd", None)]
    kept, problems = results.limit_durations(sources, [])
    assert kept == sources and problems == []


def test_a_limit_matches_the_duration_not_the_label():
    sources = [_source("24h (A)", 24.0), _source("24h (B)", 24.0), _source("72h", 72.0)]
    kept, _ = results.limit_durations(sources, [24])
    assert [source.label for source in kept] == ["24h (A)", "24h (B)"]


def test_a_duration_that_cannot_be_read_is_left_out_under_a_limit_and_said():
    kept, problems = results.limit_durations([_source("odd", None)], [24])
    assert kept == [] and "could not be read" in problems[0]


# -- the rule: default, all, or chosen -------------------------------------------

@pytest.mark.parametrize("value, expected", [
    (None, [36.0]), ([], [36.0]),           # follows the default; [] is the old "all"
    ("all", []), ("ALL", []),               # every duration, whatever the default
    ([72, "72", 36], [36.0, 72.0]),         # its own, tidied
])
def test_a_table_follows_the_default_unless_it_says_otherwise(study, value, expected):
    studies.set_default_durations(study, [36])
    spec = {"durations": value} if value is not None else {}
    assert studies.durations_for(study, spec) == expected


def test_without_a_default_every_duration_counts(study):
    assert studies.default_durations(study) == []
    assert studies.durations_for(study, {}) == []


def test_the_default_is_kept_in_the_study_file_and_cleared_by_none(study):
    studies.set_default_durations(study, ["72", 36, "nonsense"])
    study.save()
    assert studies.default_durations(studies.load_study(study.path)) == [36.0, 72.0]
    studies.set_default_durations(study, [])
    study.save()
    assert "durations" not in studies.load_study(study.path).extra


def test_the_choices_offered_are_the_durations_the_study_ran(study):
    assert studies.durations_in(study) == [36.0, 72.0]
    assert studies.durations_in(study, [{"run": "E099 RFSL", "group": GROUP}]) == [36.0, 72.0]
    assert studies.durations_in(study, [{"run": "E099 RFSL", "group": "nope"}]) == []


def test_durations_read_as_text():
    assert studies.durations_text([]) == "all durations"
    assert studies.durations_text([96, 6, 12.5]) == "6, 12.5, 96 h"


# -- tables ----------------------------------------------------------------------

def test_leaving_a_duration_out_moves_the_critical_duration_and_its_flows(study):
    every = by_aep(rt.build(study, design_spec()))
    assert every["1,000"] == ["1,000", "5,000", "4,400", "217.30", "36"]
    only_72 = by_aep(rt.build(study, design_spec(durations=[72])))
    # 36 h was critical at 1 in 1,000; without it the 72 h level and inflow are read
    assert only_72["1,000"] == ["1,000", f"{INFLOW[72][4]:,}", "4,360",
                                f"{LEVEL[72][4]:.2f}", "72"]
    assert only_72["10"] == every["10"]                  # 72 h already governed there


def test_the_study_default_reaches_a_table_that_follows_it(study):
    studies.set_default_durations(study, [36])
    table = by_aep(rt.build(study, design_spec()))
    assert table["10"][-1] == "36" and table["10"][3] == f"{LEVEL[36][2]:.2f}"
    # ... and not one that asked for every duration
    assert by_aep(rt.build(study, design_spec(durations="all")))["10"][-1] == "72"


def test_a_filtered_read_and_an_unfiltered_one_are_cached_apart(study):
    whole = rt.group_curves(study, "E099 RFSL", GROUP)
    part = rt.group_curves(study, "E099 RFSL", GROUP, [72])
    assert list(whole.comparisons["level"].frame.columns) == ["36h", "72h"]
    assert list(part.comparisons["level"].frame.columns) == ["72h"]
    assert rt.group_curves(study, "E099 RFSL", GROUP) is whole


def test_a_limit_no_run_meets_is_reported_not_drawn_empty(study):
    table = rt.build(study, design_spec(durations=[48]))
    assert any("durations chosen (48 h)" in problem for problem in table.problems)


def test_every_kind_that_reads_a_design_curve_takes_the_limit(study):
    spec = rt.new_spec(rt.PEAK_AT_AEP)
    spec["aep"] = 1000
    spec["sections"] = [{"heading": "", "rows": [
        {"label": "Near Term", "run": "E099 RFSL", "group": GROUP}]}]
    assert rt.build(study, spec).rows[0].cells[-1] == "36"
    spec["durations"] = [72]
    assert rt.build(study, spec).rows[0].cells[-1] == "72"


def test_an_old_frequent_table_s_blank_limit_now_follows_the_default(study):
    spec = rt.new_spec(rt.FREQUENT)
    spec.update(columns=["GWL 1.3"], durations=[],         # as such tables were saved
                sections=[{"heading": "", "groups": [{"run": "E099 RFSL", "group": GROUP}]}])
    assert rt.build(study, spec).rows[0].cells[-1] == "72"     # 1 in 2: 72 h governs
    studies.set_default_durations(study, [36])
    assert rt.build(study, spec).rows[0].cells[-1] == "36"


def test_a_representative_event_outside_the_limit_is_flagged(study):
    from core import events

    folder = study.folder / "runs" / "E099" / "sims_mc" / "results"
    events.save_targets(events.selection_path(folder, GROUP), [events.Target.from_dict(
        {"kind": "aep", "value": 100, "result_type": "level", "picked": 4980,
         "output_file": "CLD_mc_72h_GWL1p3_RFSL"})])
    spec = rt.new_spec(rt.REPRESENTATIVE)
    spec["sections"] = [{"heading": "", "run": "E099 RFSL", "group": GROUP, "pmf": None}]
    assert not any("outside" in p for p in rt.build(study, spec).problems)
    spec["durations"] = [36]
    table = rt.build(study, spec)
    assert table.rows[0].cells[-1] == "72h"                # still the event's own run
    assert any("72h run, which is outside" in p for p in table.problems)


# -- figures ---------------------------------------------------------------------

def level_figure(**overrides):
    spec = figures.new_spec(filename="levels", type="level", curves=[
        {"kind": figures.GROUP, "run": "E099 RFSL", "group": GROUP, "label": "GWL 1.3"}])
    spec.update(overrides)
    return spec


def test_a_figure_s_group_curve_is_the_envelope_over_its_durations(study, tmp_path):
    curve = dict(figures.build(study, level_figure(durations=[72])).series[0].points)
    assert curve[1_000.0] == LEVEL[72][4]
    studies.set_default_durations(study, [36])
    data = figures.build(study, level_figure())
    assert dict(data.series[0].points)[5.0] == LEVEL[36][1]
    assert data.to_job(tmp_path / "x.png")["durations"] == [36.0]
    assert figures.build(study, level_figure(durations="all")) \
        .to_job(tmp_path / "x.png")["durations"] == "all"


# -- the study default's select: applied once its list closes --------------------

def test_picks_made_with_the_list_open_are_applied_once_when_it_closes():
    pytest.importorskip("nicegui")
    from durationpicker import Deferred

    applied = []
    pending = Deferred([24], applied.append)
    pending.opened()
    for value in ([24, 36], [24, 36, 72], [36, 72]):
        pending.changed(value)
    assert applied == []                           # nothing rebuilt while picking
    pending.closed()
    assert applied == [[36.0, 72.0]]
    pending.opened()
    pending.closed()                               # opened and shut, nothing changed
    assert applied == [[36.0, 72.0]]


def test_a_change_with_the_list_shut_applies_at_once():
    pytest.importorskip("nicegui")
    from durationpicker import Deferred

    applied = []
    pending = Deferred([24, 36], applied.append)
    pending.changed([36])                          # a chip's cross
    pending.changed([])                            # the clear button
    assert applied == [[36.0], []]


# -- the decision, noted ---------------------------------------------------------

def test_a_change_to_the_default_is_dated_and_a_non_change_is_not(study):
    import datetime as dt

    assert studies.durations_changed(study) == ""
    studies.set_default_durations(study, [36], today=dt.date(2026, 11, 6))
    studies.set_default_durations(study, [36.0], today=dt.date(2026, 12, 1))   # same
    assert studies.durations_changed(study) == "2026-11-06"
    studies.set_default_durations(study, [], today=dt.date(2026, 12, 2))       # cleared
    assert studies.durations_changed(study) == "2026-12-02"


def test_the_reason_is_kept_in_the_study_file(study):
    studies.set_durations_note(study, "  120 h left out: outside the ARR temporal patterns  ")
    study.save()
    assert studies.durations_note(studies.load_study(study.path)) \
        == "120 h left out: outside the ARR temporal patterns"
    studies.set_durations_note(study, "")
    assert "durations_note" not in study.extra


def test_what_is_left_out_is_read_against_what_was_run(study):
    assert studies.left_out(study) == []                  # no default, nothing left out
    studies.set_default_durations(study, [36])
    assert studies.left_out(study) == [72.0]
    assert studies.left_out(study, [36, 72]) == []


def test_the_reason_in_effect_is_the_study_s_or_an_override_s_own(study):
    studies.set_durations_note(study, "study reason")
    assert studies.note_for(study, {}) == ""              # the default limits nothing
    studies.set_default_durations(study, [36])
    assert studies.note_for(study, {}) == "study reason"
    assert studies.note_for(study, {"durations": "all", "durations_note": "own"}) == "own"
    assert studies.note_for(study, {"durations": [72]}) == ""


def test_the_folded_line_says_what_is_left_out_and_an_override_s_reason():
    pytest.importorskip("nicegui")
    import durationpicker
    from report_fixtures import build_study as build

    import tempfile
    with tempfile.TemporaryDirectory() as folder:
        studies.forget_runs()
        study = build(Path(folder) / "study")
        studies.set_default_durations(study, [36])
        assert durationpicker.describe(study, {}) == "leaves out 72 h (study default)"
        long = "the 36 h storms are left out because " + "x" * 80
        assert durationpicker.describe(study, {"durations": [72], "durations_note": long}) \
            == "leaves out 36 h (own choice: the 36 h storms are left out because...)"
        assert durationpicker.describe(study, {"durations": "all"}) == "all durations"


def test_an_exported_figure_records_the_reason(study, tmp_path):
    studies.set_default_durations(study, [36])
    studies.set_durations_note(study, "120 h dropped")
    job = figures.build(study, level_figure()).to_job(tmp_path / "x.png")
    assert job["durations"] == [36.0] and job["durations_note"] == "120 h dropped"
