"""Planning downstream storm generation, and keeping it a subprocess.

The Events page's saved selection is the input to the downstream generator, so
what matters here is that a plan can be read off it without the generator being
importable at all - assembling the depths needs scipy, and the launcher
environment has none.

A saved target carries the realisation and the database it came from, but not
the duration or the warming level, so those are read out of the database's own
name. That is a convention, not a guarantee, so a name that does not carry them
has to be reported rather than guessed at: picking the wrong duration would
write a plausible storm of the wrong length.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

UI_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(UI_ROOT))

from core import downstream                                          # noqa: E402


def selection(tmp_path, targets, name='GWL1p3_representative_events.json'):
    path = tmp_path / name
    path.write_text(json.dumps({'targets': targets, 'settings': {}}))
    return path


def target(database, picked=42, **kwargs):
    base = dict(kind='aep', value=1000.0, result_type='level', rain_aep=None,
                source='', output_file='CLD_mc_24h', database=str(database),
                count=10, picked=picked, comment='')
    base.update(kwargs)
    return base


def database(tmp_path, name='CLD_mc_24h_E009_no-pbp_GWL1p3_RGN__mcdf.parquet'):
    path = tmp_path / name
    path.write_bytes(b'')                      # only its name and existence are read
    return path


def test_duration_and_warming_level_come_off_the_database_name(tmp_path):
    db = database(tmp_path)
    plan = downstream.plan(selection(tmp_path, [target(db)]))
    storm, = plan.storms
    assert (storm.duration, storm.gwl) == (24.0, 1.3)
    assert storm.filename == '000042_24h_GWL1p3.24'


def test_a_fractional_duration_survives_the_p_notation(tmp_path):
    # The regional model runs 6-120 h, all whole numbers, but the naming
    # convention carries 4p5 and the parser has to read it rather than stop at
    # the digits - a duration silently read as 4 would be a valid storm of the
    # wrong length.
    db = database(tmp_path, 'CLD_mc_4p5h_E009_GWL0p3_RGN__mcdf.parquet')
    storm, = downstream.plan(selection(tmp_path, [target(db)])).storms
    assert (storm.duration, storm.gwl) == (4.5, 0.3)
    assert storm.filename.startswith('000042_4p5h_GWL0p3')


def test_a_name_carrying_neither_is_reported_not_guessed(tmp_path):
    db = database(tmp_path, 'results__mcdf.parquet')
    storm, = downstream.plan(selection(tmp_path, [target(db)])).storms
    assert not storm.ready
    assert 'duration' in storm.problem


def test_an_override_beats_the_name(tmp_path):
    db = database(tmp_path, 'results__mcdf.parquet')
    plan = downstream.plan(selection(tmp_path, [target(db)]), duration=72, gwl=2.7)
    storm, = plan.storms
    assert storm.ready and (storm.duration, storm.gwl) == (72, 2.7)


def test_a_loading_with_no_event_picked_is_not_a_storm(tmp_path):
    db = database(tmp_path)
    plan = downstream.plan(selection(tmp_path, [target(db), target(db, picked=None)]))
    assert len(plan.storms) == 1


def test_a_database_that_is_not_on_disk_is_a_problem(tmp_path):
    plan = downstream.plan(selection(tmp_path, [target(tmp_path / 'gone_24h_GWL1p3.parquet')]))
    storm, = plan.storms
    assert not storm.ready and 'not on disk' in storm.problem


def test_selections_are_found_newest_first(tmp_path):
    import os, time
    older = selection(tmp_path, [], 'GWL0p3_representative_events.json')
    time.sleep(0.01)
    newer = selection(tmp_path, [], 'GWL1p3_representative_events.json')
    os.utime(older, (1, 1))
    found = downstream.find_selections(tmp_path)
    assert [Path(f).name for f in found][0] == newer.name


def test_the_command_names_the_generator_and_the_interpreter(tmp_path):
    argv = downstream.command('sel.json', 'cfg.json', 'model.json',
                              bryan_python='/bryan/.venv/bin/python', duration=24)
    assert argv[0] == '/bryan/.venv/bin/python'
    assert argv[2].endswith('DownstreamStorms.py')
    assert '--duration' in argv and '24' in argv


def test_the_generator_is_launched_never_imported():
    # The whole reason core/downstream.py exists as a launcher: importing the
    # generator would pull scipy into an environment that deliberately has none.
    source = (UI_ROOT / 'core' / 'downstream.py').read_text()
    assert 'DownstreamStorms' in source
    assert 'from lib.DownstreamStorms' not in source
    assert 'import lib.DownstreamStorms' not in source


def test_summarise_flags_what_needs_a_look(tmp_path):
    csv = tmp_path / 'r.csv'
    csv.write_text('filename,duration_dip,embedded_bursts\n'
                   'a.24,,No embedded bursts\n'
                   'b.24,GRV;PRO,No embedded bursts\n'
                   'c.24,,Embedded burst in the 12h sub-burst\n')
    out = downstream.summarise(str(csv))
    assert out['written'] == 3
    assert [r['filename'] for r in out['flagged']] == ['b.24', 'c.24']


# -- two groups that share a realisation ---------------------------------------------

def two_groups(tmp_path):
    """RFSL and FSL selections holding the same realisation of the same database."""
    db = database(tmp_path)
    rfsl = selection(tmp_path, [target(db)], 'RFSL_Design_GWL1p3_representative_events.json')
    fsl = selection(tmp_path, [target(db)], 'FSL_Design_GWL1p3_representative_events.json')
    return rfsl, fsl


def test_a_suffix_goes_before_the_extension(tmp_path):
    plan = downstream.plan(selection(tmp_path, [target(database(tmp_path))]), suffix='_RFSL')
    assert plan.storms[0].filename == '000042_24h_GWL1p3_RFSL.24'


def test_two_groups_sharing_a_realisation_clash(tmp_path):
    rfsl, fsl = two_groups(tmp_path)
    assert downstream.clashes(rfsl, [rfsl, fsl]) == [fsl]
    assert downstream.clashes(rfsl, [rfsl, fsl], suffix='_RFSL') == []


def test_the_suggested_suffix_is_the_word_that_tells_the_groups_apart(tmp_path):
    rfsl, fsl = two_groups(tmp_path)
    assert downstream.suggested_suffix(rfsl, [rfsl, fsl]) == '_RFSL'
    assert downstream.suggested_suffix(fsl, [rfsl, fsl]) == '_FSL'


def test_no_suffix_is_suggested_where_nothing_clashes(tmp_path):
    only = selection(tmp_path, [target(database(tmp_path))])
    assert downstream.suggested_suffix(only, [only]) == ''


def test_groups_with_the_same_leading_word_get_the_whole_name(tmp_path):
    db = database(tmp_path)
    one = selection(tmp_path, [target(db)], 'RFSL_Design_GWL1p3_representative_events.json')
    two = selection(tmp_path, [target(db)], 'RFSL_Sens_GWL1p3_representative_events.json')
    assert downstream.suggested_suffix(one, [one, two]) == '_RFSL_Design_GWL1p3'


# -- running URBS, and the dam's outflow ----------------------------------------------

def test_urbs_is_run_only_when_asked_and_never_on_a_check(tmp_path):
    base = ('sel.json', 'cfg.json', 'model.json')
    assert '--run-urbs' not in downstream.command(*base)
    assert '--run-urbs' in downstream.command(*base, run_urbs=True)
    checked = downstream.command(*base, run_urbs=True, dry_run=True)
    assert '--dry-run' in checked and '--run-urbs' not in checked


def test_the_suffix_and_the_sims_config_reach_the_generator():
    argv = downstream.command('sel.json', 'cfg.json', 'model.json', suffix='_RFSL',
                              sims_config='sims_config.json')
    assert argv[argv.index('--suffix') + 1] == '_RFSL'
    assert argv[argv.index('--sims-config') + 1] == 'sims_config.json'


def test_the_generator_runs_unbuffered_so_its_log_can_be_followed():
    argv = downstream.command('sel.json', 'cfg.json', 'model.json', bryan_python='/py')
    assert argv[:2] == ['/py', '-u'] and argv[2].endswith('DownstreamStorms.py')


def test_the_status_line_is_the_last_line_worth_showing():
    log = ("7 loadings in X.json, 7 with an event picked\n"
           "\naep 100 (level) -> realisation 5044, 72 h, GWL 1.3\n"
           "  wrote 005044_72h_GWL1p3_RFSL.72\n"
           "  antecedent storage for 005044_72h_GWL1p3_RFSL:\n"
           "  initial_lake_level: 93,501 ML from the event ADV\n\n")
    assert downstream.last_line(log) == "wrote 005044_72h_GWL1p3_RFSL.72"
    assert downstream.last_line("") == ""
    assert len(downstream.last_line("x" * 500)) == downstream.STATUS_WIDTH


def test_the_outcome_reads_the_generators_closing_lines():
    log = ("...\n\n7 storm files; what went into them: C:\a\b_downstream_storms.csv\n"
           "3 need a look - see duration_dip and embedded_bursts\n")
    assert downstream.outcome(0, log, "x.log") == (
        True, "Finished. 7 storm files. 3 need a look - see duration_dip and embedded_bursts.")
    ok, said = downstream.outcome(1, "Traceback...\nSystemExit: no such file\n", "x.log")
    assert not ok and "exit code 1" in said and "x.log" in said and "no such file" in said


def test_the_sims_config_comes_from_the_selection_where_it_records_one(tmp_path):
    from core.bryan import representative_events
    events = representative_events()
    (tmp_path / "runs").mkdir()
    recorded = tmp_path / "runs" / "CLD_RFSL_mc_sims_01.json"
    recorded.write_text("{}")
    selection = tmp_path / "RFSL_Design_GWL1p3_representative_events.json"
    selection.write_text(json.dumps(events.selection_payload(
        [], sims_config=r"runs\CLD_RFSL_mc_sims_01.json")))
    assert downstream.sims_config_for(selection, tmp_path, "open.json") == str(recorded)
    # one saved before the file recorded it, or recording one no longer on disk: the open one
    selection.write_text(json.dumps({"targets": []}))
    assert downstream.sims_config_for(selection, tmp_path, "open.json") == "open.json"
    selection.write_text(json.dumps(events.selection_payload([], sims_config="gone.json")))
    assert downstream.sims_config_for(selection, tmp_path, "open.json") == "open.json"
