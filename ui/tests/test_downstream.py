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
    assert argv[1].endswith('DownstreamStorms.py')
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
