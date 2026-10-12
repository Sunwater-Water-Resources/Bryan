"""A saved selection's database is found wherever the reader is running.

The Events page saves each target's database relative to the project folder. The
launcher runs the downstream generator in that folder, so the path resolved there and
nowhere else: the same command typed in Bryan's own folder failed on a file that was
on disk all along. ``locate_database`` is what makes the selection file enough on its
own, and these pin the three places it looks and what it does when none has the file.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

from lib import RepresentativeEvents as events                       # noqa: E402

DATABASE = 'sims_mc\\results\\CLD_mc_72h_E013_GWL1p3_RFSL__mcdf.csv'


def project(tmp_path):
    """A project folder with one database and a selection saved beside it."""
    results = tmp_path / 'project' / 'sims_mc' / 'results'
    results.mkdir(parents=True)
    database = results / 'CLD_mc_72h_E013_GWL1p3_RFSL__mcdf.csv'
    database.write_text('x')
    selection = results / 'RFSL_Design_GWL1p3_representative_events.json'
    selection.write_text('{}')
    return tmp_path / 'project', database, selection


def test_found_from_the_selection_file_when_run_somewhere_else(tmp_path, monkeypatch):
    _, database, selection = project(tmp_path)
    elsewhere = tmp_path / 'elsewhere'
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    assert Path(events.locate_database(DATABASE, selection)) == database


def test_found_in_the_current_folder_as_the_launcher_runs_it(tmp_path, monkeypatch):
    folder, database, selection = project(tmp_path)
    monkeypatch.chdir(folder)
    assert Path(events.locate_database(DATABASE, selection)) == database


def test_a_named_project_folder_is_looked_in_first(tmp_path, monkeypatch):
    # the selection has been copied out of the project, so walking up from it finds nothing
    folder, database, _ = project(tmp_path)
    stray = tmp_path / 'copied' / 'selection.json'
    stray.parent.mkdir()
    stray.write_text('{}')
    monkeypatch.chdir(stray.parent)
    assert events.locate_database(DATABASE, stray) == DATABASE
    assert Path(events.locate_database(DATABASE, stray, folder)) == database


def test_forward_slashes_resolve_too(tmp_path, monkeypatch):
    _, database, selection = project(tmp_path)
    monkeypatch.chdir(tmp_path)
    assert Path(events.locate_database(DATABASE.replace('\\', '/'), selection)) == database


def test_an_absolute_path_is_left_alone(tmp_path):
    _, database, selection = project(tmp_path)
    assert events.locate_database(str(database), selection) == str(database)


def test_a_missing_database_comes_back_as_it_was_saved(tmp_path, monkeypatch):
    _, _, selection = project(tmp_path)
    monkeypatch.chdir(tmp_path)
    missing = 'sims_mc\\results\\not_there__mcdf.csv'
    assert events.locate_database(missing, selection) == missing
    assert not os.path.isfile(events.locate_database(missing, selection))
