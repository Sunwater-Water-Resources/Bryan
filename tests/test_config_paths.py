"""Config paths written on Windows have to resolve on Linux too.

Bryan's configs are authored on Windows and record their relative paths with
backslashes. Windows treats ``\\`` and ``/`` alike, so the choice is invisible
there; on Linux a backslash is an ordinary filename character, and
``os.path.join('storm_data', 'ifd\\_ifd_files.json')`` names one file that has
never existed instead of reaching a file inside a directory that has. Every
config-relative path in Bryan goes through ``ConfigPaths.resolve`` so that the
same config runs on both.

The tests build real files under ``tmp_path`` rather than asserting on strings,
because the whole point is whether the path opens - a comparison against an
expected string would pass on the platform that was never broken.

``resolve`` deliberately does **not** touch capitalisation. Windows resolves
``ifd`` to a directory named ``IFD``; Linux does not, and a resolver that
picked a similarly-named directory would turn a config typo into a silent
substitution of one input file for another. Those are fixed in the config,
where they are also wrong on Windows - just forgiven there.

Run with Bryan's own interpreter:

    python -m pytest tests -q
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

from lib.ConfigPaths import resolve, resolve_all                    # noqa: E402


@pytest.fixture
def study(tmp_path):
    """A miniature study tree, laid out the way a real one is."""
    (tmp_path / 'storm_data' / 'IFD' / 'Regionalised').mkdir(parents=True)
    (tmp_path / 'storm_data' / 'IFD' / '_ifd_files.json').write_text('{}')
    (tmp_path / 'storm_data' / 'IFD' / 'Regionalised' / 'average_ifd_24.csv').write_text('')
    return tmp_path


def test_a_windows_path_opens_on_this_platform(study):
    path = resolve(str(study / 'storm_data'), 'IFD\\_ifd_files.json')
    assert os.path.isfile(path), path


def test_nested_windows_path_opens_too(study):
    path = resolve(str(study / 'storm_data' / 'IFD'),
                   'Regionalised\\average_ifd_24.csv')
    assert os.path.isfile(path), path


def test_a_posix_path_is_unaffected(study):
    path = resolve(str(study / 'storm_data'), 'IFD/_ifd_files.json')
    assert os.path.isfile(path), path


def test_both_separators_give_the_same_path(study):
    folder = str(study / 'storm_data')
    assert resolve(folder, 'IFD\\_ifd_files.json') == resolve(folder, 'IFD/_ifd_files.json')


def test_an_absolute_path_is_left_alone(study):
    # model_exe is an absolute Windows path, and must survive intact for the
    # platform that can run it rather than being joined onto a study folder.
    absolute = str(study / 'storm_data' / 'IFD' / '_ifd_files.json')
    assert resolve(str(study), absolute) == os.path.normpath(absolute)


def test_capitalisation_is_never_guessed(study):
    # 'ifd' is not 'IFD'. Windows forgives it; resolve does not invent a match,
    # because the alternative is quietly reading a different file.
    path = resolve(str(study / 'storm_data'), 'ifd\\_ifd_files.json')
    assert not os.path.exists(path)


def test_resolve_all_rewrites_a_filepaths_mapping(study):
    filepaths = {'rare_ifds': 'IFD\\_ifd_files.json',
                 'areal': 'IFD\\Regionalised\\average_ifd_24.csv'}
    resolve_all(str(study / 'storm_data'), filepaths)
    for key, path in filepaths.items():
        assert os.path.isfile(path), f'{key} -> {path}'


@pytest.mark.skipif(os.sep != '/', reason='posix-only behaviour')
def test_on_posix_the_backslash_becomes_a_separator(study):
    assert '\\' not in resolve(str(study / 'storm_data'), 'IFD\\_ifd_files.json')
