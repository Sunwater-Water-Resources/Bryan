"""Writing the analysis file: the folder, and the one question worth asking.

The Monte Carlo run ended by asking "the file may be open in Excel" on a
console where Excel was not open, and pressing enter wrote the file happily.
That is not a lock at all - it is ``os.makedirs('')`` raising FileNotFoundError
for an ``Output file`` with no folder part, caught by an ``except IOError``
that assumed only one cause.
"""

from __future__ import annotations

import os

import pandas as pd
import pytest

from lib import FileTools


@pytest.fixture(autouse=True)
def no_prompt(monkeypatch):
    """Nothing may block on the console unless a test says so."""
    def refuse(prompt=''):
        raise AssertionError(f'asked the user a question it should not: {prompt!r}')
    monkeypatch.setattr('builtins.input', refuse)
    monkeypatch.setattr(FileTools, 'INITIAL_DELAY', 0)


@pytest.fixture
def frame():
    return pd.DataFrame({'peak': [1.0, 2.0]}, index=[0, 1])


def test_a_name_with_no_folder_is_written_without_a_question(tmp_path, frame,
                                                             monkeypatch):
    """The reported bug: a bare Output file, blamed on Excel."""
    monkeypatch.chdir(tmp_path)
    FileTools.write_csv(frame, 'CLD_run__mcdf.csv', 'monte carlo')
    assert (tmp_path / 'CLD_run__mcdf.csv').is_file()


def test_a_missing_folder_is_created(tmp_path, frame):
    path = tmp_path / 'results' / 'deeper' / 'out__mcdf.csv'
    FileTools.write_csv(frame, str(path), 'monte carlo')
    assert path.is_file()


def test_an_existing_folder_is_left_alone(tmp_path, frame):
    (tmp_path / 'results').mkdir()
    path = tmp_path / 'results' / 'out__mcdf.csv'
    FileTools.write_csv(frame, str(path), 'monte carlo')
    assert pd.read_csv(path, index_col=0)['peak'].tolist() == [1.0, 2.0]


class _Refuses:
    """A frame whose write fails a set number of times, then succeeds."""

    def __init__(self, error, failures):
        self.error = error
        self.failures = failures
        self.attempts = 0

    def to_csv(self, path):
        self.attempts += 1
        if self.attempts <= self.failures:
            raise self.error
        open(path, 'w').write('written\n')


def test_a_transient_lock_is_retried_rather_than_asked_about(tmp_path):
    """An antivirus scanner or sync agent clears in well under a second."""
    path = tmp_path / 'out__mcdf.csv'
    path.write_text('held')                      # it exists, so a lock is possible
    frame = _Refuses(PermissionError(13, 'in use'), failures=2)

    FileTools.write_csv(frame, str(path), 'monte carlo')
    assert frame.attempts == 3
    assert path.read_text() == 'written\n'


def test_a_held_file_asks_once_and_then_writes(tmp_path, monkeypatch):
    """Excel holds a csv until it is closed, so this one is a real question."""
    path = tmp_path / 'out__mcdf.csv'
    path.write_text('held')
    frame = _Refuses(PermissionError(13, 'in use'), failures=FileTools.ATTEMPTS)

    asked = []
    monkeypatch.setattr('builtins.input', lambda prompt='': asked.append(prompt))
    FileTools.write_csv(frame, str(path), 'monte carlo method analysis')

    assert len(asked) == 1
    assert 'monte carlo method analysis' in asked[0]
    assert str(path) in asked[0]
    assert path.read_text() == 'written\n'


def test_a_windowless_run_still_gets_eoferror_from_a_held_file(tmp_path,
                                                               monkeypatch):
    """What the launcher depends on: stdin=DEVNULL raises instead of hanging."""
    path = tmp_path / 'out__mcdf.csv'
    path.write_text('held')
    frame = _Refuses(PermissionError(13, 'in use'), failures=FileTools.ATTEMPTS)

    def eof(prompt=''):
        raise EOFError('EOF when reading a line')
    monkeypatch.setattr('builtins.input', eof)
    with pytest.raises(EOFError):
        FileTools.write_csv(frame, str(path), 'monte carlo')


def test_a_bad_path_is_raised_not_blamed_on_excel(tmp_path):
    """A directory in the way is not a locked file, and pressing enter is no fix."""
    path = tmp_path / 'out__mcdf.csv'
    path.mkdir()                                  # writing over a folder
    frame = _Refuses(IsADirectoryError(21, 'is a directory'), failures=1)
    with pytest.raises(IsADirectoryError):
        FileTools.write_csv(frame, str(path), 'monte carlo')
    assert frame.attempts == 1                    # no retry, no prompt


def test_a_missing_drive_is_raised_not_blamed_on_excel(tmp_path):
    """The folder cannot be made, so the file was never going to be written."""
    frame = _Refuses(PermissionError(13, 'in use'), failures=0)
    blocker = tmp_path / 'blocker'
    blocker.write_text('not a folder')
    with pytest.raises(OSError):
        FileTools.write_csv(frame, str(blocker / 'deeper' / 'out.csv'), 'monte carlo')
    assert frame.attempts == 0                    # it never got as far as writing
