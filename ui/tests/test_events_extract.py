"""Running the representative-event extraction from the Events page.

The command used to be ``python util/RepresentativeEvents.py ...`` - right only
when typed in Bryan's folder with an interpreter that has scipy and matplotlib -
and had to be run by hand. It is now built from absolute paths and Bryan's own
interpreter, and the page runs it and shows the plots it wrote.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

from core import events


def project_at(config_path):
    return SimpleNamespace(config=SimpleNamespace(config_path=Path(config_path)))


def test_the_command_is_bryan_s_interpreter_and_absolute_paths(tmp_path):
    argv = events.extract_argv(project_at(tmp_path / "sims.json"),
                               tmp_path / "g_representative_events.json", "C:/env/python.exe")
    assert argv[0] == "C:/env/python.exe"
    assert Path(argv[2]).is_absolute() and argv[2].endswith("RepresentativeEvents.py")
    assert argv[argv.index("--config") + 1] == str(tmp_path / "sims.json")


def test_a_path_with_spaces_is_quoted_in_the_copied_command(tmp_path):
    folder = tmp_path / "Design estimate"
    command = events.extract_command(project_at(folder / "sims.json"),
                                     folder / "x_representative_events.json", "python")
    assert f'"{folder / "sims.json"}"' in command


def test_the_plots_are_the_ones_the_run_says_it_wrote(tmp_path):
    written = tmp_path / "event_1_in_100_sim00042.png"
    written.write_bytes(b"png")
    missing = tmp_path / "never_made.png"
    script = (f"print('  wrote {written.as_posix()}'); print('  wrote {missing.as_posix()}');"
              f"print('wrote {(tmp_path / 'events.xlsx').as_posix()}')")
    result = events.run_extract([sys.executable, "-c", script], tmp_path)
    assert result.ok
    assert result.plots == [Path(written.as_posix())]


def test_a_failed_run_keeps_its_output(tmp_path):
    result = events.run_extract([sys.executable, "-c", "import sys; print('boom');"
                                 " sys.exit(3)"], tmp_path)
    assert not result.ok and result.returncode == 3 and "boom" in result.output


def test_an_interpreter_that_is_not_there_is_reported_not_raised(tmp_path):
    result = events.run_extract([str(tmp_path / "no_python.exe")], tmp_path)
    assert not result.ok and "could not start" in result.output
