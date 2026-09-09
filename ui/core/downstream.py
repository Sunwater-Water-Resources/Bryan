"""Planning and launching downstream storm generation.

The Events page already picks one realisation per design loading and saves the
list as ``<group>_representative_events.json``. This turns that saved list into
a plan - which storm files would be written, for which realisation, at which
duration and warming level - and launches ``DownstreamStorms.py`` to write them.

**It launches rather than imports.** Assembling the depths drives
``lib/StormGenerator.py`` and needs scipy, and the launcher environment
deliberately has neither scipy nor matplotlib. So this is a subprocess, exactly
as ``Main.py`` is - see ``core/launcher.py``. Importing the generator here would
fail ``test_dependency_direction`` and deserve to.

No nicegui in this module: the page is presentation over what is planned here,
and the planning is what is worth testing.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .bryan import representative_events

EVENTS = representative_events()

#: A saved selection does not carry the duration or the warming level, so they
#: come out of the database's own name - the naming the runs already use.
DURATION = re.compile(r'_(\d+(?:p\d+)?)h_')
GWL = re.compile(r'GWL(\d+(?:p\d+)?)')
SELECTION_GLOB = '*representative_events.json'


def _number(text: str) -> float:
    return float(text.replace('p', '.'))


def read_from_name(pattern: re.Pattern, name: str):
    found = pattern.search(os.path.basename(str(name or '')))
    return _number(found.group(1)) if found else None


@dataclass
class PlannedStorm:
    """One storm file the run would write."""
    loading: str
    result_type: str
    realisation: int
    database: str
    duration: float | None
    gwl: float | None
    filename: str | None = None
    problem: str = ''

    @property
    def ready(self) -> bool:
        return not self.problem


@dataclass
class Plan:
    selection: str
    storms: list = field(default_factory=list)

    @property
    def ready(self):
        return [s for s in self.storms if s.ready]

    @property
    def problems(self):
        return [s for s in self.storms if not s.ready]


def find_selections(folder) -> list:
    """Every saved Events page selection under a project folder, newest first."""
    found = sorted(Path(folder).rglob(SELECTION_GLOB),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    return [str(p) for p in found]


def plan(selection_path, duration=None, gwl=None) -> Plan:
    """What a run over this selection would write, and what would stop it.

    Nothing here is fatal on its own: a loading with no event picked is simply
    not a storm, and one whose duration cannot be read is reported so it can be
    given explicitly rather than guessed at.
    """
    targets, _ = EVENTS.read_selection(selection_path)
    out = Plan(selection=str(selection_path))
    for target in targets:
        if target.picked is None:
            continue
        d = duration if duration is not None else read_from_name(DURATION, target.database)
        g = gwl if gwl is not None else read_from_name(GWL, target.database)
        problem = ''
        if not target.database or not os.path.isfile(str(target.database)):
            problem = 'the database this event came from is not on disk'
        elif d is None:
            problem = 'no duration in the database name - give one'
        elif g is None:
            problem = 'no warming level in the database name - give one'
        name = None
        if not problem:
            name = f'{int(target.picked):06d}_{d:g}h_GWL{g:g}'.replace('.', 'p') + f'.{d:g}'
        out.storms.append(PlannedStorm(
            loading=f'{target.kind} {target.value:g}', result_type=target.result_type,
            realisation=int(target.picked), database=str(target.database),
            duration=d, gwl=g, filename=name, problem=problem))
    return out


def command(selection, config, model, bryan_python=None, duration=None, gwl=None,
            rain_lag=0.0, suffix='', dry_run=False) -> list:
    """The command line the run uses. Built here so a test can read it."""
    script = Path(__file__).resolve().parents[2] / 'DownstreamStorms.py'
    argv = [str(bryan_python or sys.executable), str(script),
            '--config', str(config), '--selection', str(selection), '--model', str(model)]
    if duration is not None:
        argv += ['--duration', str(duration)]
    if gwl is not None:
        argv += ['--gwl', str(gwl)]
    if rain_lag:
        argv += ['--rain-lag', str(rain_lag)]
    if suffix:
        argv += ['--suffix', suffix]
    if dry_run:
        argv += ['--dry-run']
    return argv


def launch(argv, cwd, log_path=None):
    """Start the generator. stdin is DEVNULL for the same reason Main.py's is."""
    handle = open(log_path, 'w') if log_path else subprocess.PIPE
    return subprocess.Popen(argv, cwd=str(cwd), stdin=subprocess.DEVNULL,
                            stdout=handle, stderr=subprocess.STDOUT, text=True)


def results_path(selection) -> str:
    """Where the run writes its record of what went into each storm file."""
    return os.path.splitext(str(selection))[0] + '_downstream_storms.csv'


def summarise(results_csv) -> dict:
    """What the finished run produced, for the page to show without pandas games."""
    import csv
    if not os.path.isfile(results_csv):
        return {'written': 0, 'flagged': []}
    with open(results_csv, newline='') as handle:
        rows = list(csv.DictReader(handle))
    flagged = [r for r in rows
               if r.get('duration_dip') or
               r.get('embedded_bursts', '') not in ('', 'No embedded bursts')]
    return {'written': len(rows), 'flagged': flagged, 'rows': rows}
