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


def group_of(selection_path) -> str:
    """The group a selection was saved for, from its file name; '' where it carries none."""
    name = os.path.basename(str(selection_path))
    tail = '_' + EVENTS.SELECTION_FILE
    return name[:-len(tail)] if name.endswith(tail) else ''


def storm_names(selection_path, suffix='') -> set:
    """The storm files a selection would write, for comparing one selection with another."""
    return {storm.filename for storm in plan(selection_path, suffix=suffix).ready}


def clashes(selection_path, selections, suffix='') -> list:
    """The other selections in the project that would write a file of the same name.

    A storm file is named for its realisation, duration and warming level. Two groups that
    share a realisation at the same duration and warming level - the same Monte Carlo run
    re-routed at two full supply levels, say - therefore write the same name, and the second
    silently replaces the first: the storm, and the URBS results named after it. Each is
    compared as it would be written with no suffix of its own, since that is what a run of it
    left behind unless one was given.
    """
    mine = storm_names(selection_path, suffix)
    here = os.path.abspath(str(selection_path))
    return [other for other in selections
            if os.path.abspath(str(other)) != here and mine & storm_names(other)]


def suggested_suffix(selection_path, selections) -> str:
    """A suffix that keeps this selection's files apart from the others', or '' if none is needed.

    The leading word of the group name where that is enough to tell the clashing groups apart
    (``RFSL_Design_GWL1p3`` and ``FSL_Design_GWL1p3`` give ``_RFSL`` and ``_FSL``), and the
    whole group name where it is not.
    """
    others = clashes(selection_path, selections)
    group = group_of(selection_path)
    if not others or not group:
        return ''
    lead = group.split('_')[0]
    if all(group_of(other).split('_')[0] != lead for other in others):
        return f'_{lead}'
    return f'_{group}'


def plan(selection_path, duration=None, gwl=None, suffix='') -> Plan:
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
        if not target.database or not os.path.isfile(
                EVENTS.locate_database(target.database, selection_path)):
            problem = 'the database this event came from is not on disk'
        elif d is None:
            problem = 'no duration in the database name - give one'
        elif g is None:
            problem = 'no warming level in the database name - give one'
        name = None
        if not problem:
            name = (f'{int(target.picked):06d}_{d:g}h_GWL{g:g}'.replace('.', 'p')
                    + f'{suffix or ""}.{d:g}')
        out.storms.append(PlannedStorm(
            loading=f'{target.kind} {target.value:g}', result_type=target.result_type,
            realisation=int(target.picked), database=str(target.database),
            duration=d, gwl=g, filename=name, problem=problem))
    return out


def command(selection, config, model, bryan_python=None, duration=None, gwl=None,
            rain_lag=0.0, suffix='', dry_run=False, sims_config=None, run_urbs=False) -> list:
    """The command line the run uses. Built here so a test can read it."""
    script = Path(__file__).resolve().parents[2] / 'DownstreamStorms.py'
    argv = [str(bryan_python or sys.executable), str(script),
            '--config', str(config), '--selection', str(selection), '--model', str(model)]
    if sims_config:
        # so an event whose embedded bursts were filtered keeps its filtered pattern
        argv += ['--sims-config', str(sims_config)]
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
    elif run_urbs:
        # never with --dry-run: that writes no storm for URBS to run, and the generator refuses
        argv += ['--run-urbs']
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
