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
    method: str = 'monte carlo'

    @property
    def ready(self) -> bool:
        return not self.problem

    @property
    def event_label(self) -> str:
        """'sim 4920' for a realisation, 'event 42' for an ensemble row."""
        return (f'event {self.realisation}' if self.method == EVENTS.ENSEMBLE
                else f'sim {self.realisation}')


@dataclass
class Plan:
    selection: str
    storms: list = field(default_factory=list)
    method: str = 'monte carlo'

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


def method_of(selection_path) -> str:
    return EVENTS.read_selection_file(selection_path).method


def selection_label(selection_path) -> str:
    """How the page lists a selection: its group, marked where it is an ensemble one."""
    group = group_of(selection_path) or os.path.basename(str(selection_path))
    return f'{group}  (ensemble)' if method_of(selection_path) == EVENTS.ENSEMBLE else group


def extra_dams(model_config) -> list:
    """The model config's additional dams, as (location, ensemble_start) - read from the
    JSON, so a config with a problem still gives the page something to say."""
    try:
        with open(model_config, encoding='utf-8') as handle:
            data = json.load(handle)
    except (OSError, ValueError, TypeError):
        return []
    return [(str(entry.get('location', '')), entry.get('ensemble_start', 'fsl'))
            for entry in data.get('additional_dams') or [] if isinstance(entry, dict)]


def method_line(selection_path, model_config) -> str:
    """One sentence on how each dam starts, for the selection's kind of event."""
    dams = extra_dams(model_config) if model_config else []
    if method_of(selection_path) == EVENTS.ENSEMBLE:
        said = "Ensemble events: the main dam starts from the event's ADV"
        for location, start in dams:
            try:
                where = ('full supply' if str(start).strip().lower() == 'fsl'
                         else f'{float(start):,.0f} ML')
            except ValueError:
                where = f'"{start}" (not a volume - the run will refuse it)'
            said += f'; {location} at {where}'
        return said + (". Change that with ensemble_start in the model config's "
                       'additional_dams.' if dams else '.')
    said = "Monte Carlo events: the main dam starts from the event's ADV"
    if dams:
        said += ('; ' + ', '.join(location for location, _ in dams)
                 + " from the event's lake_z through its own lake config")
    return said + '.'


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

    The words of the group name - split at underscores - that no clashing group has, in order:
    ``RFSL_Design_GWL1p3`` against ``FSL_Design_GWL1p3`` gives ``_RFSL``, and so does
    ``CLD_enb_PMF_..._GWL1p3_RFSL`` against ``..._FSL``, where the word that differs is last.
    Only the leading word used to be tried, and a shared one fell back to the whole group name.
    That is still the answer where no word of the name is its own.
    """
    others = clashes(selection_path, selections)
    group = group_of(selection_path)
    if not others or not group:
        return ''
    theirs = {word for other in others for word in group_of(other).split('_')}
    own = [word for word in group.split('_') if word and word not in theirs]
    return '_' + '_'.join(own) if own else f'_{group}'


def sims_config_for(selection_path, project_folder, fallback=None) -> str | None:
    """The sims config to rebuild the selection's patterns from: the one the selection
    records, where it is on disk, else ``fallback`` - the open project's. A selection saved
    before the file recorded one only works with the sims config it was saved under open."""
    recorded = EVENTS.read_selection_file(selection_path).sims_config
    if recorded:
        found = EVENTS.locate_database(recorded, selection_path, project_folder)
        if os.path.isfile(found):
            return found
    return str(fallback) if fallback else None


def storm_filename(picked, duration, gwl, suffix='', method='monte carlo') -> str:
    """The storm file's name. ``DownstreamStorms.storm_name`` writes the same, and
    ``tests/test_downstream_naming.py`` holds the two together: an ensemble event, a row
    of its own database, is ``e`` and five digits so it can never take a realisation's."""
    event = f'e{int(picked):05d}' if method == EVENTS.ENSEMBLE else f'{int(picked):06d}'
    return (f'{event}_{duration:g}h_GWL{gwl:g}'.replace('.', 'p')
            # 'p' in the extension as well: '.4.5' made the run name end '.4'
            + f'{suffix or ""}.' + f'{duration:g}'.replace('.', 'p'))


def plan(selection_path, duration=None, gwl=None, suffix='') -> Plan:
    """What a run over this selection would write, and what would stop it.

    Nothing here is fatal on its own: a loading with no event picked is simply
    not a storm, and one whose duration cannot be read is reported so it can be
    given explicitly rather than guessed at.
    """
    selection = EVENTS.read_selection_file(selection_path)
    out = Plan(selection=str(selection_path), method=selection.method)
    for target in selection.targets:
        if target.picked is None:
            continue
        # an ensemble database holds every duration; the selection records the event's own
        d = (duration if duration is not None
             else target.duration if target.duration is not None
             else read_from_name(DURATION, target.database))
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
            name = storm_filename(target.picked, d, g, suffix, selection.method)
        out.storms.append(PlannedStorm(
            loading=target.title, result_type=target.result_type,
            realisation=int(target.picked), database=str(target.database),
            duration=d, gwl=g, filename=name, problem=problem, method=selection.method))
    return out


def command(selection, config, model, bryan_python=None, duration=None, gwl=None,
            rain_lag=0.0, suffix='', dry_run=False, sims_config=None, run_urbs=False) -> list:
    """The command line the run uses. Built here so a test can read it."""
    script = Path(__file__).resolve().parents[2] / 'DownstreamStorms.py'
    # -u: the output goes to a file, and block-buffered it reached the log in 8 kB lumps,
    # so the page's status line and anyone watching the log saw nothing for minutes
    argv = [str(bryan_python or sys.executable), '-u', str(script),
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
    """Start the generator. stdin is DEVNULL for the same reason Main.py's is.

    The child inherits its own handle on the log, so the page's copy is closed as soon
    as the child has it; holding it open kept the file locked after the run ended."""
    if not log_path:
        return subprocess.Popen(argv, cwd=str(cwd), stdin=subprocess.DEVNULL,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(log_path, 'w') as handle:
        return subprocess.Popen(argv, cwd=str(cwd), stdin=subprocess.DEVNULL,
                                stdout=handle, stderr=subprocess.STDOUT, text=True)


#: The longest log line the status line shows.
STATUS_WIDTH = 140


def last_line(log_text: str) -> str:
    """The last line of the generator's log worth showing while it runs: not blank, and not
    one of the indented antecedent or pre-burst table lines that follow each storm."""
    lines = [line.rstrip() for line in log_text.splitlines() if line.strip()]
    for line in reversed(lines):
        if line.lstrip().startswith(('wrote ', 'would write ', 'batch file')) or not line[:1].isspace():
            text = line.strip()
            return text if len(text) <= STATUS_WIDTH else text[:STATUS_WIDTH - 3] + '...'
    return ''


def outcome(returncode: int, log_text: str, log_name: str) -> tuple:
    """(succeeded, what to say) once the generator has stopped."""
    if returncode == 0:
        # DownstreamStorms.main's closing lines: "7 storm files; what went into them: <csv>"
        # and, where any were flagged, "3 need a look - see duration_dip and embedded_bursts"
        said = []
        for line in log_text.splitlines():
            line = line.strip()
            if re.match(r'\d+ storm files', line):
                said.append(line.split(';')[0] + '.')
            elif 'need a look' in line or re.match(r'\d+ skipped', line):
                said.append(line + '.')
        return True, " ".join(["Finished."] + said)
    return False, (f"Stopped with an error (exit code {returncode}) - see {log_name}. "
                   f"{last_line(log_text)}").strip()


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
