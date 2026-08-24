"""Exporting the critical duration analysis from the Results page.

The page computes its own envelope and critical duration - two lines of pandas
- but the exported table also carries the smoothed confidence percentiles, and
that is a degree-5 polyfit in log space with monotonic accumulation
(``util/UtilModule.smoothen_percentiles``). Reimplementing it here would be a
second copy to drift, so the export runs the real thing:

    <bryan python> <bryan root>/util/CriticalDurationAnalysis.py --sim ...

as a subprocess, the way ``core/launcher`` runs Main.py. The exported files are
then identical in format to every study output already produced, and the UI
environment still needs neither scipy nor matplotlib.

One run per result type, because that is how the analysis is shaped: the
critical duration for level is not the critical duration for inflow.
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .bryan import BRYAN_ROOT
from .grouping import strip_duration_token

SCRIPT = BRYAN_ROOT / "util" / "CriticalDurationAnalysis.py"
TIMEOUT_SECONDS = 600

# The output name becomes a filename, so it cannot carry a path.
UNSAFE_NAME = re.compile(r"[\\/:*?\"<>|]")


@dataclass(frozen=True)
class ExportJob:
    """One result type's export: the command, and what it will write."""

    key: str
    command: tuple
    csv: Path
    png: Path

    @property
    def outputs(self) -> tuple:
        return (self.csv, self.png)

    @property
    def existing(self) -> tuple:
        return tuple(path for path in self.outputs if path.is_file())


@dataclass
class ExportPlan:
    jobs: list = field(default_factory=list)
    problems: list = field(default_factory=list)

    @property
    def can_run(self) -> bool:
        return bool(self.jobs) and not self.problems

    @property
    def existing(self) -> list:
        return [path for job in self.jobs for path in job.existing]


@dataclass(frozen=True)
class ExportResult:
    key: str
    returncode: int
    output: str
    written: tuple = ()

    @property
    def ok(self) -> bool:
        return self.returncode == 0


def _number(value) -> str:
    """A duration or AEP as a command-line argument.

    Whole numbers stay whole: ``%g`` renders 2000000 as '2e+06', which parses
    but is unreadable in a log, and leaves an AEP to be matched against an
    integer index as a float.
    """
    number = float(value)
    return f"{number:.0f}" if number.is_integer() else f"{number:g}"


def default_base_name(sources) -> str:
    """The group's name with the duration taken out of it.

    Uses the same rule the grouping does - a duration token is only removed
    when it is that row's own duration - so the exported name matches what the
    Select page calls the group.
    """
    for source in sources:
        if not source.output_name:
            continue
        stripped, _ = strip_duration_token(source.output_name, source.duration)
        if stripped:
            return stripped
    return "critical_durations"


def output_name(base: str, key: str) -> str:
    return f"{base}_{key}_critical"


def default_folder(sources) -> Path | None:
    """Beside the quantile files, which is where the results already live."""
    for source in sources:
        return Path(source.path).parent
    return None


def plan(sources_by_key, keys, *, folder, base_name, python, drop_aeps=(),
         plot=True) -> ExportPlan:
    """What the export would run and write. Nothing is executed here."""
    plan = ExportPlan()

    if not keys:
        plan.problems.append("No result type has been chosen.")
    if not base_name or UNSAFE_NAME.search(base_name):
        plan.problems.append(
            "The output name is empty or holds a path separator - it becomes a "
            "filename, so it cannot contain \\ / : * ? \" < > |.")
    if not python:
        plan.problems.append("Bryan's Python interpreter has not been set - "
                             "set it on the Project page.")
    elif not Path(python).exists():
        plan.problems.append(f"Interpreter not found: {python}")
    if not SCRIPT.is_file():
        plan.problems.append(f"Script not found: {SCRIPT}")
    if folder is None:
        plan.problems.append("No output folder.")

    if plan.problems:
        return plan

    folder = Path(folder)
    for key in keys:
        sources = sources_by_key.get(key, [])
        if not sources:
            plan.problems.append(f"Nothing selected for {key}.")
            continue
        name = output_name(base_name, key)
        command = [str(python), "-u", str(SCRIPT),
                   "--result-type", sources[0].column,
                   "--output-folder", str(folder),
                   "--output-name", name]
        for source in sources:
            duration = source.duration
            if duration is None:
                plan.problems.append(
                    f"{source.label} has no duration, so it cannot go into a "
                    f"critical duration analysis.")
                continue
            command += ["--sim", _number(duration), str(source.path)]
        for aep in drop_aeps:
            command += ["--drop-aep", _number(aep)]
        if not plot:
            command.append("--no-plot")
        plan.jobs.append(ExportJob(key=key, command=tuple(command),
                                   csv=folder / f"{name}.csv",
                                   png=folder / f"{name}_durations.png"))
    return plan


def run(job: ExportJob) -> ExportResult:
    """Run one export. Blocking - call it off the UI thread."""
    environment = dict(os.environ)
    environment["PYTHONUNBUFFERED"] = "1"
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    job.csv.parent.mkdir(parents=True, exist_ok=True)
    try:
        finished = subprocess.run(
            list(job.command), capture_output=True, text=True,
            timeout=TIMEOUT_SECONDS, env=environment,
            stdin=subprocess.DEVNULL,        # never let it wait on a console
            cwd=str(job.csv.parent))
    except subprocess.TimeoutExpired:
        return ExportResult(job.key, 1,
                            f"timed out after {TIMEOUT_SECONDS} s")
    except OSError as exc:
        return ExportResult(job.key, 1, f"could not start the export: {exc}")

    output = (finished.stdout or "") + (finished.stderr or "")
    written = tuple(path for path in job.outputs if path.is_file())
    return ExportResult(job.key, finished.returncode, output, written)
