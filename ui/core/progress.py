"""Working out how far a running chunk has got.

Four sources, each answering a different question, because no single one is
enough:

- **the run log CSV** is authoritative for *finished* simulations, but
  lib/RunLog.py only appends after each one completes, so it reads zero all the
  way through a first hour-long simulation;
- **the console log** names the simulation currently running - Main.py prints a
  banner per simulation before starting it;
- **the per-simulation log file** growing tells you it is alive and doing work,
  which is what distinguishes 'grinding through realisation 4000' from 'wedged';
- **the process return code** says how the whole chunk ended (Main.py:128 exits
  1 if anything failed).

Everything here is a pure function of file contents, so it is fully testable
without launching anything.

One trap worth naming: run-log rows must be matched to simulations
**positionally**, not by name. ``RunLog.start_entry`` records ``Simulation`` =
``Output file`` (lib/RunLog.py:36) and no duration, and in
CLD_RFSL_mc_sims_01.xlsx twenty-four rows share one ``Output file``. The Nth
appended row is the Nth simulation of the chunk; nothing else is reliable.
"""

from __future__ import annotations

import csv
import re
import time
from dataclasses import dataclass
from pathlib import Path

from .logtail import Tail, tail

# --- the strings Main.py and RunLog.print_summary print -------------------
#
# These mirror literal text in Bryan. ui/tests/test_progress.py greps Main.py
# and lib/RunLog.py for them, so a change there fails the suite loudly instead
# of quietly degrading the progress display to nothing.

BANNER_RE = re.compile(
    r"^Simulation (?P<position>\d+) of (?P<total>\d+): (?P<name>.*)$", re.MULTILINE)
METHOD_RE = re.compile(
    r"^  method: (?P<method>.*?)    method config: (?P<config>.*)$", re.MULTILINE)
FAILURE_RE = re.compile(
    r"^(?P<status>FAILED[^:]*): (?P<name>.*)$", re.MULTILINE)
ALL_COMPLETED_RE = re.compile(r"^All (?P<total>\d+) simulations completed\.$", re.MULTILINE)
SOME_FAILED_RE = re.compile(
    r"^(?P<failed>\d+) of (?P<total>\d+) simulations did not complete:$", re.MULTILINE)
LOG_RENAME_RE = re.compile(r"simulation\(s\) in this batch share a ")

PENDING = "pending"
RUNNING = "running"
COMPLETED = "completed"
FAILED = "FAILED"
FAILED_MISSING_INPUT = "FAILED - missing input"
CANCELLED = "cancelled"
UNKNOWN = "unknown"

FINISHED_STATES = (COMPLETED, FAILED, FAILED_MISSING_INPUT, CANCELLED)


@dataclass(frozen=True)
class SimProgress:
    position: int                 # 1-based, as Main.py counts
    row: int | None               # frame index, when known
    name: str
    status: str = PENDING
    started: str = ""
    ended: str = ""
    error: str = ""
    log_path: Path | None = None
    log_bytes: int = 0
    log_mtime: float = 0.0

    @property
    def finished(self) -> bool:
        return self.status in FINISHED_STATES

    @property
    def failed(self) -> bool:
        return self.status.startswith("FAILED")


@dataclass(frozen=True)
class ChunkProgress:
    index: int
    total: int
    sims: tuple = ()
    finished: bool = False
    returncode: int | None = None
    console: Tail | None = None
    summary: str = ""

    @property
    def completed(self) -> int:
        return sum(1 for sim in self.sims if sim.finished)

    @property
    def failures(self) -> list:
        return [sim for sim in self.sims if sim.failed]

    @property
    def current(self):
        for sim in self.sims:
            if sim.status == RUNNING:
                return sim
        return None

    @property
    def fraction(self) -> float:
        if not self.total:
            return 0.0
        return min(1.0, self.completed / self.total)


def read_run_log(path) -> list:
    """The run log's rows, in the order Bryan appended them."""
    try:
        with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as stream:
            return list(csv.DictReader(stream))
    except (OSError, csv.Error):
        return []


def parse_console(text: str) -> dict:
    """What the console output says about progress."""
    banners = [
        {"position": int(match.group("position")),
         "total": int(match.group("total")),
         "name": match.group("name").strip()}
        for match in BANNER_RE.finditer(text)
    ]
    completed = ALL_COMPLETED_RE.search(text)
    failed_summary = SOME_FAILED_RE.search(text)
    return {
        "banners": banners,
        "last": banners[-1] if banners else None,
        "total": banners[-1]["total"] if banners else None,
        "all_completed": bool(completed),
        "summary": (
            f"All {completed.group('total')} simulations completed."
            if completed else
            f"{failed_summary.group('failed')} of {failed_summary.group('total')} "
            f"simulations did not complete."
            if failed_summary else ""
        ),
        "finished": bool(completed or failed_summary),
    }


def read_chunk_progress(chunk, *, log_paths=None, names=None,
                        returncode=None, alive=True, cancelled_positions=()):
    """Fold the four sources into one view of a chunk.

    ``chunk`` needs ``index``, ``rows``, ``run_log`` and ``console_log`` - a
    ``ChunkFiles`` from core/runwriter.py, or anything shaped like one.
    """
    rows = list(getattr(chunk, "rows", ()) or ())
    total = len(rows)
    names = list(names or [])
    log_paths = list(log_paths or [])

    entries = read_run_log(getattr(chunk, "run_log", ""))
    console_text = tail(getattr(chunk, "console_log", ""))
    console = parse_console(console_text.text)

    sims = []
    for position in range(1, total + 1):
        row = rows[position - 1] if position <= len(rows) else None
        name = names[position - 1] if position <= len(names) else ""
        log_path = log_paths[position - 1] if position <= len(log_paths) else None

        status, started, ended, error = PENDING, "", "", ""
        # Positional match - see the module docstring.
        if position <= len(entries):
            entry = entries[position - 1]
            status = (entry.get("Status") or "").strip() or UNKNOWN
            started = (entry.get("Start time") or "").strip()
            ended = (entry.get("End time") or "").strip()
            error = (entry.get("Error") or "").strip()
            if not name:
                name = (entry.get("Simulation") or "").strip()
        elif console["last"] and console["last"]["position"] == position:
            # Bryan announced it but has not logged it, so it is in flight -
            # unless the process is already gone, in which case it died there.
            status = RUNNING if alive else CANCELLED
        if position in cancelled_positions:
            status = CANCELLED

        stat_bytes, stat_mtime = 0, 0.0
        if log_path:
            try:
                stat = Path(log_path).stat()
                stat_bytes, stat_mtime = stat.st_size, stat.st_mtime
            except OSError:
                pass

        sims.append(SimProgress(
            position=position, row=row, name=name, status=status,
            started=started, ended=ended, error=error, log_path=log_path,
            log_bytes=stat_bytes, log_mtime=stat_mtime,
        ))

    finished = console["finished"] or returncode is not None
    return ChunkProgress(
        index=getattr(chunk, "index", 0),
        total=total,
        sims=tuple(sims),
        finished=bool(finished),
        returncode=returncode,
        console=console_text,
        summary=console["summary"],
    )


def explain_returncode(returncode) -> str:
    """What the exit status means, in Bryan's terms."""
    if returncode is None:
        return "still running"
    if returncode == 0:
        return "every simulation completed (Main.py exits 0)"
    if returncode == 1:
        return ("at least one simulation failed - Main.py:128 exits 1 when "
                "anything did not complete. Check the Status column of the "
                "run log.")
    if returncode < 0:
        return f"killed by signal {-returncode} (cancelled, or the OS stopped it)"
    return f"exited {returncode} before finishing the list - it crashed early"


def explain_error(error: str) -> str:
    """Turn a run-log Error into something actionable.

    ``EOFError`` is the one worth explaining: the UI launches Bryan with
    ``stdin=DEVNULL``, so the ``input()`` prompts in MCScheme.store_simulations
    (lib/MCScheme.py:162) and lib/EnbScheme.py:60 raise instead of hanging a
    windowless process forever. Bryan catches that per row and carries on.

    There used to be a third, in lib/URBSmodel.py, when the URBS executable was
    not where the model config said. That one now warns and carries on instead,
    so a missing executable no longer shows up here - it fails later, and only
    for a simulation that actually runs the model.
    """
    if not error:
        return ""
    if error.startswith("EOFError"):
        return (
            "Bryan asked a question on the console, which a background run "
            "cannot answer, so the simulation failed instead of hanging. The "
            "cause is an output CSV still open in Excel "
            "(MCScheme/EnbScheme.store_simulations)."
        )
    if error.startswith(("FileNotFoundError", "NotADirectoryError")):
        return ("An input file was not where the sims list said. Paths resolve "
                "against the project folder.")
    return ""


def elapsed_text(seconds: float) -> str:
    if seconds < 0:
        return "-"
    hours, rest = divmod(int(seconds), 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {secs:02d}s"
    return f"{secs}s"


def stalled_for(sim: SimProgress, now=None) -> float:
    """Seconds since this simulation's log last grew. 0 when unknown."""
    if not sim.log_mtime:
        return 0.0
    return max(0.0, (now or time.time()) - sim.log_mtime)
