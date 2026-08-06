"""Spawning Bryan, watching it, and stopping it.

Bryan has no callable seam: ``MonteCarloSimulator(sim, filepaths, test_runs)``
does the entire run in ``__init__``, which also reassigns ``sys.stdout``
globally to a ``Logger``. So the UI runs Bryan as a subprocess - never
in-process - and derives progress by reading files.

Three decisions here are load-bearing.

**cwd = the project folder.** Nothing in Main.py chdirs, so every path *inside*
the sims list resolves against the process CWD. The launcher .bat does
``cd /D "%~dp0"`` first; so does this.

**stdout to a file, not a pipe.** No reader thread to keep fed, nothing
deadlocks on a full pipe buffer, and the capture survives the UI being closed.
``Logger`` keeps ``sys.__stdout__`` and writes to both, so the file gets
everything: Main.py's banners, each simulator's output, tracebacks and the
closing summary. ``ConsoleTitle`` only emits escape codes when
``sys.__stdout__.isatty()``, and a file is not a tty, so nothing corrupts it.

**stdin = DEVNULL.** Bryan calls ``input()`` in three places - MCScheme.store_simulations
and lib/EnbScheme.py:60 when an output CSV is open in Excel, lib/URBSmodel.py:41
when the URBS exe is missing. A windowless background process would hang there
forever. With DEVNULL they raise EOFError at once, Main.py catches it per row,
and the run log records a failure. core/progress.explain_error spells that out
when it happens.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

from . import runstate
from .runstate import (CANCELLED, COMPLETED, DIED, FAILED, QUEUED, RUNNING,
                       ChunkRecord, RunRecord)

try:
    import psutil
except ImportError:      # pragma: no cover - psutil is a hard dependency
    psutil = None

TERMINATE_GRACE_SECONDS = 5.0


class LaunchError(Exception):
    """Bryan could not be started."""


def _popen_kwargs(show_console: bool = False) -> dict:
    """Platform flags that make a whole process tree killable."""
    if os.name == "nt":
        flags = subprocess.CREATE_NEW_PROCESS_GROUP
        if not show_console:
            # URBS is launched through cmd.exe with shell=True, and some
            # console software dislikes having no console at all, so this is
            # behind a setting.
            flags |= getattr(subprocess, "CREATE_NO_WINDOW", 0)
        return {"creationflags": flags}
    return {"start_new_session": True}      # own process group, so killpg works


class RunManager:
    """Owns the Bryan processes for the app's lifetime.

    Shaped after ``swe2d_ui/runner.py``'s RunManager - start / poll / cancel /
    registry - but the transport is a subprocess and a log file rather than a
    multiprocessing Queue, because Bryan has no in-process seam.
    """

    def __init__(self, bryan_python=None, bryan_main=None, show_console=False):
        self.bryan_python = str(bryan_python) if bryan_python else sys.executable
        self.bryan_main = str(bryan_main) if bryan_main else ""
        self.show_console = show_console
        self.runs: dict = {}                  # run_id -> RunRecord
        self._processes: dict = {}            # (run_id, chunk index) -> Popen
        self._handles: dict = {}              # (run_id, chunk index) -> file
        self._lock = threading.Lock()

    # -- starting ---------------------------------------------------------

    def submit(self, record: RunRecord) -> RunRecord:
        """Register a run and start as many chunks as the parallelism allows."""
        if not self.bryan_main:
            raise LaunchError(
                "Bryan's Main.py has not been configured - set it on the "
                "Project page before launching."
            )
        if not Path(self.bryan_main).is_file():
            raise LaunchError(f"Main.py not found: {self.bryan_main}")
        if not Path(self.bryan_python).exists():
            raise LaunchError(f"Python interpreter not found: {self.bryan_python}")

        with self._lock:
            self.runs[record.run_id] = record
            record.save()
        self.pump()
        return record

    def pump(self) -> None:
        """Start queued chunks while there is room. Safe to call any time."""
        with self._lock:
            for record in self.runs.values():
                if record.stop_requested:
                    continue
                live = sum(1 for chunk in record.chunks if chunk.status == RUNNING)
                for chunk in record.chunks:
                    if live >= record.max_parallel:
                        break
                    if chunk.status == QUEUED:
                        self._start_chunk(record, chunk)
                        live += 1
                record.save()

    def _start_chunk(self, record: RunRecord, chunk: ChunkRecord) -> None:
        project = Path(record.project_folder)
        config = Path(chunk.config)
        try:
            relative = os.path.relpath(config, project)
        except ValueError:                     # different drive on Windows
            relative = str(config)

        console = Path(chunk.console_log)
        console.parent.mkdir(parents=True, exist_ok=True)
        handle = console.open("ab")

        environment = dict(os.environ)
        # -u alone is not enough when stdout is a file: without these the child
        # block-buffers 8 KB and the live view lags by minutes.
        environment["PYTHONUNBUFFERED"] = "1"
        environment.setdefault("PYTHONIOENCODING", "utf-8")

        command = [self.bryan_python, "-u", self.bryan_main, relative]
        try:
            process = subprocess.Popen(
                command,
                cwd=str(project),
                stdin=subprocess.DEVNULL,
                stdout=handle,
                stderr=subprocess.STDOUT,
                env=environment,
                **_popen_kwargs(self.show_console),
            )
        except OSError as exc:
            handle.close()
            chunk.status = FAILED
            chunk.note = f"could not start Bryan: {exc}"
            chunk.ended = time.time()
            return

        chunk.status = RUNNING
        chunk.pid = process.pid
        chunk.started = time.time()
        chunk.note = ""
        try:
            chunk.create_time = psutil.Process(process.pid).create_time() if psutil else None
        except Exception:
            chunk.create_time = None

        self._processes[(record.run_id, chunk.index)] = process
        self._handles[(record.run_id, chunk.index)] = handle

    # -- watching ---------------------------------------------------------

    def poll(self) -> None:
        """Reap finished chunks and start whatever the free slots allow."""
        with self._lock:
            for record in self.runs.values():
                changed = False
                for chunk in record.chunks:
                    if chunk.status != RUNNING:
                        continue
                    process = self._processes.get((record.run_id, chunk.index))
                    if process is None:
                        if not self._external_alive(chunk):
                            self._finish_unknown(chunk)
                            changed = True
                        continue
                    code = process.poll()
                    if code is None:
                        continue
                    chunk.returncode = code
                    chunk.ended = time.time()
                    chunk.status = COMPLETED if code == 0 else FAILED
                    if chunk.status == FAILED and code not in (0, 1):
                        chunk.note = f"Bryan exited {code} before finishing the list"
                    self._close(record.run_id, chunk.index)
                    changed = True
                if changed:
                    record.save()
        self.pump()

    def _external_alive(self, chunk: ChunkRecord) -> bool:
        """Is the recorded pid still THIS chunk's process?

        The create-time token is what makes this safe: a bare pid check would
        happily reattach to whatever reused the number.
        """
        if not chunk.pid or psutil is None:
            return False
        try:
            process = psutil.Process(chunk.pid)
            if chunk.create_time is None:
                return process.is_running()
            return (process.is_running()
                    and abs(process.create_time() - chunk.create_time) < 1.0)
        except Exception:
            return False

    def _finish_unknown(self, chunk: ChunkRecord) -> None:
        """A chunk whose process is gone and which we never reaped.

        Reconcile from evidence rather than guessing: the console log's closing
        summary is what Main.py/RunLog.print_summary write at the very end.
        """
        from .progress import parse_console
        from .logtail import tail

        chunk.ended = time.time()
        console = parse_console(tail(chunk.console_log).text)
        if console["all_completed"]:
            chunk.status = COMPLETED
            chunk.returncode = 0
        elif console["finished"]:
            chunk.status = FAILED
            chunk.returncode = 1
            chunk.note = console["summary"]
        else:
            chunk.status = DIED
            chunk.note = ("the process is gone but Bryan never printed its "
                          "closing summary - it was killed, or it crashed")

    def reattach(self, run_root) -> list:
        """Re-find runs from disk after a UI restart."""
        found = []
        for record in runstate.discover(run_root):
            with self._lock:
                existing = self.runs.get(record.run_id)
                if existing is not None:
                    found.append(existing)
                    continue
                for chunk in record.chunks:
                    if chunk.status == RUNNING and not self._external_alive(chunk):
                        self._finish_unknown(chunk)
                self.runs[record.run_id] = record
                record.save()
            found.append(record)
        return found

    # -- stopping ---------------------------------------------------------

    def request_stop(self, run_id: str) -> None:
        """Let the running chunk finish; start nothing more.

        The safe default - no signal reaches Bryan, so nothing is left half
        written.
        """
        with self._lock:
            record = self.runs.get(run_id)
            if record is None:
                return
            record.stop_requested = True
            for chunk in record.chunks:
                if chunk.status == QUEUED:
                    chunk.status = CANCELLED
                    chunk.note = "not started - stop was requested"
            record.save()

    def cancel(self, run_id: str, chunk_index=None) -> None:
        """Kill now, including the URBS/RORB processes underneath.

        Killing only the Bryan process orphans ``cmd.exe`` and ``urbs32.exe``,
        which keep running and keep writing. The in-flight simulation never
        reaches ``RunLog.write_entry``, so it gets **no run-log row at all** -
        which is why its status comes from here, not from the log.
        """
        with self._lock:
            record = self.runs.get(run_id)
            if record is None:
                return
            record.stop_requested = True
            for chunk in record.chunks:
                if chunk_index is not None and chunk.index != chunk_index:
                    continue
                if chunk.status == QUEUED:
                    chunk.status = CANCELLED
                    chunk.note = "not started - the run was cancelled"
                    continue
                if chunk.status != RUNNING:
                    continue
                process = self._processes.get((record.run_id, chunk.index))
                pid = process.pid if process else chunk.pid
                if pid:
                    kill_tree(pid)
                if process is not None:
                    try:
                        process.wait(timeout=TERMINATE_GRACE_SECONDS)
                    except subprocess.TimeoutExpired:
                        pass
                    chunk.returncode = process.returncode
                chunk.status = CANCELLED
                chunk.ended = time.time()
                chunk.note = ("cancelled - the simulation that was running has "
                              "no run-log entry, and its working folder and log "
                              "are left part-written")
                self._close(record.run_id, chunk.index)
            record.save()

    def shutdown(self) -> None:
        for run_id in list(self.runs):
            self.cancel(run_id)

    def _close(self, run_id, chunk_index) -> None:
        handle = self._handles.pop((run_id, chunk_index), None)
        if handle is not None:
            try:
                handle.close()
            except OSError:
                pass
        self._processes.pop((run_id, chunk_index), None)

    # -- queries ----------------------------------------------------------

    def live_runs(self) -> list:
        return [record for record in self.runs.values() if record.is_live]

    def returncode_of(self, run_id, chunk_index):
        record = self.runs.get(run_id)
        chunk = record.chunk(chunk_index) if record else None
        return chunk.returncode if chunk else None

    def is_alive(self, run_id, chunk_index) -> bool:
        record = self.runs.get(run_id)
        chunk = record.chunk(chunk_index) if record else None
        if chunk is None or chunk.status != RUNNING:
            return False
        process = self._processes.get((run_id, chunk_index))
        if process is not None:
            return process.poll() is None
        return self._external_alive(chunk)


def kill_tree(pid: int, grace: float = TERMINATE_GRACE_SECONDS) -> None:
    """Kill a process and everything it started."""
    if psutil is not None:
        try:
            parent = psutil.Process(pid)
            children = parent.children(recursive=True)
            for target in [*children, parent]:
                try:
                    target.terminate()
                except Exception:
                    pass
            gone, alive = psutil.wait_procs([*children, parent], timeout=grace)
            for target in alive:
                try:
                    target.kill()
                except Exception:
                    pass
            return
        except Exception:
            pass

    # Fallbacks, in case psutil is unavailable or the process vanished mid-walk.
    if os.name == "nt":
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(pid)],
                       capture_output=True, check=False)
        return
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
        time.sleep(min(grace, 1.0))
        os.killpg(os.getpgid(pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError, OSError):
        pass
