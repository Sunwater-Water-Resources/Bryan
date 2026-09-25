"""A child process run to the end, or until it is cancelled or runs out of time.

The Lake record's steps and the Lake levels bands run Bryan's util scripts for
minutes at a time. ``subprocess.run`` offers no way to stop one early, so they
run through here: the output goes to a temporary file rather than a pipe (a
long log cannot then fill the pipe and stall the script), and the process is
checked every ``poll`` seconds for the cancel event and the time limit.
"""

from __future__ import annotations

import subprocess
import tempfile
import time
from dataclasses import dataclass


@dataclass
class Finished:
    returncode: int
    output: str
    cancelled: bool = False
    timed_out: bool = False


def _stop(process) -> None:
    """Stop the process and anything it started."""
    try:
        import psutil
        for child in psutil.Process(process.pid).children(recursive=True):
            try:
                child.kill()
            except psutil.Error:
                pass
    except Exception:                    # noqa: BLE001 - psutil missing, or the pid gone
        pass
    try:
        process.kill()
    except OSError:
        pass
    process.wait()


def run(argv, *, cwd=None, env=None, timeout=None, cancel=None, poll=0.2) -> Finished:
    """Run ``argv``; ``cancel`` is a threading.Event that stops it when set."""
    with tempfile.TemporaryFile() as sink:
        try:
            process = subprocess.Popen([str(part) for part in argv], stdout=sink,
                                       stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                       cwd=None if cwd is None else str(cwd), env=env)
        except OSError as exc:
            return Finished(1, f"could not start {argv[0]}: {exc}")
        started = time.monotonic()
        cancelled = timed_out = False
        while True:
            try:
                process.wait(timeout=poll)
                break
            except subprocess.TimeoutExpired:
                pass
            if cancel is not None and cancel.is_set():
                _stop(process)
                cancelled = True
                break
            if timeout is not None and time.monotonic() - started > timeout:
                _stop(process)
                timed_out = True
                break
        sink.seek(0)
        output = sink.read().decode("utf-8", errors="replace")
    if cancelled:
        output += "\nCancelled."
    if timed_out:
        output += f"\nStopped: it ran past {timeout:.0f} s."
    return Finished(process.returncode if not (cancelled or timed_out) else 1,
                    output, cancelled, timed_out)
