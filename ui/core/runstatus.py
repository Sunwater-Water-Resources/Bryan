"""Each of a study's runs summed up, for the runs panel: how many rows are in each state.

Judging a row's state stats its results and every input it reads, and the panel
is drawn on every page, so the sums are kept: for ``MAX_AGE`` seconds, or until
the sims list or its config changes. The open run's sums are replaced whenever
the Select page judges its rows, so what the panel says of it is never older
than what the table says.
"""

from __future__ import annotations

import threading
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from . import completion
from . import study as studies
from .config import load_sims_config

MAX_AGE = 300.0

# Worst first: what the panel leads with.
ORDER = (completion.NEEDS_PRIOR, completion.STALE, completion.INCOMPLETE,
         completion.NOT_RUN, completion.UNKNOWN, completion.UP_TO_DATE)


@dataclass(frozen=True)
class RunStatus:
    name: str
    config_path: Path | None
    counts: dict = field(default_factory=dict)     # state -> rows
    problem: str = ""                              # why it could not be judged

    @property
    def total(self) -> int:
        return sum(self.counts.values())

    def in_order(self) -> list:
        """(state, rows) for each state present, worst first."""
        return [(state, self.counts[state]) for state in ORDER if self.counts.get(state)]

    def describe(self) -> str:
        if self.problem:
            return self.problem
        if not self.total:
            return "no rows"
        return ", ".join(f"{count} {state}" for state, count in self.in_order())


_CACHE: dict = {}          # resolved config path -> (stamp, taken, RunStatus)
_LOCK = threading.Lock()


def _key(config_path) -> str:
    return str(Path(config_path).resolve())


def _stamp(config_path) -> tuple | None:
    try:
        config = load_sims_config(config_path)
        return (config.config_path.stat().st_mtime_ns,
                config.sims_list_path.stat().st_mtime_ns)
    except Exception:                               # noqa: BLE001 - judged below
        return None


def cached(study, name, *, max_age=MAX_AGE, now=time.monotonic) -> RunStatus | None:
    """The kept sums for a run, or None when there are none young enough."""
    path = study.run_config_path(name)
    if path is None:
        return None
    with _LOCK:
        entry = _CACHE.get(_key(path))
    if entry is None or now() - entry[1] > max_age:
        return None
    return RunStatus(name, entry[2].config_path, entry[2].counts, entry[2].problem)


def judge(study, name, *, now=time.monotonic) -> RunStatus:
    """Sum a run's row states afresh, and keep them."""
    path = study.run_config_path(name)
    if path is None:
        return RunStatus(name, None, problem="not in the study")
    if not path.is_file():
        return RunStatus(name, path, problem="sims_config.json not found")
    try:
        run = studies.open_run(name, path)
    except studies.StudyError as exc:
        return RunStatus(name, path, problem=str(exc))
    states = completion.assess_frame(run.frame, run.config)
    status = RunStatus(name, path, dict(Counter(state.state for state in states.values())))
    with _LOCK:
        _CACHE[_key(path)] = (_stamp(path), now(), status)
    return status


def status_of(study, name, *, max_age=MAX_AGE, now=time.monotonic) -> RunStatus:
    """The kept sums when young enough and the files unchanged, else afresh."""
    kept = cached(study, name, max_age=max_age, now=now)
    if kept is not None and kept.config_path is not None:
        with _LOCK:
            stamp = _CACHE.get(_key(kept.config_path), (None,))[0]
        if stamp is not None and stamp == _stamp(kept.config_path):
            return kept
    return judge(study, name, now=now)


def note(config_path, completions, *, now=time.monotonic) -> None:
    """Keep sums the launcher has just worked out for an open run anyway."""
    counts = dict(Counter(state.state for state in completions.values()))
    with _LOCK:
        _CACHE[_key(config_path)] = (_stamp(config_path), now(),
                                     RunStatus("", Path(config_path), counts))


def forget() -> None:
    with _LOCK:
        _CACHE.clear()
