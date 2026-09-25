"""Whether the results a table or figure is drawn from are out of date.

A report table read from a group whose results are older than an input they
were made from - a rating curve edited after the run - looks exactly like a
good one. So every place that draws on a group's results says, beside what it
draws, when that group is not up to date: "E013 FSL RFSL_Design_GWL1p3 is
stale: CLD_E013.sq changed after the run". Copy for Word says it again.

The states are the Select page's (core/completion.py), judged row by row over
the group's durations. A duration not run yet is said too: the group's
envelope and critical duration are then taken without it.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from . import completion
from . import study as studies

ATTENTION = (completion.STALE, completion.INCOMPLETE, completion.NOT_RUN,
             completion.NEEDS_PRIOR)


@dataclass(frozen=True)
class GroupState:
    run: str                        # the study's name for the run, or the sims list's
    group: str
    counts: dict = field(default_factory=dict)       # state -> durations
    changed: tuple = ()             # inputs newer than the results, newest first

    @property
    def total(self) -> int:
        return sum(self.counts.values())

    @property
    def is_current(self) -> bool:
        return not any(self.counts.get(state) for state in ATTENTION)

    def messages(self) -> list:
        name = f"{self.run} {self.group}".strip()
        said = []
        stale = self.counts.get(completion.STALE, 0)
        if stale:
            what = (f"{self.changed[0].name} changed after the run" if self.changed
                    else "an input changed after the run")
            said.append(f"{name} is stale: {what}{self._of(stale)}.")
        incomplete = self.counts.get(completion.INCOMPLETE, 0)
        if incomplete:
            said.append(f"{name} is incomplete{self._of(incomplete)}: a run stopped "
                        f"early.")
        missing = self.counts.get(completion.NOT_RUN, 0) + \
            self.counts.get(completion.NEEDS_PRIOR, 0)
        if missing:
            said.append(f"{name}: {missing} of {self.total} durations not run yet, so "
                        f"what is drawn from it leaves them out.")
        return said

    def _of(self, count) -> str:
        return "" if count == self.total else f" ({count} of {self.total} durations)"


def of_rows(run_name, group, frame, config, rows) -> GroupState:
    """Judge the rows of one group - ``rows`` are frame indices."""
    states = completion.assess_frame(frame, config, rows=list(rows))
    changed = sorted({state.newest_input for state in states.values()
                      if state.state == completion.STALE and state.newest_input},
                     key=lambda path: -_mtime(path))
    return GroupState(run_name, group,
                      dict(Counter(state.state for state in states.values())),
                      tuple(changed))


def _mtime(path: Path) -> float:
    try:
        return path.stat().st_mtime
    except OSError:
        return 0.0


def of_group(study, run_name, group) -> GroupState | None:
    """A study run's group judged, or None when the run cannot be read (the
    table or figure says so itself) or holds no such group."""
    if not run_name or not group:
        return None
    try:
        run = study.open_run(run_name)
    except studies.StudyError:
        return None
    rows = run.rows_in(group)
    if not rows:
        return None
    return of_rows(run_name, group, run.frame, run.config, rows)


def of_open_group(project, run_name, group) -> list:
    """What to say about a group of the open sims list (the Results and Ensemble
    pages' groups are its group keys)."""
    keys = project.group_keys()
    rows = [index for index in project.frame.index if keys[index] == group]
    if not rows:
        return []
    state = of_rows(run_name, group, project.frame, project.config, rows)
    return [] if state.is_current else state.messages()


def for_sources(study, sources) -> list:
    """What to say about the {run, group} sources a table or figure reads, once
    per group, in the order they are first read."""
    seen, said = set(), []
    for source in sources:
        key = (source.get("run") or "", source.get("group") or "")
        if key in seen:
            continue
        seen.add(key)
        state = of_group(study, *key)
        if state is not None and not state.is_current:
            said += state.messages()
    return said
