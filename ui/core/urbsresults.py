"""The regional model's runs, read back for the Downstream page's Results tab.

Where everything is comes from the regional model config, through the same
``lib.ConfigPaths`` helpers ``DownstreamStorms.py`` ran URBS with, so the page looks in
the folder the results were written to and not one it guessed.

**A run's status is read off what is on disk and what the page knows**, never off a
clock. URBS's results present is *finished*. Otherwise: while the generator the page
started is still going, an event whose batch file is there is the one URBS is *running*
now, and the rest are *waiting*. Once it has stopped, a batch file left behind is either a
run that failed or one written to run by hand: if the page ran URBS this session it is
*failed*, if not it is *not run*. A launcher restarted mid-run knows neither, and says so.

pandas and the standard library only; no nicegui.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .bryan import config_paths

FINISHED = "finished"
RUNNING = "running"
WAITING = "waiting"
FAILED = "failed"
NOT_RUN = "not run"
NOT_WRITTEN = "not written"
UNKNOWN = "not finished"

#: Shown beside each status.
EXPLAIN = {
    FINISHED: "URBS's results are in the results folder",
    RUNNING: "URBS is running this event now",
    WAITING: "the generator has not reached this event yet",
    FAILED: "URBS stopped without results - the batch file is left to run by hand; "
            "see urbserr.log in the results folder",
    NOT_RUN: "the storm and its batch file are written; URBS was not asked to run it",
    NOT_WRITTEN: "no storm file yet",
    UNKNOWN: "the launcher was restarted during the run - Reload once URBS has finished",
}


class ModelConfigError(Exception):
    """A regional model config that cannot be read."""


@dataclass(frozen=True)
class Folders:
    storms: Path
    results: Path


def folders(model_config) -> Folders:
    """Where the regional model's storms and results are, as DownstreamStorms.py put them."""
    path = Path(model_config)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ModelConfigError(f"{path.name} could not be read ({exc})") from exc
    paths = config_paths()
    try:
        storms = paths.resolve(str(path.parent), data["storms_folder"])
        results = paths.regional_results_folder(str(path), data)
    except KeyError as exc:
        raise ModelConfigError(f"{path.name} has no {exc.args[0]}") from exc
    return Folders(Path(storms), Path(results))


def finished(results_folder, run_name) -> bool:
    """URBS has written this run's results: the flows, which every run writes."""
    folder = Path(results_folder)
    return (folder / f"{run_name}.q").is_file() or (folder / f"{run_name}_i.csv").is_file()


def run_status(found: Folders, run_name: str, storm_file: str, *,
               generator_running: bool = False, ran_urbs: bool | None = None) -> str:
    """One run's status. ``ran_urbs`` is True or False for a run this page started and
    has seen finish, None where the page does not know (a restarted launcher)."""
    if finished(found.results, run_name):
        return FINISHED
    batch = (found.results / f"_run_{run_name}.bat").is_file()
    if generator_running:
        return RUNNING if batch and ran_urbs else WAITING
    if not (found.storms / storm_file).is_file():
        return NOT_WRITTEN
    if not batch:
        return NOT_RUN
    if ran_urbs is None:
        return UNKNOWN
    return FAILED if ran_urbs else NOT_RUN


def storms_record(selection) -> dict:
    """What went into each storm, from ``<selection>_downstream_storms.csv``, by run name."""
    path = Path(os.path.splitext(str(selection))[0] + "_downstream_storms.csv")
    if not path.is_file():
        return {}
    try:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception:                               # noqa: BLE001 - a record, not a must
        return {}
    if "filename" not in frame:
        return {}
    names = (frame["run_name"] if "run_name" in frame
             else frame["filename"].map(lambda name: os.path.splitext(name)[0]))
    return {name: row.to_dict() for name, (_, row) in zip(names, frame.iterrows())}


def _number(text, digits) -> str:
    try:
        return f"{float(text):,.{digits}f}"
    except (TypeError, ValueError):
        return ""


def runs_rows(plan, found: Folders, record: dict, *, generator_running=False,
              ran_urbs=None) -> list:
    """The Runs table: one row per storm the plan names."""
    rows = []
    for storm in plan.storms:
        if not storm.filename:
            continue
        run = os.path.splitext(storm.filename)[0]
        status = run_status(found, run, storm.filename,
                            generator_running=generator_running, ran_urbs=ran_urbs)
        said = record.get(run, {})
        flags = [text for text in (said.get("embedded_bursts", ""), said.get("duration_dip", ""))
                 if text and text != "No embedded bursts"]
        rows.append({"loading": storm.loading, "event": storm.event_label, "run": run,
                     "status": status, "why": EXPLAIN[status],
                     "IL (mm)": _number(said.get("initial_loss"), 1),
                     "CL (mm/h)": _number(said.get("continuing_loss"), 1),
                     "mean depth (mm)": _number(said.get("mean_depth_mm"), 1),
                     "pattern": said.get("temporal_pattern", ""),
                     "starting storage": said.get("antecedent", ""),
                     "flags": "; ".join(flags)})
    return rows


def counts(rows) -> str:
    """'5 finished, 1 running, 1 waiting', in a fixed order, leaving out what is none."""
    order = (FINISHED, RUNNING, WAITING, FAILED, NOT_RUN, UNKNOWN, NOT_WRITTEN)
    tally = {status: sum(row["status"] == status for row in rows) for status in order}
    return ", ".join(f"{tally[s]} {s}" for s in order if tally[s]) or "no runs"
