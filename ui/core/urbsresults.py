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


# -- reading a run's results ---------------------------------------------------------

_READ: dict = {}


def _cached(path: Path, reader):
    """``reader(path)``, kept per file version: a regional run's files are read on every
    Reload and redraw, and they only change when URBS writes them again."""
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path), stat.st_mtime, stat.st_size, reader.__name__)
    if key not in _READ:
        _READ[key] = reader(path)
    return _READ[key]


def forget_cached() -> None:
    _READ.clear()


def _parse_levels(path: Path) -> pd.DataFrame:
    """URBS's ``<run>.h``: fixed-width text, one column per level location.

    The names are on the ninth line, cut to eight characters by URBS, and the data start
    after the ``Time  metres`` line and the rule under it, ending at ``*END``. A name that
    appears twice gets ``_1`` on its second, as ``UrbsModel.get_hydrographs`` does."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    seen: dict = {}
    columns = []
    for name in lines[8].split():
        columns.append(name if name not in seen else f"{name}_{seen[name]}")
        seen[name] = seen.get(name, 0) + 1
    start = next((i for i, line in enumerate(lines)
                  if line.split()[:2] == ["Time", "metres"]), None)
    if start is None:
        raise ValueError(f"{path.name} has no 'Time  metres' line")
    rows = []
    for line in lines[start + 2:]:
        if line.startswith("*END"):
            break
        parts = line.split()
        if len(parts) == len(columns) + 1:
            rows.append([float(part) for part in parts])
    return pd.DataFrame(rows, columns=["Time"] + columns).set_index("Time")


def read_levels(results_folder, run_name) -> pd.DataFrame | None:
    """A run's level hydrographs by location (eight characters), or None if not there."""
    try:
        return _cached(Path(results_folder) / f"{run_name}.h", _parse_levels)
    except (OSError, ValueError, IndexError):
        return None


# -- the dams, against their reference levels ------------------------------------------

#: The level the comparison is drawn against, where a dam has one.
CREST = "DCL"


def dam_setup(model_config) -> tuple:
    """(dams, levels, problems) from the regional model config, read as JSON.

    ``dams`` is every dam the model sets a starting level for - the main one (max_keys
    level) and the additional ones - in that order. ``levels`` is ``dam_levels`` with its
    notes dropped. Problems are reported, not raised: Bryan refuses a bad dam_levels when
    it opens the model, and the page still shows what it can."""
    try:
        data = json.loads(Path(model_config).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return [], {}, []
    dams = []
    main = (data.get("max_keys") or {}).get("level")
    if main:
        dams.append(str(main))
    dams += [str(entry.get("location")) for entry in data.get("additional_dams") or []
             if isinstance(entry, dict) and entry.get("location")]
    levels, problems = {}, []
    for location, named in (data.get("dam_levels") or {}).items():
        if str(location).startswith("_"):
            continue
        if not isinstance(named, dict):
            problems.append(f"dam_levels for {location} is not a mapping of level names")
            continue
        levels[str(location)] = {str(key): value for key, value in named.items()
                                 if not str(key).startswith("_")}
        if str(location) not in dams:
            dams.append(str(location))
    return dams, levels, problems


def level_at(value, gwl) -> float | None:
    """A reference level at a run's warming level: the number, or that GWL's entry of a
    per-GWL map. None where it is not given."""
    if isinstance(value, dict):
        if gwl is None:
            return None
        for key, level in value.items():
            try:
                if abs(float(key) - float(gwl)) < 1e-9:
                    return float(level)
            except (TypeError, ValueError):
                continue
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


@dataclass
class DamPeak:
    dam: str
    peak: float | None = None
    time: float | None = None
    margins: dict | None = None     # level name -> peak less level; None where not given
    start: float | None = None      # the level the run started at

    @property
    def found(self) -> bool:
        return self.peak is not None


def dam_peak(levels, dam, references, gwl) -> DamPeak:
    """One dam's peak in one run, and its margin to each reference level."""
    column = dam[:8]
    margins = {name: None for name in references}
    if levels is None or column not in levels:
        return DamPeak(dam, margins=margins)
    series = levels[column].dropna()
    if series.empty:
        return DamPeak(dam, margins=margins)
    peak = float(series.max())
    for name, value in references.items():
        level = level_at(value, gwl)
        margins[name] = None if level is None else peak - level
    return DamPeak(dam, peak, float(series.idxmax()), margins, float(series.iloc[0]))


def margin_text(margin) -> str:
    return "not given" if margin is None else f"{margin:+.2f} m"


def dams_rows(plan, results_folder, dams, levels) -> list:
    """The Dams table: per finished event, each dam's start, peak, time of peak and margin
    to each reference level - the peak less the level, positive above it. A level the config
    does not give, or not for the run's GWL, says 'not given', never blank."""
    # every dam gets every level name any dam has, so a missing one shows as 'not given'
    names = list(dict.fromkeys(name for dam in dams for name in levels.get(dam, {})))
    rows = []
    for storm in plan.storms:
        if not storm.filename:
            continue
        run = os.path.splitext(storm.filename)[0]
        frame = read_levels(results_folder, run)
        if frame is None:
            continue
        row = {"loading": storm.loading, "run": run}
        for dam in dams:
            references = {name: levels.get(dam, {}).get(name) for name in names}
            peak = dam_peak(frame, dam, references, storm.gwl)
            row[f"{dam} start"] = "" if peak.start is None else f"{peak.start:.2f}"
            row[f"{dam} peak"] = "" if not peak.found else f"{peak.peak:.2f}"
            row[f"{dam} at (h)"] = "" if not peak.found else f"{peak.time:g}"
            for name, margin in peak.margins.items():
                row[f"{dam} to {name}"] = margin_text(margin) if peak.found else ""
        rows.append(row)
    return rows


def reaching_crest(rows, dams) -> int:
    """How many events have every dam that has a crest level at or above it."""
    with_crest = [dam for dam in dams if any(f"{dam} to {CREST}" in row for row in rows)]
    if not with_crest:
        return 0
    count = 0
    for row in rows:
        texts = [row.get(f"{dam} to {CREST}", "") for dam in with_crest]
        if all(text.endswith(" m") and float(text.split()[0]) >= 0 for text in texts):
            count += 1
    return count


def crest_series(frame, dam, references, gwl):
    """The dam's level less its crest level through the run, or None without a crest."""
    crest = level_at(references.get(CREST), gwl)
    column = dam[:8]
    if frame is None or crest is None or column not in frame:
        return None
    return frame[column] - crest


def record_runs(selection) -> list:
    """(run name, gwl) for every storm a selection's last run wrote, from its record."""
    out = []
    for run, row in storms_record(selection).items():
        try:
            gwl = float(row.get("gwl", ""))
        except ValueError:
            gwl = None
        out.append((run, gwl))
    return out


def margin_points(selections, results_folder, dams, levels, label_of) -> dict:
    """Across every selection in the project: per group, a point per finished run - the first
    two dams' margins to their crest levels, with the run's name. Empty without two crests."""
    if len(dams) < 2:
        return {}
    first, second = dams[:2]
    points: dict = {}
    for selection in selections:
        for run, gwl in record_runs(selection):
            frame = read_levels(results_folder, run)
            if frame is None:
                continue
            x = dam_peak(frame, first, {CREST: levels.get(first, {}).get(CREST)}, gwl)
            y = dam_peak(frame, second, {CREST: levels.get(second, {}).get(CREST)}, gwl)
            if x.margins.get(CREST) is None or y.margins.get(CREST) is None:
                continue
            points.setdefault(label_of(selection), []).append(
                (x.margins[CREST], y.margins[CREST], run))
    return points
