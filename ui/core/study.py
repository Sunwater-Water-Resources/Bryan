"""The study file: one dam study's runs and the report tables drawn from them.

A sims_config.json is one run - E001, say, at one full supply level. A design
flood report is written across several: the RFSL and FSL sims lists are separate
configs, the adopted run sits beside the one it replaced, and a sensitivity comes
from a third. The study file is the level above: it names those runs once, and
says which group of which run each report table is drawn from.

It is a plain JSON file the analyst keeps at the top of the study, and every path
in it is stored relative to that file wherever it can be (``portable``), so a
study copied from ``C:\\PythonProjects`` to an external disk still opens. Nothing
here imports nicegui; ``pages/report.py`` is the view.

    {
      "bryan_study": 1,
      "name": "Juniper Creek Dam design flood hydrology",
      "runs": [{"name": "E001 RFSL",
                "sims_config": "03_DESIGN/runs/E001/JCD_RFSL_mc_sims_01.json"}],
      "tables": [{"id": "near-term-rfsl", "kind": "design_floods",
                  "title": "Table 26: Near-Term design hydrology results",
                  "source": {"run": "E001 RFSL", "group": "RFSL_Design_GWL1p3"}}]
    }

A table names its run by the run's **name**, never by path, so pointing the
study at a re-run (E002 for E001) is one edit to ``runs`` rather than one per
table.
"""

from __future__ import annotations

import copy
import os
import re
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

from . import grouping
from .config import ConfigError, load_sims_config
from .paths import atomic_write_json, clean_path_text, normalise_sep, read_json
from .simslist import read_sims_list

FORMAT_KEY = "bryan_study"
FORMAT_VERSION = 1
DEFAULT_NAME = "bryan_study.json"

SLUG_RE = re.compile(r"[^a-z0-9]+")


class StudyError(Exception):
    """A study file that cannot be used."""


# -- paths ---------------------------------------------------------------------

def portable(base: Path, text) -> str:
    """A path as it should be stored: relative when it sits under ``base``.

    Forward slashes either way, so the file reads the same on Linux.
    """
    text = clean_path_text(text or "")
    if not text:
        return ""
    path = Path(normalise_sep(text))
    if not path.is_absolute():
        return Path(normalise_sep(text)).as_posix()
    try:
        if path.resolve().is_relative_to(Path(base).resolve()):
            return Path(os.path.relpath(path.resolve(), Path(base).resolve())).as_posix()
    except (OSError, ValueError):
        pass
    return text


def resolve(base: Path, text) -> Path | None:
    text = clean_path_text(text or "")
    if not text:
        return None
    path = Path(normalise_sep(text))
    return path if path.is_absolute() else (Path(base) / path).resolve()


def guess_run_name(path) -> str:
    """A run's name from its sims_config.json path, as the report names runs.

    'runs/E001/JCD_RFSL_mc_sims_01.json' -> 'E001 RFSL'; a PMF list gets 'PMF'.
    """
    text = clean_path_text(path or "")
    if not text:
        return ""
    local = Path(normalise_sep(text))
    parent = local.parent.name
    stem = local.stem.upper()
    state = next((token for token in ("RFSL", "FSL") if token in stem.split("_")), "")
    prefix = "PMF" if stem.startswith("PMF") else ""
    return " ".join(part for part in (parent, prefix, state) if part) or local.stem


def slug(text: str, limit: int = 40) -> str:
    """A table id from its title: short, because a caption is a sentence."""
    text = SLUG_RE.sub("-", str(text).lower()).strip("-")
    if len(text) > limit:
        text = text[:limit].rsplit("-", 1)[0]
    return text or "table"


# -- runs ----------------------------------------------------------------------

@dataclass
class RunData:
    """One run opened for reading: its config, its sims list and its groups."""

    name: str
    config: object                 # core.config.SimsConfig
    frame: object                  # the sims list, pandas
    keys: object                   # group key per row

    @property
    def project_folder(self) -> Path:
        return self.config.project_folder

    def groups(self) -> list:
        return grouping.groups_in_order(self.frame)

    def rows_in(self, group: str) -> list:
        return list(self.keys[self.keys == group].index)


# Sims lists are read off disks that may be slow (an external drive), and a page
# builds every table on every refresh, so a run is read once per change to its
# sims list or config.
_RUN_CACHE: dict = {}


def open_run(name: str, config_path: Path) -> RunData:
    """Read a run, reusing the last read while neither file has changed."""
    config_path = Path(config_path)
    try:
        config = load_sims_config(config_path)
    except ConfigError as exc:
        raise StudyError(str(exc)) from exc
    try:
        stamp = (config.config_path.stat().st_mtime_ns,
                 config.sims_list_path.stat().st_mtime_ns)
    except OSError as exc:
        raise StudyError(f"{name}: the sims list could not be read "
                         f"({config.sims_list_path})") from exc
    key = str(config.config_path)
    cached = _RUN_CACHE.get(key)
    if cached is not None and cached[0] == stamp:
        data = cached[1]
        return RunData(name=name, config=data.config, frame=data.frame, keys=data.keys)
    try:
        frame = read_sims_list(config.sims_list_path).frame
    except Exception as exc:                     # noqa: BLE001 - report, never crash
        raise StudyError(f"{name}: the sims list could not be read ({exc})") from exc
    data = RunData(name=name, config=config, frame=frame,
                   keys=grouping.add_group_keys(frame))
    _RUN_CACHE[key] = (stamp, data)
    return data


def forget_runs() -> None:
    _RUN_CACHE.clear()


# -- the study -----------------------------------------------------------------

@dataclass
class Study:
    """An open study file."""

    path: Path
    name: str = ""
    runs: list = field(default_factory=list)      # [{"name", "sims_config"}]
    tables: list = field(default_factory=list)    # table specs, see reporttables
    extra: dict = field(default_factory=dict)     # keys this version does not use

    @property
    def folder(self) -> Path:
        return self.path.parent

    # -- runs --------------------------------------------------------------

    def run_names(self) -> list:
        return [run["name"] for run in self.runs]

    def run_entry(self, name: str) -> dict | None:
        return next((run for run in self.runs if run["name"] == name), None)

    def run_config_path(self, name: str) -> Path | None:
        entry = self.run_entry(name)
        return resolve(self.folder, entry["sims_config"]) if entry else None

    def run_named_for(self, config_path) -> str | None:
        """The study's name for a sims_config.json, or None when it is not a run."""
        target = Path(config_path).resolve()
        for entry in self.runs:
            path = resolve(self.folder, entry["sims_config"])
            if path is not None and path.resolve() == target:
                return entry["name"]
        return None

    def open_run(self, name: str) -> RunData:
        path = self.run_config_path(name)
        if path is None:
            raise StudyError(f"the study has no run called {name!r}")
        return open_run(name, path)

    def add_run(self, name: str, sims_config) -> dict:
        name = str(name).strip()
        if not name:
            raise StudyError("a run needs a name")
        if self.run_entry(name) is not None:
            raise StudyError(f"there is already a run called {name!r}")
        path = resolve(self.folder, sims_config)
        if path is None:
            raise StudyError("give the run's sims_config.json")
        if path.is_dir():
            configs = sorted(item.name for item in path.glob("*.json"))[:6]
            raise StudyError(f"{path} is a folder - give the sims_config.json in it"
                             + (f" ({', '.join(configs)})" if configs else ""))
        if not path.is_file():
            relative = not Path(normalise_sep(sims_config)).is_absolute()
            raise StudyError(
                f"sims config not found: {path}"
                + (f" - a relative path is read from the study file's folder, "
                   f"{self.folder}" if relative else ""))
        entry = {"name": name, "sims_config": portable(self.folder, str(path))}
        self.runs.append(entry)
        return entry

    def rename_run(self, old: str, new: str) -> None:
        """Rename a run and every table that reads it, so none is orphaned."""
        new = str(new).strip()
        entry = self.run_entry(old)
        if entry is None or not new or new == old:
            return
        if self.run_entry(new) is not None:
            raise StudyError(f"there is already a run called {new!r}")
        entry["name"] = new
        for source in self.all_sources():
            if source.get("run") == old:
                source["run"] = new

    def remove_run(self, name: str) -> list:
        """Drop a run. Returns what still names it - tables and PMF groups."""
        self.runs = [run for run in self.runs if run["name"] != name]
        left = [table.get("title") or table.get("id")
                for table in self.tables
                if any(source.get("run") == name for source in table_sources(table))]
        left += [f"PMF {entry.get('label') or ''}".strip()
                 for entry in pmf_entries(self)
                 if any((entry.get(key) or {}).get("run") == name
                        for key in ("ensemble", "mc"))]
        return left

    def all_sources(self) -> list:
        """Every {run, group} the study reads - the tables' and the PMF page's."""
        found = [source for table in self.tables for source in table_sources(table)]
        for entry in pmf_entries(self):
            found += [entry[key] for key in ("ensemble", "mc")
                      if isinstance(entry.get(key), dict)]
        for spec in self.extra.get("figures") or []:        # core/figures.py owns these
            found += [curve for curve in (spec.get("curves") or [])
                      if isinstance(curve, dict) and "run" in curve]
        return found

    # -- tables ------------------------------------------------------------

    def table(self, table_id: str) -> dict | None:
        return next((table for table in self.tables if table.get("id") == table_id), None)

    def unique_id(self, title: str) -> str:
        base = slug(title)
        taken = {table.get("id") for table in self.tables}
        candidate, n = base, 2
        while candidate in taken:
            candidate, n = f"{base}-{n}", n + 1
        return candidate

    def put_table(self, spec: dict) -> dict:
        """Add a table, or replace the one with the same id."""
        spec = copy.deepcopy(spec)
        if not spec.get("id"):
            spec["id"] = self.unique_id(spec.get("title") or spec.get("kind", "table"))
        for position, table in enumerate(self.tables):
            if table.get("id") == spec["id"]:
                self.tables[position] = spec
                return spec
        self.tables.append(spec)
        return spec

    def remove_table(self, table_id: str) -> None:
        self.tables = [table for table in self.tables if table.get("id") != table_id]

    def move_table(self, table_id: str, step: int) -> None:
        ids = [table.get("id") for table in self.tables]
        if table_id not in ids:
            return
        here = ids.index(table_id)
        there = max(0, min(len(ids) - 1, here + step))
        self.tables.insert(there, self.tables.pop(here))

    # -- disk --------------------------------------------------------------

    def to_json(self) -> dict:
        data = dict(self.extra)
        data.update({FORMAT_KEY: FORMAT_VERSION, "name": self.name,
                     "runs": [dict(run) for run in self.runs],
                     "tables": copy.deepcopy(self.tables)})
        return data

    def save(self) -> Path:
        atomic_write_json(self.path, self.to_json())
        return self.path


def table_sources(table: dict) -> list:
    """Every {run, group} a table reads, wherever in the spec it sits."""
    found = []
    for key in ("source", "pmf"):                    # a design table and its PMF row
        if isinstance(table.get(key), dict):
            found.append(table[key])
    for section in table.get("sections") or []:
        if "run" in section:                         # a representative-events section
            found.append(section)
        if isinstance(section.get("pmf"), dict):
            found.append(section["pmf"])
        for row in (section.get("rows") or []) + (section.get("groups") or []):
            if isinstance(row, dict) and "run" in row:      # rows, or a grid's columns
                found.append(row)
    return found


# -- the storm durations a critical duration is taken over -------------------------
#
# A table's or figure's design curve is the envelope over its group's storm
# durations, and the critical duration is whichever owns it. A study often reports
# over fewer durations than it ran - Callide's tables were made from 6-96 h, and a
# 120 h run governs the 1 in 2 level where it is allowed to - so the study keeps a
# **default** set (``"durations"`` at the top of the file) and each table and figure
# may override it with its own ``"durations"``:
#
#   absent, None or []   follow the study's default
#   "all"                every duration, whatever the default says
#   [6, 12, 24, ...]     those durations (hours)
#
# A default that is empty is every duration, so a study that never sets one reads
# exactly as before. An empty list on a table follows the default rather than
# meaning "all", because that is what the Frequent levels table wrote for "blank,
# all" before there was a default to follow.
#
# The PMF page and the ensemble peak table are not read through this: an ensemble's
# durations are the PMF's own, with their own rules (core/ensemble.py).

DURATIONS_KEY = "durations"
ALL_DURATIONS = "all"


def parse_hours(value) -> list:
    """Hours as a sorted list of distinct floats; anything not a number is dropped."""
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    found = set()
    for item in value or []:
        try:
            hours = float(item)
        except (TypeError, ValueError):
            continue
        if hours == hours and hours > 0:
            found.add(hours)
    return sorted(found)


def default_durations(study: Study) -> list:
    """The study's default durations (hours); empty is every duration."""
    return parse_hours(study.extra.get(DURATIONS_KEY))


def set_default_durations(study: Study, hours, *, today=None) -> None:
    """Store the default, and the date it last changed when it does change.

    The date is kept when the default is cleared: "every duration since 6 November"
    is as much a decision as leaving one out.
    """
    hours = parse_hours(hours)
    if hours == default_durations(study):
        return
    if hours:
        study.extra[DURATIONS_KEY] = hours
    else:
        study.extra.pop(DURATIONS_KEY, None)
    study.extra[DURATIONS_CHANGED_KEY] = (today or date.today()).isoformat()


# Why the default leaves what it leaves out, and when it last changed - a decision
# the report rests on, kept with it. A table or figure that overrides the default
# keeps its own reason under the same note key.
DURATIONS_NOTE_KEY = "durations_note"
DURATIONS_CHANGED_KEY = "durations_changed"


def durations_note(study: Study) -> str:
    return str(study.extra.get(DURATIONS_NOTE_KEY) or "").strip()


def set_durations_note(study: Study, text) -> None:
    text = str(text or "").strip()
    if text:
        study.extra[DURATIONS_NOTE_KEY] = text
    else:
        study.extra.pop(DURATIONS_NOTE_KEY, None)


def durations_changed(study: Study) -> str:
    """The date (ISO) the default last changed; '' if it never has."""
    return str(study.extra.get(DURATIONS_CHANGED_KEY) or "")


def left_out(study: Study, hours=None) -> list:
    """The durations the study ran that ``hours`` (default: the study's default) drops.

    Empty when there is no limit - nothing is left out.
    """
    hours = default_durations(study) if hours is None else parse_hours(hours)
    if not hours:
        return []
    return [each for each in durations_in(study) if each not in set(hours)]


def note_for(study: Study, spec: dict) -> str:
    """The reason in effect for a table's or figure's durations, or ''.

    Its own note when it overrides the default; the study's when it follows a
    default that limits anything.
    """
    if follows_default(spec):
        return durations_note(study) if default_durations(study) else ""
    return str((spec or {}).get(DURATIONS_NOTE_KEY) or "").strip()


def follows_default(spec: dict) -> bool:
    value = (spec or {}).get(DURATIONS_KEY)
    return value is None or (not isinstance(value, str) and not parse_hours(value))


def durations_for(study: Study, spec: dict) -> list:
    """The durations (hours) a table or figure takes its critical duration over.

    Empty means every duration.
    """
    value = (spec or {}).get(DURATIONS_KEY)
    if isinstance(value, str) and value.strip().lower() == ALL_DURATIONS:
        return []
    if follows_default(spec):
        return default_durations(study)
    return parse_hours(value)


def durations_text(hours) -> str:
    """'6, 12, 24, 48 h', or 'all durations' when there is no limit."""
    hours = parse_hours(hours)
    if not hours:
        return "all durations"
    return ", ".join(f"{hour:g}" for hour in hours) + " h"


def durations_in(study: Study, sources=None) -> list:
    """The storm durations (hours) run by the groups ``sources`` names.

    None means every row of every run in the study - what the default is chosen
    from. A run that cannot be opened contributes nothing.
    """
    from .results import duration_of             # local: results is the heavier module

    found = set()
    if sources is None:
        wanted = {name: None for name in study.run_names()}
    else:
        wanted = {}
        for source in sources:
            if source.get("run") and source.get("group"):
                wanted.setdefault(source["run"], set()).add(source["group"])
    for run_name, groups in wanted.items():
        try:
            run = study.open_run(run_name)
        except StudyError:
            continue
        rows = (list(run.frame.index) if groups is None
                else [index for group in groups for index in run.rows_in(group)])
        for index in rows:
            row = run.frame.loc[index]
            output = str(row.get("Output file") or "").replace("\\", "/").rsplit("/", 1)[-1]
            duration = duration_of(row, output)
            if duration is not None:
                found.add(float(duration))
    return sorted(found)


def pmf_entries(study: Study) -> list:
    """The PMF page's groups, as stored (core/ensemble.py owns the section)."""
    section = study.extra.get("pmf")
    return [entry for entry in (section or {}).get("groups") or [] if isinstance(entry, dict)]


def new_study(path, name: str = "") -> Study:
    path = Path(clean_path_text(path)).expanduser()
    if path.suffix.lower() != ".json":
        path = path / DEFAULT_NAME
    if path.exists():
        raise StudyError(f"{path} already exists - open it instead")
    study = Study(path=path.resolve() if path.parent.exists() else path,
                  name=name or path.parent.name)
    study.save()
    return study


def is_study_file(path) -> bool:
    """Whether a JSON file is a study file rather than a sims config: it carries the
    study format key. Unreadable is False - the caller's own open reports why."""
    try:
        data = read_json(Path(path))
    except Exception:                                # noqa: BLE001 - unreadable is "no"
        return False
    return isinstance(data, dict) and FORMAT_KEY in data


def load_study(path) -> Study:
    path = Path(clean_path_text(path)).expanduser()
    if path.is_dir():
        path = path / DEFAULT_NAME
    if not path.is_file():
        raise StudyError(f"study file not found: {path}")
    data = read_json(path, default=None)
    if not isinstance(data, dict):
        raise StudyError(f"{path.name} is not a JSON object - a trailing comma after "
                         f"the last entry is the usual cause")
    if FORMAT_KEY not in data:
        raise StudyError(f"{path.name} is not a Bryan study file (no {FORMAT_KEY!r} key)"
                         f" - a sims_config.json is opened on the Project page")
    version = data.get(FORMAT_KEY)
    if not isinstance(version, int) or version > FORMAT_VERSION:
        raise StudyError(f"{path.name} is study format {version!r}; this Bryan reads "
                         f"up to {FORMAT_VERSION}")
    runs = []
    for entry in data.get("runs") or []:
        if isinstance(entry, dict) and entry.get("name") and entry.get("sims_config"):
            runs.append({"name": str(entry["name"]),
                         "sims_config": str(entry["sims_config"])})
    tables = [table for table in data.get("tables") or [] if isinstance(table, dict)]
    extra = {key: value for key, value in data.items()
             if key not in (FORMAT_KEY, "name", "runs", "tables")}
    return Study(path=path.resolve(), name=str(data.get("name") or ""),
                 runs=runs, tables=tables, extra=extra)
