"""The sims_config.json chain, resolved exactly as Main.py resolves it.

Main.py uses TWO different base folders, which is easy to get wrong:

- ``project_folder`` = ``dirname(config)``, unless the config's own
  ``project_folder`` key says otherwise (Main.py:29-32). The simulation list is
  joined onto this (Main.py:35).
- ``filepaths`` (model/storm/climate config) resolve against
  ``dirname(simulation_file)`` - always the config's own folder, never the
  overridden project folder (Main.py:45-50).

And nothing in Main.py chdirs, so every path *inside* the sims list resolves
against the **process CWD**. The launcher .bat does ``cd /D "%~dp0"`` first,
which is why core/launcher.py runs Bryan with ``cwd=project_folder``.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path

FILEPATH_KEYS = ("model_config", "storm_config", "climate_config")


class ConfigError(Exception):
    """A sims_config.json that cannot be used."""


@dataclass(frozen=True)
class SimsConfig:
    """A loaded sims_config.json, with every path resolved as Bryan resolves it."""

    config_path: Path          # the ORIGINAL sims_config.json, absolute
    project_folder: Path       # what simulation_list is joined onto
    config_folder: Path        # what filepaths resolve against
    simulation_list_rel: str   # the raw json string, verbatim
    sims_list_path: Path
    filepaths: dict            # key -> absolute Path
    test_runs: int
    raw: dict = field(default_factory=dict)   # verbatim, for the run copies

    @property
    def run_root(self) -> Path:
        """Where the UI puts its run folders."""
        return self.project_folder / "_ui_runs"

    @property
    def master_run_log(self) -> Path:
        """The run log Bryan writes beside the master sims list.

        Mirrors lib/RunLog.py:29, which builds the path from the raw
        ``simulation_list`` string and resolves it against the CWD - which for
        a batch-file run is the project folder.
        """
        return run_log_path_for(self.simulation_list_rel, self.project_folder)


def run_log_path_for(simulation_list_value: str, cwd) -> Path:
    """The run log path for a given ``simulation_list`` value.

    Reproduces lib/RunLog.py:29 including its quirk - it is a substring
    ``.replace('.xlsx', '_log.csv')`` on the raw string, not a suffix swap. A
    value with '.xlsx' anywhere but the end would be mangled; core/paths.py
    constrains run ids so the UI can never generate one.
    """
    replaced = str(simulation_list_value).replace(".xlsx", "_log.csv")
    from .paths import normalise_sep  # local: keeps this module import-light
    local = Path(normalise_sep(replaced))
    return local if local.is_absolute() else Path(cwd) / local


def load_sims_config(path) -> SimsConfig:
    """Load a sims_config.json the way Main.py does, and say so when it can't.

    Main.py has no validation at all - a missing key is a bare KeyError deep in
    the run. Everything that would blow up there is checked here instead.
    """
    config_path = Path(path).expanduser()
    if not config_path.is_file():
        raise ConfigError(f"config file not found: {config_path}")
    config_path = config_path.resolve()

    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        # Manual/SubDocs/config/json_files.md exists mostly to explain this one.
        raise ConfigError(
            f"{config_path.name} is not valid JSON: {exc}. A trailing comma "
            f"after the last entry is the usual cause."
        ) from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{config_path.name} must hold a JSON object")

    if "simulation_list" not in raw:
        raise ConfigError(
            f"{config_path.name} has no 'simulation_list' key - Main.py:35 "
            f"needs it to find the sims list"
        )

    # Main.py:29-32
    config_folder = config_path.parent
    project_folder = config_folder
    override = raw.get("project_folder")
    if override and str(override).strip().lower() != "default":
        from .paths import normalise_sep
        candidate = Path(normalise_sep(str(override)))
        project_folder = candidate if candidate.is_absolute() else config_folder / candidate
    project_folder = project_folder.resolve()

    # Main.py:35
    simulation_list_rel = str(raw["simulation_list"])
    from .paths import normalise_sep
    sims_list_path = project_folder / normalise_sep(simulation_list_rel)

    # Main.py:44-50 - relative to the CONFIG folder, then normalised absolute.
    filepaths = {}
    for key, value in (raw.get("filepaths") or {}).items():
        local = Path(normalise_sep(str(value)))
        resolved = local if local.is_absolute() else config_folder / local
        filepaths[key] = Path(os.path.normpath(resolved))

    # Main.py:55-60. Note the manual calls this key 'test run'; the code reads
    # 'test_runs'. Accept the documented spelling so a config written from the
    # manual is not silently ignored.
    test_runs = raw.get("test_runs", raw.get("test run", 0))
    try:
        test_runs = int(test_runs)
    except (TypeError, ValueError):
        raise ConfigError(f"'test_runs' must be a whole number, got {test_runs!r}")

    return SimsConfig(
        config_path=config_path,
        project_folder=project_folder,
        config_folder=config_folder,
        simulation_list_rel=simulation_list_rel,
        sims_list_path=sims_list_path,
        filepaths=filepaths,
        test_runs=test_runs,
        raw=raw,
    )


def config_issues(config: SimsConfig) -> list[str]:
    """Problems worth showing before anything is launched. Empty means fine."""
    issues = []
    if not config.sims_list_path.is_file():
        issues.append(
            f"simulation list not found: {config.sims_list_path} "
            f"(from 'simulation_list': {config.simulation_list_rel!r})"
        )
    for key in FILEPATH_KEYS:
        if key not in config.filepaths:
            issues.append(f"'filepaths' has no {key!r} entry")
        elif not config.filepaths[key].is_file():
            issues.append(f"{key} not found: {config.filepaths[key]}")
    for key, value in config.filepaths.items():
        if key not in FILEPATH_KEYS and not value.exists():
            issues.append(f"{key} not found: {value}")
    return issues
