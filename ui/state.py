"""The app-lifetime singleton: the open project, the selection, the runs.

Mirrors ``swe2d_ui/state.py``. Everything durable lives on disk - the sims list
is the source of truth for what exists, and ``_ui_runs/<id>/run.json`` for what
is running - so this holds only what a session needs to remember.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field

import pandas as pd

from core import completion, grouping, preflight, runplan, runstate, runwriter
from core.config import SimsConfig, config_issues, load_sims_config
from core.launcher import RunManager
from core.outputs import log_path_for
from core.paths import cell_text
from core.simslist import SimsList, read_sims_list
from settings import UiSettings


@dataclass
class Project:
    """An open sims_config.json and the sims list it names."""

    config: SimsConfig
    sims: SimsList
    issues: list = field(default_factory=list)

    @classmethod
    def open(cls, path) -> "Project":
        config = load_sims_config(path)
        sims = read_sims_list(config.sims_list_path)
        return cls(config=config, sims=sims, issues=config_issues(config))

    def reload(self) -> "Project":
        return Project.open(self.config.config_path)

    @property
    def name(self) -> str:
        return self.config.config_path.name

    @property
    def frame(self) -> pd.DataFrame:
        return self.sims.frame

    def group_keys(self) -> pd.Series:
        return grouping.add_group_keys(self.frame)

    def groups(self) -> list:
        return grouping.groups_in_order(self.frame)

    def group_report(self):
        return grouping.derived_group_report(self.frame)

    def label(self, index) -> str:
        return cell_text(self.frame.loc[index].get("Output file")) or f"row {index + 2}"

    def log_paths(self) -> dict:
        return {index: log_path_for(self.frame.loc[index], self.config.project_folder)
                for index in self.frame.index}


class AppState:
    def __init__(self) -> None:
        self.settings = UiSettings.load()
        self.project: Project | None = None
        self.selected: set = set()
        self.completions: dict = {}
        self.manager = RunManager(
            bryan_python=self.settings.bryan_python,
            bryan_main=self.settings.bryan_main,
            show_console=self.settings.show_consoles,
        )
        self.active_run_id: str | None = None
        self._lock = threading.Lock()

    # -- project ----------------------------------------------------------

    def open_project(self, path) -> Project:
        project = Project.open(path)
        with self._lock:
            self.project = project
            self.selected = {
                index for index in project.frame.index
                if project.frame.loc[index].get("Include") == "yes"
            }
            self.completions = {}
        self.settings.remember(project.config.config_path)
        self.manager.reattach(project.config.run_root)
        return project

    def reload_project(self) -> Project | None:
        if self.project is None:
            return None
        keep = {self.project.label(index) for index in self.selected}
        project = self.project.reload()
        with self._lock:
            self.project = project
            self.selected = {index for index in project.frame.index
                             if project.label(index) in keep}
            self.completions = {}
        return project

    def apply_settings(self) -> None:
        self.settings.save()
        self.manager.bryan_python = self.settings.bryan_python
        self.manager.bryan_main = self.settings.bryan_main
        self.manager.show_console = self.settings.show_consoles

    # -- selection --------------------------------------------------------

    def set_selected(self, rows) -> None:
        with self._lock:
            self.selected = set(rows)

    def toggle(self, index, on: bool) -> None:
        with self._lock:
            if on:
                self.selected.add(index)
            else:
                self.selected.discard(index)

    def selected_in_order(self) -> list:
        if self.project is None:
            return []
        return [index for index in self.project.frame.index if index in self.selected]

    # -- completion -------------------------------------------------------

    def refresh_completion(self, *, check_truncation=False) -> dict:
        if self.project is None:
            return {}
        self.completions = completion.assess_frame(
            self.project.frame, self.project.config,
            check_truncation=check_truncation)
        return self.completions

    def completion_of(self, index):
        return self.completions.get(index)

    def rows_needing_rerun(self) -> list:
        return [index for index, state in self.completions.items()
                if state.state in completion.RERUN_WORTHY]

    def already_done(self, rows) -> list:
        return [index for index in rows
                if (self.completions.get(index) is not None
                    and self.completions[index].needs_confirm_to_overwrite)]

    # -- runs -------------------------------------------------------------

    def plan(self, rows=None, n_chunks=None):
        rows = self.selected_in_order() if rows is None else list(rows)
        return runplan.plan_run(
            self.project.sims, rows,
            n_chunks=n_chunks or self.settings.max_parallel,
            keep_groups_together=self.settings.keep_groups_together,
            history=runplan.observed_minutes(self.project.config),
        )

    def preflight(self, rows=None) -> list:
        rows = self.selected_in_order() if rows is None else list(rows)
        return preflight.check(self.project.sims, self.project.config, rows,
                               completions=self.completions or None)

    def launch(self, plan, *, test_runs=None):
        project = self.project
        folder = runwriter.write_run(project.sims, project.config, plan,
                                     test_runs=test_runs)
        names = {index: project.label(index) for index in project.frame.index}
        record = runstate.from_run_folder(
            folder, project.config,
            max_parallel=self.settings.max_parallel,
            names=names, log_paths=project.log_paths())
        self.manager.submit(record)
        self.active_run_id = record.run_id
        return record, folder

    def runs(self) -> list:
        if self.project is None:
            return []
        self.manager.poll()
        return sorted(self.manager.runs.values(),
                      key=lambda record: record.created, reverse=True)

    def active_run(self):
        if self.active_run_id:
            return self.manager.runs.get(self.active_run_id)
        live = self.manager.live_runs()
        return live[0] if live else None


STATE = AppState()
