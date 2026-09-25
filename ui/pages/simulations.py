"""Simulations: which sims list the run pages work on.

A sims_config.json and its simulation list are one **run** - E013 RFSL, say. The
Select, Run, Results, Ensemble, Events and Lake levels pages all work on the one
open here. With a study open, its runs are listed first, each opened with one
click; any other sims_config.json is opened below and can be added to the study.
"""

from __future__ import annotations

from pathlib import Path

from nicegui import ui

from core import study as studies
from core.config import ConfigError
from layout import page_frame, severity_banner
from state import STATE


def simulations_page() -> None:
    with page_frame("Simulations"):
        if STATE.study is not None:
            _study_runs_card()
        _open_card()
        if STATE.project is not None:
            if STATE.study is not None:
                _add_to_study_card()
            _summary_card()


def _open_run(path) -> None:
    """Open a sims list and go to choosing its simulations."""
    try:
        project = STATE.open_project(path)
    except (ConfigError, FileNotFoundError, OSError) as exc:
        ui.notify(str(exc), type="negative", timeout=0, close_button=True)
        return
    ui.notify(f"Opened {project.name} - {len(project.frame)} simulations",
              type="positive")
    ui.navigate.to("/select")


# -- the study's runs --------------------------------------------------------------

def _study_runs_card() -> None:
    study = STATE.study
    current = STATE.project.config.config_path if STATE.project is not None else None
    with ui.card().classes("w-full").mark("study-runs-card"):
        ui.label("The study's runs").classes("text-lg font-bold")
        ui.label(f"From {study.name or study.path.name}. Open one to choose and run its "
                 f"simulations.").classes("text-sm text-body")
        if not study.runs:
            ui.label("No runs yet. Open a sims_config.json below and add it to the study."
                     ).classes("text-sm text-muted")
            return
        with ui.element("div").classes("w-full grid gap-x-4 gap-y-1 items-center") \
                .style("grid-template-columns: 12rem 1fr 7rem"):
            for entry in study.runs:
                path = study.run_config_path(entry["name"])
                is_open = current is not None and path is not None \
                    and path.resolve() == Path(current).resolve()
                ui.label(entry["name"]).classes(
                    "text-sm " + ("font-bold text-ink" if is_open else "text-body"))
                ui.label(entry["sims_config"]).classes("mono text-xs text-muted")
                if is_open:
                    ui.label("open").classes("text-xs text-positive") \
                        .mark(f"run-open-{entry['name']}")
                elif path is None or not path.is_file():
                    ui.label("not found").classes("text-xs text-negative") \
                        .tooltip(str(path))
                else:
                    ui.button("Open", on_click=lambda _, p=path: _open_run(p)) \
                        .props("flat dense no-caps").mark(f"open-run-{entry['name']}")


def _add_to_study_card() -> None:
    """Offered when the open sims list is not one of the study's runs."""
    study, project = STATE.study, STATE.project
    config = project.config.config_path
    if study.run_named_for(config) is not None:
        return
    with ui.card().classes("w-full").mark("add-to-study-card"):
        ui.label("This sims list is not in the study").classes("font-bold")
        ui.label("Add it to name it as a run, for the report tables, PMF and figures to "
                 "read.").classes("text-sm text-body")
        with ui.row().classes("items-end gap-2"):
            name = ui.input("Run name", value=studies.guess_run_name(config)) \
                .classes("w-64").props("dense").mark("add-to-study-name")

            def add() -> None:
                try:
                    study.add_run(name.value, config)
                    study.save()
                except (studies.StudyError, OSError) as exc:
                    ui.notify(str(exc), type="negative", multi_line=True)
                    return
                ui.notify(f"Added {name.value} to the study", type="positive")
                ui.navigate.reload()

            ui.button("Add to study", icon="playlist_add", on_click=add) \
                .mark("add-to-study")


def _open_card() -> None:
    with ui.card().classes("w-full").mark("open-sims-card"):
        ui.label("Open another sims_config.json" if STATE.study is not None
                 else "Open a sims_config.json").classes("text-lg font-bold")
        ui.label("Point at the sims_config.json a batch file would pass to "
                 "Main.py.").classes("text-sm text-body")

        path_input = ui.input("sims_config.json",
                              value=str(STATE.project.config.config_path)
                              if STATE.project else "").classes("w-full") \
            .mark("sims-config-path")

        def do_open(path=None) -> None:
            target = path or path_input.value
            if not target:
                ui.notify("Give a path to a sims_config.json", type="warning")
                return
            _open_run(target)

        with ui.row().classes("items-center gap-2"):
            ui.button("Open", on_click=lambda: do_open()).props("color=primary").mark("open-sims")
            if STATE.project is not None:
                ui.button("Reload from disk",
                          on_click=lambda: (STATE.reload_project(),
                                            ui.notify("Reloaded"),
                                            ui.navigate.to("/select"))
                          ).props("flat")

        recent = [entry for entry in STATE.settings.recent_configs
                  if Path(entry).is_file()]
        if recent:
            ui.label("Recent").classes("text-sm font-bold mt-2")
            for entry in recent:
                ui.link(entry, "#").on(
                    "click", lambda _, target=entry: do_open(target)
                ).classes("text-sm")


def _summary_card() -> None:
    project = STATE.project
    sims = project.sims
    with ui.card().classes("w-full"):
        ui.label("This sims list").classes("text-lg font-bold")
        rows = [
            ("Project folder", str(project.config.project_folder)),
            ("Simulation list", f"{sims.path.name}  ({len(sims.frame)} rows, "
                                f"sheet {sims.sheet_name!r})"),
            ("Run log", project.config.master_run_log.name),
        ]
        for key, value in project.config.filepaths.items():
            rows.append((key, str(value)))
        with ui.grid(columns=2).classes("gap-x-6 gap-y-1 text-sm"):
            for key, value in rows:
                ui.label(key).classes("font-bold")
                ui.label(value)

        if sims.other_sheets:
            ui.label(f"Other sheets, not read by Bryan and not copied into run "
                     f"folders: {', '.join(sims.other_sheets)}"
                     ).classes("text-xs text-muted")

        for issue in project.issues:
            severity_banner("warn", issue)

        if sims.audit.is_damaged:
            severity_banner(
                "block",
                f"{sims.path.name}: {sims.audit.describe()}",
                "Until it is fixed, the affected rows cannot be selected: they "
                "would run with no output name and no input files.",
            )
