"""Opening a project, and telling Bryan where it lives."""

from __future__ import annotations

from pathlib import Path

from nicegui import ui

from core.config import ConfigError
from layout import page_frame, severity_banner
from state import STATE


def project_page() -> None:
    with page_frame("Project"):
        _open_card()
        _bryan_card()
        if STATE.project is not None:
            _summary_card()


def _open_card() -> None:
    with ui.card().classes("w-full"):
        ui.label("Open a project").classes("text-lg font-bold")
        ui.label("Point at the sims_config.json a batch file would pass to "
                 "Main.py.").classes("text-sm text-gray-600")

        path_input = ui.input("sims_config.json",
                              value=str(STATE.project.config.config_path)
                              if STATE.project else "").classes("w-full")

        def do_open(path=None) -> None:
            target = path or path_input.value
            if not target:
                ui.notify("Give a path to a sims_config.json", type="warning")
                return
            try:
                project = STATE.open_project(target)
            except (ConfigError, FileNotFoundError, OSError) as exc:
                ui.notify(str(exc), type="negative", timeout=0, close_button=True)
                return
            path_input.value = str(project.config.config_path)
            ui.notify(f"Opened {project.name} - {len(project.frame)} simulations",
                      type="positive")
            ui.navigate.to("/select")

        with ui.row().classes("items-center gap-2"):
            ui.button("Open", on_click=lambda: do_open()).props("color=primary")
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


def _bryan_card() -> None:
    settings = STATE.settings
    with ui.card().classes("w-full"):
        ui.label("Bryan").classes("text-lg font-bold")
        ui.label("The UI shells out to Bryan rather than importing it, so the "
                 "environment that reproduces study results is left alone. "
                 "These are the two values the per-model .bat files set as "
                 "VENV_PY and PYFILE.").classes("text-sm text-gray-600")

        python_input = ui.input("Python interpreter",
                                value=settings.bryan_python).classes("w-full")
        main_input = ui.input("Main.py", value=settings.bryan_main).classes("w-full")

        with ui.row().classes("items-center gap-4 flex-wrap"):
            parallel = ui.number("Run at once", value=settings.max_parallel,
                                 min=1, max=16, step=1).classes("w-32")
            poll = ui.number("Refresh (s)", value=settings.poll_seconds,
                             min=0.5, max=30, step=0.5).classes("w-32")
            groups = ui.switch("Keep groups together",
                               value=settings.keep_groups_together)
            consoles = ui.switch("Show Bryan console windows",
                                 value=settings.show_consoles)

        ui.label("Running several at once rarely helps: reservoir routing "
                 "takes seconds, and a Monte Carlo simulation already drives "
                 "thousands of model runs one after another. More processes "
                 "mostly multiply memory and storm files. The per-chunk time "
                 "estimate on the Run page is the thing to judge it by."
                 ).classes("text-xs text-gray-500")

        def save() -> None:
            settings.bryan_python = python_input.value.strip()
            settings.bryan_main = main_input.value.strip()
            settings.max_parallel = int(parallel.value or 1)
            settings.poll_seconds = float(poll.value or 2.0)
            settings.keep_groups_together = bool(groups.value)
            settings.show_consoles = bool(consoles.value)
            STATE.apply_settings()
            problems = settings.problems()
            if problems:
                ui.notify("; ".join(problems), type="warning", timeout=0,
                          close_button=True)
            else:
                ui.notify("Saved", type="positive")

        ui.button("Save", on_click=save).props("color=primary")

        for problem in settings.problems():
            severity_banner("warn", problem)


def _summary_card() -> None:
    project = STATE.project
    sims = project.sims
    with ui.card().classes("w-full"):
        ui.label("This project").classes("text-lg font-bold")
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
                     ).classes("text-xs text-gray-500")

        for issue in project.issues:
            severity_banner("warn", issue)

        if sims.audit.is_damaged:
            severity_banner(
                "block",
                f"{sims.path.name}: {sims.audit.describe()}",
                "Until it is fixed, the affected rows cannot be selected: they "
                "would run with no output name and no input files.",
            )
