"""The settings panel, behind the gear in the menu bar: this computer's settings.

What is kept here is the user's, not the study's or the run's - where Bryan and
its interpreter are, how runs are launched, where the AWAP grids are - saved to
``~/.bryan_ui.json``. **Check setup** asks Bryan's interpreter what it can
import, looks for Main.py and the model executable, and writes the answer to a
file that can be sent on (core/setupcheck.py).
"""

from __future__ import annotations

from nicegui import run, ui

import settings as settings_module
from core import setupcheck
from core.paths import clean_path_text
from state import STATE

STATUS_LOOK = {setupcheck.OK: ("check_circle", "text-positive"),
               setupcheck.WARN: ("warning", "text-warning"),
               setupcheck.FAIL: ("cancel", "text-negative")}


def open_settings() -> None:
    settings = STATE.settings
    with ui.dialog() as dialog, \
            ui.card().classes("w-[44rem] max-w-full").mark("settings-panel"):
        with ui.row().classes("w-full items-center justify-between no-wrap"):
            ui.label("Settings").classes("text-lg font-bold")
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        ui.label(f"This computer's, kept in {settings_module.SETTINGS_PATH}. The interpreter "
                 f"and Main.py are the pair the model's batch files set as VENV_PY and "
                 f"PYFILE: the launcher runs Bryan with them rather than importing it, so "
                 f"the environment that reproduces study results is left alone."
                 ).classes("text-sm text-body")

        python_input = ui.input("Bryan's Python interpreter", value=settings.bryan_python) \
            .classes("w-full").props("dense").mark("settings-python")
        main_input = ui.input("Main.py", value=settings.bryan_main) \
            .classes("w-full").props("dense").mark("settings-main")
        awap_input = ui.input("Folder of daily AWAP / AWRA-L grids (optional)",
                              value=settings.awap_folder) \
            .classes("w-full").props("dense").mark("settings-awap")

        with ui.row().classes("items-center gap-4 flex-wrap"):
            parallel = ui.number("Run at once", value=settings.max_parallel,
                                 min=1, max=16, step=1).classes("w-28").props("dense")
            poll = ui.number("Refresh (s)", value=settings.poll_seconds,
                             min=0.5, max=30, step=0.5).classes("w-28").props("dense")
            groups = ui.switch("Keep groups together", value=settings.keep_groups_together)
            consoles = ui.switch("Show Bryan console windows", value=settings.show_consoles)
        ui.label("Running several at once rarely helps: reservoir routing takes seconds, "
                 "and a Monte Carlo simulation already drives thousands of model runs one "
                 "after another. More processes mostly multiply memory and storm files. "
                 "The per-chunk time estimate on the Run page is the thing to judge it by."
                 ).classes("text-xs text-muted")

        problems_box = ui.column().classes("w-full gap-1")

        def show_problems() -> None:
            problems_box.clear()
            with problems_box:
                for problem in settings.problems():
                    with ui.row().classes("items-center gap-2 no-wrap"):
                        ui.icon("warning").classes("text-warning")
                        ui.label(problem).classes("text-sm")

        def save() -> None:
            settings.bryan_python = clean_path_text(python_input.value)
            settings.bryan_main = clean_path_text(main_input.value)
            settings.awap_folder = clean_path_text(awap_input.value)
            settings.max_parallel = int(parallel.value or 1)
            settings.poll_seconds = float(poll.value or 2.0)
            settings.keep_groups_together = bool(groups.value)
            settings.show_consoles = bool(consoles.value)
            STATE.apply_settings()
            show_problems()
            if settings.problems():
                ui.notify("Saved, with problems - see the panel", type="warning")
            else:
                ui.notify("Saved", type="positive")

        results = ui.column().classes("w-full gap-1").mark("setup-results")

        async def check() -> None:
            results.clear()
            with results:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Asking Bryan's interpreter what it can import...") \
                        .classes("text-sm text-body")
            project = STATE.project
            checks = await run.io_bound(
                setupcheck.run_checks,
                bryan_python=clean_path_text(python_input.value),
                bryan_main=clean_path_text(main_input.value),
                awap_folder=clean_path_text(awap_input.value),
                sims_config=project.config if project is not None else None)
            try:
                report = setupcheck.write_report(checks, settings_module.SETTINGS_PATH.parent)
            except OSError as exc:
                report, report_problem = None, str(exc)
            results.clear()
            with results:
                _draw_checks(checks)
                if report is not None:
                    ui.label(f"Written to {report} - send it on when asking for help."
                             ).classes("text-xs text-muted").mark("setup-report")
                else:
                    ui.label(f"The report could not be written: {report_problem}"
                             ).classes("text-xs text-negative")

        with ui.row().classes("w-full items-center gap-2"):
            ui.button("Save", on_click=save).props("color=primary").mark("save-settings")
            ui.button("Check setup", icon="fact_check", on_click=check) \
                .props("outline").mark("check-setup") \
                .tooltip("Uses the values above, saved or not")
        show_problems()
    dialog.open()


def _draw_checks(checks) -> None:
    verdict = {setupcheck.OK: "Ready to run.",
               setupcheck.WARN: "Ready to run, with warnings.",
               setupcheck.FAIL: "Not ready - see the crosses."}[setupcheck.worst(checks)]
    ui.label(verdict).classes("font-bold").mark("setup-verdict")
    area = None
    with ui.element("div").classes("w-full grid gap-x-3 gap-y-1 items-start") \
            .style("grid-template-columns: 1.25rem 8rem 1fr"):
        for item in checks:
            if item.area != area:
                area = item.area
                ui.label(area).classes("text-xs font-bold text-muted uppercase mt-2") \
                    .style("grid-column: 1 / -1")
            icon, colour = STATUS_LOOK[item.status]
            ui.icon(icon).classes(f"{colour} text-lg")
            ui.label(item.name).classes("text-sm")
            ui.label(item.detail).classes("text-sm text-body break-all") \
                .mark(f"setup-{item.area}-{item.name}")
