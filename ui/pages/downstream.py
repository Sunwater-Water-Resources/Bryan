"""Generating the regional model's storm files from a chosen event set.

The Events page picks the events; this page turns a saved selection into storm
files for the downstream regional model. It shows the plan first - which
realisation, at which duration and warming level, into which filename - because
a duration read wrongly from a database name produces a perfectly valid storm of
the wrong length, and that is far cheaper to catch here than in a hydrograph.

**The generator runs as a subprocess.** Assembling the depths needs scipy and
this environment has none; see ``core/downstream.py``.
"""
from __future__ import annotations

from pathlib import Path

from nicegui import ui

from core import downstream
from layout import page_frame, require_project
from state import STATE


def downstream_page() -> None:
    with page_frame("Downstream storms"):
        project = require_project()
        if project is None:
            return
        folder = project.config.project_folder

        selections = downstream.find_selections(folder)
        if not selections:
            ui.label("No saved event selections under this project.").classes("text-lg")
            ui.markdown(
                "Pick events on the **Events** page and save them; this page reads "
                "the `<group>_representative_events.json` it writes.")
            return

        # The two config paths are remembered against this project (see
        # UiSettings.downstream_for), so they survive leaving the page and reloading the
        # project. The duration and warming level are deliberately not: blank means "read
        # it from the database name", and a stale remembered override is worse than
        # retyping one.
        remembered = STATE.settings.downstream_for(project.config.config_path)
        state = {"selection": selections[0], "duration": None, "gwl": None,
                 "config": remembered["config"], "model": remembered["model"],
                 "process": None}

        def remember() -> None:
            STATE.settings.remember_downstream(project.config.config_path,
                                               state["config"], state["model"])
        plan_area = ui.column().classes("w-full gap-2")

        def refresh() -> None:
            plan_area.clear()
            plan = downstream.plan(state["selection"], state["duration"], state["gwl"])
            with plan_area:
                if not plan.storms:
                    ui.label("No loading in this selection has an event picked.")
                    return
                rows = [{"loading": s.loading, "type": s.result_type,
                         "realisation": s.realisation,
                         "duration": "" if s.duration is None else f"{s.duration:g} h",
                         "gwl": "" if s.gwl is None else f"{s.gwl:g}",
                         "file": s.filename or "-", "problem": s.problem}
                        for s in plan.storms]
                ui.table(columns=[{"name": k, "label": k.title(), "field": k}
                                  for k in rows[0]], rows=rows).classes("w-full")
                ready, problems = len(plan.ready), len(plan.problems)
                if problems:
                    ui.label(f"{problems} of {len(plan.storms)} cannot be written yet - "
                             f"give the duration or warming level below.").classes("text-orange-700")
                ui.label(f"{ready} storm file{'s' if ready != 1 else ''} would be written.")

        with ui.row().classes("w-full items-end gap-4"):
            ui.select(selections, value=selections[0], label="Event selection",
                      on_change=lambda e: (state.update(selection=e.value), refresh())
                      ).classes("grow")
            ui.number(label="Duration (h)", value=None, format="%g",
                      on_change=lambda e: (state.update(duration=e.value), refresh())
                      ).props("clearable").classes("w-32")
            ui.number(label="GWL (°C)", value=None, format="%g",
                      on_change=lambda e: (state.update(gwl=e.value), refresh())
                      ).props("clearable").classes("w-32")

        with ui.row().classes("w-full items-end gap-4"):
            ui.input(label="Downstream storm config", value=state["config"],
                     on_change=lambda e: (state.update(config=e.value), remember())
                     ).classes("grow")
            ui.input(label="Regional model URBS config", value=state["model"],
                     on_change=lambda e: (state.update(model=e.value), remember())
                     ).classes("grow")

        def generate(dry_run: bool) -> None:
            if not state["config"] or not state["model"]:
                ui.notify("Name the downstream storm config and the model config first.",
                          type="warning")
                return
            argv = downstream.command(
                state["selection"], state["config"], state["model"],
                bryan_python=STATE.settings.bryan_python or None,
                duration=state["duration"], gwl=state["gwl"], dry_run=dry_run)
            log = Path(state["selection"]).with_suffix(".downstream.log")
            state["process"] = downstream.launch(argv, cwd=folder, log_path=log)
            ui.notify(f"{'Checking' if dry_run else 'Generating'} - output in {log.name}")

        with ui.row().classes("gap-3"):
            ui.button("Check without writing", on_click=lambda: generate(True)).props("outline")
            ui.button("Generate storm files", on_click=lambda: generate(False))

        refresh()
