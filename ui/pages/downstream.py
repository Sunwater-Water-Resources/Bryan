"""Generating the regional model's storm files from a chosen event set.

The Events page picks the events; this page turns a saved selection into storm
files for the downstream regional model. It shows the plan first - which
realisation, at which duration and warming level, into which filename - because
a duration read wrongly from a database name produces a perfectly valid storm of
the wrong length, and that is far cheaper to catch here than in a hydrograph.

It then carries the event set the rest of the way to a hydraulic model: URBS can be
run on each storm as it is written, and the dam's column in the inflow files URBS
writes can be replaced by the outflow of a dam failure model.

**The generator runs as a subprocess.** Assembling the depths needs scipy and
this environment has none; see ``core/downstream.py``.
"""
from __future__ import annotations

from pathlib import Path

from nicegui import ui

from core import downstream
from layout import page_frame, require_project
from state import STATE
from widgets import OUTPUT_FOLDER, path_input


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
        boundaries = STATE.settings.boundaries_for(project.config.config_path)
        state = {"selection": selections[0], "duration": None, "gwl": None,
                 "config": remembered["config"], "model": remembered["model"],
                 # The suffix follows the selection until it is typed in: two groups that
                 # share a realisation need one each, and which one depends on the group.
                 "suffix": downstream.suggested_suffix(selections[0], selections),
                 "suffix_typed": False, "run_urbs": False,
                 "dam_outflow": boundaries["dam_outflow"],
                 "boundaries_folder": boundaries["boundaries_folder"],
                 "name_contains": boundaries["name_contains"], "gate_scenario": "",
                 "process": None}

        def remember() -> None:
            STATE.settings.remember_downstream(project.config.config_path,
                                               state["config"], state["model"])

        def remember_boundaries() -> None:
            STATE.settings.remember_boundaries(
                project.config.config_path, dam_outflow=state["dam_outflow"],
                boundaries_folder=state["boundaries_folder"],
                name_contains=state["name_contains"])

        plan_area = ui.column().classes("w-full gap-2")

        def refresh() -> None:
            plan_area.clear()
            plan = downstream.plan(state["selection"], state["duration"], state["gwl"],
                                   suffix=state["suffix"])
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
                             f"give the duration or warming level below.").classes("text-attention")
                ui.label(f"{ready} storm file{'s' if ready != 1 else ''} would be written.")
                clashing = downstream.clashes(state["selection"], selections, state["suffix"])
                if clashing:
                    names = ", ".join(downstream.group_of(c) or Path(c).name for c in clashing)
                    ui.label(f"These file names are also what {names} would write, so one "
                             f"set replaces the other - storms and URBS results. Give this "
                             f"selection a suffix.").classes("text-attention")

        def choose(selection) -> None:
            state["selection"] = selection
            if not state["suffix_typed"]:
                state["suffix"] = downstream.suggested_suffix(selection, selections)
                suffix_box.value = state["suffix"]
            refresh()

        def typed_suffix(text) -> None:
            state["suffix"] = (text or "").strip()
            refresh()

        with ui.row().classes("w-full items-end gap-4"):
            ui.select(selections, value=selections[0], label="Event selection",
                      on_change=lambda e: choose(e.value)).classes("grow")
            ui.number(label="Duration (h)", value=None, format="%g",
                      on_change=lambda e: (state.update(duration=e.value), refresh())
                      ).props("clearable").classes("w-32")
            ui.number(label="GWL (°C)", value=None, format="%g",
                      on_change=lambda e: (state.update(gwl=e.value), refresh())
                      ).props("clearable").classes("w-32")
            suffix_box = ui.input(label="File name suffix", value=state["suffix"],
                                  on_change=lambda e: typed_suffix(e.value)).classes("w-40")
            # typing is what stops the suffix following the selection, not the page setting it
            suffix_box.on("keydown", lambda: state.update(suffix_typed=True))

        with ui.row().classes("w-full items-end gap-4"):
            # Bryan runs with the project folder as its working folder, so a
            # relative path is read from there.
            for label, key in (("Downstream storm config", "config"),
                               ("Regional model URBS config", "model")):
                path_input(label, state[key], base=folder, base_name="the project folder",
                           suffixes=(".json",), classes="grow",
                           on_change=lambda text, k=key: (state.update({k: text}),
                                                          remember()))

        def generate(dry_run: bool) -> None:
            if not state["config"] or not state["model"]:
                ui.notify("Name the downstream storm config and the model config first.",
                          type="warning")
                return
            argv = downstream.command(
                state["selection"], state["config"], state["model"],
                bryan_python=STATE.settings.bryan_python or None,
                duration=state["duration"], gwl=state["gwl"], dry_run=dry_run,
                suffix=state["suffix"], run_urbs=state["run_urbs"],
                sims_config=project.config.config_path)
            log = Path(state["selection"]).with_suffix(".downstream.log")
            state["process"] = downstream.launch(argv, cwd=folder, log_path=log)
            doing = "Checking" if dry_run else (
                "Generating and running URBS" if state["run_urbs"] else "Generating")
            ui.notify(f"{doing} - output in {log.name}")

        with ui.row().classes("items-center gap-3"):
            ui.button("Check without writing", on_click=lambda: generate(True)).props("outline")
            ui.button("Generate storm files", on_click=lambda: generate(False))
            ui.switch("Run URBS on each storm",
                      on_change=lambda e: state.update(run_urbs=bool(e.value)))
        ui.label("Without it a batch file is written beside each storm, to run by hand. "
                 "With it the results land in the regional model's folder, and the inflow "
                 "file for a hydraulic model is <storm name>_i.csv.").classes("text-sm")

        # -- the dam's outflow, from a dam failure model -------------------------------
        ui.separator()
        ui.label("Dam outflow from a dam failure model").classes("text-lg")
        ui.markdown(
            "URBS's inflow file carries the regional model's own routing of the dam. Where "
            "the dam is routed in a dam failure model instead, name that model's simulation "
            "list here. A copy of each event's inflow file is written with the dam's column "
            "replaced. The dam runs have to reach the end of the regional model's run.")

        with ui.row().classes("w-full items-end gap-4"):
            path_input("Dam failure simulation list", state["dam_outflow"], base=folder,
                       base_name="the project folder", suffixes=(".xlsx",), classes="grow",
                       on_change=lambda text: (state.update(dam_outflow=text),
                                               remember_boundaries()))
            path_input("Write the inflow files to", state["boundaries_folder"], base=folder,
                       base_name="the project folder", expect=OUTPUT_FOLDER, classes="grow",
                       on_change=lambda text: (state.update(boundaries_folder=text),
                                               remember_boundaries()))
        with ui.row().classes("w-full items-end gap-4"):
            ui.input(label="Gate scenario (blank: gates working)", value="",
                     on_change=lambda e: state.update(gate_scenario=(e.value or "").strip())
                     ).classes("w-72")
            ui.input(label="Only simulations whose name contains", value=state["name_contains"],
                     on_change=lambda e: (state.update(name_contains=(e.value or "").strip()),
                                          remember_boundaries())).classes("w-72")

        def write_boundaries() -> None:
            if not state["model"]:
                ui.notify("Name the regional model URBS config first.", type="warning")
                return
            if not state["dam_outflow"] or not state["boundaries_folder"]:
                ui.notify("Name the dam failure simulation list and the folder to write to.",
                          type="warning")
                return
            argv = downstream.boundaries_command(
                state["selection"], state["model"], state["dam_outflow"],
                state["boundaries_folder"],
                bryan_python=STATE.settings.bryan_python or None,
                duration=state["duration"], gwl=state["gwl"], suffix=state["suffix"],
                gate_scenario=state["gate_scenario"], name_contains=state["name_contains"])
            log = Path(state["selection"]).with_suffix(".boundaries.log")
            state["process"] = downstream.launch(argv, cwd=folder, log_path=log)
            ui.notify(f"Writing the inflow files - output in {log.name}")

        ui.button("Write hydraulic model inflow files", on_click=write_boundaries)

        refresh()
