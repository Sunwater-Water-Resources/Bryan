"""Taking a chosen event set downstream: the regional model's storms, and what came of them.

The Events page picks the events; this page turns a saved selection into storm files for
the downstream regional model, runs URBS on them if asked, and reads the runs back.

**Storms** shows the plan first - which event, at which duration and warming level, into
which filename - because a duration read wrongly from a database name produces a perfectly
valid storm of the wrong length, and that is far cheaper to catch here than in a hydrograph.
**Results** reads the regional runs back: which have finished, and what went into each.

The selection and the regional model config sit above both tabs because both tabs are
about them. What the page shows follows the selection's method: an ensemble selection
records each event's duration, so it has no duration box to fill in.

Where the dam is routed in a dam failure model, that project puts the dam's outflow into
the inflow files URBS writes, after Bryan; it is not done here, so that one tool writes
the boundaries.

**The generator runs as a subprocess.** Assembling the depths needs scipy and this
environment has none; see ``core/downstream.py``.
"""
from __future__ import annotations

import os
from pathlib import Path

from nicegui import ui

import address
from clipboard import copy_buttons, table_from_rows
from core import damchart, downstream, logtail, urbsresults
from layout import page_frame, require_project, severity_banner
from state import STATE
from theme import house_echart
from widgets import Running, path_input

TABS = ("Storms", "Results")


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
                "Pick events on the **Events** page; this page reads the "
                "`<group>_representative_events.json` it saves.")
            return
        _DownstreamView(project, selections).build()


class _DownstreamView:
    def __init__(self, project, selections) -> None:
        self.project = project
        self.folder = project.config.project_folder
        self.selections = selections
        # The two config paths are remembered against this project (see
        # UiSettings.downstream_for), so they survive leaving the page and reloading the
        # project. The duration and warming level are deliberately not: blank means "read
        # it from the database name", and a stale remembered override is worse than
        # retyping one.
        remembered = STATE.settings.downstream_for(project.config.config_path)
        # The selection the address names - the page reopens on the one last looked at -
        # else the newest.
        named = address.param("selection")
        chosen = next((path for path in selections if downstream.group_of(path) == named),
                      selections[0])
        self.state = {"selection": chosen, "duration": None, "gwl": None,
                      "config": remembered["config"], "model": remembered["model"],
                      # The suffix is the one the selection's last run used, read off the
                      # names it wrote, or blank for one never run - never a guess. Where
                      # two groups share a realisation the plan says so and suggests one.
                      "suffix": self._recorded_suffix(chosen),
                      "run_urbs": False, "process": None,
                      # whether the last run this page started ran URBS: None until one has
                      "ran_urbs": None}
        self.results_built = False

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        state = self.state
        options = {path: downstream.selection_label(path) for path in self.selections}
        with ui.row().classes("w-full items-end gap-4"):
            ui.select(options, value=state["selection"], label="Event selection",
                      on_change=lambda e: self._choose(e.value)).classes("grow") \
                .mark("selection")
            # Bryan runs with the project folder as its working folder, so a relative
            # path is read from there.
            path_input("Regional model URBS config", state["model"], base=self.folder,
                       base_name="the project folder", suffixes=(".json",), classes="grow",
                       on_change=lambda text: self._set_path("model", text))
        self.method_label = ui.label("").classes("text-xs text-muted").mark("method-line")

        tab = address.param("tab", TABS[0])
        tab = tab if tab in TABS else TABS[0]
        with ui.tabs().classes("w-full").mark("downstream-tabs") as tabs:
            for name in TABS:
                ui.tab(name)
        with ui.tab_panels(tabs, value=tab).classes("w-full"):
            with ui.tab_panel("Storms").classes("p-0 gap-3"):
                self._storms_tab()
            with ui.tab_panel("Results").classes("p-0 gap-3"):
                self.results_box = ui.column().classes("w-full gap-3")

        def on_tab(event) -> None:
            address.keep(tab=event.value if event.value != TABS[0] else "")
            if event.value == "Results":
                self.refresh_results()

        tabs.on_value_change(on_tab)
        self._method_changed()
        self.refresh_plan()
        if tab == "Results":
            self.refresh_results()

    def _storms_tab(self) -> None:
        state = self.state
        # The inputs in a card of their own: under the tabs they read as the plan table's
        # header row and were easy to miss.
        with ui.card().classes("w-full gap-2").mark("storm-inputs"):
            ui.label("Storm inputs").classes("font-bold")
            with ui.row().classes("w-full items-end gap-4"):
                path_input("Downstream storm config", state["config"], base=self.folder,
                           base_name="the project folder", suffixes=(".json",),
                           classes="grow",
                           on_change=lambda text: self._set_path("config", text))
            with ui.row().classes("w-full items-end gap-4"):
                self.suffix_box = ui.input(
                    label="File name suffix", value=state["suffix"],
                    placeholder="e.g. _RFSL",
                    on_change=lambda e: self._typed_suffix(e.value)).classes("w-48") \
                    .props("outlined dense").mark("suffix")
                self.duration_box = ui.number(
                    label="Duration (h)", value=None, format="%g",
                    on_change=lambda e: (state.update(duration=e.value), self.refresh_plan())
                ).props("clearable outlined dense").classes("w-36").tooltip(
                    "Blank reads it from the database name.")
                ui.number(label="GWL (°C)", value=None, format="%g",
                          on_change=lambda e: (state.update(gwl=e.value), self.refresh_plan())
                          ).props("clearable outlined dense").classes("w-36") \
                    .tooltip("Blank reads it from the database name.")
            ui.label("The suffix is added to every storm and run name. Duration and GWL are "
                     "read from the database name when left blank."
                     ).classes("text-xs text-muted")
        self.plan_area = ui.column().classes("w-full gap-2")
        with ui.row().classes("items-center gap-3"):
            ui.button("Check without writing",
                      on_click=lambda: self._generate(True)).props("outline")
            ui.button("Generate storm files", on_click=lambda: self._generate(False))
            ui.switch("Run URBS on each storm",
                      on_change=lambda e: state.update(run_urbs=bool(e.value)))
        self.spinner = ui.row().classes("items-center")
        self.status = ui.label("").classes("text-sm text-muted").mark("downstream-status")
        ui.label("Without Run URBS a batch file is written beside each storm, to run by hand. "
                 "Results land in the results folder named in the regional model's config. "
                 "Keep a group's suffix the same from run to run: a dam failure project "
                 "finds the runs by name.").classes("text-xs text-muted")

    # -- the Storms tab ----------------------------------------------------

    def refresh_plan(self) -> None:
        state = self.state
        self.plan_area.clear()
        plan = downstream.plan(state["selection"], state["duration"], state["gwl"],
                               suffix=state["suffix"])
        with self.plan_area:
            if not plan.storms:
                ui.label("No loading in this selection has an event picked.")
                return
            rows = [{"loading": s.loading, "type": s.result_type,
                     "event": s.event_label,
                     "duration": "" if s.duration is None else f"{s.duration:g} h",
                     "gwl": "" if s.gwl is None else f"{s.gwl:g}",
                     "file": s.filename or "-", "problem": s.problem}
                    for s in plan.storms]
            ui.table(columns=[{"name": k, "label": k.title(), "field": k}
                              for k in rows[0]], rows=rows).classes("w-full")
            ready, problems = len(plan.ready), len(plan.problems)
            if problems:
                ui.label(f"{problems} of {len(plan.storms)} cannot be written yet - "
                         f"give the duration or warming level above.").classes("text-attention")
            ui.label(f"{ready} storm file{'s' if ready != 1 else ''} would be written.")
            clashing = downstream.clashes(state["selection"], self.selections, state["suffix"])
            if clashing:
                names = ", ".join(downstream.group_of(c) or Path(c).name for c in clashing)
                suggested = downstream.suggested_suffix(state["selection"], self.selections)
                with ui.row().classes("items-center gap-3"):
                    ui.label(f"These file names are also what {names} would write, so one "
                             f"set replaces the other - storms and URBS results. Give this "
                             f"selection a suffix ({suggested} would keep them apart)."
                             ).classes("text-attention")
                    ui.button(f"Use {suggested}",
                              on_click=lambda: setattr(self.suffix_box, "value", suggested)) \
                        .props("outline dense no-caps").mark("use-suffix")

    @staticmethod
    def _recorded_suffix(selection) -> str:
        return urbsresults.recorded_suffix(downstream.plan(selection), selection)

    def _choose(self, selection) -> None:
        self.state["selection"] = selection
        address.keep(selection=downstream.group_of(selection))
        # the suffix this selection was last run with, so its results are found and a rerun
        # keeps its names; setting the box runs _typed_suffix, which redraws the plan
        self.state["suffix"] = self._recorded_suffix(selection)
        self.suffix_box.value = self.state["suffix"]
        self._method_changed()
        self.refresh_plan()
        if self.results_built:
            self.refresh_results()

    def _method_changed(self) -> None:
        state = self.state
        ensemble = downstream.method_of(state["selection"]) == downstream.EVENTS.ENSEMBLE
        # an ensemble selection records each event's own duration
        self.duration_box.set_visibility(not ensemble)
        if ensemble and state["duration"] is not None:
            state["duration"] = None
            self.duration_box.value = None
        self.method_label.text = downstream.method_line(state["selection"],
                                                        self._absolute(state["model"]))

    def _typed_suffix(self, text) -> None:
        self.state["suffix"] = (text or "").strip()
        self.refresh_plan()

    def _set_path(self, key, text) -> None:
        self.state[key] = text
        STATE.settings.remember_downstream(self.project.config.config_path,
                                           self.state["config"], self.state["model"])
        if key == "model":
            self._method_changed()
            if self.results_built:
                self.refresh_results()

    def _absolute(self, text) -> str | None:
        if not text:
            return None
        path = Path(text)
        return str(path if path.is_absolute() else Path(self.folder) / path)

    def _generate(self, dry_run: bool) -> None:
        state = self.state
        if not state["config"] or not state["model"]:
            ui.notify("Name the downstream storm config and the regional model config first.",
                      type="warning")
            return
        if state["process"] is not None and state["process"].poll() is None:
            ui.notify("The last run is still going.", type="warning")
            return
        argv = downstream.command(
            state["selection"], state["config"], state["model"],
            bryan_python=STATE.settings.bryan_python or None,
            duration=state["duration"], gwl=state["gwl"], dry_run=dry_run,
            suffix=state["suffix"], run_urbs=state["run_urbs"],
            sims_config=downstream.sims_config_for(state["selection"], self.folder,
                                                   self.project.config.config_path))
        log = Path(state["selection"]).with_suffix(".downstream.log")
        state["process"] = downstream.launch(argv, cwd=self.folder, log_path=log)
        if not dry_run:
            state["ran_urbs"] = bool(state["run_urbs"])
        doing = "Checking" if dry_run else (
            "Generating and running URBS" if state["run_urbs"] else "Generating")
        self._follow(doing, log)

    def _follow(self, doing: str, log: Path) -> None:
        """The one thing on this page that updates by itself: while the generator runs,
        how long it has been going and the last thing it wrote to its log."""
        state = self.state
        running = Running(self.spinner, doing, mark="downstream-running")
        self.status.text = f"Output in {log.name}"
        self.status.classes(replace="text-sm text-muted")

        def tick() -> None:
            process = state["process"]
            text = logtail.tail(log).text
            if running.cancel.is_set() and process.poll() is None:
                process.terminate()
            code = process.poll()
            if code is None:
                self.status.text = downstream.last_line(text) or f"Output in {log.name}"
                return
            timer.cancel()
            running.done()
            ok, said = downstream.outcome(code, text, log.name)
            self.status.text = said
            self.status.classes(replace="text-sm " + ("text-body" if ok else "text-negative"))
            if self.results_built:          # once, at the end - the tab never polls
                self.refresh_results()

        timer = ui.timer(1.0, tick)

    # -- the Results tab ---------------------------------------------------

    def refresh_results(self) -> None:
        """Read the runs back. Only on opening the tab, Reload, or a run ending: results
        change when URBS finishes an event, and the status line says when that is."""
        state = self.state
        self.results_built = True
        self.results_box.clear()
        with self.results_box:
            model = self._absolute(state["model"])
            if not model:
                ui.label("Name the regional model URBS config above to read its results."
                         ).classes("text-muted text-sm")
                return
            try:
                found = urbsresults.folders(model)
            except urbsresults.ModelConfigError as exc:
                severity_banner("warn", str(exc))
                return
            # the runs as the last run named them, whatever is in the suffix box now
            plan = urbsresults.as_run(
                downstream.plan(state["selection"], state["duration"], state["gwl"],
                                suffix=state["suffix"]), state["selection"])
            process = state["process"]
            rows = urbsresults.runs_rows(
                plan, found, urbsresults.storms_record(state["selection"]),
                generator_running=process is not None and process.poll() is None,
                ran_urbs=state["ran_urbs"])
            with ui.row().classes("items-center gap-3"):
                ui.button("Reload", icon="refresh", on_click=self.refresh_results) \
                    .props("flat dense").mark("results-reload")
                ui.label(urbsresults.counts(rows)).classes("text-sm text-body") \
                    .mark("results-counts")
                ui.label(f"results: {found.results}").classes("text-xs text-muted")
            if not rows:
                ui.label("No storm in this selection can be named yet - see Storms."
                         ).classes("text-muted text-sm")
                return
            columns = [{"name": key, "label": key, "field": key, "align": "left"}
                       for key in rows[0] if key != "why"]
            table = ui.table(columns=columns, rows=rows, row_key="run") \
                .classes("w-full").props("dense flat bordered wrap-cells").mark("runs-table")
            table.add_slot("body-cell-status", r'''
                <q-td :props="props"><span :title="props.row.why">{{ props.value }}</span></q-td>
            ''')
            ui.label("Hover a status for what it means. A run's starting storage and losses "
                     "are what the generator recorded when it wrote the storm."
                     ).classes("text-xs text-muted")
            self._dams(model, found, plan)
            self._print_locations(model, found, plan)

    def _print_locations(self, model, found, plan) -> None:
        """The flow at the * print locations - the hydraulic model's boundaries - every
        finished event overlaid, and whether any run stopped before its flood had passed."""
        flows = {}
        for storm in plan.storms:
            if storm.filename:
                run = os.path.splitext(storm.filename)[0]
                frame = urbsresults.read_flows(found.results, run)
                if frame is not None:
                    flows[run] = frame
        if not flows:
            return
        ui.separator()
        ui.label("Print locations").classes("text-lg").mark("locations-heading")
        cut = {run: urbsresults.end_of_run(frame) for run, frame in flows.items()}
        cut = {run: where for run, where in cut.items() if where}
        if cut:
            said = "; ".join(f"{run} at {', '.join(where[:4])}"
                             + (f" and {len(where) - 4} more" if len(where) > 4 else "")
                             for run, where in cut.items())
            severity_banner(
                "warn", f"{len(cut)} of {len(flows)} runs stop before the flood has passed - "
                        f"{said}.",
                "Still above a tenth of its peak at the end, or peaking at the very end. "
                "Lengthen simulation_periods in the regional model config for those durations.")
        else:
            ui.label(f"All {len(flows)} runs have receded at every print location by the end "
                     "of the run.").classes("text-sm text-body").mark("end-of-run-ok")
        locations = list(next(iter(flows.values())).columns)
        first = urbsresults.default_location(model, locations)

        def draw(chosen) -> None:
            charts.clear()
            with charts:
                for location in chosen or []:
                    series = {run: frame[location] for run, frame in flows.items()
                              if location in frame}
                    house_echart(damchart.flows_chart(series, location)) \
                        .classes("w-full h-80").mark(f"flows-{location}")

        ui.select(locations, multiple=True, value=[first] if first else [],
                  label="Print locations (* in the vec)",
                  on_change=lambda e: draw(e.value)).classes("w-full").props("use-chips") \
            .mark("locations")
        charts = ui.column().classes("w-full")
        draw([first] if first else [])

        # the peaks at every * location, for the report
        loading_of = {os.path.splitext(s.filename)[0]: s.loading
                      for s in plan.storms if s.filename}
        shown = {"what": urbsresults.PEAK}

        def peaks() -> list:
            return urbsresults.peaks_rows(flows, loading_of, shown["what"])

        def draw_peaks() -> None:
            peaks_box.clear()
            rows = peaks()
            with peaks_box:
                ui.table(columns=[{"name": key, "label": key, "field": key,
                                   "align": "left" if key in ("loading", "run") else "right"}
                                  for key in rows[0]], rows=rows, row_key="run") \
                    .classes("w-full").props("dense flat bordered").mark("peaks-table")

        with ui.row().classes("items-center gap-3"):
            ui.label("Peaks at every print location").classes("font-bold")
            ui.toggle({urbsresults.PEAK: "Peak flow (m³/s)", urbsresults.TIME: "Time of peak (h)"},
                      value=urbsresults.PEAK,
                      on_change=lambda e: (shown.update(what=e.value), draw_peaks())
                      ).props("no-caps dense").mark("peaks-show")
            copy_buttons(lambda: table_from_rows(list(peaks()[0]), peaks()), mark="peaks")
        peaks_box = ui.column().classes("w-full")
        draw_peaks()

    def _dams(self, model, found, plan) -> None:
        """Every dam against its reference levels: shown where the model has more than one
        dam or the config gives dam_levels, and nowhere else."""
        dams, levels, problems = urbsresults.dam_setup(model)
        if len(dams) < 2 and not levels:
            return
        ui.separator()
        ui.label("Dams").classes("text-lg").mark("dams-heading")
        for problem in problems:
            severity_banner("warn", problem)
        if not levels:
            ui.label("Give each dam's levels as dam_levels in the regional model config - "
                     'e.g. {"CALLIDE": {"DCL": 219.13, "PMF": {"1.3": 221.50}}} - to compare '
                     "its peaks against them.").classes("text-sm text-muted")
        rows = urbsresults.dams_rows(plan, found.results, dams, levels)
        if not rows:
            ui.label("No run of this selection has finished yet.").classes("text-sm text-muted")
            return
        reached = urbsresults.reaching_crest(rows, dams)
        with_crest = [dam for dam in dams if urbsresults.CREST in levels.get(dam, {})]
        if len(with_crest) > 1:
            ui.label(f"{reached} of {len(rows)} events take every dam to its crest level."
                     ).classes("text-sm text-body").mark("dams-reached")
        ui.table(columns=[{"name": key, "label": key, "field": key, "align": "left"}
                          for key in rows[0]], rows=rows, row_key="run") \
            .classes("w-full").props("dense flat bordered").mark("dams-table")
        copy_buttons(lambda: table_from_rows(list(rows[0]), rows), mark="dams")

        gwl_of = {os.path.splitext(s.filename)[0]: s.gwl for s in plan.storms if s.filename}
        runs = [row["run"] for row in rows]
        chart_box = ui.column().classes("w-full")
        shown = {"run": runs[0], "absolute": False}

        def draw() -> None:
            run = shown["run"]
            chart_box.clear()
            frame = urbsresults.read_levels(found.results, run)
            series = {dam: urbsresults.crest_series(frame, dam, levels.get(dam, {}), gwl_of.get(run))
                      for dam in with_crest}
            options = damchart.crest_chart(series, levels, gwl_of.get(run), title=run,
                                           absolute=shown["absolute"])
            with chart_box:
                if options["series"]:
                    house_echart(options).classes("w-full h-96").mark("crest-chart")
                else:
                    ui.label("No dam with a crest level in this run.").classes("text-sm text-muted")

        def choose(key, value) -> None:
            shown[key] = value
            draw()

        if with_crest:
            with ui.row().classes("items-center gap-4"):
                ui.select(runs, value=runs[0], label="Event",
                          on_change=lambda e: choose("run", e.value)).classes("w-80")
                ui.switch("Levels in m AHD", value=False,
                          on_change=lambda e: choose("absolute", e.value)) \
                    .mark("crest-absolute") \
                    .tooltip("Off: height above each dam's own crest, one zero line. "
                             "On: m AHD, each dam against its own crest line.")
            draw()
        if len(with_crest) > 1:
            points = urbsresults.margin_points(self.selections, found.results, with_crest,
                                               levels, downstream.selection_label)
            if points:
                house_echart(damchart.margin_scatter(points, *with_crest[:2])) \
                    .classes("w-full h-96").mark("margin-scatter")
        ui.label(
            "Read with care. Each loading is one event, so this shows what the chosen events "
            "did, not a probability that the dams reach their crests together. The other dams' "
            "starting storage is an assumption: the same lake_z through each dam's own curve "
            "for a Monte Carlo event, full supply for an ensemble one. And the rain on the other "
            "dams' catchments is conditional on the event over the main dam."
        ).classes("text-xs text-muted").mark("dams-caveats")
