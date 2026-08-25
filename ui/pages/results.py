"""Viewing results, on two tabs.

**Durations** compares storm durations inside one group: the maximum envelope,
who owns it, and whether the range that was run brackets the critical duration
- so the warnings are as much the point as the plot.

**Groups** compares groups with each other, one line each. That line is the
envelope, because the envelope *is* the design quantile, so the question it
answers is what a warmer climate or a raised full supply level does to the
flood - with the change from a baseline drawn underneath, in metres for level
and percent for flows.

Nothing here polls: results only change when a run finishes, so the page has a
Reload button instead of a timer. Refreshing still rebuilds the warnings, the
table and the file list in place, so the expansions around them are built once
and only their contents are cleared - rebuilding an expansion closes it, which
is what went wrong with the Run page's console output.

The one thing tabs cost: ECharts sizes itself when it is created, and a chart
created inside a hidden ``tab_panel`` measures zero and stays blank. So the
Groups tab draws nothing until it is first shown, and asks its charts to resize
when it is - see ``_GroupView.activate``.
"""

from __future__ import annotations

from pathlib import Path

from nicegui import run as nicerun
from nicegui import ui

from core import critexport, grouping, overlay, results, resultchart
from core.paths import cell_text
from layout import page_frame, require_project, severity_banner
from state import STATE

# The order the result types are offered in; volume windows follow.
TYPE_ORDER = ("inflow", "level", "outflow")

# How many groups the Groups tab ticks for you. Above this it ticks two - a
# baseline and one comparison - rather than drawing a dozen curves nobody
# asked for. Every group is still listed with its own checkbox either way.
INITIAL_GROUPS = 8


def results_page() -> None:
    with page_frame("Results"):
        project = require_project()
        if project is None:
            return
        available = scan(project)
        if not available:
            _empty_state()
            return

        durations = _ResultsView(project, available)
        groups = _GroupView(project, available)
        with ui.tabs().classes("w-full").mark("result-tabs") as tabs:
            ui.tab("Durations")
            ui.tab("Groups")
        with ui.tab_panels(tabs, value="Durations").classes("w-full"):
            with ui.tab_panel("Durations").classes("p-0"):
                durations.build()
            with ui.tab_panel("Groups").classes("p-0"):
                groups.build()
        tabs.on_value_change(
            lambda event: groups.activate(durations.key)
            if event.value == "Groups" else None)


def scan(project) -> dict:
    """Which groups have anything plottable. Probes files, runs nothing."""
    frame = project.frame
    folder = project.config.project_folder
    keys = project.group_keys()
    available = {}
    for key in grouping.groups_in_order(frame):
        rows = [index for index in frame.index if keys[index] == key]
        found = results.sources_for_rows(frame, folder, rows)
        if found:
            available[key] = found
    return available


def _empty_state() -> None:
    with ui.card().classes("w-full items-center p-8"):
        ui.icon("insights").classes("text-5xl text-gray-400")
        ui.label("No analysed results found for this sims list."
                 ).classes("text-gray-500")
        ui.label("Rows write their quantile files when 'Analyse results' "
                 "is yes. The ensemble method does its own critical "
                 "duration analysis instead, so it does not appear here."
                 ).classes("text-xs text-gray-500 max-w-lg text-center")
        ui.button("Choose simulations", on_click=lambda: ui.navigate.to("/select"))


def _ordered_keys(keys) -> list:
    """Result types in the page's order, volume windows after the three peaks."""
    keys = list(keys)
    lead = [key for key in TYPE_ORDER if key in keys]
    return lead + sorted(key for key in keys if key not in TYPE_ORDER)


def _reload() -> None:
    STATE.reload_project()
    ui.navigate.to("/results")


def _set_options(chart, options) -> None:
    """Replace a chart's options.

    ``ui.echart.options`` is a read-only property - it hands back the dict
    inside ``_props``, so a new one has to be written into that same dict
    rather than assigned over the property.
    """
    chart.options.clear()
    chart.options.update(options)
    chart.update()


class _ResultsView:
    def __init__(self, project, available) -> None:
        self.project = project
        self.available = available   # group key -> {result key: [CurveSource]}
        self.group = next(iter(available), None)
        self.key = None
        self.selected = set()        # labels
        self.show_envelope = True
        self.show_markup = True
        self.noise_floor = None      # None -> the default for the result type
        self.aep_from = None
        self.aep_to = None
        self._range_aeps = None      # only redraw the range selects when it changes

        self.chart = None
        self.critical_chart = None
        self.curve_box = None
        self.type_box = None
        self.range_box = None
        self.floor_box = None
        self._floor_kind = None
        self.warning_box = None
        self.table_box = None
        self.files_box = None

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        self._controls()
        self._chart_card()
        self._warnings_card()
        self._table_card()
        self._load_group(self.group)

    def _controls(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.select(list(self.available), value=self.group, label="Group",
                          on_change=lambda event: self._load_group(event.value)
                          ).classes("min-w-96")
                ui.button("Reload", icon="refresh", on_click=_reload
                          ).props("flat dense")
                ui.button("Export critical durations", icon="download",
                          on_click=self._export_dialog).props("outline dense")
            self.type_box = ui.row().classes("items-center gap-4 flex-wrap")
            ui.separator()
            ui.label("Durations").classes("text-xs text-gray-500")
            self.curve_box = ui.row().classes("items-center gap-x-4 gap-y-1 flex-wrap")
            with ui.row().classes("items-center gap-2"):
                ui.button("All", on_click=lambda: self._all(True)).props("flat dense")
                ui.button("None", on_click=lambda: self._all(False)).props("flat dense")
            ui.separator()
            with ui.row().classes("items-center gap-6 flex-wrap"):
                ui.checkbox("Maximum envelope", value=self.show_envelope,
                            on_change=self._on_envelope)
                ui.checkbox("Mark up critical durations", value=self.show_markup,
                            on_change=self._on_markup)
                self.range_box = ui.row().classes("items-center gap-2")
                self.floor_box = ui.row().classes("items-center gap-2")

    def _chart_card(self) -> None:
        with ui.card().classes("w-full"):
            self.chart = ui.echart({"series": []}).classes("w-full h-96") \
                .mark("duration-chart")
            with ui.expansion("Critical duration against AEP").classes("w-full"):
                ui.label(
                    "For lake level this normally falls from left to right: the "
                    "long durations win while the storage is filling, the short "
                    "ones on the rare tail once the dam is acting as a "
                    "conveyance. Peak inflow follows the catchment, not the "
                    "storage, so it does not do this."
                ).classes("text-xs text-gray-500")
                self.critical_chart = ui.echart({"series": []}) \
                    .classes("w-full h-64").mark("critical-duration-chart")
            # Built once, not per refresh: rebuilding an expansion closes it,
            # which is what the Run page's console output had to be taught.
            with ui.expansion("Files read").classes("w-full"):
                self.files_box = ui.column().classes("w-full gap-0")

    def _warnings_card(self) -> None:
        self.warning_box = ui.column().classes("w-full gap-2")

    def _table_card(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Critical durations").classes("font-bold")
            ui.label("'margin %' is how far the critical duration beat the "
                     "runner-up. A few tenths of a percent is sampling noise, "
                     "not a crossover.").classes("text-xs text-gray-500")
            self.table_box = ui.column().classes("w-full")

    # -- state ------------------------------------------------------------

    def _load_group(self, group) -> None:
        self.group = group
        found = self.available.get(group, {})
        keys = _ordered_keys(found)
        if self.key not in keys:
            self.key = keys[0] if keys else None
        self.selected = {source.label for source in found.get(self.key, [])}
        self.aep_from = self.aep_to = self._range_aeps = None
        self._draw_types()
        self._draw_curves()
        self.refresh()

    def _sources(self) -> list:
        found = self.available.get(self.group, {}).get(self.key, [])
        return [source for source in found if source.label in self.selected]

    def _draw_types(self) -> None:
        found = self.available.get(self.group, {})
        self.type_box.clear()
        with self.type_box:
            keys = _ordered_keys(found)
            options = {key: results.y_axis(key)[0].split(" (")[0] for key in keys}
            ui.toggle(options, value=self.key, on_change=self._on_type
                      ).props("no-caps dense")

    def _draw_curves(self) -> None:
        sources = self.available.get(self.group, {}).get(self.key, [])
        self.curve_box.clear()
        with self.curve_box:
            if not sources:
                ui.label("nothing to plot").classes("text-gray-500 text-sm")
            for source in sources:
                checkbox = ui.checkbox(source.label,
                                       value=source.label in self.selected)
                checkbox.mark(f"curve-{source.label}")
                checkbox.on_value_change(
                    lambda event, label=source.label: self._toggle(label, event.value))
                with checkbox:
                    ui.tooltip(str(source.path))

    def _draw_range(self, aeps) -> None:
        """Trimming the ends, the way util's drop_aeps does for its plots.

        Only rebuilt when the AEP set itself changes: redrawing these on every
        refresh would mean clearing the very select whose event is running.
        """
        if aeps == self._range_aeps:
            return
        self._range_aeps = list(aeps)
        self.range_box.clear()
        options = {aep: f"1 in {results.format_aep(aep)}" for aep in aeps}
        with self.range_box:
            ui.select(options, value=self.aep_from if self.aep_from in options
                      else aeps[0], label="From",
                      on_change=self._on_from).classes("min-w-36").props("dense")
            ui.select(options, value=self.aep_to if self.aep_to in options
                      else aeps[-1], label="To",
                      on_change=self._on_to).classes("min-w-36").props("dense")

    def _toggle(self, label, on) -> None:
        if on:
            self.selected.add(label)
        else:
            self.selected.discard(label)
        self.refresh()

    def _all(self, on) -> None:
        sources = self.available.get(self.group, {}).get(self.key, [])
        self.selected = {source.label for source in sources} if on else set()
        self._draw_curves()
        self.refresh()

    def _on_type(self, event) -> None:
        self.key = event.value
        sources = self.available.get(self.group, {}).get(self.key, [])
        self.selected = {source.label for source in sources}
        self.aep_from = self.aep_to = self._range_aeps = None
        self.noise_floor = self._floor_kind = None
        self._draw_curves()
        self.refresh()

    def _on_envelope(self, event) -> None:
        self.show_envelope = event.value
        self.refresh()

    def _on_markup(self, event) -> None:
        self.show_markup = event.value
        self.refresh()

    def _on_from(self, event) -> None:
        self.aep_from = event.value
        self.refresh()

    def _on_to(self, event) -> None:
        self.aep_to = event.value
        self.refresh()

    # -- export -----------------------------------------------------------

    def _export_dialog(self) -> None:
        """Everything the export will write, before it writes any of it."""
        sources_by_key = self.available.get(self.group, {})
        chosen = {self.key} if self.key else set()
        folder = critexport.default_folder(self._sources())
        base = critexport.default_base_name(self._sources())

        with ui.dialog() as dialog, ui.card().classes("min-w-[36rem]"):
            ui.label("Export critical durations").classes("text-lg font-bold")
            ui.label("Runs util/CriticalDurationAnalysis.py with Bryan's own "
                     "interpreter, so the files match the ones the study "
                     "post-processing already produces."
                     ).classes("text-xs text-gray-500")

            folder_input = ui.input("Output folder", value=str(folder or "")
                                    ).classes("w-full").props("dense")
            name_input = ui.input("Base name", value=base
                                  ).classes("w-full").props("dense")

            ui.label("Result types").classes("text-xs text-gray-500 mt-2")
            with ui.row().classes("gap-x-4 flex-wrap"):
                for key in _ordered_keys(sources_by_key):
                    box = ui.checkbox(key, value=key in chosen)
                    box.on_value_change(
                        lambda event, k=key: (chosen.add(k) if event.value
                                              else chosen.discard(k)))
            trimmed = self._trimmed_aeps()
            drop = ui.checkbox(
                f"Leave the {len(trimmed)} trimmed AEP(s) off the plot",
                value=bool(trimmed))
            if not trimmed:
                drop.set_visibility(False)

            preview = ui.column().classes("w-full gap-0")
            report = ui.column().classes("w-full gap-1")

            def refresh_preview() -> None:
                preview.clear()
                built = critexport.plan(
                    sources_by_key, sorted(chosen),
                    folder=Path(folder_input.value) if folder_input.value else None,
                    base_name=name_input.value.strip(),
                    python=STATE.settings.bryan_python,
                    drop_aeps=trimmed if drop.value else (),
                )
                with preview:
                    for problem in built.problems:
                        severity_banner("block", problem)
                    for job in built.jobs:
                        for path in job.outputs:
                            exists = " - exists, will be overwritten" if path.is_file() else ""
                            ui.label(f"{path}{exists}").classes(
                                "text-xs " + ("text-orange-700" if exists
                                              else "text-gray-500"))
                run_button.set_enabled(built.can_run)
                return built

            with ui.row().classes("justify-end gap-2 w-full"):
                ui.button("Close", on_click=dialog.close).props("flat")
                run_button = ui.button("Export")

            async def export() -> None:
                built = refresh_preview()
                if not built.can_run:
                    return
                run_button.set_enabled(False)
                report.clear()
                for job in built.jobs:
                    with report:
                        ui.label(f"{job.key}...").classes("text-xs text-gray-500")
                    result = await nicerun.io_bound(critexport.run, job)
                    report.clear()
                    with report:
                        if result.ok:
                            severity_banner(
                                "info", f"{job.key}: wrote "
                                        f"{', '.join(p.name for p in result.written)}")
                        else:
                            severity_banner(
                                "block", f"{job.key}: the export failed "
                                         f"(exit {result.returncode}).",
                                result.output[-1500:])
                run_button.set_enabled(True)
                ui.notify("Export finished")

            run_button.on_click(export)
            for element in (folder_input, name_input):
                element.on("blur", refresh_preview)
            drop.on_value_change(refresh_preview)
            refresh_preview()
        dialog.open()

    def _trimmed_aeps(self) -> tuple:
        """The AEPs the From/To range currently leaves off the plot."""
        sources = self._sources()
        if not sources:
            return ()
        full = results.compare(sources)
        if full.is_empty or self._range_aeps is None:
            return ()
        aeps = list(full.frame.index)
        low = self.aep_from if self.aep_from in aeps else aeps[0]
        high = self.aep_to if self.aep_to in aeps else aeps[-1]
        if low > high:
            low, high = high, low
        return tuple(aep for aep in aeps if aep < low or aep > high)

    # -- drawing ----------------------------------------------------------

    def refresh(self) -> None:
        comparison = results.compare(self._sources())
        comparison = self._trim(comparison)
        analysis = results.analyse(comparison, noise_floor=self.noise_floor)
        self._draw_floor(analysis)

        title = f"{self.group} - {results.y_axis(self.key)[0] if self.key else ''}"
        _set_options(self.chart, resultchart.duration_chart(
            comparison, analysis, self.key, show_envelope=self.show_envelope,
            show_markup=self.show_markup, title=title))
        _set_options(self.critical_chart,
                     resultchart.critical_duration_chart(comparison, analysis))

        self._draw_warnings(comparison, analysis)
        self._draw_table(comparison, analysis)
        self._draw_files(comparison)

    def _trim(self, comparison):
        if comparison.is_empty:
            return comparison
        aeps = list(comparison.frame.index)
        self._draw_range(aeps)
        low = self.aep_from if self.aep_from in aeps else aeps[0]
        high = self.aep_to if self.aep_to in aeps else aeps[-1]
        if low > high:
            low, high = high, low
        kept = comparison.frame.loc[low:high]
        return results.Comparison(frame=kept, durations=comparison.durations,
                                  problems=comparison.problems)

    def _draw_floor(self, analysis) -> None:
        """The noise floor, in the units of the result type.

        Rebuilt only when those units change - a level plot is judged in metres
        and a flow plot in percent, so the control cannot be one shared number.
        """
        if analysis.margin_kind == self._floor_kind:
            return
        self._floor_kind = analysis.margin_kind
        self.noise_floor = analysis.noise_floor
        self.floor_box.clear()
        with self.floor_box:
            unit = "m" if analysis.margin_kind == results.ABSOLUTE else "%"
            field = ui.number(f"Noise floor ({unit})", value=analysis.noise_floor,
                              min=0, step=0.01, format="%.3f",
                              on_change=self._on_floor).classes("w-40").props("dense")
            with field:
                ui.tooltip("A crossover whose new duration never gets this far "
                           "clear of the next one is reported as sampling "
                           "noise rather than marked up.")

    def _on_floor(self, event) -> None:
        if event.value is None:
            return
        self.noise_floor = float(event.value)
        self.refresh()

    def _draw_warnings(self, comparison, analysis) -> None:
        self.warning_box.clear()
        with self.warning_box:
            for label, why in comparison.problems:
                severity_banner("warn", f"{label}: {why}")
            for text in analysis.warnings:
                severity_banner("info", text)
            if analysis.never_critical:
                severity_banner(
                    "info",
                    f"Never critical at any AEP shown: "
                    f"{', '.join(analysis.never_critical)}.",
                    "Worth knowing before running them again.")

    def _draw_table(self, comparison, analysis) -> None:
        self.table_box.clear()
        table = results.table(comparison, analysis)
        with self.table_box:
            if table.empty:
                ui.label("nothing selected").classes("text-gray-500 text-sm")
                return
            columns = [{"name": "aep", "label": "AEP (1 in X)", "field": "aep",
                        "align": "left"}]
            columns += [{"name": str(name), "label": str(name), "field": str(name)}
                        for name in table.columns]
            rows = []
            for aep, values in table.iterrows():
                row = {"aep": results.format_aep(aep)}
                for name, value in values.items():
                    row[str(name)] = (
                        results.format_margin(value, analysis.margin_kind)
                        if name == analysis.margin_label else _cell(name, value))
                rows.append(row)
            ui.table(columns=columns, rows=rows, row_key="aep"
                     ).classes("w-full").props("dense flat bordered")

    def _draw_files(self, comparison) -> None:
        self.files_box.clear()
        sources = [source for source in self._sources()
                   if source.label in comparison.frame.columns]
        with self.files_box:
            for source in sources:
                ui.label(f"{source.label}: {source.path}"
                         ).classes("text-xs text-gray-500")



def _initial_groups(offered) -> list:
    offered = list(offered)
    return offered if len(offered) <= INITIAL_GROUPS else offered[:2]


class _GroupView:
    """Overlaying one envelope per group.

    A group contributes one line however many durations it ran, because the
    envelope over its durations *is* the design quantile
    (``util/MaxQuantiles.py:96``). The maximum *across* groups is not a
    quantity at all - they are scenarios, not alternatives - so what goes
    underneath is the change from a baseline instead of an envelope.
    """

    def __init__(self, project, available) -> None:
        self.project = project
        self.available = available
        self.key = None
        self.selected = []           # group keys, in the order they were scanned
        self.baseline = None         # a group key
        self.aep_from = None
        self.aep_to = None
        self._range_aeps = None
        self._baseline_options = None
        self._filter = ""
        self._started = False

        self.chart = None
        self.delta_chart = None
        self.critical_chart = None
        self.type_box = None
        self.group_box = None
        self.baseline_box = None
        self.range_box = None
        self.delta_box = None
        self.note_box = None
        self.table_box = None
        self.files_box = None

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        self._controls()
        self._chart_card()
        self._notes_card()
        self._table_card()

    def activate(self, preferred_key=None) -> None:
        """First sight of the tab: choose what to show, then draw it.

        Deferred because a chart built inside a hidden tab panel measures zero
        and renders blank - and for the same reason the charts are asked to
        resize now that the panel has a size.
        """
        if not self._started:
            self._started = True
            keys = _ordered_keys(overlay.result_keys(self.available))
            self.key = preferred_key if preferred_key in keys else (
                keys[0] if keys else None)
            self.selected = _initial_groups(
                overlay.groups_with(self.available, self.key))
            self.baseline = self.selected[0] if self.selected else None
            self._draw_types()
            self._draw_groups()
            self._draw_baseline()
        self.refresh()
        self._resize()

    def _resize(self) -> None:
        for chart in (self.chart, self.delta_chart, self.critical_chart):
            if chart is None:
                continue
            try:
                chart.run_chart_method("resize")
            except Exception:        # noqa: BLE001 - never worth failing a page over
                pass

    def _controls(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.label("One line per group: the maximum over that group's "
                         "storm durations, which is the design quantile."
                         ).classes("text-xs text-gray-500")
                ui.button("Reload", icon="refresh", on_click=_reload
                          ).props("flat dense")
            self.type_box = ui.row().classes("items-center gap-4 flex-wrap")
            ui.separator()
            with ui.row().classes("items-center gap-2 flex-wrap"):
                ui.label("Groups to overlay").classes("text-xs text-gray-500")
                ui.input(placeholder="filter", on_change=self._on_filter
                         ).props("dense clearable").classes("w-64")
                ui.button("All", on_click=lambda: self._all(True)).props("flat dense")
                ui.button("None", on_click=lambda: self._all(False)).props("flat dense")
            self.group_box = ui.column().classes("w-full gap-0")
            ui.separator()
            with ui.row().classes("items-center gap-6 flex-wrap"):
                self.baseline_box = ui.row().classes("items-center gap-2")
                self.range_box = ui.row().classes("items-center gap-2")

    def _chart_card(self) -> None:
        with ui.card().classes("w-full"):
            self.chart = ui.echart({"series": []}).classes("w-full h-96") \
                .mark("overlay-chart")
            self.delta_box = ui.column().classes("w-full gap-0")
            with self.delta_box:
                ui.label(
                    "Change from the baseline. Level is compared in metres and "
                    "flows and volumes in percent - a percentage of a level on "
                    "an arbitrary datum says nothing."
                ).classes("text-xs text-gray-500")
                self.delta_chart = ui.echart({"series": []}) \
                    .classes("w-full h-64").mark("delta-chart")
            # nothing to compare until a baseline and a second group exist
            self.delta_box.set_visibility(False)
            with ui.expansion("Critical duration against AEP").classes("w-full"):
                ui.label(
                    "Whether the critical duration itself moves between groups "
                    "- a warmer climate or a raised full supply level can shift "
                    "it, and the single-group view cannot show that."
                ).classes("text-xs text-gray-500")
                self.critical_chart = ui.echart({"series": []}) \
                    .classes("w-full h-64").mark("critical-overlay-chart")
            with ui.expansion("Files read").classes("w-full"):
                self.files_box = ui.column().classes("w-full gap-0")

    def _notes_card(self) -> None:
        self.note_box = ui.column().classes("w-full gap-2")

    def _table_card(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Envelopes by group").classes("font-bold")
            ui.label("The design quantile each group produced, and the change "
                     "from the baseline.").classes("text-xs text-gray-500")
            self.table_box = ui.column().classes("w-full")

    # -- state ------------------------------------------------------------

    def _labels(self) -> dict:
        """Trimmed over every group that offers this result type.

        Over the *offered* groups rather than the selected ones, so ticking a
        group on or off never renames the rest of them mid-comparison.
        """
        return overlay.distinguishing_labels(
            overlay.groups_with(self.available, self.key))

    def _baseline_label(self) -> str:
        return self._labels().get(self.baseline, "")

    def _draw_types(self) -> None:
        keys = _ordered_keys(overlay.result_keys(self.available))
        self.type_box.clear()
        with self.type_box:
            options = {key: results.y_axis(key)[0].split(" (")[0] for key in keys}
            ui.toggle(options, value=self.key, on_change=self._on_type
                      ).props("no-caps dense").mark("group-types")

    def _draw_groups(self) -> None:
        offered = overlay.groups_with(self.available, self.key)
        labels = self._labels()
        needle = self._filter.lower()
        self.group_box.clear()
        with self.group_box:
            if not offered:
                ui.label("no group has this result type"
                         ).classes("text-gray-500 text-sm")
                return
            shown = [group for group in offered
                     if needle in f"{labels[group]} {group}".lower()]
            if not shown:
                ui.label(f"nothing matches '{self._filter}'"
                         ).classes("text-gray-500 text-sm")
            for group in shown:
                count = len(self.available[group][self.key])
                with ui.row().classes("items-center gap-2"):
                    box = ui.checkbox(labels[group], value=group in self.selected)
                    box.mark(f"group-{labels[group]}")
                    box.on_value_change(
                        lambda event, key=group: self._toggle(key, event.value))
                    with box:
                        ui.tooltip(group)
                    ui.label(f"{count} duration{'' if count == 1 else 's'}"
                             ).classes("text-xs text-gray-500")
            hidden = len(offered) - len(shown)
            if hidden:
                ui.label(f"{hidden} group(s) hidden by the filter"
                         ).classes("text-xs text-gray-500")
            if len(offered) == 1:
                ui.label("Only one group has results of this type, so there is "
                         "nothing to overlay it against yet."
                         ).classes("text-xs text-gray-500")

    def _draw_baseline(self) -> None:
        """Only rebuilt when the option set changes - see ``_draw_range``."""
        labels = self._labels()
        options = {group: labels[group] for group in self.selected}
        if list(options) == self._baseline_options:
            return
        self._baseline_options = list(options)
        self.baseline_box.clear()
        if not options:
            return
        with self.baseline_box:
            ui.select(options, value=self.baseline, label="Baseline",
                      on_change=self._on_baseline).classes("min-w-48").props("dense")

    def _draw_range(self, aeps) -> None:
        if aeps == self._range_aeps:
            return
        self._range_aeps = list(aeps)
        self.range_box.clear()
        options = {aep: f"1 in {results.format_aep(aep)}" for aep in aeps}
        with self.range_box:
            ui.select(options, value=self.aep_from if self.aep_from in options
                      else aeps[0], label="From",
                      on_change=self._on_from).classes("min-w-36").props("dense")
            ui.select(options, value=self.aep_to if self.aep_to in options
                      else aeps[-1], label="To",
                      on_change=self._on_to).classes("min-w-36").props("dense")

    def _toggle(self, group, on) -> None:
        chosen = set(self.selected)
        if on:
            chosen.add(group)
        else:
            chosen.discard(group)
        self._set_selection(chosen)

    def _all(self, on) -> None:
        offered = overlay.groups_with(self.available, self.key)
        self._set_selection(set(offered) if on else set())
        self._draw_groups()

    def _set_selection(self, chosen) -> None:
        offered = overlay.groups_with(self.available, self.key)
        self.selected = [group for group in offered if group in chosen]
        if self.baseline not in self.selected:
            self.baseline = self.selected[0] if self.selected else None
        self._draw_baseline()
        self.refresh()

    def _on_type(self, event) -> None:
        self.key = event.value
        offered = overlay.groups_with(self.available, self.key)
        kept = [group for group in self.selected if group in offered]
        self.selected = kept or _initial_groups(offered)
        if self.baseline not in self.selected:
            self.baseline = self.selected[0] if self.selected else None
        self.aep_from = self.aep_to = self._range_aeps = None
        self._baseline_options = None
        self._draw_groups()
        self._draw_baseline()
        self.refresh()

    def _on_filter(self, event) -> None:
        self._filter = event.value or ""
        self._draw_groups()

    def _on_baseline(self, event) -> None:
        self.baseline = event.value
        self.refresh()

    def _on_from(self, event) -> None:
        self.aep_from = event.value
        self.refresh()

    def _on_to(self, event) -> None:
        self.aep_to = event.value
        self.refresh()

    # -- drawing ----------------------------------------------------------

    def refresh(self) -> None:
        built = overlay.build(self.available, self.selected, self.key,
                              self.project.frame, labels=self._labels())
        built = self._trim(built)
        change = overlay.deltas(built, self._baseline_label())
        critical, skipped = overlay.critical_frame(built)

        title = results.y_axis(self.key)[0] if self.key else ""
        _set_options(self.chart, resultchart.overlay_chart(built, title=title))
        _set_options(self.delta_chart, resultchart.delta_chart(built, change))
        _set_options(self.critical_chart, resultchart.critical_overlay_chart(
            built, self._clip(critical)))
        self.delta_box.set_visibility(not change.is_empty)

        self._draw_notes(built, change, skipped)
        self._draw_table(built, change)
        self._draw_files(built)

    def _trim(self, built):
        if built.is_empty:
            return built
        self._draw_range(list(built.frame.index))
        return overlay.Overlay(frame=self._clip(built.frame), curves=built.curves,
                               key=built.key, problems=built.problems,
                               notes=built.notes)

    def _clip(self, frame):
        """The AEP range the From/To selects are showing."""
        if frame is None or frame.empty:
            return frame
        aeps = list(frame.index)
        low = self.aep_from if self.aep_from in aeps else aeps[0]
        high = self.aep_to if self.aep_to in aeps else aeps[-1]
        if low > high:
            low, high = high, low
        return frame.loc[low:high]

    def _draw_notes(self, built, change, skipped) -> None:
        self.note_box.clear()
        with self.note_box:
            for who, why in built.problems:
                severity_banner("warn", f"{who}: {why}")
            for note in built.notes:
                severity_banner("info", note)
            if change.undefined:
                severity_banner(
                    "info",
                    f"{change.baseline} is zero or missing at 1 in "
                    f"{', '.join(results.format_aep(aep) for aep in change.undefined)}"
                    f", so no change is shown there.",
                    "A dam that does not spill has no outflow quantile at the "
                    "frequent end.")
            if skipped:
                severity_banner(
                    "info",
                    f"Critical duration not plotted for {', '.join(skipped)}.",
                    "Their storm durations could not be read off the sims list.")

    def _draw_table(self, built, change) -> None:
        self.table_box.clear()
        table = overlay.table(built, change)
        with self.table_box:
            if table.empty:
                ui.label("nothing selected").classes("text-gray-500 text-sm")
                return
            changed = set() if change.is_empty else {
                f"{column} {change.label}" for column in change.frame.columns}
            columns = [{"name": "aep", "label": "AEP (1 in X)", "field": "aep",
                        "align": "left"}]
            columns += [{"name": str(name), "label": str(name), "field": str(name)}
                        for name in table.columns]
            rows = []
            for aep, values in table.iterrows():
                row = {"aep": results.format_aep(aep)}
                for name, value in values.items():
                    row[str(name)] = (results.format_margin(value, change.kind)
                                      if name in changed else _cell(name, value))
                rows.append(row)
            ui.table(columns=columns, rows=rows, row_key="aep"
                     ).classes("w-full").props("dense flat bordered")

    def _draw_files(self, built) -> None:
        self.files_box.clear()
        with self.files_box:
            for curve in built.curves:
                ui.label(curve.key).classes("text-xs font-medium")
                for source in curve.sources:
                    ui.label(f"    {curve.label} {source.label}: {source.path}"
                             ).classes("text-xs text-gray-500")


def _cell(name, value) -> str:
    if isinstance(value, str):
        return value
    try:
        number = float(value)
    except (TypeError, ValueError):
        return cell_text(value)
    if number != number:                       # NaN - the duration never got here
        return "-"
    return f"{number:,.6g}"
