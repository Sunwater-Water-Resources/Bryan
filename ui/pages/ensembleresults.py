"""Ensemble: the median-pattern design curves of an ensemble run, and its patterns.

The ensemble counterpart of the Results page's Durations tab, which reads the
Monte Carlo quantile tables. Pick an ensemble group and a result type:

- **the durations** - each duration's median pattern against AEP, with the
  envelope over them, marked up by which duration is critical, exactly as the
  Results page draws the Monte Carlo curves;
- **the critical durations** - per AEP, the critical duration, its median and
  which pattern gave it, the margin over the runner-up (metres for level), and
  the single highest event; checked against Bryan's own ``_critical.csv``;
- **the patterns at one AEP** - the PMF page's box plot, one box per duration.

Reads the open sims list, not the study. See ``core/enbresults.py``.
"""

from __future__ import annotations

import math

from nicegui import ui

from core import enbresults, ensemble, pmfchart, resultchart, results
from core.results import format_aep
from layout import page_frame, require_project, severity_banner
from clipboard import copy_buttons, table_from_rows
from state import STATE
from theme import house_echart

RESULT_OPTIONS = {"level": "Level", "inflow": "Inflow", "outflow": "Outflow"}


def ensemble_page() -> None:
    with page_frame("Ensemble"):
        project = require_project()
        if project is None:
            return
        available = enbresults.sources_by_group(project)
        if not available:
            _empty_state()
            return
        _EnsembleView(available).build()


def _empty_state() -> None:
    with ui.card().classes("w-full items-center p-8"):
        ui.icon("stacked_bar_chart").classes("text-5xl text-muted")
        ui.label("No ensemble results found for this sims list.").classes("text-muted")
        ui.label("An ensemble row writes its database once it has run. Monte Carlo "
                 "results are on the Results page.").classes(
            "text-xs text-muted max-w-lg text-center")
        ui.button("Choose simulations", on_click=lambda: ui.navigate.to("/select"))


def _reload() -> None:
    STATE.reload_project()
    ui.navigate.to("/ensemble")


def _set_options(chart, options) -> None:
    """``ui.echart.options`` is read-only: write into the dict it hands back."""
    chart.options.clear()
    chart.options.update(options)
    chart.update()


def _cell(value, result) -> str:
    if isinstance(value, str):
        return value
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "-"
    if math.isnan(number):
        return "-"
    return f"{number:,.2f}" if result == "level" else f"{number:,.0f}"


class _EnsembleView:
    def __init__(self, available) -> None:
        self.available = available       # group key -> [enbresults.Source]
        self.group = next(iter(available))
        self.result = "level"
        self.show_envelope = True
        self.show_markup = True
        self.aep = None
        self.database = None
        self.error = ""

        self.chart = None
        self.critical_chart = None
        self.files_box = None
        self.warning_box = None
        self.margin_note = None
        self.table_box = None
        self.aep_box = None
        self.box_chart = None
        self.pattern_box = None
        self.shown_table = None      # what Copy for Word copies: the table as drawn

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        self._controls()
        self._chart_card()
        self.warning_box = ui.column().classes("w-full gap-2")
        self._table_card()
        self._aep_card()
        self._load(self.group)

    def _controls(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.select(list(self.available), value=self.group, label="Group",
                          on_change=lambda event: self._load(event.value)) \
                    .classes("min-w-96").mark("ensemble-group")
                ui.toggle(RESULT_OPTIONS, value=self.result, on_change=self._on_result) \
                    .props("no-caps dense").mark("ensemble-result")
                ui.button("Reload", icon="refresh", on_click=_reload).props("flat dense")
            with ui.row().classes("items-center gap-6 flex-wrap"):
                ui.checkbox("Maximum envelope", value=self.show_envelope,
                            on_change=self._on_envelope)
                ui.checkbox("Mark up critical durations", value=self.show_markup,
                            on_change=self._on_markup)
            ui.label("Each line is one duration's median pattern at each AEP - Bryan's "
                     "own pick, the sixth of ten patterns in ascending order. The "
                     "critical duration is the one whose median is largest, as in "
                     "Bryan's _critical.csv. Each result type has its own.") \
                .classes("text-xs text-muted")

    def _chart_card(self) -> None:
        with ui.card().classes("w-full"):
            self.chart = house_echart({"series": []}).classes("w-full h-96") \
                .mark("ensemble-chart")
            with ui.expansion("Critical duration against AEP").classes("w-full"):
                self.critical_chart = house_echart({"series": []}) \
                    .classes("w-full h-64").mark("ensemble-critical-chart")
            with ui.expansion("Files read").classes("w-full"):
                self.files_box = ui.column().classes("w-full gap-0")

    def _table_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Critical durations").classes("font-bold")
                with ui.row().classes("gap-1"):
                    copy_buttons(lambda: self.shown_table, mark="ensemble")
            self.margin_note = ui.label("").classes("text-xs text-muted")
            self.table_box = ui.column().classes("w-full")

    def _aep_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("The patterns at one AEP").classes("font-bold")
                self.aep_box = ui.row().classes("items-center")
            ui.label("One box per duration over the temporal patterns, with every "
                     "pattern as a dot. The middle line is Bryan's median pick; the "
                     "diamond is the highest event.").classes("text-xs text-muted")
            self.box_chart = house_echart({"series": []}).classes("w-full h-80") \
                .mark("ensemble-box")
            self.pattern_box = ui.column().classes("w-full")

    # -- state ------------------------------------------------------------

    def _load(self, group) -> None:
        self.group = group
        try:
            self.database = enbresults.load(self.available[group])
            self.error = ""
        except enbresults.EnsembleError as exc:
            self.database, self.error = None, str(exc)
        aeps = self.database.aeps if self.database is not None else []
        if self.aep not in aeps:
            self.aep = aeps[-1] if aeps else None     # the rarest - usually the PMP
        self._draw_files()
        self.refresh()

    def _on_result(self, event) -> None:
        self.result = event.value
        self.refresh()

    def _on_envelope(self, event) -> None:
        self.show_envelope = event.value
        self.refresh()

    def _on_markup(self, event) -> None:
        self.show_markup = event.value
        self.refresh()

    def _on_aep(self, event) -> None:
        self.aep = event.value
        self._draw_aep()

    # -- drawing ----------------------------------------------------------

    def refresh(self) -> None:
        self.warning_box.clear()
        if self.database is None:
            with self.warning_box:
                severity_banner("block", self.error or "Nothing to show.")
            for chart in (self.chart, self.critical_chart, self.box_chart):
                _set_options(chart, {"series": []})
            self.table_box.clear()
            self.pattern_box.clear()
            self.aep_box.clear()
            return
        frame = self.database.frame
        found = enbresults.medians(frame, self.result)
        analysis = results.analyse(found.comparison)
        title = f"{self.group} - {results.y_axis(self.result)[0]}, median pattern"
        _set_options(self.chart, resultchart.duration_chart(
            found.comparison, analysis, self.result, show_envelope=self.show_envelope,
            show_markup=self.show_markup, title=title))
        _set_options(self.critical_chart,
                     resultchart.critical_duration_chart(found.comparison, analysis))
        self._draw_warnings(analysis, enbresults.check(self.database.sources, analysis,
                                                       self.result))
        self._draw_table(found, analysis, enbresults.highest(frame, self.result))
        self._draw_aep()

    def _draw_warnings(self, analysis, check) -> None:
        with self.warning_box:
            for note in self.database.notes:
                severity_banner("warn", note)
            if check is not None and check.compared:
                if check.agrees:
                    severity_banner("info", f"Agrees with Bryan's {check.path.name} at "
                                            f"all {check.compared} AEPs.")
                else:
                    hint = "\n".join(check.differences[:8])
                    if check.older:
                        hint += ("\nThe csv is older than the database, so it predates "
                                 "the last run - re-analyse to bring it up to date.")
                    severity_banner("warn", f"Differs from Bryan's {check.path.name} at "
                                            f"{len(check.differences)} of "
                                            f"{check.compared} AEPs.", hint)
            for text in analysis.warnings:
                severity_banner("info", text)
            if analysis.never_critical:
                severity_banner("info", f"Never critical at any AEP: "
                                        f"{', '.join(analysis.never_critical)}.")

    def _draw_table(self, found, analysis, top) -> None:
        self.table_box.clear()
        self.margin_note.set_text(results.margin_note(analysis.margin_kind))
        table = enbresults.critical_table(found, analysis, top)
        self.shown_table = None
        with self.table_box:
            if table.empty:
                ui.label("nothing to show").classes("text-muted text-sm")
                return
            # The highest event's duration goes in its own cell, and its pattern is
            # left to the box plot's table, so this fits beside a dozen durations.
            extra = ("highest duration", "highest pattern")
            shown = [name for name in table.columns if name not in extra]
            columns = [{"name": "aep", "label": "AEP (1 in X)", "field": "aep",
                        "align": "left"}]
            columns += [{"name": str(name), "label": str(name), "field": str(name)}
                        for name in shown]
            rows = []
            for aep, values in table.iterrows():
                row = {"aep": format_aep(aep)}
                for name in shown:
                    value = values[name]
                    if name == analysis.margin_label:
                        row[str(name)] = results.format_margin(value, analysis.margin_kind)
                    elif name == "highest" and isinstance(values.get(extra[0]), str):
                        row[name] = f"{_cell(value, self.result)} ({values[extra[0]]})"
                    else:
                        row[str(name)] = _cell(value, self.result)
                rows.append(row)
            self.shown_table = table_from_rows(
                columns, rows, title=f"{self.group} - {RESULT_OPTIONS.get(self.result, self.result)}")
            ui.table(columns=columns, rows=rows, row_key="aep") \
                .classes("w-full").props("dense flat bordered").mark("ensemble-table")

    def _draw_aep(self) -> None:
        self.aep_box.clear()
        self.pattern_box.clear()
        if self.database is None or self.aep is None:
            _set_options(self.box_chart, {"series": []})
            return
        with self.aep_box:
            ui.select({aep: f"1 in {format_aep(aep)}" for aep in self.database.aeps},
                      value=self.aep, label="AEP", on_change=self._on_aep) \
                .classes("w-48").props("dense").mark("ensemble-aep")
        local = enbresults.at_aep(self.database.frame, self.aep)
        stats = ensemble.by_duration(local, self.result)
        top = ensemble.pick(local, ensemble.HIGHEST, self.result)
        _set_options(self.box_chart, pmfchart.box_chart(
            stats, ensemble.pattern_points(local, self.result), self.result, highest=top))
        with self.pattern_box:
            if stats.empty:
                return
            columns = [{"name": k, "label": l, "field": k, "align": "right"}
                       for k, l in (("duration", "Duration (h)"), ("n", "Patterns"),
                                    ("median", "Median"), ("median_pattern", "Median pattern"),
                                    ("max", "Highest"), ("max_pattern", "Highest pattern"))]
            rows = [{"duration": f"{hours:g}", "n": int(row["n"]),
                     "median": _cell(row["median"], self.result),
                     "median_pattern": row["median_pattern"],
                     "max": _cell(row["max"], self.result),
                     "max_pattern": row["max_pattern"]}
                    for hours, row in stats.iterrows()]
            ui.table(columns=columns, rows=rows, row_key="duration") \
                .props("dense flat bordered").classes("w-full").mark("ensemble-patterns")

    def _draw_files(self) -> None:
        self.files_box.clear()
        with self.files_box:
            for source in self.available.get(self.group, []):
                ui.label(str(source.path)).classes("text-xs text-muted")
