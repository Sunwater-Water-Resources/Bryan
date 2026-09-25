"""Figures: the report's frequency-curve figures, drawn from the study's runs.

Each figure lives in the study file and names its curves by run and group, so a
re-run moves every figure with it. The preview is the data the PNG will be drawn
from; the label of each curve is edited beside it and belongs to this figure
only, because one curve is 'URBS' against the FFA and 'GWL 1.3 °C' among the
horizons. **Export PNG** runs util/ReportFigure.py under Bryan's interpreter, in
the style of PlotFrequencyCurves_v03.py, and leaves the job (``<name>.json``)
beside the PNG as a record of what was drawn. See ``core/figures.py``.
"""

from __future__ import annotations

import copy
from pathlib import Path

from nicegui import app, run, ui

from core import figurechart, figures, staleness, study as studies
from layout import confirm, no_study, page_frame, severity_banner
from state import STATE
from theme import house_echart
from widgets import OUTPUT_FOLDER, path_input

TYPES = {"level": "Lake level", "inflow": "Peak inflow", "outflow": "Peak outflow",
         "inflowVol24h": "24 h inflow volume", "inflowVol48h": "48 h inflow volume",
         "inflowVol72h": "72 h inflow volume"}
KIND_LABELS = {figures.GROUP: "A group's design curve", figures.FILE: "A curve from a file",
               figures.FFA: "A flood frequency analysis (RMC Bestfit export)"}


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


def figures_page() -> None:
    with page_frame("Figures"):
        _FiguresView().build()


def _number(value):
    try:
        return None if value in (None, "") else float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None


class _FiguresView:
    def __init__(self) -> None:
        self.study = STATE.study
        self.body = None

    def build(self) -> None:
        if self.study is None:
            no_study("Figures are kept in the study file.")
            return
        with ui.row().classes("w-full items-center justify-between"):
            ui.label(f"Figures are kept in {self.study.path}").classes("text-xs text-muted")
            with ui.row().classes("gap-2"):
                ui.button("Add figure", icon="add", on_click=self._new) \
                    .props("outline").mark("add-figure")
                ui.button("Export all", icon="collections", on_click=self._export_all) \
                    .props("outline").mark("export-all")
        self.body = ui.column().classes("w-full gap-4")
        self.redraw()

    def save(self) -> None:
        try:
            self.study.save()
        except OSError as exc:
            ui.notify(f"Could not save {self.study.path.name}: {exc}", type="warning")

    def redraw(self) -> None:
        self.body.clear()
        with self.body:
            specs = figures.figures(self.study)
            if not specs:
                ui.label("No figures yet. Add one: a lake level, inflow or outflow "
                         "frequency plot of any groups in the study's runs, with FFA "
                         "results and reference levels if wanted.").classes("text-sm text-muted")
            for spec in specs:
                _FigureCard(self, spec).build()

    def store(self, spec) -> dict:
        stored = figures.put(self.study, spec)
        self.save()
        self.redraw()
        return stored

    def _new(self) -> None:
        runs = self.study.run_names()
        spec = figures.new_spec(filename="new_figure")
        if runs:
            spec["curves"].append({"kind": figures.GROUP, "run": runs[0], "group": "",
                                   "label": ""})
        _FigureEditor(self, spec, is_new=True).open()

    async def _export_all(self) -> None:
        done, failed = 0, []
        for spec in figures.figures(self.study):
            result = await _off_thread(_export, self.study, spec)
            if result.ok:
                done += 1
            else:
                failed.append(f"{spec.get('filename')}: {result.output.strip()[-200:]}")
        if failed:
            ui.notify(f"{done} exported; failed: " + "; ".join(failed), type="warning",
                      multi_line=True)
        else:
            ui.notify(f"All {done} figures exported")
        self.redraw()


def _export(study, spec):
    data = figures.build(study, spec)
    return figures.export(data, figures.output_path(study, spec), STATE.settings.bryan_python)


class _FigureCard:
    def __init__(self, view: _FiguresView, spec: dict) -> None:
        self.view = view
        self.study = view.study
        self.spec = spec
        self.figure_id = spec.get("id")

    def build(self) -> None:
        with ui.card().classes("w-full").mark(f"figure-{self.figure_id}"):
            with ui.row().classes("w-full items-start justify-between no-wrap"):
                with ui.column().classes("gap-0"):
                    ui.label(self.spec.get("filename") or self.figure_id).classes("font-bold")
                    ui.label(f"{TYPES.get(self.spec.get('type'), self.spec.get('type'))} - "
                             f"{len(self.spec.get('curves') or [])} curves - "
                             f"{figures.output_path(self.study, self.spec)}") \
                        .classes("text-xs text-muted")
                with ui.row().classes("gap-1 no-wrap"):
                    ui.button("Export PNG", icon="image", on_click=self._export) \
                        .props("dense no-caps").mark(f"export-{self.figure_id}")
                    ui.button(icon="edit", on_click=self._edit) \
                        .props("flat dense round").tooltip("Edit").mark(f"edit-{self.figure_id}")
                    ui.button(icon="content_paste_go", on_click=self._duplicate) \
                        .props("flat dense round").tooltip(
                            "Duplicate - swap a curve for a sensitivity figure")
                    ui.button(icon="delete", on_click=self._ask_delete) \
                        .props("flat dense round").tooltip("Remove") \
                        .mark(f"remove-{self.figure_id}")
            with ui.row().classes("w-full items-start gap-4 no-wrap"):
                self.chart_box = ui.column().classes("grow min-w-0")
                with self.chart_box:
                    ui.spinner()
                with ui.column().classes("w-72 shrink-0 gap-1"):
                    ui.label("Labels in this figure").classes("text-sm text-muted")
                    for index, curve in enumerate(self.spec.get("curves") or []):
                        ui.input(value=curve.get("label", ""),
                                 placeholder=curve.get("group") or Path(
                                     str(curve.get("file", ""))).stem) \
                            .props("dense").classes("w-full") \
                            .mark(f"label-{self.figure_id}-{index}") \
                            .on("blur", lambda e, i=index: self._relabel(i, e.sender.value))
                    png = figures.output_path(self.study, self.spec)
                    if png.is_file():
                        ui.label("Last export").classes("text-sm text-muted pt-2")
                        ui.image(png).classes("w-full").props("fit=contain") \
                            .mark(f"last-export-{self.figure_id}")
        ui.timer(0.01, self._fill, once=True)

    async def _fill(self) -> None:
        data = await _off_thread(figures.build, self.study, self.spec)
        sources = [curve for curve in self.spec.get("curves") or []
                   if curve.get("kind", figures.GROUP) == figures.GROUP]
        stale = await _off_thread(staleness.for_sources, self.study, sources) or []
        if self.chart_box.is_deleted:
            return
        self.chart_box.clear()
        with self.chart_box:
            if stale:
                severity_banner("warn", "\n".join(stale),
                                "Re-run the group before exporting this figure.")
            if data.problems:
                severity_banner("warn", "\n".join(data.problems[:6]))
            house_echart(figurechart.preview(data)).classes("w-full h-96") \
                .mark(f"preview-{self.figure_id}")

    def _relabel(self, index, value) -> None:
        curves = self.spec.get("curves") or []
        if index >= len(curves) or curves[index].get("label", "") == value:
            return
        spec = copy.deepcopy(self.spec)
        spec["curves"][index]["label"] = value
        self.view.store(spec)

    async def _export(self) -> None:
        result = await _off_thread(_export, self.study, self.spec)
        if result.ok:
            ui.notify(f"Wrote {result.png}")
        else:
            ui.notify(f"Export failed: {result.output.strip()[-400:]}", type="negative",
                      multi_line=True)
        self.view.redraw()

    def _edit(self) -> None:
        _FigureEditor(self.view, copy.deepcopy(self.spec)).open()

    def _duplicate(self) -> None:
        spec = copy.deepcopy(self.spec)
        spec["id"] = ""
        spec["filename"] = f"{spec.get('filename') or 'figure'}_copy"
        _FigureEditor(self.view, spec, is_new=True).open()

    def _ask_delete(self) -> None:
        confirm(f"Remove the figure {self.spec.get('filename') or self.figure_id}?",
                "Its exported PNG is left where it is.", self._delete)

    def _delete(self) -> None:
        figures.remove(self.study, self.figure_id)
        self.view.save()
        self.view.redraw()


class _FigureEditor:
    """A dialog over a working copy of one figure; nothing is saved until Save."""

    def __init__(self, view: _FiguresView, spec: dict, *, is_new=False) -> None:
        self.view = view
        self.study = view.study
        self.spec = spec
        self.is_new = is_new

    def open(self) -> None:
        with ui.dialog() as self.dialog, ui.card().classes("min-w-[56rem] max-w-[68rem]"):
            ui.label("New figure" if self.is_new else "Edit figure").classes("text-lg font-bold")
            self.form = ui.column().classes("w-full gap-2")
            self.draw()
            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=self.dialog.close).props("flat")
                ui.button("Save", on_click=self.save).mark("save-figure")
        self.dialog.open()

    def save(self) -> None:
        if not str(self.spec.get("filename") or "").strip():
            ui.notify("Give the figure a file name", type="warning")
            return
        self.view.store(self.spec)
        self.dialog.close()

    def _path_base(self) -> dict:
        return {"base": self.study.folder, "base_name": "the study folder"}

    def draw(self) -> None:
        spec = self.spec
        self.form.clear()
        with self.form:
            with ui.row().classes("w-full gap-2 no-wrap"):
                ui.input("File name (no extension)", value=spec.get("filename", "")) \
                    .classes("w-64").props("dense").mark("figure-filename") \
                    .on_value_change(lambda e: spec.update(filename=e.value))
                path_input("Folder (blank: figures/ beside the study)",
                           spec.get("folder", ""), **self._path_base(),
                           expect=OUTPUT_FOLDER, classes="grow",
                           on_change=lambda text: spec.update(folder=text))
                ui.select(TYPES, value=spec.get("type", "level"), label="Result",
                          on_change=lambda e: spec.update(type=e.value)) \
                    .classes("w-48").props("dense").mark("figure-type")
            with ui.row().classes("w-full gap-2 no-wrap items-center"):
                for key, label, width in (("min_aep", "From 1 in", "w-28"),
                                          ("max_aep", "to 1 in", "w-32"),
                                          ("aep_of_pmp", "AEP of PMP line (1 in)", "w-44"),
                                          ("y_min", "y from", "w-24"),
                                          ("y_max", "y to", "w-24")):
                    value = spec.get(key)
                    ui.input(label, value="" if value is None else f"{value:g}") \
                        .classes(width).props("dense") \
                        .on_value_change(lambda e, k=key: spec.update({k: _number(e.value)}))
                ui.input("Title", value=spec.get("title", "")).classes("grow").props("dense") \
                    .on_value_change(lambda e: spec.update(title=e.value))
                ui.checkbox("on the figure", value=spec.get("show_title", False),
                            on_change=lambda e: spec.update(show_title=e.value))
            ui.textarea("Reference levels, one 'label = level' or 'label = level = colour' "
                        "per line", value=figures.reference_levels_text(
                            spec.get("reference_levels"))) \
                .classes("w-full").props("dense autogrow").mark("figure-levels") \
                .on_value_change(lambda e: spec.update(
                    reference_levels=figures.parse_reference_levels(e.value)))
            ui.label("Curves, in legend order").classes("text-sm font-bold pt-2")
            curves = spec.setdefault("curves", [])
            for index, curve in enumerate(curves):
                self._curve_row(curves, index, curve)
            with ui.row().classes("gap-2"):
                for kind, label in KIND_LABELS.items():
                    ui.button(label, icon="add",
                              on_click=lambda _, k=kind: self._add(curves, k)) \
                        .props("flat dense no-caps").mark(f"add-curve-{kind}")

    def _curve_row(self, curves, index, curve) -> None:
        kind = curve.get("kind", figures.GROUP)
        with ui.card().classes("w-full").props("flat bordered"):
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                ui.label(f"{index + 1}.").classes("text-muted")
                ui.input("Label", value=curve.get("label", "")).classes("w-48") \
                    .props("dense").mark(f"curve-label-{index}") \
                    .on_value_change(lambda e: curve.update(label=e.value))
                if kind == figures.GROUP:
                    self._group_picker(curve, index)
                elif kind == figures.FILE:
                    path_input("File (AEP in the first column)", curve.get("file", ""),
                               **self._path_base(), suffixes=(".csv",), classes="grow",
                               on_change=lambda text: curve.update(file=text))
                    ui.input("Column (blank: the result)", value=curve.get("column", "")) \
                        .classes("w-40").props("dense") \
                        .on_value_change(lambda e: curve.update(column=e.value))
                else:
                    path_input("RMC Bestfit export (.csv)", curve.get("file", ""),
                               **self._path_base(), suffixes=(".csv",), classes="grow",
                               on_change=lambda text: curve.update(file=text))
                    ui.select(list(figures.POSTERIORS), value=curve.get("posterior", "Both"),
                              label="Posterior",
                              on_change=lambda e: curve.update(posterior=e.value)) \
                        .classes("w-32").props("dense")
                ui.button(icon="arrow_upward", on_click=lambda: self._move(curves, index, -1)) \
                    .props("flat dense round")
                ui.button(icon="close", on_click=lambda: self._drop(curves, index)) \
                    .props("flat dense round").tooltip("Remove the curve")

    def _group_picker(self, curve, index) -> None:
        runs = self.study.run_names()
        groups = []
        if curve.get("run"):
            try:
                groups = self.study.open_run(curve["run"]).groups()
            except studies.StudyError:
                groups = []

        def on_run(event) -> None:
            curve["run"], curve["group"] = event.value, ""
            self.draw()

        ui.select(runs, value=curve.get("run") if curve.get("run") in runs else None,
                  label="Run", on_change=on_run).classes("w-40").props("dense") \
            .mark(f"curve-run-{index}")
        ui.select(groups, value=curve.get("group") if curve.get("group") in groups else None,
                  label="Group", with_input=True,
                  on_change=lambda e: curve.update(group=e.value)) \
            .classes("grow").props("dense").mark(f"curve-group-{index}")

    def _add(self, curves, kind) -> None:
        curve = {"kind": kind, "label": ""}
        if kind == figures.GROUP:
            runs = self.study.run_names()
            curve.update(run=runs[0] if runs else "", group="")
        elif kind == figures.FFA:
            curve.update(file="", posterior="Both")
        else:
            curve.update(file="", column="")
        curves.append(curve)
        self.draw()

    def _move(self, curves, index, step) -> None:
        there = index + step
        if 0 <= there < len(curves):
            curves.insert(there, curves.pop(index))
        self.draw()

    def _drop(self, curves, index) -> None:
        if 0 <= index < len(curves):
            curves.pop(index)
        self.draw()

