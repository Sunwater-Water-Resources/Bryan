"""Report: the design flood report's result tables, straight from the runs.

A **study file** names the runs a report is written from - one sims_config.json
each, RFSL and FSL, the adopted run and the one it replaced - and holds one entry
per report table saying which group of which run fills it. Each table is shown
here as it will paste, with **Copy for Word** (a formatted table) and **Copy as
text** (tab-separated, for Excel). When a rating is revised, re-run the groups
and copy the tables again; nothing is typed.

The study file is saved on every change. See ``core/study.py`` for its format and
``core/reporttables.py`` for what each kind of table computes.
"""

from __future__ import annotations

import copy
from pathlib import Path

from nicegui import app, run, ui

from core import reporttables, study as studies, wordtable
from layout import page_frame, severity_banner
from state import STATE


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


def report_page() -> None:
    with page_frame("Report"):
        _ReportView().build()


class _ReportView:
    def __init__(self) -> None:
        self.study = STATE.study
        self.built: dict = {}
        self.body = None

    # -- layout ------------------------------------------------------------

    def build(self) -> None:
        self.body = ui.column().classes("w-full gap-4")
        self.redraw()

    def redraw(self) -> None:
        self.body.clear()
        with self.body:
            if self.study is None:
                self._open_card()
                return
            self._study_card()
            self._runs_card()
            self._tables()

    def save(self) -> None:
        try:
            self.study.save()
        except OSError as exc:
            ui.notify(f"Could not save {self.study.path.name}: {exc}", type="warning")

    # -- opening a study ---------------------------------------------------

    def _open_card(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Study file").classes("text-lg font-bold")
            ui.label("A study file names the runs a report is written from and says "
                     "which group of which run fills each report table. Keep it at "
                     "the top of the study folder; the paths in it are stored "
                     "relative to it, so the study can be moved."
                     ).classes("text-sm text-body")
            default = STATE.settings.last_study or ""
            if not default and STATE.project is not None:
                default = str(STATE.project.config.config_path.parent
                              / studies.DEFAULT_NAME)
            path = ui.input("Study file (or its folder)", value=default) \
                .classes("w-full").props("dense").mark("study-path")
            name = ui.input("Name, for a new study", value="").classes("w-full") \
                .props("dense")
            with ui.row().classes("gap-2"):
                ui.button("Open", icon="folder_open",
                          on_click=lambda: self._open(path.value)).mark("open-study")
                ui.button("New study", icon="add",
                          on_click=lambda: self._open(path.value, create=True,
                                                      name=name.value)) \
                    .props("outline").mark("new-study")
            recent = [entry for entry in STATE.settings.recent_studies if entry]
            if recent:
                ui.label("Recent").classes("text-sm text-muted pt-2")
                for entry in recent:
                    ui.button(entry, on_click=lambda _, e=entry: self._open(e)) \
                        .props("flat dense no-caps align=left").classes("mono text-sm")

    def _open(self, path, *, create=False, name="") -> None:
        text = str(path or "").strip().strip('"')
        if not text:
            ui.notify("Give the study file's path", type="warning")
            return
        try:
            self.study = STATE.open_study(text, create=create, name=name)
        except (studies.StudyError, OSError) as exc:
            ui.notify(str(exc), type="negative", multi_line=True)
            return
        self.built = {}
        self.redraw()

    def _close(self) -> None:
        STATE.close_study()
        self.study = None
        self.redraw()

    # -- the study ---------------------------------------------------------

    def _study_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between no-wrap"):
                with ui.column().classes("gap-0"):
                    name = ui.input("Study", value=self.study.name) \
                        .props("dense borderless").classes("text-lg font-bold")
                    name.on("blur", lambda: self._rename_study(name.value))
                    ui.label(str(self.study.path)).classes("mono text-xs text-muted")
                with ui.row().classes("gap-2"):
                    ui.button("Refresh", icon="refresh", on_click=self._refresh) \
                        .props("outline").mark("refresh-tables") \
                        .tooltip("Re-read the results - after a re-run")
                    ui.button("Close", icon="close", on_click=self._close).props("flat")

    def _rename_study(self, value) -> None:
        if value != self.study.name:
            self.study.name = value
            self.save()

    def _refresh(self) -> None:
        studies.forget_runs()
        reporttables.forget_cached()
        self.built = {}
        self.redraw()

    # -- runs --------------------------------------------------------------

    def _runs_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Runs").classes("text-lg font-bold")
                ui.button("Add run", icon="add", on_click=self._add_run_dialog) \
                    .props("outline dense").mark("add-run")
            if not self.study.runs:
                ui.label("No runs yet. Add each sims_config.json the report draws on - "
                         "the RFSL and FSL lists are separate runs."
                         ).classes("text-sm text-muted")
                return
            for entry in list(self.study.runs):
                try:
                    count = len(self.study.open_run(entry["name"]).groups())
                    groups = f"{count} group{'' if count == 1 else 's'}"
                    colour = "text-body"
                except studies.StudyError as exc:
                    groups, colour = f"cannot open: {exc}", "text-negative"
                with ui.row().classes("w-full items-center gap-3 no-wrap") \
                        .mark(f"run-row-{entry['name']}"):
                    name = ui.input(value=entry["name"]).props("dense").classes("w-48") \
                        .tooltip("Rename - the tables that read this run follow it")
                    name.on("blur", lambda _, old=entry["name"], box=name:
                            self._rename_run(old, box.value))
                    ui.label(entry["sims_config"]).classes("mono text-xs grow")
                    ui.label(groups).classes(f"text-xs {colour}")
                    ui.button(icon="delete",
                              on_click=lambda _, n=entry["name"]: self._remove_run(n)) \
                        .props("flat dense round").tooltip("Remove from the study")

    def _rename_run(self, old, new) -> None:
        if not new or new == old:
            return
        try:
            self.study.rename_run(old, new)
        except studies.StudyError as exc:
            ui.notify(str(exc), type="negative")
            return
        self.save()
        self.redraw()

    def _remove_run(self, name) -> None:
        orphans = self.study.remove_run(name)
        self.save()
        if orphans:
            ui.notify(f"Removed {name}; these tables still name it: "
                      f"{', '.join(orphans)}", type="warning", multi_line=True)
        self.redraw()

    def _add_run_dialog(self) -> None:
        suggested = ""
        if STATE.project is not None:
            suggested = str(STATE.project.config.config_path)
        with ui.dialog() as dialog, ui.card().classes("min-w-[36rem]"):
            ui.label("Add a run").classes("text-lg font-bold")
            path = ui.input("sims_config.json", value=suggested).classes("w-full") \
                .props("dense").mark("run-path")
            name = ui.input("Name", value=_guess_run_name(suggested)).classes("w-full") \
                .props("dense").mark("run-name")
            path.on("blur", lambda: name.value or name.set_value(_guess_run_name(path.value)))

            def add() -> None:
                try:
                    self.study.add_run(name.value, path.value)
                except studies.StudyError as exc:
                    ui.notify(str(exc), type="negative")
                    return
                self.save()
                dialog.close()
                self.redraw()

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=dialog.close).props("flat")
                ui.button("Add", on_click=add).mark("confirm-add-run")
        dialog.open()

    # -- tables ------------------------------------------------------------

    def _tables(self) -> None:
        with ui.row().classes("w-full items-center justify-between"):
            ui.label("Tables").classes("text-lg font-bold")
            with ui.button("Add table", icon="add").props("outline").mark("add-table"):
                with ui.menu():
                    for kind in reporttables.KINDS.values():
                        with ui.menu_item(on_click=lambda _, k=kind.key: self._new_table(k)):
                            with ui.column().classes("gap-0"):
                                ui.label(kind.label)
                                ui.label(kind.help).classes("text-xs text-muted")
        if not self.study.tables:
            ui.label("No tables yet.").classes("text-sm text-muted")
            return
        for spec in self.study.tables:
            self._table_card(spec)

    def _table_card(self, spec: dict) -> None:
        table_id = spec.get("id")
        kind = reporttables.KINDS.get(spec.get("kind"))
        with ui.card().classes("w-full").mark(f"table-{table_id}"):
            with ui.row().classes("w-full items-start justify-between no-wrap"):
                with ui.column().classes("gap-0"):
                    ui.label(spec.get("title") or "(untitled)").classes("font-bold")
                    ui.label(_describe(spec, kind)).classes("text-xs text-muted")
                with ui.row().classes("gap-1 no-wrap"):
                    ui.button("Copy for Word", icon="content_copy",
                              on_click=lambda: self._copy(table_id, rich=True)) \
                        .props("dense no-caps").mark(f"copy-word-{table_id}")
                    ui.button(icon="notes",
                              on_click=lambda: self._copy(table_id, rich=False)) \
                        .props("flat dense round").mark(f"copy-text-{table_id}") \
                        .tooltip("Copy as text (tab-separated, for Excel)")
                    ui.button(icon="edit", on_click=lambda: self._edit(table_id)) \
                        .props("flat dense round").mark(f"edit-{table_id}").tooltip("Edit")
                    ui.button(icon="content_paste_go",
                              on_click=lambda: self._duplicate(table_id)) \
                        .props("flat dense round").tooltip(
                            "Duplicate - for the same table from another group")
                    ui.button(icon="arrow_upward",
                              on_click=lambda: self._move(table_id, -1)) \
                        .props("flat dense round").tooltip("Move up")
                    ui.button(icon="arrow_downward",
                              on_click=lambda: self._move(table_id, 1)) \
                        .props("flat dense round").tooltip("Move down")
                    ui.button(icon="delete", on_click=lambda: self._delete(table_id)) \
                        .props("flat dense round").tooltip("Remove")
            preview = ui.column().classes("w-full gap-1")
            with preview:
                ui.spinner(size="sm")
            ui.timer(0.01, lambda: self._fill(table_id, preview), once=True)

    async def _fill(self, table_id, box) -> None:
        spec = self.study.table(table_id)
        if spec is None:
            return
        built = await _off_thread(reporttables.build, self.study, spec)
        self.built[table_id] = built
        if box.is_deleted:
            return
        box.clear()
        with box:
            if built.problems:
                shown = built.problems[:6]
                more = len(built.problems) - len(shown)
                severity_banner("warn", "\n".join(shown)
                                + (f"\n... and {more} more" if more > 0 else ""))
            ui.html(wordtable.to_html(built, self.study.extra.get("word")),
                    sanitize=False).classes("overflow-x-auto").mark(f"preview-{table_id}")

    async def _copy(self, table_id, *, rich: bool) -> None:
        spec = self.study.table(table_id)
        built = self.built.get(table_id)
        if built is None and spec is not None:
            built = await _off_thread(reporttables.build, self.study, spec)
        if built is None:
            return
        text = wordtable.to_text(built)
        if not rich:
            ui.clipboard.write(text)
            ui.notify("Copied as text")
            return
        fragment = wordtable.to_html(built, self.study.extra.get("word"))
        outcome = await ui.run_javascript(
            wordtable.clipboard_script(wordtable.clipboard_document(fragment), text),
            timeout=5.0)
        if outcome == "html":
            ui.notify("Copied - paste into Word")
        elif outcome == "text":
            ui.notify("This browser would only take text; copied as text",
                      type="warning")
        else:
            ui.notify(f"Could not copy: {outcome}", type="negative")

    def _new_table(self, kind: str) -> None:
        spec = reporttables.new_spec(kind)
        spec["title"] = reporttables.KINDS[kind].label
        if not reporttables.KINDS[kind].multi and self.study.runs:
            spec["source"]["run"] = self.study.runs[0]["name"]
        _TableEditor(self, spec, is_new=True).open()

    def _edit(self, table_id) -> None:
        spec = self.study.table(table_id)
        if spec is not None:
            _TableEditor(self, copy.deepcopy(spec)).open()

    def _duplicate(self, table_id) -> None:
        spec = copy.deepcopy(self.study.table(table_id) or {})
        if not spec:
            return
        spec.pop("id", None)
        spec["title"] = f"{spec.get('title', '')} (copy)"
        _TableEditor(self, spec, is_new=True).open()

    def _move(self, table_id, step) -> None:
        self.study.move_table(table_id, step)
        self.save()
        self.redraw()

    def _delete(self, table_id) -> None:
        spec = self.study.table(table_id) or {}
        with ui.dialog() as dialog, ui.card():
            ui.label(f"Remove '{spec.get('title') or table_id}' from the study?")
            ui.label("The results are not touched.").classes("text-sm text-muted")
            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=dialog.close).props("flat")

                def remove() -> None:
                    self.study.remove_table(table_id)
                    self.built.pop(table_id, None)
                    self.save()
                    dialog.close()
                    self.redraw()

                ui.button("Remove", on_click=remove).props("color=negative")
        dialog.open()

    def store(self, spec: dict) -> None:
        stored = self.study.put_table(spec)
        self.built.pop(stored["id"], None)
        self.save()
        self.redraw()


def _guess_run_name(path) -> str:
    """'runs/E013/CLD_RFSL_mc_sims_01.json' -> 'E013 RFSL'."""
    text = str(path or "").strip().strip('"')
    if not text:
        return ""
    local = Path(text)
    parent = local.parent.name
    stem = local.stem.upper()
    state = next((token for token in ("RFSL", "FSL") if token in stem.split("_")), "")
    prefix = "PMF " if stem.startswith("PMF") else ""
    return " ".join(part for part in (parent, prefix.strip(), state) if part) or local.stem


def _describe(spec, kind) -> str:
    label = kind.label if kind else f"unknown kind {spec.get('kind')!r}"
    sources = studies.table_sources(spec)
    if not sources:
        return label
    if len(sources) == 1:
        source = sources[0]
        return f"{label} - {source.get('group') or '(no group)'} in {source.get('run') or '(no run)'}"
    runs = sorted({source.get("run") for source in sources if source.get("run")})
    return f"{label} - {len(sources)} groups from {', '.join(runs)}"


# -- editing one table ---------------------------------------------------------

class _TableEditor:
    """A dialog over a working copy of one spec; nothing is saved until Save."""

    def __init__(self, view: _ReportView, spec: dict, *, is_new=False) -> None:
        self.view = view
        self.study = view.study
        self.spec = spec
        self.is_new = is_new
        self.kind = reporttables.KINDS[spec["kind"]]
        self.dialog = None

    def open(self) -> None:
        with ui.dialog() as self.dialog, \
                ui.card().classes("min-w-[52rem] max-w-[64rem]"):
            ui.label(f"{'New' if self.is_new else 'Edit'} table - {self.kind.label}"
                     ).classes("text-lg font-bold")
            ui.label(self.kind.help).classes("text-xs text-muted")
            ui.input("Title (for your reference - Word's caption stays yours)",
                     value=self.spec.get("title", "")) \
                .classes("w-full").props("dense").mark("table-title") \
                .on_value_change(lambda e: self.spec.update(title=e.value))
            self.form = ui.column().classes("w-full gap-2")
            self.draw()
            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=self.dialog.close).props("flat")
                ui.button("Save", on_click=self.save).mark("save-table")
        self.dialog.open()

    def save(self) -> None:
        self.view.store(self.spec)
        self.dialog.close()

    def draw(self) -> None:
        self.form.clear()
        with self.form:
            if self.kind.key == reporttables.DESIGN_FLOODS:
                self._design_floods()
            else:
                self._multi()
            ui.input("Footnote, pasted under the table", value=self.spec.get("footnote", "")) \
                .classes("w-full").props("dense") \
                .on_value_change(lambda e: self.spec.update(footnote=e.value))

    # -- one source --------------------------------------------------------

    def _groups(self, run_name) -> list:
        if not run_name:
            return []
        try:
            return self.study.open_run(run_name).groups()
        except studies.StudyError:
            return []

    def _source_picker(self, holder: dict, *, mark="") -> None:
        runs = self.study.run_names()
        if not runs:
            ui.label("Add a run to the study first.").classes("text-sm text-muted")
            return

        def on_run(event) -> None:
            holder["run"] = event.value
            holder["group"] = ""
            self.draw()

        groups = self._groups(holder.get("run"))
        with ui.row().classes("w-full gap-2 no-wrap"):
            ui.select(runs, value=holder.get("run") if holder.get("run") in runs else None,
                      label="Run", on_change=on_run).classes("w-48").props("dense") \
                .mark(f"run-{mark}" if mark else "run")
            ui.select(groups, value=holder.get("group") if holder.get("group") in groups else None,
                      label="Group", with_input=True,
                      on_change=lambda e: holder.update(group=e.value)) \
                .classes("grow").props("dense").mark(f"group-{mark}" if mark else "group")

    def _design_floods(self) -> None:
        self.spec.setdefault("source", {"run": "", "group": ""})
        self._source_picker(self.spec["source"], mark="source")
        with ui.row().classes("w-full gap-2 no-wrap"):
            ui.input("AEP of the PMP (1 in x)",
                     value=f"{self.spec.get('pmp_aep') or ''}") \
                .classes("w-48").props("dense") \
                .on_value_change(lambda e: self.spec.update(
                    pmp_aep=(reporttables.parse_aeps(e.value) or [None])[0]))
            ui.input("Its label", value=self.spec.get("pmp_label", "")) \
                .classes("w-32").props("dense") \
                .on_value_change(lambda e: self.spec.update(pmp_label=e.value))
            ui.checkbox("Thousands separators", value=self.spec.get("thousands", True),
                        on_change=lambda e: self.spec.update(thousands=e.value))
        ui.textarea("AEPs (1 in x), in order", value=reporttables.aeps_text(
            self.spec.get("aeps") or reporttables.DEFAULT_AEPS)) \
            .classes("w-full").props("dense autogrow") \
            .on_value_change(lambda e: self.spec.update(
                aeps=reporttables.parse_aeps(e.value)))

    # -- sections of rows ----------------------------------------------------

    def _multi(self) -> None:
        ui.input("First column heading", value=self.spec.get("first_column", "")) \
            .classes("w-64").props("dense") \
            .on_value_change(lambda e: self.spec.update(first_column=e.value))
        key = self.kind.key
        if key == reporttables.FLOOD_LEVELS:
            with ui.row().classes("w-full gap-4 no-wrap items-start"):
                ui.textarea("Levels, one 'label = level (m AHD)' per line",
                            value=reporttables.levels_text(self.spec.get("levels"))) \
                    .classes("grow").props("dense autogrow").mark("levels") \
                    .on_value_change(lambda e: self.spec.update(
                        levels=reporttables.parse_levels(e.value)))
                with ui.column().classes("gap-1"):
                    ui.number("Round AEPs to", value=self.spec.get("round_to", 10),
                              min=1, step=1) \
                        .props("dense").classes("w-32") \
                        .on_value_change(lambda e: self.spec.update(round_to=e.value))
                    ui.checkbox("Critical duration column",
                                value=self.spec.get("duration", True),
                                on_change=lambda e: self.spec.update(duration=e.value))
        if key == reporttables.PEAK_AT_AEP:
            ui.input("AEP (1 in x)", value=f"{self.spec.get('aep') or ''}") \
                .classes("w-48").props("dense") \
                .on_value_change(lambda e: self.spec.update(
                    aep=(reporttables.parse_aeps(e.value) or [None])[0]))
        if key in (reporttables.PEAK_AT_AEP, reporttables.ENSEMBLE_PEAK):
            ui.radio({reporttables.AT_LEVEL: "Inflow and outflow of the storm that "
                                             "gives the peak level",
                      reporttables.AT_OWN: "Each result's own maximum (what the old "
                                           "scripts did)"},
                     value=self.spec.get("flows_at", reporttables.AT_LEVEL),
                     on_change=lambda e: self.spec.update(flows_at=e.value)) \
                .props("dense").mark("flows-at")

        sections = self.spec.setdefault("sections", [])
        for s_index, section in enumerate(sections):
            with ui.card().classes("w-full").props("flat bordered"):
                with ui.row().classes("w-full items-center gap-2 no-wrap"):
                    ui.input("Section heading (blank for none)",
                             value=section.get("heading", "")) \
                        .classes("grow").props("dense") \
                        .on_value_change(lambda e, s=section: s.update(heading=e.value))
                    ui.button(icon="delete", on_click=lambda _, i=s_index: self._drop(
                        sections, i)).props("flat dense round").tooltip("Remove the section")
                rows = section.setdefault("rows", [])
                for r_index, row in enumerate(rows):
                    self._row_editor(rows, r_index, row, f"{s_index}-{r_index}")
                with ui.row().classes("gap-2"):
                    ui.button("Row from a group", icon="add",
                              on_click=lambda _, r=rows: self._append(
                                  r, {"label": "", "run": self._default_run(),
                                      "group": ""})) \
                        .props("flat dense no-caps").mark(f"add-row-{s_index}")
                    ui.button("Row of fixed values", icon="add",
                              on_click=lambda _, r=rows: self._append(
                                  r, {"label": "", "values": []})) \
                        .props("flat dense no-caps")
        ui.button("Add section", icon="add",
                  on_click=lambda: self._append(sections, {"heading": "", "rows": []})) \
            .props("outline dense no-caps").mark("add-section")

    def _row_editor(self, rows, index, row, mark) -> None:
        with ui.row().classes("w-full items-center gap-2 no-wrap"):
            ui.input("Row label", value=row.get("label", "")) \
                .classes("w-48").props("dense").mark(f"row-label-{mark}") \
                .on_value_change(lambda e, r=row: r.update(label=e.value))
            if "values" in row:
                ui.input("Cells, separated by |", value=" | ".join(row.get("values") or [])) \
                    .classes("grow").props("dense") \
                    .on_value_change(lambda e, r=row: r.update(
                        values=reporttables.parse_values(e.value)))
            else:
                with ui.element("div").classes("grow"):
                    self._source_picker(row, mark=mark)
            ui.button(icon="close", on_click=lambda: self._drop(rows, index)) \
                .props("flat dense round").tooltip("Remove the row")

    def _default_run(self) -> str:
        runs = self.study.run_names()
        return runs[0] if runs else ""

    def _append(self, target: list, item: dict) -> None:
        target.append(item)
        self.draw()

    def _drop(self, target: list, index: int) -> None:
        if 0 <= index < len(target):
            target.pop(index)
        self.draw()
