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

from nicegui import app, run, ui

from core import reporttables, staleness, study as studies, wordtable
from layout import confirm, no_study, page_frame, severity_banner
import address
from clipboard import copy_table
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
        self.stale: dict = {}           # table id -> what is out of date in its results
        self.body = None
        self.cards: dict = {}           # table id -> its card's parts
        self._went_to = False           # the address's table shown, once
        self.open: set = (STATE.settings.open_tables_for(self.study.path)
                          if self.study is not None else set())

    # -- layout ------------------------------------------------------------

    def build(self) -> None:
        self.body = ui.column().classes("w-full gap-4")
        self.redraw()

    def redraw(self) -> None:
        self.body.clear()
        self.cards = {}
        with self.body:
            if self.study is None:
                no_study("The report tables are kept in the study file.")
                return
            self._study_card()
            self._runs_card()
            self._tables()

    def save(self) -> None:
        try:
            self.study.save()
        except OSError as exc:
            ui.notify(f"Could not save {self.study.path.name}: {exc}", type="warning")

    # -- the study ---------------------------------------------------------

    def _study_card(self) -> None:
        """The study, named - it is opened, renamed and closed on the Study page."""
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between no-wrap"):
                with ui.column().classes("gap-0"):
                    ui.label(self.study.name or self.study.path.stem) \
                        .classes("text-lg font-bold").mark("report-study")
                    ui.label(str(self.study.path)).classes("mono text-xs text-muted")
                with ui.row().classes("gap-2"):
                    ui.button("Refresh", icon="refresh", on_click=self._refresh) \
                        .props("outline").mark("refresh-tables") \
                        .tooltip("Re-read the results - after a re-run")
                    ui.button("Change study", icon="swap_horiz",
                              on_click=lambda: ui.navigate.to("/study")).props("flat")

    def _refresh(self) -> None:
        studies.forget_runs()
        reporttables.forget_cached()
        self.built, self.stale = {}, {}
        self.redraw()

    # -- runs --------------------------------------------------------------

    def _runs_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Runs").classes("text-lg font-bold")
                ui.button("Add runs on the Simulations page", icon="playlist_add",
                          on_click=lambda: ui.navigate.to("/")) \
                    .props("flat dense no-caps").mark("add-run")
            if not self.study.runs:
                ui.label("No runs yet. Open each sims_config.json the report draws on "
                         "on the Simulations page and add it to the study - the RFSL "
                         "and FSL lists are separate runs."
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
                              on_click=lambda _, n=entry["name"]: self._ask_remove_run(n)) \
                        .props("flat dense round").tooltip("Remove from the study") \
                        .mark(f"remove-run-{entry['name']}")

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

    def _ask_remove_run(self, name) -> None:
        readers = [spec.get("title") or spec.get("id") for spec in self.study.tables
                   if any(source.get("run") == name
                          for source in studies.table_sources(spec))]
        detail = ("Its sims list and results are not touched. "
                  + (f"These tables read it, and would be left naming a run that is not "
                     f"there: {', '.join(readers)}." if readers else ""))
        confirm(f"Remove the run {name} from the study?", detail,
                lambda: self._remove_run(name))

    def _remove_run(self, name) -> None:
        orphans = self.study.remove_run(name)
        self.save()
        if orphans:
            ui.notify(f"Removed {name}; these tables still name it: "
                      f"{', '.join(orphans)}", type="warning", multi_line=True)
        self.redraw()

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
        # Tables open folded: one line each, built only when opened - most of the
        # page's loading time was spent reading runs for tables nobody was reading.
        self.open = self.open & {spec.get("id") for spec in self.study.tables}
        self._contents()
        for spec in self.study.tables:
            self._table_card(spec)
        named = address.param("table")           # /report?table=... opens and shows it
        if named in self.cards and not self._went_to:
            self._went_to = True
            ui.timer(0.3, lambda: self._go_to(named), once=True)

    def _contents(self) -> None:
        with ui.card().classes("w-full gap-1").mark("table-contents"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label("Contents").classes("font-bold")
                with ui.row().classes("gap-1"):
                    ui.button("Open all", icon="unfold_more",
                              on_click=lambda: self._open_all(True)) \
                        .props("flat dense no-caps").mark("open-all-tables")
                    ui.button("Fold all", icon="unfold_less",
                              on_click=lambda: self._open_all(False)) \
                        .props("flat dense no-caps").mark("fold-all-tables")
            for number, spec in enumerate(self.study.tables, start=1):
                table_id = spec.get("id")
                # A label, not a link: a link's own navigation would fight the scroll.
                ui.label(f"{number}. {spec.get('title') or '(untitled)'}") \
                    .classes("text-sm text-primary cursor-pointer hover:underline") \
                    .on("click", lambda _, t=table_id: self._go_to(t)) \
                    .mark(f"contents-{table_id}")

    def _table_card(self, spec: dict) -> None:
        table_id = spec.get("id")
        kind = reporttables.KINDS.get(spec.get("kind"))
        is_open = table_id in self.open
        with ui.card().classes("w-full gap-2").mark(f"table-{table_id}") as card:
            with ui.row().classes("w-full items-center justify-between no-wrap"):
                with ui.row().classes("items-center gap-2 no-wrap grow cursor-pointer") \
                        .on("click", lambda: self._set_open(table_id, table_id not in self.open)) \
                        .mark(f"fold-{table_id}"):
                    chevron = ui.icon("expand_more" if is_open else "chevron_right") \
                        .classes("text-2xl text-muted").mark(f"chevron-{table_id}")
                    with ui.column().classes("gap-0"):
                        ui.label(spec.get("title") or "(untitled)").classes("font-bold")
                        ui.label(_describe(spec, kind)).classes("text-xs text-muted")
                status = ui.row().classes("items-center no-wrap").mark(f"status-{table_id}")
                with ui.row().classes("gap-1 no-wrap"):
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
                        .props("flat dense round").tooltip("Remove") \
                        .mark(f"remove-{table_id}")
            body = ui.column().classes("w-full gap-2")
            body.set_visibility(is_open)
        self.cards[table_id] = {"card": card, "chevron": chevron, "status": status,
                                "body": body, "drawn": False}
        self._draw_status(table_id)
        if is_open:
            self._draw_body(table_id)

    def _draw_status(self, table_id) -> None:
        parts = self.cards.get(table_id)
        if parts is None or parts["status"].is_deleted:
            return
        parts["status"].clear()
        built = self.built.get(table_id)
        with parts["status"]:
            if built is None:
                ui.label("not built yet" if table_id not in self.open else "building...") \
                    .classes("text-xs text-muted")
                return
            if self.stale.get(table_id):
                ui.chip("results out of date", icon="update") \
                    .props("color=warning text-color=white dense square") \
                    .tooltip("\n".join(self.stale[table_id]))
            if built.problems:
                count = len(built.problems)
                ui.chip(f"{count} problem{'' if count == 1 else 's'} to read",
                        icon="warning").props("color=warning text-color=white dense square")
            else:
                ui.chip("built", icon="check_circle") \
                    .props("color=positive text-color=white dense square")

    def _draw_body(self, table_id) -> None:
        """The Copy buttons and the preview, made the first time a table is opened."""
        parts = self.cards[table_id]
        if parts["drawn"]:
            return
        parts["drawn"] = True
        with parts["body"]:
            with ui.row().classes("gap-1 no-wrap"):
                ui.button("Copy for Word", icon="content_copy",
                          on_click=lambda: self._copy(table_id, rich=True)) \
                    .props("dense no-caps").mark(f"copy-word-{table_id}")
                ui.button("Copy as text", icon="notes",
                          on_click=lambda: self._copy(table_id, rich=False)) \
                    .props("flat dense no-caps").mark(f"copy-text-{table_id}") \
                    .tooltip("Tab-separated, for Excel")
            preview = ui.column().classes("w-full gap-1")
            with preview:
                ui.spinner(size="sm")
        self._draw_status(table_id)
        ui.timer(0.01, lambda: self._fill(table_id, preview), once=True)

    def _set_open(self, table_id, on: bool, *, remember=True) -> None:
        parts = self.cards.get(table_id)
        if parts is None:
            return
        if on:
            self.open.add(table_id)
        else:
            self.open.discard(table_id)
        parts["chevron"].name = "expand_more" if on else "chevron_right"
        parts["body"].set_visibility(on)
        if on:
            self._draw_body(table_id)
        if remember:
            STATE.settings.remember_open_tables(self.study.path, self.open)
            # The address names the table last opened, so it can be sent.
            if on:
                address.keep(table=table_id)
            elif address.param("table") == table_id:
                address.keep(table="")

    def _open_all(self, on: bool) -> None:
        for table_id in list(self.cards):
            self._set_open(table_id, on, remember=False)
        STATE.settings.remember_open_tables(self.study.path, self.open)

    def _go_to(self, table_id) -> None:
        self._set_open(table_id, True)
        card = self.cards[table_id]["card"]
        ui.run_javascript(f"document.getElementById('c{card.id}')"
                          f".scrollIntoView({{behavior: 'smooth', block: 'start'}})")

    async def _fill(self, table_id, box) -> None:
        spec = self.study.table(table_id)
        if spec is None:
            return
        built = await _off_thread(reporttables.build, self.study, spec)
        self.built[table_id] = built
        self.stale[table_id] = await _off_thread(
            staleness.for_sources, self.study, studies.table_sources(spec)) or []
        self._draw_status(table_id)
        if box.is_deleted:
            return
        box.clear()
        with box:
            if self.stale[table_id]:
                severity_banner("warn", "\n".join(self.stale[table_id]),
                                "Re-run the group, then Refresh, before copying this "
                                "table into the report.")
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
        await copy_table(built, rich=rich, warnings=self.stale.get(table_id, []))

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

        def remove() -> None:
            self.study.remove_table(table_id)
            self.built.pop(table_id, None)
            self.save()
            self.redraw()

        confirm(f"Remove '{spec.get('title') or table_id}' from the study?",
                "The results are not touched.", remove)

    def store(self, spec: dict) -> None:
        stored = self.study.put_table(spec)
        self.built.pop(stored["id"], None)
        self.open.add(stored["id"])                 # show what was just made or changed
        STATE.settings.remember_open_tables(self.study.path, self.open)
        self.save()
        self.redraw()


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
            elif self.kind.key == reporttables.REPRESENTATIVE:
                self._representative()
            elif self.kind.key == reporttables.FREQUENT:
                self._frequent()
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

        # Table 1's two extra rows, each placed in AEP order among the rest
        dcf = self.spec.get("dcf")
        with ui.row().classes("w-full items-center gap-2 no-wrap"):
            ui.checkbox("Dam crest flood row", value=bool(dcf),
                        on_change=lambda e: self._toggle(
                            "dcf", {"level": 219.13, "label": "DCF"} if e.value else None)) \
                .mark("dcf-row")
            if dcf:
                ui.number("at level (m AHD)", value=dcf.get("level"), format="%.2f") \
                    .classes("w-40").props("dense") \
                    .on_value_change(lambda e: dcf.update(level=e.value))
                ui.input("Label", value=dcf.get("label", "DCF")).classes("w-24") \
                    .props("dense").on_value_change(lambda e: dcf.update(label=e.value))
                ui.label("AEP from the critical duration's realisations"
                         ).classes("text-xs text-muted")
        pmf = self.spec.get("pmf")
        ui.checkbox("PMF row, at the notional AEP adopted on the PMF page", value=bool(pmf),
                    on_change=lambda e: self._toggle(
                        "pmf", {"run": self._default_run(), "group": "", "label": "PMF"}
                        if e.value else None)).mark("pmf-row")
        if pmf:
            self._source_picker(pmf, mark="pmf")

    def _toggle(self, key, value) -> None:
        self.spec[key] = value
        self.draw()

    # -- representative events ---------------------------------------------

    def _representative(self) -> None:
        spec = self.spec
        ui.label("One section per group: its events are the ones saved on the Events "
                 "page for that group.").classes("text-xs text-muted")
        with ui.row().classes("w-full gap-4 no-wrap items-start"):
            ui.textarea("Trigger names for level loadings, one 'label = level' per line",
                        value=reporttables.levels_text(spec.get("triggers"))) \
                .classes("grow").props("dense autogrow").mark("triggers") \
                .on_value_change(lambda e: spec.update(
                    triggers=reporttables.parse_levels(e.value)))
            with ui.column().classes("gap-1"):
                ui.input("AEP of the PMP (1 in x)", value=f"{spec.get('pmp_aep') or ''}") \
                    .classes("w-48").props("dense") \
                    .on_value_change(lambda e: spec.update(
                        pmp_aep=(reporttables.parse_aeps(e.value) or [None])[0]))
                ui.radio({reporttables.MCDF: "Level AEPs from the realisations",
                          reporttables.ENVELOPE: "Level AEPs off the design curve"},
                         value=spec.get("method", reporttables.MCDF),
                         on_change=lambda e: spec.update(method=e.value)).props("dense")
        sections = spec.setdefault("sections", [])
        for index, section in enumerate(sections):
            with ui.card().classes("w-full").props("flat bordered"):
                with ui.row().classes("w-full items-center gap-2 no-wrap"):
                    ui.input("Section heading", value=section.get("heading", "")) \
                        .classes("grow").props("dense") \
                        .on_value_change(lambda e, s=section: s.update(heading=e.value))
                    ui.button(icon="delete", on_click=lambda _, i=index: self._drop(
                        sections, i)).props("flat dense round")
                self._source_picker(section, mark=f"events-{index}")
                pmf = section.get("pmf")
                ui.checkbox("PMF row, from the highest event of an ensemble group",
                            value=bool(pmf),
                            on_change=lambda e, s=section: self._set_pmf(s, e.value))
                if pmf:
                    self._source_picker(pmf, mark=f"events-pmf-{index}")
        ui.button("Add section", icon="add", on_click=lambda: self._append(
            sections, {"heading": "", "run": self._default_run(), "group": "",
                       "pmf": None})).props("outline dense no-caps").mark("add-section")

    # -- frequent levels ----------------------------------------------------

    def _frequent(self) -> None:
        spec = self.spec
        with ui.row().classes("w-full gap-4 no-wrap items-start"):
            ui.textarea("Frequencies, one 'label = AEP (1 in x)' per line",
                        value=reporttables.levels_text(
                            [{"label": f.get("label", ""), "level": f.get("aep")}
                             for f in spec.get("frequencies") or []])) \
                .classes("grow").props("dense autogrow").mark("frequencies") \
                .on_value_change(lambda e: spec.update(frequencies=[
                    {"label": item["label"], "aep": item["level"]}
                    for item in reporttables.parse_levels(e.value)]))
            with ui.column().classes("gap-1 w-80"):
                ui.input("Column headings, separated by |",
                         value=" | ".join(spec.get("columns") or [])) \
                    .classes("w-full").props("dense") \
                    .on("blur", lambda e: self._set_columns(e.sender.value))
                ui.input("Durations to consider (h), blank for all",
                         value=" ".join(f"{h:g}" for h in spec.get("durations") or [])) \
                    .classes("w-full").props("dense") \
                    .on_value_change(lambda e: spec.update(
                        durations=[float(token) for token in str(e.value).replace(",", " ").split()
                                   if token.replace(".", "", 1).isdigit()]))
        columns = spec.get("columns") or []
        sections = spec.setdefault("sections", [])
        for index, section in enumerate(sections):
            groups = section.setdefault("groups", [])
            while len(groups) < len(columns):
                groups.append({"run": self._default_run(), "group": ""})
            with ui.card().classes("w-full").props("flat bordered"):
                with ui.row().classes("w-full items-center gap-2 no-wrap"):
                    ui.input("Section heading", value=section.get("heading", "")) \
                        .classes("grow").props("dense") \
                        .on_value_change(lambda e, s=section: s.update(heading=e.value))
                    ui.button(icon="delete", on_click=lambda _, i=index: self._drop(
                        sections, i)).props("flat dense round")
                for position, column in enumerate(columns):
                    with ui.row().classes("w-full items-center gap-2 no-wrap"):
                        ui.label(column).classes("w-28 text-sm")
                        with ui.element("div").classes("grow"):
                            self._source_picker(groups[position],
                                                mark=f"frequent-{index}-{position}")
        ui.button("Add section", icon="add", on_click=lambda: self._append(
            sections, {"heading": "", "groups": []})).props("outline dense no-caps") \
            .mark("add-section")

    def _set_columns(self, text) -> None:
        columns = [part.strip() for part in str(text).split("|") if part.strip()]
        if columns != self.spec.get("columns"):
            self.spec["columns"] = columns
            self.draw()

    def _set_pmf(self, section, on) -> None:
        section["pmf"] = ({"run": self._default_run(), "group": "", "label": "PMF"}
                          if on else None)
        self.draw()

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
            ui.radio({reporttables.MCDF: "From the critical duration's realisations "
                                         "(the mcdf) - the finer answer",
                      reporttables.ENVELOPE: "Off the design curve (the envelope, between "
                                             "standard AEPs)"},
                     value=self.spec.get("method", reporttables.MCDF),
                     on_change=lambda e: self.spec.update(method=e.value)) \
                .props("dense").mark("level-method")
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
