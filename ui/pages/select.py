"""Choosing what to run: group filter, completion status, pre-flight, launch."""

from __future__ import annotations

from nicegui import ui

from core import completion, preflight
from core.launcher import LaunchError
from core.paths import cell_text
from core.runplan import describe
from layout import STATE_COLOUR, page_frame, require_project, severity_banner
from state import STATE

ALL_GROUPS = "(all groups)"
ALL_STATES = "(any state)"

# Columns worth showing before the rest; everything else follows.
LEAD_COLUMNS = ("Output file", "Duration", "Method", "GWL", "Output suffix",
                "Run models", "Analyse results")


def select_page() -> None:
    with page_frame("Select"):
        project = require_project()
        if project is None:
            return
        _SelectView(project).build()


class _SelectView:
    def __init__(self, project) -> None:
        self.project = project
        self.group_filter = ALL_GROUPS
        self.state_filter = ALL_STATES
        self.search = ""
        self.table = None
        self.summary = None
        self.issues_box = None
        self.groups = project.group_keys()

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        if not STATE.completions:
            STATE.refresh_completion()
        self._banners()
        self._filters()
        self._table()
        self._actions()
        self.refresh()

    def _banners(self) -> None:
        audit = self.project.sims.audit
        if audit.is_damaged:
            severity_banner(
                "block",
                f"{self.project.sims.path.name}: {audit.describe()}",
                "The affected rows are shown greyed out and cannot be "
                "selected.",
            )
        report = self.project.group_report()
        if not report.is_useful:
            severity_banner("warn", f"Grouping: {report.degenerate}")

    def _filters(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.select([ALL_GROUPS] + self.project.groups(),
                          value=ALL_GROUPS, label="Group",
                          on_change=self._on_group).classes("min-w-72")
                ui.select([ALL_STATES, completion.NOT_RUN, completion.STALE,
                           completion.INCOMPLETE, completion.UP_TO_DATE,
                           completion.NEEDS_PRIOR],
                          value=ALL_STATES, label="Status",
                          on_change=self._on_state).classes("min-w-48")
                ui.input("Search output name",
                         on_change=self._on_search).classes("min-w-64").props("clearable")
                ui.button("Re-check results", icon="refresh",
                          on_click=self._recheck).props("flat dense")

            with ui.row().classes("items-center gap-2 flex-wrap"):
                ui.button("Select filtered", on_click=lambda: self._bulk("filtered")
                          ).props("flat dense")
                ui.button("Select not up to date",
                          on_click=lambda: self._bulk("rerun")).props("flat dense")
                ui.button("Deselect completed",
                          on_click=lambda: self._bulk("drop-done")).props("flat dense")
                ui.button("Clear", on_click=lambda: self._bulk("none")).props("flat dense")
                report = self.project.group_report()
                ui.label(
                    f"{report.group_count} group(s)"
                    + ("" if report.from_column
                       else " - derived from the output names, no Group column")
                ).classes("text-xs text-gray-500")

    def _table(self) -> None:
        frame = self.project.frame
        ordered = [c for c in LEAD_COLUMNS if c in frame.columns]
        ordered += [c for c in frame.columns if c not in ordered]
        self.display_columns = ordered[:9]

        columns = [
            {"name": "run", "label": "Run", "field": "run", "align": "center"},
            {"name": "status", "label": "Status", "field": "status"},
            {"name": "group", "label": "Group", "field": "group",
             "sortable": True, "classes": "max-w-xs truncate"},
        ] + [
            {"name": name, "label": name, "field": name, "sortable": True}
            for name in self.display_columns
        ]

        self.table = ui.table(columns=columns, rows=[], row_key="index"
                              ).classes("w-full").props("dense flat bordered")
        self.table.add_slot("body-cell-run", r"""
            <q-td :props="props">
              <q-checkbox :model-value="props.row.run" :disable="props.row.blocked"
                          @update:model-value="() => $parent.$emit('toggle', props.row)" />
            </q-td>
        """)
        # The status is only useful with its reason attached - 'stale' means
        # nothing until you know which input changed.
        self.table.add_slot("body-cell-status", r"""
            <q-td :props="props">
              <q-badge :color="props.row.status_colour" :label="props.row.status" />
              <q-tooltip v-if="props.row.status_detail" class="text-body2"
                         style="max-width: 32rem">
                {{ props.row.status_detail }}
              </q-tooltip>
            </q-td>
        """)
        self.table.on("toggle", self._on_toggle)

    def _actions(self) -> None:
        with ui.card().classes("w-full"):
            self.summary = ui.label().classes("text-sm")
            self.issues_box = ui.column().classes("w-full gap-2")
            with ui.row().classes("items-center gap-2"):
                ui.button("Check and run", icon="play_arrow",
                          on_click=self._confirm).props("color=primary")
                ui.button("Open the Run page",
                          on_click=lambda: ui.navigate.to("/run")).props("flat")

    # -- events -----------------------------------------------------------

    def _on_group(self, event) -> None:
        self.group_filter = event.value
        self.refresh()

    def _on_state(self, event) -> None:
        self.state_filter = event.value
        self.refresh()

    def _on_search(self, event) -> None:
        self.search = (event.value or "").strip().lower()
        self.refresh()

    def _on_toggle(self, event) -> None:
        row = event.args
        index = row.get("index")
        if index is None:
            return
        STATE.toggle(index, not row.get("run"))
        self.refresh()

    def _recheck(self) -> None:
        STATE.refresh_completion(check_truncation=True)
        ui.notify("Re-checked the results on disk")
        self.refresh()

    def _bulk(self, what: str) -> None:
        visible = set(self._visible_rows())
        if what == "filtered":
            STATE.set_selected(set(STATE.selected) | visible)
        elif what == "rerun":
            STATE.set_selected(set(STATE.rows_needing_rerun()) & visible or
                               set(STATE.rows_needing_rerun()))
        elif what == "drop-done":
            STATE.set_selected(set(STATE.selected) - set(STATE.already_done(STATE.selected)))
        else:
            STATE.set_selected(set())
        self.refresh()

    # -- rendering --------------------------------------------------------

    def _visible_rows(self) -> list:
        frame = self.project.frame
        rows = []
        for index in frame.index:
            if self.group_filter != ALL_GROUPS and self.groups.get(index) != self.group_filter:
                continue
            state = STATE.completion_of(index)
            if self.state_filter != ALL_STATES:
                if state is None or state.state != self.state_filter:
                    continue
            if self.search and self.search not in self.project.label(index).lower():
                continue
            rows.append(index)
        return rows

    def refresh(self) -> None:
        frame = self.project.frame
        audit = self.project.sims.audit
        rows = []
        for index in self._visible_rows():
            state = STATE.completion_of(index)
            record = {
                "index": int(index),
                "run": index in STATE.selected,
                "blocked": audit.affects(index),
                "status": state.state if state else "-",
                "status_colour": STATE_COLOUR.get(
                    state.state if state else "", "grey-5")[0],
                "status_detail": state.detail if state else "",
                "group": self.groups.get(index, ""),
            }
            for name in self.display_columns:
                record[name] = cell_text(frame.loc[index].get(name))
            rows.append(record)

        self.table.rows = rows
        self.table.update()

        selected = STATE.selected_in_order()
        done = STATE.already_done(selected)
        text = (f"{len(selected)} of {len(frame)} selected; "
                f"{len(rows)} shown")
        if done:
            text += f"  -  {len(done)} would overwrite existing results"
        self.summary.text = text

        self.issues_box.clear()
        if selected:
            with self.issues_box:
                for issue in STATE.preflight(selected)[:8]:
                    severity_banner(issue.severity, issue.message, issue.fix_hint)

    # -- launching --------------------------------------------------------

    def _confirm(self) -> None:
        selected = STATE.selected_in_order()
        if not selected:
            ui.notify("Nothing selected", type="warning")
            return

        issues = STATE.preflight(selected)
        blocking = preflight.blocking(issues)
        plan = STATE.plan(selected)
        blocking += [h for h in plan.hazards if h.blocks]

        with ui.dialog() as dialog, ui.card().classes("min-w-[40rem]"):
            ui.label("Run these simulations?").classes("text-lg font-bold")

            if blocking:
                ui.label(f"{len(blocking)} problem(s) must be fixed first:"
                         ).classes("text-negative font-bold")
                for issue in blocking[:6]:
                    severity_banner("block", issue.message,
                                    getattr(issue, "fix_hint", ""))
                ui.button("Close", on_click=dialog.close).props("flat")
                dialog.open()
                return

            ui.label(describe(plan, self.project.frame)).classes(
                "text-sm whitespace-pre-wrap font-mono")
            ui.label(f"Estimated wall clock: about "
                     f"{plan.estimated_minutes:.0f} minutes across "
                     f"{len(plan.chunks)} process(es).").classes("text-sm")

            for hazard in plan.hazards:
                if hazard.severity != "block":
                    severity_banner(hazard.severity, hazard.message,
                                    hazard.fix_hint)

            done = STATE.already_done(selected)
            overwrite = None
            if done:
                names = ", ".join(self.project.label(i) for i in done[:6])
                more = f" and {len(done) - 6} more" if len(done) > 6 else ""
                severity_banner(
                    "warn",
                    f"{len(done)} of these have up-to-date results that would "
                    f"be overwritten: {names}{more}.",
                    "Bryan rewrites the results files and deletes the model "
                    "working folder on entry. There is no undo.",
                )
                overwrite = ui.checkbox("Yes, overwrite those results")

            test_mode = ui.checkbox(
                "Test mode - stop each Monte Carlo simulation after a few "
                "realisations")

            def go() -> None:
                if overwrite is not None and not overwrite.value:
                    ui.notify("Tick the overwrite box, or deselect those rows",
                              type="warning")
                    return
                dialog.close()
                try:
                    STATE.launch(plan, test_runs=4 if test_mode.value else None)
                except LaunchError as exc:
                    ui.notify(str(exc), type="negative", timeout=0,
                              close_button=True)
                    return
                ui.notify(f"Started {len(plan.chunks)} process(es)",
                          type="positive")
                ui.navigate.to("/run")

            with ui.row().classes("justify-end gap-2 w-full"):
                ui.button("Cancel", on_click=dialog.close).props("flat")
                ui.button("Run", on_click=go).props("color=primary")

        dialog.open()
