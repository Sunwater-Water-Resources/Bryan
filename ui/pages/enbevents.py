"""The Events page for an ensemble group: the PMF and design AEPs, chosen by rule.

The Monte Carlo view ranks realisations; an ensemble has none to rank, so this view has
no ranking controls at all. Each loading is the PMF (the highest event) or a design AEP
(the median pattern at the critical duration), and its card shows the per-duration spread
the rule chose from. The selection is saved as the page goes, as the Monte Carlo view's is,
with ``"method": "ensemble"`` so the Downstream page treats it as one.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import pandas as pd
from nicegui import ui

import address
from core import enbevents, ensemble, events, results
from layout import severity_banner
from state import STATE

KIND_LABELS = {enbevents.PMF: "PMF", enbevents.AEP: "Design AEP"}
TYPE_LABELS = {"level": "Lake level", "inflow": "Peak inflow", "outflow": "Peak outflow"}
CRITICAL = "critical"
TYPING_PAUSE_MS = 500
CHOSEN_ROW_JS = "row => row.chosen ? 'bg-water-soft text-weight-medium' : ''"


def group_options(monte_carlo, ensemble) -> dict:
    """Every group the page can open, labelled so an ensemble one is told apart."""
    options = {key: key for key in monte_carlo}
    options.update({key: f"{key}  (ensemble)" for key in ensemble if key not in monte_carlo})
    return options


def switch_group(group) -> None:
    """A group of the other kind is a different view: open the page on it afresh."""
    ui.navigate.to(f"/events?group={group}")


class EnsembleEventsView:
    def __init__(self, project, available, options) -> None:
        self.project = project
        self.available = available             # group -> [enbevents.Source]
        self.options = options
        named = address.param("group")
        self.group = named if named in available else next(iter(available))
        self.result_type = "level"
        self.targets: list = []
        self.choices: list = []
        self.frame = None
        self.problem = ""
        self.folder = None
        self.open_cards: set[int] = {0}
        self.target_box = self.detail_box = self.summary_box = self.saved_label = None

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.select(self.options, value=self.group, label="Group",
                          on_change=lambda e: self._group(e.value)).classes("min-w-96")
                ui.toggle(TYPE_LABELS, value=self.result_type,
                          on_change=lambda e: self._type(e.value)).props("no-caps dense") \
                    .mark("result-type")
                ui.button("Reload", icon="refresh", on_click=self._reload).props("flat dense")
            ui.label("An ensemble ran every pattern of every duration, so there is nothing to "
                     "rank. The PMF is the highest event; a design AEP is "
                     "the median pattern at the critical duration."
                     ).classes("text-xs text-muted").mark("ensemble-rules")
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-2"):
                ui.label("Loadings").classes("font-bold")
                ui.space()
                ui.button("Add loading", icon="add", on_click=self._add).props("flat dense")
                self.saved_label = ui.label("").classes("text-xs text-muted").mark("saved-state")
                self.save_button = ui.button("Save", icon="save", on_click=self._autosave) \
                    .props("dense no-caps").mark("save-selection")
                ui.button("Export csv", icon="download", on_click=self._export
                          ).props("flat dense")
            self.target_box = ui.column().classes("w-full gap-1")
        self.detail_box = ui.column().classes("w-full gap-2")
        with ui.card().classes("w-full"):
            ui.label("Representative events").classes("font-bold")
            ui.label("'hydrograph' is the column to pull out of the stored inflow, level and "
                     "outflow files.").classes("text-xs text-muted")
            self.summary_box = ui.column().classes("w-full")
        self._load_group()

    # -- state ------------------------------------------------------------

    def _sources(self) -> list:
        return self.available.get(self.group, [])

    def _load_group(self) -> None:
        address.keep(group=self.group)
        sources = self._sources()
        self.folder = sources[0].path.parent if sources else None
        self.problem = ""
        try:
            self.frame = enbevents.load(sources)
        except enbevents.EnsembleEventError as exc:
            self.frame, self.problem = None, str(exc)
        saved = events.load_selection(events.selection_path(self.folder, self.group)) \
            if self.folder is not None else None
        if saved is not None and saved.settings.get("result_type") in TYPE_LABELS:
            self.result_type = saved.settings["result_type"]
        if saved is not None and saved.targets and saved.method == events.EVENTS.ENSEMBLE:
            self.targets = saved.targets
        elif self.frame is not None:
            self.targets = enbevents.default_targets(self.frame, self.result_type)
            adopted = self._adopted_aep()
            for target in self.targets:
                if target.kind == enbevents.PMF and adopted:
                    target.notional_aep = adopted
        else:
            self.targets = []
        self.refresh(redraw_targets=True)
        self._show_saved_state()

    @staticmethod
    def _adopted_aep():
        """The notional AEP of the PMF adopted on the PMF page, where a study is open."""
        if STATE.study is None:
            return None
        try:
            return ensemble.settings(STATE.study).get("adopted_aep") or None
        except Exception:                            # noqa: BLE001 - a default, not a must
            return None

    def refresh(self, redraw_targets: bool = False) -> None:
        if self.frame is None:
            self.choices = []
        else:
            for target in self.targets:
                target.result_type = self.result_type
            self.choices = [enbevents.choose(self.frame, target, self._sources())
                            for target in self.targets]
        if redraw_targets:
            self._draw_targets()
        self._draw_details()
        self._draw_summary()

    # -- drawing ----------------------------------------------------------

    def _draw_targets(self) -> None:
        self.target_box.clear()
        aep_options = ({aep: f"1 in {results.format_aep(aep)}"
                        for aep in enbevents.aeps(self.frame)} if self.frame is not None else {})
        duration_options = {CRITICAL: "critical"}
        if self.frame is not None:
            duration_options.update({f"{d:g}": f"{d:g} h" for d in enbevents.durations(self.frame)})
        with self.target_box:
            if not self.targets:
                ui.label("No loadings yet - add one.").classes("text-muted text-sm")
            for position, target in enumerate(self.targets):
                with ui.row().classes("items-center gap-2 flex-wrap"):
                    ui.input("Name (optional)", value=target.name, placeholder="e.g. PMF",
                             on_change=lambda e, i=position: self._edit(i, "name", e.value)
                             ).classes("w-56").props(f"dense debounce={TYPING_PAUSE_MS}") \
                        .mark(f"loading-name-{position}")
                    ui.toggle(KIND_LABELS, value=target.kind,
                              on_change=lambda e, i=position: self._edit(i, "kind", e.value)
                              ).props("no-caps dense").mark(f"loading-kind-{position}")
                    if target.kind == enbevents.PMF:
                        ui.number("Notional AEP (1 in X)", value=target.notional_aep,
                                  format="%.0f",
                                  on_change=lambda e, i=position: self._edit(
                                      i, "notional_aep", e.value)
                                  ).classes("w-48").props(f"clearable debounce={TYPING_PAUSE_MS}") \
                            .tooltip("A label only. The concurrent rainfall downstream is set "
                                     "by the event's own rainfall AEP - the AEP of the PMP.")
                    else:
                        value = float(target.value) if target.value in aep_options else None
                        ui.select(aep_options, value=value, label="AEP",
                                  on_change=lambda e, i=position: self._edit(i, "value", e.value)
                                  ).classes("w-44").mark(f"loading-value-{position}")
                    source = str(target.source or "") or CRITICAL
                    ui.select(duration_options,
                              value=source if source in duration_options else CRITICAL,
                              label="Duration",
                              on_change=lambda e, i=position: self._edit(
                                  i, "source", "" if e.value == CRITICAL else e.value)
                              ).classes("w-32") \
                        .tooltip("Critical takes the rule's own duration; a duration "
                                 "overrides it.")
                    up = ui.button(icon="arrow_upward",
                                   on_click=lambda _, i=position: self._move(i, -1)
                                   ).props("flat dense round").tooltip("Move up")
                    up.set_enabled(position > 0)
                    down = ui.button(icon="arrow_downward",
                                     on_click=lambda _, i=position: self._move(i, 1)
                                     ).props("flat dense round").tooltip("Move down")
                    down.set_enabled(position < len(self.targets) - 1)
                    ui.button(icon="delete", on_click=lambda _, i=position: self._remove(i)
                              ).props("flat dense round")

    def _draw_details(self) -> None:
        self.detail_box.clear()
        with self.detail_box:
            if self.problem:
                severity_banner("warn", self.problem)
                return
            for position, choice in enumerate(self.choices):
                self._draw_choice(position, choice)

    def _draw_choice(self, position, choice) -> None:
        headline = choice.target.title
        if choice.found:
            headline += f"  -  {choice.duration:g} h  -  event {choice.event}"
        with ui.expansion(headline, value=position in self.open_cards,
                          on_value_change=lambda e, i=position: self._toggled(i, e.value)
                          ).classes("w-full").mark(f"target-{position}"):
            rule = enbevents.RULE_TEXT.get(choice.rule, "")
            if rule:
                override = str(choice.target.source or "")
                ui.label(f"Rule: {rule}" + (f", held at {override} h" if override else "")
                         + ".").classes("text-xs text-muted").mark(f"rule-{position}")
            if choice.problem:
                severity_banner("warn", choice.problem)
                return
            values = choice.values
            said = [f"event {choice.event}", f"{choice.duration:g} h", choice.pattern]
            for name, unit, digits in (("level", "m AHD", 3), ("inflow", "m3/s", 0),
                                       ("outflow", "m3/s", 0)):
                if name in values:
                    said.append(f"{name} {values[name]:,.{digits}f} {unit}")
            if "rain_aep" in values:
                said.append(f"rain AEP 1 in {results.format_aep(values['rain_aep'])}")
            if "ADV" in values:
                said.append(f"ADV {values['ADV']:,.0f} ML")
            ui.label("Chosen: " + ", ".join(said)).classes("text-sm text-body") \
                .mark(f"chosen-{position}")
            rows = enbevents.duration_rows(choice)
            if rows:
                ui.label("Click a heading to sort; click again for highest first."
                         ).classes("text-xs text-muted")
                table = ui.table(columns=enbevents.DURATION_COLUMNS, rows=rows,
                                 row_key="duration").classes("w-full") \
                    .props("dense flat bordered").mark(f"durations-{position}")
                table.props(f":row-class='{CHOSEN_ROW_JS}'")

    def _draw_summary(self) -> None:
        self.summary_box.clear()
        rows = enbevents.summary_rows(self.choices)
        with self.summary_box:
            if not rows:
                ui.label("Nothing chosen yet.").classes("text-muted text-sm")
                return
            ui.table(columns=[{"name": key, "label": key, "field": key} for key in rows[0]],
                     rows=rows, row_key="loading").classes("w-full") \
                .props("dense flat bordered").mark("summary-table")

    # -- events -----------------------------------------------------------

    def _group(self, group) -> None:
        if group in self.available:
            self.group = group
            self.open_cards = {0}
            self._load_group()
        else:
            switch_group(group)

    def _type(self, value) -> None:
        self.result_type = value
        self.refresh()
        self._autosave()

    def _edit(self, position, name, value) -> None:
        if position >= len(self.targets):
            return
        target = self.targets[position]
        if name == "value" and value in (None, ""):
            return
        if name in ("value", "notional_aep"):
            value = None if value in (None, "") else float(value)
        if getattr(target, name) == value:
            return
        if name == "name":
            target.name = str(value or "")
            self._draw_details()
            self._draw_summary()
            self._autosave()
            return
        setattr(target, name, value)
        if name == "kind":
            target.picked = None
            if value == enbevents.AEP and self.frame is not None and not target.value:
                target.value = enbevents.aeps(self.frame)[-1]
            self.refresh(redraw_targets=True)
        else:
            self.refresh()
        self._autosave()

    def _add(self) -> None:
        value = enbevents.aeps(self.frame)[-1] if self.frame is not None else 0.0
        self.targets.append(events.Target(kind=enbevents.AEP, value=value,
                                          result_type=self.result_type))
        self.open_cards.add(len(self.targets) - 1)
        self.refresh(redraw_targets=True)
        self._autosave()

    def _remove(self, position) -> None:
        if position < len(self.targets):
            self.targets.pop(position)
        self.open_cards = {index if index < position else index - 1
                           for index in self.open_cards if index != position}
        self.refresh(redraw_targets=True)
        self._autosave()

    def _move(self, position, step) -> None:
        other = position + step
        if not (0 <= position < len(self.targets) and 0 <= other < len(self.targets)):
            return
        self.targets[position], self.targets[other] = self.targets[other], self.targets[position]
        swap = {position: other, other: position}
        self.open_cards = {swap.get(index, index) for index in self.open_cards}
        self.refresh(redraw_targets=True)
        self._autosave()

    def _toggled(self, position, is_open) -> None:
        (self.open_cards.add if is_open else self.open_cards.discard)(position)

    def _reload(self) -> None:
        enbevents.forget_cached()
        STATE.reload_project()
        ui.navigate.to(f"/events?group={self.group}")

    def _autosave(self) -> None:
        if self.folder is None:
            self._saved("not saved - there is no database folder to save beside")
            return
        relative = lambda path: events.project_relative(self.project, path)  # noqa: E731
        stored = [enbevents.stored_target(choice, relative) for choice in self.choices]
        path = events.selection_path(self.folder, self.group)
        try:
            events.save_targets(path, stored, {"result_type": self.result_type},
                                method=events.EVENTS.ENSEMBLE,
                                sims_config=relative(self.project.config.config_path))
        except OSError as exc:
            self._saved(f"NOT saved: {exc}")
            ui.notify(f"Could not save {Path(path).name}: {exc}", type="warning")
            return
        self._saved(f"saved {time.strftime('%H:%M:%S')} to {Path(path).name}")

    def _show_saved_state(self) -> None:
        """Say whether this group's selection is on disk; until it is, there is a Save button.
        A group left at its defaults has no change to save it, and the Downstream page lists
        only saved selections."""
        path = events.selection_path(self.folder, self.group) if self.folder is not None else None
        if path is not None and Path(path).is_file():
            self._saved(f"saved to {Path(path).name}")
        else:
            self._saved("not saved yet - change anything, or Save")

    def _saved(self, text) -> None:
        if self.saved_label is not None:
            self.saved_label.set_text(text)
        if getattr(self, "save_button", None) is not None:
            self.save_button.set_visibility(not text.startswith("saved"))

    def _export(self) -> None:
        rows = enbevents.summary_rows(self.choices)
        if not rows or self.folder is None:
            ui.notify("Nothing to export", type="warning")
            return
        path = Path(events.selection_path(self.folder, self.group)).with_suffix(".csv")
        pd.DataFrame(rows).to_csv(path, index=False)
        ui.notify(f"Wrote {path}")
