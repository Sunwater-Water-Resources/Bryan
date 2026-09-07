"""Choosing a representative event for each design flood loading.

You give it loadings - a design AEP, or a lake level off the frequency curve -
and it ranks the realisations of the Monte Carlo database against each one. The
ranking is closeness to the loading in standard normal variate space, on both
axes at once: the flood should be as rare as the loading asks, and the rainfall
that produced it should be about as rare as the flood. What the table adds is
everything that would make an event indefensible even so - an embedded burst,
a pre-burst or an antecedent lake level far off the median.

Nothing is excluded quietly. Flags are shown against every candidate and the
filters that actually drop events are opt-in, because the choice between the
closest match and the cleanest storm is the user's, and a tool that made it
silently would not be trusted twice.

The page holds the chosen list, saves it beside the results, and exports it.
Extracting the hydrographs for the chosen events is the util script's job.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from nicegui import ui

from core import events, eventchart, results
from layout import page_frame, require_project, severity_banner
from state import STATE

TYPE_LABELS = {"level": "Lake level", "inflow": "Peak inflow",
               "outflow": "Peak outflow"}

KIND_LABELS = {"aep": "AEP", "level": "Lake level"}

# What "closest" means. Neutrality is usually secondary to reaching the loading
# - an event exists to put the lake at a level - so both are offered and the
# one not ranked on is still reported in the table.
ORDER_LABELS = {events.EVENTS.DELTA_Z: "\u0394z (AEP neutral)",
                events.EVENTS.RESULT: "Closest result"}

# Enough to see the shape of the cloud around the target without turning the
# table into a second database.
DEFAULT_COUNT = 10

# How long a number field waits after the last keystroke before it reports the
# value. Re-ranking an mcdf on every digit is wasted work, and "1000" would be
# evaluated as 1, 10, 100 and 1000 on the way past.
TYPING_PAUSE_MS = 500


def events_page() -> None:
    with page_frame("Events"):
        project = require_project()
        if project is None:
            return
        available = events.sources_by_group(project)
        if not available:
            _empty_state()
            return
        _EventsView(project, available).build()


def _empty_state() -> None:
    with ui.card().classes("w-full items-center p-8 gap-2"):
        ui.icon("scatter_plot", size="3rem").classes("text-gray-400")
        ui.label("No Monte Carlo database found for this sims list."
                 ).classes("text-gray-500")
        ui.label(
            "A representative event is one realisation, so this page needs the "
            "mcdf a monte carlo or reservoir routing row writes - not just the "
            "quantile tables. An ensemble run has no realisations to choose "
            "between: it ran every combination by design."
        ).classes("text-xs text-gray-500 max-w-lg text-center")
        ui.button("Choose simulations", on_click=lambda: ui.navigate.to("/select"))


class _EventsView:
    def __init__(self, project, available) -> None:
        self.project = project
        self.available = available            # group -> [EventSource]
        self.group = next(iter(available), None)
        self.result_type = "level"
        self.order = events.EVENTS.DELTA_Z
        # How close counts as the same result, per result type and in the units
        # the field asks for (millimetres of lake level, m3/s of flow).
        self.bands = dict(events.DEFAULT_BANDS)
        self.targets: list = []
        self.filters = events.Filters()
        self.outcomes: list = []
        self.folder = None
        self._curve = None                    # the level envelope, per group
        self._curves: dict = {}               # design envelopes by result type
        # Which loading cards are open. Every redraw rebuilds the expansions,
        # so without this the first one springs open and the one being worked
        # on shuts every time an event is picked.
        self.open_cards: set[int] = {0}

        self.band_input = None
        self.band_units = None
        self.target_box = None
        self.detail_box = None
        self.summary_box = None
        self.command_box = None
        self.status = None

    # -- build ------------------------------------------------------------

    def build(self) -> None:
        self._controls()
        self._targets_card()
        self._details_card()
        self._summary_card()
        self._load_group(self.group)

    def _controls(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.select(list(self.available), value=self.group, label="Group",
                          on_change=lambda event: self._load_group(event.value)
                          ).classes("min-w-96")
                ui.toggle(TYPE_LABELS, value=self.result_type,
                          on_change=self._on_type).props("no-caps dense") \
                    .mark("result-type")
                ui.toggle(ORDER_LABELS, value=self.order,
                          on_change=self._on_order).props("no-caps dense") \
                    .mark("rank-order") \
                    .tooltip("\u0394z ranks on reaching the loading and being "
                             "AEP neutral about it, at once. Closest result "
                             "ranks on the loading itself - the level in "
                             "metres, or the design value at that AEP - and "
                             "leaves neutrality to be read off the table.")
                self.band_input = ui.number(
                    "Same result within", value=self.bands[self.result_type],
                    format="%g", step=5,
                    on_change=lambda e: self._on_band(e.value)) \
                    .classes("w-44").props(f"debounce={TYPING_PAUSE_MS}") \
                    .mark("result-band") \
                    .tooltip("Only used by 'Closest result'. Events this close "
                             "to the loading count as reaching it, and are then "
                             "ordered by AEP neutrality - a lake level is not "
                             "meaningful to the millimetre, and without a band "
                             "the sort is decided by noise.")
                self.band_units = ui.label(
                    events.band_units(self.result_type)[0]
                ).classes("text-xs text-gray-500")
                ui.button("Reload", icon="refresh", on_click=self._reload
                          ).props("flat dense")
            ui.separator()
            with ui.row().classes("items-center gap-4 flex-wrap"):
                ui.number("AEP of the PMP (1 in X)",
                          value=self.filters.aep_of_pmp, format="%.0f",
                          on_change=lambda e: self._on_filter("aep_of_pmp", e.value)
                          ).classes("w-56") \
                    .props(f"clearable debounce={TYPING_PAUSE_MS}") \
                    .tooltip("Rainfall rarer than this cannot be sampled, so "
                             "events near it are the edge of the scheme. Read "
                             "from the IFD files config the storm config points at.")
                ui.number("Furthest Δz to offer", value=self.filters.max_delta_z,
                          format="%.2f", step=0.1,
                          on_change=lambda e: self._on_filter("max_delta_z", e.value)
                          ).classes("w-48") \
                    .props(f"clearable debounce={TYPING_PAUSE_MS}")
                ui.checkbox("Drop events with an embedded burst",
                            value=self.filters.exclude_embedded,
                            on_change=lambda e: self._on_filter("exclude_embedded",
                                                                e.value))
                ui.checkbox("Drop anything flagged",
                            value=self.filters.exclude_flagged,
                            on_change=lambda e: self._on_filter("exclude_flagged",
                                                                e.value))
            self.status = ui.label().classes("text-xs text-gray-500")

    def _targets_card(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-2"):
                ui.label("Loadings").classes("font-bold")
                ui.space()
                ui.button("Add loading", icon="add", on_click=self._add
                          ).props("flat dense")
                ui.button("Save", icon="save", on_click=self._save).props("flat dense")
                ui.button("Export csv", icon="download", on_click=self._export
                          ).props("flat dense")
            self.target_box = ui.column().classes("w-full gap-1")

    def _details_card(self) -> None:
        self.detail_box = ui.column().classes("w-full gap-2")

    def _summary_card(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Representative events").classes("font-bold")
            ui.label("'hydrograph' is the column to pull out of the stored "
                     "inflow, level and outflow files."
                     ).classes("text-xs text-gray-500")
            self.summary_box = ui.column().classes("w-full")
            with ui.expansion("Extract the hydrographs and plot them").classes("w-full"):
                ui.label(
                    "Save first, then run this with Bryan's interpreter. It "
                    "writes a three-panel plot per event - the rebuilt "
                    "hyetograph, the inflow and outflow, and the lake level - "
                    "and a workbook of the series."
                ).classes("text-xs text-gray-500")
                self.command_box = ui.column().classes("w-full")

    # -- state ------------------------------------------------------------

    def _sources(self) -> list:
        return self.available.get(self.group, [])

    def _load_group(self, group) -> None:
        self.group = group
        self._curve = None
        self._curves = {}
        sources = self._sources()
        self.folder = events.default_folder(sources)

        saved, settings = ([], {})
        if self.folder is not None:
            saved, settings = events.load_targets(
                events.selection_path(self.folder, group))
        # Before the defaults are built, so they are made for the right type.
        if settings.get("result_type") in TYPE_LABELS:
            self.result_type = settings["result_type"]
        if settings.get("order") in ORDER_LABELS:
            self.order = settings["order"]
        saved_bands = settings.get("bands")
        if isinstance(saved_bands, dict):
            self.bands = {**events.DEFAULT_BANDS,
                          **{key: float(value) for key, value in saved_bands.items()
                             if key in events.DEFAULT_BANDS}}
        if self.band_input is not None:
            self.band_input.value = self.bands.get(self.result_type, 0.0)
        if self.band_units is not None:
            self.band_units.set_text(events.band_units(self.result_type)[0])
        self.targets = saved or [
            events.Target(kind="aep", value=aep, result_type=self.result_type,
                          count=DEFAULT_COUNT)
            for aep in (100, 1000, 10_000)
        ]
        # A different group is a different set of loadings, so start again
        # with the first one open.
        self.open_cards = {0}
        self.filters = events.Filters(
            aep_of_pmp=settings.get("aep_of_pmp") or self._pmp_aep(),
            max_delta_z=settings.get("max_delta_z"),
            exclude_embedded=bool(settings.get("exclude_embedded")),
            exclude_flagged=bool(settings.get("exclude_flagged")),
        )
        self.refresh(redraw_targets=True)

    def _pmp_aep(self):
        """From the storm config chain, as a Monte Carlo run gets it."""
        storm = self.project.config.filepaths.get("storm_config")
        return events.EVENTS.pmp_aep_from_storm_config(storm) if storm else None

    def _curve_for_levels(self):
        if self._curve is None:
            self._curve = events.level_curve(self.project, self._sources())
        return self._curve

    def refresh(self, redraw_targets: bool = False) -> None:
        """Re-rank and redraw.

        `redraw_targets` rebuilds the loading rows, which is what adding or
        removing one needs - and what editing a value must *not* do: clearing
        the column destroys the input being typed into, so the field would
        take the first digit and lose focus.
        """
        sources = self._sources()
        self.outcomes = [
            events.evaluate(self.project, sources, target, self.filters,
                            curve=self._curve_for_levels(), order=self.order,
                            curves=self._curves, band=self._band())
            for target in self.targets
        ]
        if redraw_targets:
            self._draw_targets()
        self._draw_details()
        self._draw_summary()
        self._draw_command()
        if self.status is not None:
            files = ", ".join(sorted({source.path.name for source in sources}))
            self.status.set_text(f"{len(sources)} database(s): {files}")

    # -- drawing ----------------------------------------------------------

    def _draw_targets(self) -> None:
        self.target_box.clear()
        labels = [source.label for source in self._sources()]
        with self.target_box:
            if not self.targets:
                ui.label("No loadings yet - add one.").classes("text-gray-500 text-sm")
            for position, target in enumerate(self.targets):
                with ui.row().classes("items-center gap-2 flex-wrap"):
                    ui.toggle(KIND_LABELS, value=target.kind,
                              on_change=lambda e, i=position: self._edit(i, "kind", e.value)
                              ).props("no-caps dense")
                    ui.number("Value", value=target.value, format="%g",
                              on_change=lambda e, i=position: self._edit(i, "value", e.value)
                              ).classes("w-32").props(f"debounce={TYPING_PAUSE_MS}") \
                        .mark(f"loading-value-{position}")
                    ui.select([""] + labels, value=target.source, label="From",
                              on_change=lambda e, i=position: self._edit(i, "source", e.value)
                              ).classes("w-40") \
                        .tooltip("Blank uses the duration that is critical at "
                                 "this loading's AEP.")
                    ui.number("Rain AEP", value=target.rain_aep, format="%g",
                              on_change=lambda e, i=position: self._edit(i, "rain_aep", e.value)
                              ).classes("w-32") \
                        .props(f"clearable debounce={TYPING_PAUSE_MS}") \
                        .tooltip("Blank judges the rainfall against the loading "
                                 "itself, which is the AEP-neutral case.")
                    ui.number("Show", value=target.count, format="%d",
                              on_change=lambda e, i=position: self._edit(i, "count", e.value)
                              ).classes("w-24").props(f"debounce={TYPING_PAUSE_MS}")
                    ui.button(icon="delete", on_click=lambda _, i=position: self._remove(i)
                              ).props("flat dense round")

    def _draw_details(self) -> None:
        self.detail_box.clear()
        with self.detail_box:
            for position, outcome in enumerate(self.outcomes):
                self._draw_outcome(position, outcome)

    def _draw_outcome(self, position, outcome) -> None:
        target = outcome.target
        picked = outcome.picked
        headline = f"{target.label}"
        if outcome.aep:
            headline += f"  -  1 in {results.format_aep(outcome.aep)}"
        if outcome.source is not None:
            headline += f"  -  from {outcome.source.label}"
        if picked is not None:
            headline += f"  -  sim {int(picked.name)}"

        with ui.expansion(headline, value=position in self.open_cards,
                          on_value_change=lambda e, i=position:
                              self._card_toggled(i, e.value)).classes("w-full") \
                .mark(f"target-{position}"):
            for note in outcome.notes:
                ui.label(note).classes("text-xs text-gray-500")
            if outcome.problem:
                severity_banner("warn", outcome.problem)
                return
            if outcome.ranking is not None and outcome.ranking.excluded:
                ui.label("Left out: " + ", ".join(
                    f"{count} {reason}"
                    for reason, count in outcome.ranking.excluded.items())
                ).classes("text-xs text-gray-500")

            rows = events.candidate_rows(outcome)
            if not rows:
                severity_banner("warn", "No candidate events survive the filters.")
                return
            self._candidate_table(position, rows)
            chart = ui.echart(eventchart.neutrality_chart(
                outcome, self.result_type)).classes("w-full h-96")
            chart.mark(f"neutrality-{position}")

    def _candidate_table(self, position, rows) -> None:
        columns = [{"name": "pick", "label": "", "field": "pick", "align": "center"}] + [
            {"name": name, "label": label, "field": name, "sortable": True}
            for name, label in events.CANDIDATE_COLUMNS
        ]
        table = ui.table(columns=columns, rows=rows, row_key="sim") \
            .classes("w-full").props("dense flat bordered") \
            .mark(f"candidates-{position}")
        table.add_slot("body-cell-pick", r"""
            <q-td :props="props">
              <q-radio :model-value="props.row.picked" :val="true"
                       @update:model-value="() => $parent.$emit('pick', props.row)" />
            </q-td>
        """)
        table.on("pick", lambda event, i=position: self._pick(i, event.args))

    def _card_toggled(self, position, is_open) -> None:
        """Remember an open card, and nothing else - no redraw.

        Redrawing here would destroy the expansion that raised the event.
        """
        if is_open:
            self.open_cards.add(position)
        else:
            self.open_cards.discard(position)

    def _draw_command(self) -> None:
        if self.command_box is None:
            return
        self.command_box.clear()
        if self.folder is None:
            return
        path = events.selection_path(self.folder, self.group)
        with self.command_box:
            ui.code(events.extract_command(self.project, path)).classes("w-full text-xs")

    def _draw_summary(self) -> None:
        self.summary_box.clear()
        rows = events.summary_rows(self.outcomes)
        with self.summary_box:
            if not rows:
                ui.label("Nothing chosen yet.").classes("text-gray-500 text-sm")
                return
            columns = [{"name": name, "label": name, "field": name}
                       for name in rows[0]]
            ui.table(columns=columns, rows=rows, row_key="loading") \
                .classes("w-full").props("dense flat bordered") \
                .mark("summary-table")

    # -- events -----------------------------------------------------------

    def _band(self) -> float:
        """The band in the result's own units - metres, not millimetres."""
        return events.band_in_result_units(self.result_type,
                                           self.bands.get(self.result_type, 0.0))

    def _on_type(self, event) -> None:
        self.result_type = event.value
        for target in self.targets:
            target.result_type = self.result_type
        # The band is per result type: 20 mm of lake level is not 20 m3/s.
        if self.band_input is not None:
            self.band_input.value = self.bands.get(self.result_type, 0.0)
        if self.band_units is not None:
            self.band_units.set_text(events.band_units(self.result_type)[0])
        self.refresh()

    def _on_band(self, value) -> None:
        if value in (None, ""):
            return                            # mid-edit, as the loading rows are
        self.bands[self.result_type] = max(float(value), 0.0)
        self.refresh()

    def _on_order(self, event) -> None:
        self.order = event.value
        self.refresh()

    def _on_filter(self, name, value) -> None:
        if name in ("aep_of_pmp", "max_delta_z"):
            value = float(value) if value not in (None, "") else None
        self.filters = events.Filters(**{
            **{field: getattr(self.filters, field)
               for field in self.filters.__dataclass_fields__},
            name: value,
        })
        self.refresh()

    def _edit(self, position, name, value) -> None:
        if position >= len(self.targets):
            return
        target = self.targets[position]
        if name in ("value", "count") and value in (None, ""):
            # An empty box is mid-edit, not a loading of zero: the row is not
            # redrawn while it is being typed into, so leave the target alone
            # and take the value when there is one.
            return
        if name in ("value", "rain_aep"):
            value = float(value) if value not in (None, "") else None
        if name == "count":
            value = max(int(value), 1)
        if getattr(target, name) == value:
            return
        setattr(target, name, value)
        if name in ("kind", "value", "source"):
            # A different loading is a different event - do not carry the pick.
            target.picked = None
        self.refresh()

    def _add(self) -> None:
        self.targets.append(events.Target(kind="aep", value=100,
                                          result_type=self.result_type,
                                          count=DEFAULT_COUNT))
        # Open the one just added: it is what the user is about to work on.
        self.open_cards.add(len(self.targets) - 1)
        self.refresh(redraw_targets=True)

    def _remove(self, position) -> None:
        if position < len(self.targets):
            self.targets.pop(position)
        # The cards are keyed by position, so the ones above a deletion move
        # down with their loading.
        self.open_cards = {index if index < position else index - 1
                           for index in self.open_cards if index != position}
        self.refresh(redraw_targets=True)

    def _pick(self, position, row) -> None:
        if position < len(self.targets) and isinstance(row, dict):
            self.targets[position].picked = int(row.get("sim"))
        self.refresh()

    def _reload(self) -> None:
        events.forget_cached()
        STATE.reload_project()
        ui.navigate.to("/events")

    def _settings(self) -> dict:
        return {"result_type": self.result_type,
                "order": self.order,
                "bands": dict(self.bands),
                "aep_of_pmp": self.filters.aep_of_pmp,
                "max_delta_z": self.filters.max_delta_z,
                "exclude_embedded": self.filters.exclude_embedded,
                "exclude_flagged": self.filters.exclude_flagged}

    def _save(self) -> None:
        if self.folder is None:
            ui.notify("Nowhere to save - no database folder", type="warning")
            return
        # Save what is on screen, pick included, so reopening shows the same
        # list - and where each event came from, which is what lets the util
        # script find the database and the stored hydrographs without being
        # told a second time.
        for target, outcome in zip(self.targets, self.outcomes):
            if outcome.source is not None:
                target.output_file = outcome.source.output_name
                target.database = events.project_relative(self.project,
                                                          outcome.source.path)
            if target.picked is None and outcome.picked_id is not None:
                target.picked = int(outcome.picked_id)
        path = events.selection_path(self.folder, self.group)
        events.save_targets(path, self.targets, self._settings())
        ui.notify(f"Saved {path}")

    def _export(self) -> None:
        rows = events.summary_rows(self.outcomes)
        if not rows or self.folder is None:
            ui.notify("Nothing to export", type="warning")
            return
        path = Path(events.selection_path(self.folder, self.group)).with_suffix(".csv")
        pd.DataFrame(rows).to_csv(path, index=False)
        ui.notify(f"Wrote {path}")
