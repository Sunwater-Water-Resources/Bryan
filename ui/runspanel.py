"""The runs panel down the left of every page: the study's runs, one click apart.

Swapping E013 RFSL for E013 FSL on the Results page means clicking the other
run here: it is opened and the page redrawn, without leaving it. Each run shows
how many of its rows are in each state (core/runstatus.py), worst first; sums
that are not kept yet are worked out after the page is drawn, so a slow disk
never holds up a page.

Open, the panel takes 300 px down the left; folded, a strip. It takes its room
rather than floating over the page, so it never hides the page's left edge.
Which it is is remembered, and the run name in the menu bar toggles it.
"""

from __future__ import annotations

from pathlib import Path

from nicegui import run, ui

from core import runstatus
from core.config import ConfigError
from core.palette import SURFACE, WATER, WATER_SOFT
from layout import STATE_COLOUR
from state import STATE

WIDTH, STRIP = 300, 44


class RunsPanel:
    def __init__(self) -> None:
        self.drawer = None
        self.slots: dict = {}          # run name -> the element its status is drawn in

    # -- open or folded ------------------------------------------------------

    def _fold(self) -> None:
        if STATE.settings.runs_panel_open:
            self.drawer.props(remove="mini")
        else:
            self.drawer.props("mini")

    def toggle(self) -> None:
        STATE.settings.runs_panel_open = not STATE.settings.runs_panel_open
        STATE.settings.save()
        self._fold()

    # -- drawing ---------------------------------------------------------------

    def draw(self) -> None:
        self.drawer = ui.left_drawer(value=True, fixed=True, bordered=True) \
            .props(f"behavior=desktop width={WIDTH} mini-width={STRIP}") \
            .style(f"background:{SURFACE}").classes("p-0").mark("runs-panel")
        self._fold()
        with self.drawer:
            # Quasar hides these by drawer state; a plain div, since a column's own
            # display:flex would win over its display:none.
            with ui.element("div").classes("q-mini-drawer-only w-full text-center pt-2"):
                ui.button(icon="view_list", on_click=self.toggle) \
                    .props("flat round dense color=primary").mark("unfold-runs") \
                    .tooltip("Show the runs")
            with ui.element("div").classes("q-mini-drawer-hide w-full"), \
                    ui.column().classes("w-full gap-1 p-3"):
                with ui.row().classes("w-full items-center justify-between no-wrap"):
                    ui.label("Runs").classes("text-lg font-bold")
                    with ui.row().classes("gap-0 no-wrap"):
                        if STATE.study is not None and STATE.study.runs:
                            ui.button(icon="refresh", on_click=self._refresh) \
                                .props("flat round dense").mark("refresh-runs") \
                                .tooltip("Judge every run's rows again")
                        ui.button(icon="chevron_left", on_click=self.toggle) \
                            .props("flat round dense").mark("fold-runs").tooltip("Fold")
                if STATE.study is None:
                    self._no_study()
                else:
                    self._study()
                self._others()

    def _no_study(self) -> None:
        ui.label("No study is open.").classes("text-sm text-muted")
        ui.link("Open one on the Study page", "/study").classes("text-sm")
        project = STATE.project
        if project is not None:
            self._run_entry(project.name, project.config.config_path, is_open=True,
                            status=None)

    def _study(self) -> None:
        study = STATE.study
        ui.label(study.name or study.path.stem).classes("text-sm font-bold text-ink") \
            .mark("panel-study")
        with ui.link(target="/study").classes("text-sm no-underline"):
            with ui.row().classes("items-center gap-1 no-wrap"):
                ui.icon("tune").classes("text-base")
                ui.label("Dam inputs")
        current = _resolved(STATE.project.config.config_path) if STATE.project else None
        missing = []
        with ui.column().classes("w-full gap-1 mt-2"):
            if not study.runs:
                ui.label("The study has no runs yet.").classes("text-sm text-muted")
            for entry in study.runs:
                name = entry["name"]
                path = study.run_config_path(name)
                status = runstatus.cached(study, name)
                if status is None:
                    missing.append(name)
                self._run_entry(name, path, is_open=path is not None
                                and _resolved(path) == current, status=status)
        if missing:
            ui.timer(0.05, lambda: self._fill(missing), once=True)

    def _run_entry(self, name, path, *, is_open, status) -> None:
        look = (f"background:{WATER_SOFT}; border-left:3px solid {WATER}" if is_open
                else "border-left:3px solid transparent")
        found = path is not None and Path(path).is_file()
        entry = ui.column().classes("w-full gap-0 px-2 py-1 rounded-sm" +
                                    (" cursor-pointer hover:bg-gray-100" if found
                                     and not is_open else "")) \
            .style(look).mark(f"panel-run-{name}")
        if found and not is_open:
            entry.on("click", lambda _, p=path: self._load(p))
        with entry:
            ui.label(name).classes("text-sm " + ("font-bold text-ink" if is_open
                                                 else "text-body"))
            self.slots[name] = ui.row().classes("items-center gap-2 no-wrap")
            if not found:
                with self.slots[name]:
                    ui.label("sims_config.json not found").classes("text-xs text-negative")
            elif status is not None:
                self._draw_status(name, status)
            elif STATE.study is not None:
                with self.slots[name]:
                    ui.spinner(size="xs")

    def _draw_status(self, name, status) -> None:
        slot = self.slots.get(name)
        if slot is None:
            return
        slot.clear()
        with slot:
            if status.problem or not status.total:
                ui.label(status.describe()).classes("text-xs text-muted")
                return
            for state, count in status.in_order():
                colour, icon = STATE_COLOUR.get(state, ("muted", "help_outline"))
                with ui.row().classes("items-center gap-0 no-wrap") \
                        .tooltip(f"{count} {state}"):
                    ui.icon(icon).classes(f"text-{colour} text-sm")
                    ui.label(str(count)).classes("text-xs text-body")

    async def _fill(self, names) -> None:
        study = STATE.study
        for name in names:
            if study is None or study is not STATE.study:
                return
            status = await run.io_bound(runstatus.status_of, study, name)
            self._draw_status(name, status)

    def _refresh(self) -> None:
        runstatus.forget()
        ui.navigate.reload()

    def _others(self) -> None:
        with ui.column().classes("w-full gap-1 mt-3"):
            project, study = STATE.project, STATE.study
            if project is not None and study is not None \
                    and study.run_named_for(project.config.config_path) is None:
                ui.label(f"Open, not in the study: {project.name}") \
                    .classes("text-xs text-body")
                ui.link("Add it to the study", "/").classes("text-sm") \
                    .mark("panel-add-to-study")
            ui.link("Open another sims_config.json...", "/").classes("text-sm") \
                .mark("panel-open-other")

    # -- loading a run ---------------------------------------------------------

    def _load(self, path) -> None:
        try:
            project = STATE.open_project(path)
        except (ConfigError, FileNotFoundError, OSError) as exc:
            ui.notify(str(exc), type="negative", timeout=0, close_button=True)
            return
        ui.notify(f"Opened {project.name}", type="positive")
        ui.navigate.reload()


def _resolved(path) -> Path:
    try:
        return Path(path).resolve()
    except OSError:
        return Path(path)


def runs_panel() -> RunsPanel:
    panel = RunsPanel()
    panel.draw()
    return panel
