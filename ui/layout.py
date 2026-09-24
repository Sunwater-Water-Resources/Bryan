"""Shared page chrome: the header, the open-project chip, and status colours."""

from __future__ import annotations

from contextlib import contextmanager

from nicegui import ui

import theme
from core import completion
from core.palette import BRAND_CYAN, INK, ON_INK_MUTED
from state import STATE

NAV = [
    ("Project", "/"),
    ("Select", "/select"),
    ("Run", "/run"),
    ("Results", "/results"),
    ("Events", "/events"),
    ("Lake levels", "/lake-levels"),
    ("Downstream", "/downstream"),
    ("Report", "/report"),
    ("PMF", "/pmf"),
    ("Figures", "/figures"),
    ("Lake record", "/lake-record"),
    ("Edit", "/edit"),
    ("History", "/history"),
]

# The house status colours, each with an icon and its own label, so no state is
# told apart by colour alone. Stale and incomplete share the attention ochre: both
# mean the result on disk should not be taken as it stands, and the icon and the
# chip's text say which. White on each is at least 5:1.
STATE_COLOUR = {
    completion.NOT_RUN: ("muted", "radio_button_unchecked"),
    completion.UP_TO_DATE: ("positive", "check_circle"),
    completion.STALE: ("warning", "update"),
    completion.INCOMPLETE: ("warning", "error_outline"),
    completion.NEEDS_PRIOR: ("negative", "report_problem"),
    completion.UNKNOWN: ("muted", "help_outline"),
}


@contextmanager
def page_frame(title: str):
    theme.apply()
    # The same bar Judith's window carries: the name, what the tool is for, and on
    # the right what is open.
    with ui.header().classes("items-center justify-between") \
            .style(f"background:{INK}; padding:8px 18px"):
        with ui.row().classes("items-center gap-1 no-wrap"):
            ui.label("Bryan").style("color:#fff; font-size:17px; font-weight:600")
            ui.label("Design flood hydrology").style(
                f"color:{ON_INK_MUTED}; font-size:13px; padding:0 18px 0 10px")
            for label, target in NAV:
                button = ui.button(label, on_click=lambda _, t=target: ui.navigate.to(t)) \
                    .props("flat color=white dense no-caps")
                if title.startswith(label):
                    button.style(f"border-bottom:2px solid {BRAND_CYAN}; border-radius:0")
        project = STATE.project
        ui.label(project.name if project else "no project open") \
            .classes("mono").style(f"color:{ON_INK_MUTED}; font-size:12px")
    with ui.column().classes("w-full max-w-7xl mx-auto p-4 gap-4"):
        yield


def require_project():
    """The empty state for pages that need a project open."""
    if STATE.project is not None:
        return STATE.project
    with ui.card().classes("w-full items-center p-8"):
        ui.icon("folder_open").classes("text-5xl text-muted")
        ui.label("Open a sims_config.json first.").classes("text-muted")
        ui.button("Go to Project", on_click=lambda: ui.navigate.to("/"))
    return None


def status_chip(state) -> None:
    """A completion state as a coloured chip with its reason on hover."""
    if state is None:
        ui.label("-").classes("text-muted")
        return
    colour, icon = STATE_COLOUR.get(state.state, STATE_COLOUR[completion.UNKNOWN])
    with ui.element("div").classes("inline-flex"):
        chip = ui.chip(state.state, icon=icon).props(
            f"color={colour} text-color=white dense square")
        if state.detail:
            with chip:
                ui.tooltip(state.detail).classes("max-w-md")


def severity_banner(severity: str, message: str, hint: str = "") -> None:
    kinds = {"block": ("negative", "error"),
             "warn": ("warning", "warning"),
             "info": ("info", "info")}
    colour, icon = kinds.get(severity, kinds["info"])
    with ui.card().classes("w-full").props("flat bordered"):
        with ui.row().classes("items-start no-wrap gap-3"):
            ui.icon(icon).classes(f"text-{colour} text-2xl")
            with ui.column().classes("gap-1"):
                ui.label(message).classes("whitespace-pre-wrap")
                if hint:
                    ui.label(hint).classes("text-sm text-body")
