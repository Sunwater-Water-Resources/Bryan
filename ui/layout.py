"""Shared page chrome: the header, the open-project chip, and status colours."""

from __future__ import annotations

from contextlib import contextmanager

from nicegui import ui

from core import completion
from state import STATE

NAV = [
    ("Project", "/"),
    ("Select", "/select"),
    ("Run", "/run"),
    ("Edit", "/edit"),
    ("History", "/history"),
]

STATE_COLOUR = {
    completion.NOT_RUN: ("grey-6", "radio_button_unchecked"),
    completion.UP_TO_DATE: ("positive", "check_circle"),
    completion.STALE: ("warning", "update"),
    completion.INCOMPLETE: ("orange-9", "error_outline"),
    completion.NEEDS_PRIOR: ("negative", "report_problem"),
    completion.UNKNOWN: ("grey-5", "help_outline"),
}


@contextmanager
def page_frame(title: str):
    with ui.header().classes("items-center justify-between bg-primary"):
        with ui.row().classes("items-center gap-1"):
            ui.label("Bryan").classes("text-lg font-bold q-mr-md")
            for label, target in NAV:
                ui.button(label, on_click=lambda _, t=target: ui.navigate.to(t)) \
                    .props("flat color=white dense no-caps") \
                    .classes("font-bold" if label == title else "")
        project = STATE.project
        ui.label(project.name if project else "no project open") \
            .classes("text-white text-sm")
    with ui.column().classes("w-full max-w-7xl mx-auto p-4 gap-4"):
        yield


def require_project():
    """The empty state for pages that need a project open."""
    if STATE.project is not None:
        return STATE.project
    with ui.card().classes("w-full items-center p-8"):
        ui.icon("folder_open").classes("text-5xl text-gray-400")
        ui.label("Open a sims_config.json first.").classes("text-gray-500")
        ui.button("Go to Project", on_click=lambda: ui.navigate.to("/"))
    return None


def status_chip(state) -> None:
    """A completion state as a coloured chip with its reason on hover."""
    if state is None:
        ui.label("-").classes("text-gray-400")
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
                    ui.label(hint).classes("text-sm text-gray-600")
