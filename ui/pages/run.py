"""Watching a run: per-simulation progress, live logs, and stopping it."""

from __future__ import annotations

import time
from pathlib import Path

from nicegui import ui

from core import progress, runstate
from core.logtail import tail
from layout import page_frame, require_project, severity_banner
from state import STATE

STATUS_COLOUR = {
    progress.COMPLETED: "positive",
    progress.FAILED: "negative",
    progress.FAILED_MISSING_INPUT: "negative",
    progress.RUNNING: "primary",
    progress.CANCELLED: "grey-7",
    progress.PENDING: "grey-5",
    progress.UNKNOWN: "grey-5",
}


class FakeChunk:
    """The shape core/progress.read_chunk_progress wants, from a ChunkRecord."""

    def __init__(self, record) -> None:
        self.index = record.index
        self.rows = tuple(record.rows)
        self.run_log = Path(record.run_log)
        self.console_log = Path(record.console_log)


def run_page() -> None:
    with page_frame("Run"):
        project = require_project()
        if project is None:
            return

        STATE.manager.reattach(project.config.run_root)
        record = STATE.active_run()
        if record is None:
            _empty_state()
            return

        _header(record)
        body = ui.column().classes("w-full gap-4")

        def render() -> None:
            STATE.manager.poll()
            body.clear()
            with body:
                for chunk in record.chunks:
                    _chunk_card(record, chunk)

        render()
        ui.timer(STATE.settings.poll_seconds, render)


def _empty_state() -> None:
    with ui.card().classes("w-full items-center p-8"):
        ui.icon("play_circle_outline").classes("text-5xl text-gray-400")
        ui.label("Nothing running.").classes("text-gray-500")
        ui.button("Choose simulations", on_click=lambda: ui.navigate.to("/select"))


def _header(record) -> None:
    with ui.card().classes("w-full"):
        with ui.row().classes("items-center justify-between w-full"):
            with ui.column().classes("gap-0"):
                ui.label(f"Run {record.run_id}").classes("text-lg font-bold")
                ui.label(str(record.folder)).classes("text-xs text-gray-500")
            with ui.row().classes("gap-2"):
                ui.button("Stop after the current chunk",
                          on_click=lambda: _request_stop(record)).props("flat")
                ui.button("Stop now", icon="stop",
                          on_click=lambda: _confirm_cancel(record)
                          ).props("color=negative flat")

        if record.stop_requested:
            severity_banner("info", "A stop has been requested - nothing "
                                    "further will be started.")


def _request_stop(record) -> None:
    STATE.manager.request_stop(record.run_id)
    ui.notify("Queued chunks cancelled; what is running will finish.")


def _confirm_cancel(record) -> None:
    with ui.dialog() as dialog, ui.card():
        ui.label("Stop now?").classes("text-lg font-bold")
        severity_banner(
            "warn",
            "This kills Bryan and the URBS/RORB processes underneath it.",
            "The simulation in flight never reaches the run log, so it will "
            "have no entry there at all. Its model working folder and its log "
            "file are left part-written.",
        )
        with ui.row().classes("justify-end gap-2 w-full"):
            ui.button("Keep running", on_click=dialog.close).props("flat")

            def kill() -> None:
                dialog.close()
                STATE.manager.cancel(record.run_id)
                ui.notify("Stopped", type="warning")

            ui.button("Stop now", on_click=kill).props("color=negative")
    dialog.open()


def _chunk_card(record, chunk) -> None:
    result = progress.read_chunk_progress(
        FakeChunk(chunk),
        log_paths=[Path(p) if p else None for p in chunk.log_paths],
        names=list(chunk.names),
        returncode=chunk.returncode,
        alive=STATE.manager.is_alive(record.run_id, chunk.index),
    )

    with ui.card().classes("w-full"):
        with ui.row().classes("items-center justify-between w-full"):
            ui.label(f"Chunk {chunk.index}").classes("font-bold")
            with ui.row().classes("items-center gap-3"):
                ui.label(f"{result.completed}/{result.total}").classes("text-sm")
                if chunk.started:
                    ended = chunk.ended or time.time()
                    ui.label(progress.elapsed_text(ended - chunk.started)
                             ).classes("text-sm text-gray-600")
                ui.chip(chunk.status).props(
                    f"color={_chunk_colour(chunk.status)} text-color=white dense square")

        ui.linear_progress(value=result.fraction, show_value=False).props("rounded")

        if chunk.note:
            ui.label(chunk.note).classes("text-sm text-gray-600")
        if chunk.returncode is not None:
            ui.label(progress.explain_returncode(chunk.returncode)
                     ).classes("text-xs text-gray-500")

        current = result.current
        if current is not None:
            stalled = progress.stalled_for(current)
            detail = f"running {current.position}/{result.total}: {current.name}"
            if current.log_bytes:
                detail += (f"  -  log {current.log_bytes:,} bytes, last written "
                           f"{progress.elapsed_text(stalled)} ago")
            ui.label(detail).classes("text-sm")

        _sim_table(result)

        with ui.expansion("Console output").classes("w-full"):
            text = tail(chunk.console_log).text or "(nothing yet)"
            ui.code(text[-8000:]).classes("w-full text-xs")


def _chunk_colour(status) -> str:
    return {
        runstate.COMPLETED: "positive",
        runstate.RUNNING: "primary",
        runstate.QUEUED: "grey-6",
        runstate.CANCELLED: "grey-7",
        runstate.FAILED: "negative",
        runstate.DIED: "negative",
    }.get(status, "grey-5")


def _sim_table(result) -> None:
    columns = [
        {"name": "position", "label": "#", "field": "position"},
        {"name": "name", "label": "Simulation", "field": "name"},
        {"name": "status", "label": "Status", "field": "status"},
        {"name": "started", "label": "Started", "field": "started"},
        {"name": "ended", "label": "Ended", "field": "ended"},
        {"name": "error", "label": "Error", "field": "error",
         "classes": "max-w-md truncate"},
    ]
    rows = [
        {"position": sim.position, "name": sim.name, "status": sim.status,
         "started": sim.started, "ended": sim.ended, "error": sim.error}
        for sim in result.sims
    ]
    ui.table(columns=columns, rows=rows, row_key="position"
             ).classes("w-full").props("dense flat")

    for sim in result.sims:
        explanation = progress.explain_error(sim.error)
        if explanation:
            severity_banner("warn", f"{sim.name}: {sim.error}", explanation)
