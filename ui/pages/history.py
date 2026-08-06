"""Past runs: what was launched, what it did, and tidying up afterwards."""

from __future__ import annotations

from datetime import datetime

from nicegui import ui

from core import completion, runstate
from layout import page_frame, require_project
from state import STATE


def history_page() -> None:
    with page_frame("History"):
        project = require_project()
        if project is None:
            return

        STATE.manager.reattach(project.config.run_root)
        _runs_card(project)
        _log_card(project)
        _maintenance_card(project)


def _runs_card(project) -> None:
    records = STATE.runs()
    with ui.card().classes("w-full"):
        ui.label("Runs launched from here").classes("text-lg font-bold")
        if not records:
            ui.label("None yet.").classes("text-gray-500 text-sm")
            return

        columns = [
            {"name": "run_id", "label": "Run", "field": "run_id"},
            {"name": "created", "label": "Started", "field": "created"},
            {"name": "status", "label": "Status", "field": "status"},
            {"name": "chunks", "label": "Chunks", "field": "chunks"},
            {"name": "sims", "label": "Simulations", "field": "sims"},
            {"name": "folder", "label": "Folder", "field": "folder"},
        ]
        rows = [{
            "run_id": record.run_id,
            "created": datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M"),
            "status": record.status,
            "chunks": len(record.chunks),
            "sims": sum(len(chunk.rows) for chunk in record.chunks),
            "folder": str(record.folder),
        } for record in records]

        table = ui.table(columns=columns, rows=rows, row_key="run_id"
                         ).classes("w-full").props("dense flat bordered")

        def open_run(event) -> None:
            STATE.active_run_id = event.args[1]["run_id"]
            ui.navigate.to("/run")

        table.on("rowClick", open_run)
        ui.label("Every run folder holds the exact sims list and config that "
                 "ran, plus launch.bat and launch.sh to repeat it without the "
                 "UI.").classes("text-xs text-gray-500")


def _log_card(project) -> None:
    """When each output last ran, from Bryan's own run logs.

    Provenance only. 'Simulation' in a run log is the Output file with no
    duration (lib/RunLog.py:36), so where several rows share a name - twenty
    four of them in the Callide list - it cannot say which one ran.
    """
    history = completion.last_run_from_logs(project.config)
    with ui.card().classes("w-full"):
        ui.label("From the run logs").classes("text-lg font-bold")
        ui.label("The last recorded run of each output name, across the "
                 "master log and every run folder.").classes("text-sm text-gray-600")
        if not history:
            ui.label("No run logs found yet.").classes("text-gray-500 text-sm")
            return

        columns = [
            {"name": "name", "label": "Output file", "field": "name", "sortable": True},
            {"name": "ended", "label": "Finished", "field": "ended", "sortable": True},
            {"name": "status", "label": "Status", "field": "status"},
            {"name": "computer", "label": "Computer", "field": "computer"},
        ]
        rows = [{"name": name, "ended": entry[0], "status": entry[1],
                 "computer": entry[2]}
                for name, entry in sorted(history.items())]
        ui.table(columns=columns, rows=rows, row_key="name"
                 ).classes("w-full").props("dense flat bordered")
        ui.label("A name here is not proof a particular row ran - several rows "
                 "can share one output name. The Select page decides from the "
                 "result files on disk.").classes("text-xs text-gray-500")


def _maintenance_card(project) -> None:
    with ui.card().classes("w-full"):
        ui.label("Tidy up").classes("text-lg font-bold")
        ui.label("Each launch leaves a folder under _ui_runs. Delete the "
                 "finished ones once they are no longer interesting - the "
                 "results themselves live elsewhere and are untouched."
                 ).classes("text-sm text-gray-600")
        with ui.row().classes("items-center gap-3"):
            days = ui.number("Older than (days)", value=30, min=1, max=3650
                             ).classes("w-40")

            def prune() -> None:
                removed = runstate.prune(project.config.run_root,
                                         older_than_days=float(days.value or 30))
                if removed:
                    for run_id in removed:
                        STATE.manager.runs.pop(run_id, None)
                    ui.notify(f"Deleted {len(removed)} run folder(s)",
                              type="positive")
                else:
                    ui.notify("Nothing old enough to delete")

            ui.button("Delete old run folders", on_click=prune).props("flat")
