"""Editing the sims list - reviewed, then saved to a new file."""

from __future__ import annotations

from nicegui import ui

from core import editing
from core.columns import METHODS
from core.paths import cell_text
from layout import page_frame, require_project, severity_banner

EDITABLE_HINT = (
    "Editing here never touches the master workbook. Review the changes, then "
    "save to a new file."
)


def edit_page() -> None:
    with page_frame("Edit"):
        project = require_project()
        if project is None:
            return
        _EditView(project).build()


class _EditView:
    def __init__(self, project) -> None:
        self.project = project
        self.edited = project.frame.copy()
        self.grid = None
        self.changes_box = None

    def build(self) -> None:
        severity_banner("info", EDITABLE_HINT, editing.SAVE_NOTE)
        self._toolbar()
        self._grid()
        self._changes()

    def _toolbar(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-2 flex-wrap"):
                ui.button("Add a Group column", icon="label",
                          on_click=self._add_groups).props("flat")
                ui.button("Duplicate the selected rows", icon="content_copy",
                          on_click=self._duplicate).props("flat")
                ui.button("Revert", icon="undo",
                          on_click=self._revert).props("flat")
                ui.space()
                ui.button("Save as...", icon="save",
                          on_click=self._save_dialog).props("color=primary")
            ui.label("'Include' is what Bryan reads from the file. Which rows "
                     "the UI will run is a separate thing, set on the Select "
                     "page - a run copy always writes Include = yes."
                     ).classes("text-xs text-gray-500")

    def _grid(self) -> None:
        columns = [
            {"headerName": name, "field": name, "editable": True,
             "sortable": True, "filter": True, "resizable": True,
             "cellEditor": ("agSelectCellEditor" if name in ("Include", "Method")
                            else None),
             "cellEditorParams": (
                 {"values": ["yes", "no"]} if name == "Include"
                 else {"values": list(METHODS)} if name == "Method"
                 else None),
             }
            for name in self.edited.columns
        ]
        self.grid = ui.aggrid({
            "columnDefs": columns,
            "rowData": self._row_data(),
            "rowSelection": "multiple",
            "defaultColDef": {"minWidth": 110},
            "stopEditingWhenCellsLoseFocus": True,
        }).classes("w-full h-96")
        self.grid.on("cellValueChanged", self._on_cell)

    def _row_data(self) -> list:
        return [
            {name: cell_text(self.edited.loc[index, name])
             for name in self.edited.columns}
            | {"_index": int(index)}
            for index in self.edited.index
        ]

    def _on_cell(self, event) -> None:
        data = event.args.get("data", {})
        field = event.args.get("colId")
        index = data.get("_index")
        if index is None or field is None:
            return
        value = data.get(field)
        # Not a plain .loc assignment: a column that is blank throughout reads
        # as float64, and pandas 3 refuses to take a string into it.
        editing.set_cell(self.edited, index, field,
                         None if value == "" else value)
        self._changes_refresh()

    def _add_groups(self) -> None:
        self.edited = editing.add_group_column(self.edited)
        ui.notify("Group column added - save to a new file to keep it")
        self._rebuild()

    async def _duplicate(self) -> None:
        rows = await self.grid.get_selected_rows()
        indices = [row["_index"] for row in rows if "_index" in row]
        if not indices:
            ui.notify("Select some rows in the grid first", type="warning")
            return
        try:
            self.edited = editing.duplicate_rows(self.edited, indices)
        except editing.EditError as exc:
            ui.notify(str(exc), type="warning")
            return
        ui.notify(f"Copied {len(indices)} row(s) to the end")
        self._rebuild()

    def _revert(self) -> None:
        self.edited = self.project.frame.copy()
        ui.notify("Reverted to the file on disk")
        self._rebuild()

    def _rebuild(self) -> None:
        self.grid.options["columnDefs"] = [
            {"headerName": name, "field": name, "editable": True,
             "sortable": True, "filter": True, "resizable": True}
            for name in self.edited.columns
        ]
        self.grid.options["rowData"] = self._row_data()
        self.grid.update()
        self._changes_refresh()

    def _changes(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Changes").classes("text-lg font-bold")
            self.changes_box = ui.column().classes("w-full gap-1")
        self._changes_refresh()

    def _changes_refresh(self) -> None:
        self.changes_box.clear()
        changes = editing.diff_frames(self.project.frame, self.edited)
        added = len(self.edited) - len(self.project.frame)
        with self.changes_box:
            if not changes and added <= 0:
                ui.label("Nothing changed.").classes("text-gray-500 text-sm")
                return
            if added > 0:
                ui.label(f"{added} row(s) added").classes("text-sm font-bold")
            for change in changes[:40]:
                ui.label(f"row {change.row + 2}  {change.column}: "
                         f"{change.before or '(blank)'} -> "
                         f"{change.after or '(blank)'}").classes("text-sm font-mono")
            if len(changes) > 40:
                ui.label(f"... and {len(changes) - 40} more").classes("text-sm")

    def _save_dialog(self) -> None:
        master = self.project.sims.path
        with ui.dialog() as dialog, ui.card().classes("min-w-[32rem]"):
            ui.label("Save a copy").classes("text-lg font-bold")
            ui.label(f"Beside {master.parent}").classes("text-xs text-gray-500")
            name = ui.input("Filename",
                            value=editing.default_save_name(master)).classes("w-full")
            replace = ui.checkbox("Replace it if it already exists")
            severity_banner("warn", editing.SAVE_NOTE)

            def save() -> None:
                destination = master.parent / name.value.strip()
                try:
                    editing.save_as(self.project.sims, self.edited, destination,
                                    overwrite=replace.value)
                except editing.EditError as exc:
                    ui.notify(str(exc), type="negative", timeout=0,
                              close_button=True)
                    return
                dialog.close()
                ui.notify(f"Wrote {destination.name}", type="positive")

            with ui.row().classes("justify-end gap-2 w-full"):
                ui.button("Cancel", on_click=dialog.close).props("flat")
                ui.button("Save", on_click=save).props("color=primary")
        dialog.open()
