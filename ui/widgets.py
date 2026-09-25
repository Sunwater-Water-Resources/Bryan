"""Path boxes that check themselves, and the Browse dialog behind them.

Every path the launcher asks for is typed into one of these. Under the box is
the full path the text resolves to and whether it is what the box wants - a
tick, *not found*, *is a folder* - and, for a relative path, what it is
relative to (core/pathcheck.py). **Browse...** lists this computer's folders
and the files of the types the box expects; a file picked under the box's base
folder is put in as a relative path, as the study keeps paths.

A box commits its value - calls ``on_commit`` - when it loses focus, on Enter,
and when a file is picked, never on every keystroke: several boxes re-read a
record when their path changes. ``on_change`` is for a box whose value is only
held until a Save button, and is called on every change.
"""

from __future__ import annotations

from pathlib import Path

from nicegui import ui

from core import pathcheck
from core.pathcheck import (EITHER, FILE, FOLDER, OUTPUT,  # noqa: F401 - for pages
                            OUTPUT_FOLDER)

LOOK = {pathcheck.OK: ("check_circle", "text-positive"),
        pathcheck.WRITTEN: ("edit_note", "text-body"),
        pathcheck.WRONG_TYPE: ("warning", "text-warning"),
        pathcheck.MISSING: ("cancel", "text-negative")}


def _status_line(holder, result) -> None:
    holder.clear()
    if result.status == pathcheck.BLANK:
        return
    icon, colour = LOOK[result.status]
    with holder:
        with ui.row().classes("items-start gap-1 no-wrap"):
            ui.icon(icon).classes(f"{colour} text-sm").style("margin-top:1px")
            with ui.row().classes("gap-x-2 gap-y-0"):
                ui.label(result.message).classes(f"text-xs break-all "
                                                 f"{'text-body' if result.is_ok else colour}")
                if result.relative_note:
                    ui.label(f"({result.relative_note})").classes("text-xs text-muted")


class PathBox:
    """One path: the box, its Browse button, and the line saying what it finds."""

    def __init__(self, label, value="", *, base=None, base_name="", expect=FILE,
                 suffixes=(), on_commit=None, on_change=None, mark="", places=(),
                 classes="w-full"):
        self.base, self.base_name = (Path(base) if base else None), base_name
        self.expect, self.suffixes = expect, tuple(suffixes)
        self.on_change = on_change
        self.on_commit, self.places, self.label = on_commit, places, label
        self.committed = value or ""
        with ui.column().classes(f"{classes} gap-0") as self.element:
            with ui.row().classes("w-full items-center gap-1 no-wrap"):
                self.input = ui.input(label, value=value or "").classes("grow").props("dense")
                self.browse = ui.button(icon="folder_open", on_click=self._browse) \
                    .props("flat dense round").tooltip("Browse...")
            self.status = ui.column().classes("w-full gap-0")
        if mark:
            self.input.mark(mark)
            self.browse.mark(f"{mark}-browse")
            self.status.mark(f"{mark}-status")
        self.input.on_value_change(lambda _: self._changed())
        self.input.on("blur", lambda _: self.commit())
        self.input.on("keydown.enter", lambda _: self.commit())
        self.refresh()

    def _changed(self) -> None:
        self.refresh()
        if self.on_change is not None:
            self.on_change(self.value)

    @property
    def value(self) -> str:
        return self.input.value or ""

    def check(self) -> pathcheck.PathCheck:
        return pathcheck.check(self.value, base=self.base, base_name=self.base_name,
                               expect=self.expect, suffixes=self.suffixes)

    def refresh(self) -> None:
        _status_line(self.status, self.check())

    def commit(self) -> None:
        if self.value == self.committed:
            return
        self.committed = self.value
        if self.on_commit is not None:
            self.on_commit(self.value)

    def set_value(self, text) -> None:
        self.input.value = text
        self.refresh()
        self.commit()

    def _browse(self) -> None:
        browse(title=self.label, start=pathcheck.start_folder(self.value, self.base),
               expect={EITHER: FILE, OUTPUT_FOLDER: FOLDER}.get(self.expect, self.expect),
               suffixes=self.suffixes,
               places=_places(self.base, self.base_name, self.places),
               name=Path(pathcheck.clean_path_text(self.value)).name
               if self.expect == OUTPUT else "",
               on_pick=lambda path: self.set_value(pathcheck.stored_text(path, self.base)))


class PathList:
    """Several paths, one per line - a record's gauge exports - each one checked."""

    def __init__(self, label, values=(), *, base=None, base_name="", suffixes=(),
                 on_commit=None, mark="", places=(), classes="w-full"):
        self.base, self.base_name, self.suffixes = (Path(base) if base else None), \
            base_name, tuple(suffixes)
        self.on_commit, self.places, self.label = on_commit, places, label
        self.committed = self.lines(list(values))
        with ui.column().classes(f"{classes} gap-0") as self.element:
            with ui.row().classes("w-full items-start gap-1 no-wrap"):
                self.input = ui.textarea(label, value="\n".join(self.committed)) \
                    .classes("grow").props("dense autogrow")
                self.browse = ui.button(icon="folder_open", on_click=self._browse) \
                    .props("flat dense round").tooltip("Browse... - adds a line")
            self.status = ui.column().classes("w-full gap-0")
        if mark:
            self.input.mark(mark)
            self.browse.mark(f"{mark}-browse")
            self.status.mark(f"{mark}-status")
        self.input.on_value_change(lambda _: self.refresh())
        self.input.on("blur", lambda _: self.commit())
        self.refresh()

    @staticmethod
    def lines(values) -> list:
        return [line.strip() for line in values if str(line).strip()]

    @property
    def values(self) -> list:
        return self.lines((self.input.value or "").splitlines())

    def refresh(self) -> None:
        self.status.clear()
        for text in self.values:
            with self.status:
                holder = ui.column().classes("w-full gap-0")
            _status_line(holder, pathcheck.check(text, base=self.base,
                                                 base_name=self.base_name,
                                                 suffixes=self.suffixes))

    def commit(self) -> None:
        if self.values == self.committed:
            return
        self.committed = self.values
        if self.on_commit is not None:
            self.on_commit(self.values)

    def _browse(self) -> None:
        last = self.values[-1] if self.values else ""

        def add(path) -> None:
            self.input.value = "\n".join(self.values + [pathcheck.stored_text(path, self.base)])
            self.refresh()
            self.commit()

        browse(title=self.label, start=pathcheck.start_folder(last, self.base), expect=FILE,
               suffixes=self.suffixes, places=_places(self.base, self.base_name, self.places),
               on_pick=add)


def path_input(label, value="", **kwargs) -> PathBox:
    return PathBox(label, value, **kwargs)


def path_list(label, values=(), **kwargs) -> PathList:
    return PathList(label, values, **kwargs)


def _places(base, base_name, extra) -> list:
    places = []
    if base is not None:
        places.append(((base_name or "base folder").removeprefix("the ").capitalize(),
                       Path(base)))
    places += [(name, Path(path)) for name, path in extra if path]
    places.append(("Home", Path.home()))
    seen, unique = set(), []
    for name, path in places:
        if str(path) not in seen:
            seen.add(str(path))
            unique.append((name, path))
    return unique


# -- the Browse dialog -------------------------------------------------------------

def browse(*, title, start, expect=FILE, suffixes=(), places=(), on_pick, name="") -> None:
    """List folders and files on this computer, and hand the chosen one to ``on_pick``."""
    state = {"folder": Path(start), "chosen": None, "show_all": False}
    with ui.dialog() as dialog, \
            ui.card().classes("w-[52rem] max-w-full gap-2").mark("browse-dialog"):
        ui.label(f"Choose: {title}").classes("text-lg font-bold")
        with ui.row().classes("w-full items-center gap-1 no-wrap"):
            ui.button(icon="arrow_upward", on_click=lambda: go(state["folder"].parent)) \
                .props("flat dense round").tooltip("Up a folder").mark("browse-up")
            where = ui.input(value=str(state["folder"])).classes("grow").props("dense") \
                .mark("browse-where")
            where.on("keydown.enter", lambda: go(Path(where.value or "")))
        with ui.row().classes("w-full items-center gap-1"):
            for drive in pathcheck.drives():
                ui.button(str(drive).rstrip("\\"), on_click=lambda _, d=drive: go(d)) \
                    .props("flat dense no-caps")
            for label, path in places:
                ui.button(label, icon="bookmark", on_click=lambda _, p=path: go(p)) \
                    .props("flat dense no-caps").mark(f"browse-place-{label}")
        entries = ui.column().classes("w-full gap-0 overflow-y-auto border rounded-sm") \
            .style("height: 22rem").mark("browse-entries")
        with ui.row().classes("w-full items-center justify-between no-wrap"):
            with ui.row().classes("items-center gap-2"):
                if suffixes and expect != FOLDER:
                    ui.switch(f"Show all files, not only {', '.join(suffixes)}",
                              on_change=lambda e: (state.update(show_all=e.value),
                                                   go(state["folder"]))).props("dense")
                note = ui.label("").classes("text-xs text-muted")
            if expect == OUTPUT:
                file_name = ui.input("File name", value=name).classes("w-64") \
                    .props("dense").mark("browse-name")
        chosen_label = ui.label("").classes("text-sm text-body break-all") \
            .mark("browse-chosen")
        with ui.row().classes("w-full justify-end gap-2"):
            ui.button("Cancel", on_click=dialog.close).props("flat")
            pick_label = {FOLDER: "Choose this folder", OUTPUT: "Write here"}.get(expect,
                                                                                 "Choose")
            ui.button(pick_label, on_click=lambda: pick()).props("color=primary") \
                .mark("browse-pick")

    def go(folder) -> None:
        folder = Path(folder)
        try:
            if not folder.is_dir():
                ui.notify(f"Not a folder: {folder}", type="warning")
                return
        except OSError as exc:
            ui.notify(f"{folder} cannot be read: {exc}", type="warning")
            return
        state["folder"], state["chosen"] = folder, None
        where.value = str(folder)
        draw()

    def choose(path) -> None:
        state["chosen"] = path
        chosen_label.text = str(path)
        draw()

    def pick(path=None) -> None:
        if expect == FOLDER:
            target = path or state["folder"]
        elif expect == OUTPUT:
            text = (file_name.value or "").strip()
            if path is not None:
                target = path
            elif text:
                target = state["folder"] / text
            else:
                ui.notify("Give the file a name", type="warning")
                return
        else:
            target = path or state["chosen"]
            if target is None:
                ui.notify("Choose a file", type="warning")
                return
        dialog.close()
        on_pick(Path(target))

    def draw() -> None:
        found = pathcheck.listing(state["folder"], suffixes=() if expect == FOLDER
                                  else suffixes, show_all=state["show_all"])
        entries.clear()
        with entries:
            if found.problem:
                ui.label(found.problem).classes("text-sm text-negative p-2")
            for folder in found.folders:
                with ui.row().classes("w-full items-center gap-2 px-2 py-1 cursor-pointer "
                                      "hover:bg-gray-100 no-wrap") \
                        .on("click", lambda _, f=folder: go(f)) \
                        .mark(f"browse-folder-{folder.name}"):
                    ui.icon("folder").classes("text-attention")
                    ui.label(folder.name).classes("text-sm")
            if expect != FOLDER:
                for file in found.files:
                    selected = state["chosen"] == file
                    row = ui.row().classes("w-full items-center gap-2 px-2 py-1 cursor-pointer "
                                           "no-wrap " + ("bg-water_soft" if selected
                                                         else "hover:bg-gray-100")) \
                        .mark(f"browse-file-{file.name}")
                    row.on("click", lambda _, f=file: choose(f))
                    row.on("dblclick", lambda _, f=file: pick(f))
                    with row:
                        ui.icon("description").classes("text-muted")
                        ui.label(file.name).classes("text-sm")
            if not found.folders and not found.files and not found.problem:
                ui.label("Nothing here of the types this box wants.").classes(
                    "text-sm text-muted p-2")
        parts = []
        if found.hidden_files:
            parts.append(f"{found.hidden_files} other file{'s' if found.hidden_files > 1 else ''}"
                         " not shown")
        if found.truncated:
            parts.append(f"only the first {pathcheck.MAX_LISTED} shown")
        note.text = "; ".join(parts)
        if expect == FOLDER:
            chosen_label.text = str(state["folder"])

    draw()
    dialog.open()
