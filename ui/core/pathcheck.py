"""Whether a typed path is there, and what it is read relative to.

Every path box in the launcher shows, under it, the full path the text resolves
to and whether that is what the box wants: a file that is there, a folder that
is there, or an output that will be written. A relative path says what it is
relative to - the study's folder, the sims_config.json's - because that is the
thing people get wrong. The browser behind each box's Browse button lists
folders here too.
"""

from __future__ import annotations

import os
import string
from dataclasses import dataclass
from pathlib import Path

from .paths import clean_path_text, normalise_sep

FILE, FOLDER, OUTPUT, EITHER = "file", "folder", "output", "file or folder"
OUTPUT_FOLDER = "output folder"

# Statuses, worst last: what the box's line says.
BLANK, OK, WRITTEN, WRONG_TYPE, MISSING = "blank", "ok", "written", "wrong type", "missing"

MAX_LISTED = 1000


@dataclass(frozen=True)
class PathCheck:
    text: str
    resolved: Path | None
    status: str
    message: str          # what the box's line says, the full path included
    relative_note: str    # "relative to the study folder (C:\\...)", or ""

    @property
    def is_ok(self) -> bool:
        return self.status in (OK, WRITTEN)


def resolve(text, base=None) -> Path | None:
    """The path ``text`` means: absolute as given, relative against ``base``."""
    text = clean_path_text(text or "")
    if not text:
        return None
    path = Path(normalise_sep(text))
    if not path.is_absolute() and base is not None:
        path = Path(base) / path
    return path


def _suffix_ok(path: Path, suffixes) -> bool:
    return not suffixes or path.suffix.lower() in {s.lower() for s in suffixes}


def check(text, *, base=None, base_name="", expect=FILE, suffixes=()) -> PathCheck:
    """Judge one path. ``base_name`` names ``base`` for the reader: "the study folder"."""
    cleaned = clean_path_text(text or "")
    path = resolve(cleaned, base)
    if path is None:
        return PathCheck("", None, BLANK, "", "")
    note = ""
    if not Path(normalise_sep(cleaned)).is_absolute():
        # The full path is in the message, so the base is named, not repeated.
        note = f"relative to {base_name or 'the launcher working folder'}"
    try:
        is_file, is_dir = path.is_file(), path.is_dir()
    except OSError as exc:                   # a share that does not answer
        return PathCheck(cleaned, path, MISSING, f"cannot be read: {path} ({exc})", note)
    types = ", ".join(suffixes)
    if expect == OUTPUT_FOLDER:
        if is_dir:
            return PathCheck(cleaned, path, WRITTEN, f"{path} - there; written into", note)
        if is_file:
            return PathCheck(cleaned, path, WRONG_TYPE, f"is a file, not a folder: {path}", note)
        return PathCheck(cleaned, path, WRITTEN, f"{path} - will be made", note)
    if expect == FOLDER:
        if is_dir:
            return PathCheck(cleaned, path, OK, str(path), note)
        if is_file:
            return PathCheck(cleaned, path, WRONG_TYPE, f"is a file, not a folder: {path}", note)
        return PathCheck(cleaned, path, MISSING, f"not found: {path}", note)
    if expect == EITHER and (is_dir or is_file):
        return PathCheck(cleaned, path, OK, str(path), note)
    if is_dir:
        return PathCheck(cleaned, path, WRONG_TYPE, f"is a folder, not a file: {path}", note)
    if expect == OUTPUT:
        if is_file:
            return PathCheck(cleaned, path, WRITTEN, f"{path} - there; it will be replaced", note)
        parent = "" if path.parent.is_dir() else f"; its folder {path.parent} will be made"
        return PathCheck(cleaned, path, WRITTEN, f"{path} - will be written{parent}", note)
    if not is_file:
        return PathCheck(cleaned, path, MISSING, f"not found: {path}", note)
    if not _suffix_ok(path, suffixes):
        return PathCheck(cleaned, path, WRONG_TYPE,
                         f"{path} - expected {types}, not {path.suffix or 'no extension'}", note)
    return PathCheck(cleaned, path, OK, str(path), note)


def stored_text(path, base=None) -> str:
    """A picked absolute path as a box shows it: relative when under ``base``."""
    path = Path(path)
    if base is not None:
        try:
            if path.resolve().is_relative_to(Path(base).resolve()):
                return Path(os.path.relpath(path.resolve(), Path(base).resolve())).as_posix()
        except (OSError, ValueError):
            pass
    return str(path)


# -- the browser -------------------------------------------------------------------

@dataclass(frozen=True)
class Listing:
    folder: Path
    folders: list
    files: list           # the files the box wants (all files when ``suffixes`` is empty)
    hidden_files: int     # files left out by type
    truncated: bool
    problem: str = ""


def listing(folder, *, suffixes=(), show_all=False) -> Listing:
    """A folder's sub-folders and files, by name, for the Browse dialog."""
    folder = Path(folder)
    folders, files, hidden = [], [], 0
    try:
        entries = sorted(os.scandir(folder), key=lambda entry: entry.name.lower())
    except OSError as exc:
        return Listing(folder, [], [], 0, False, f"cannot be listed: {exc.strerror or exc}")
    for entry in entries:
        if entry.name.startswith(("~$", ".")):
            continue
        try:
            if entry.is_dir():
                folders.append(Path(entry.path))
            elif show_all or _suffix_ok(Path(entry.name), suffixes):
                files.append(Path(entry.path))
            else:
                hidden += 1
        except OSError:
            continue
    truncated = len(folders) + len(files) > MAX_LISTED
    if truncated:
        folders = folders[:MAX_LISTED]
        files = files[:max(0, MAX_LISTED - len(folders))]
    return Listing(folder, folders, files, hidden, truncated)


def drives() -> list:
    """The drive roots on Windows (C:\\, F:\\ ...); '/' elsewhere."""
    if os.name != "nt":
        return [Path("/")]
    if hasattr(os, "listdrives"):
        try:
            return [Path(drive) for drive in os.listdrives()]
        except OSError:
            pass
    return [Path(f"{letter}:\\") for letter in string.ascii_uppercase
            if Path(f"{letter}:\\").exists()]


def start_folder(text, base=None) -> Path:
    """Where Browse opens: the box's own folder where it has one, else ``base``, else home."""
    path = resolve(text, base)
    for candidate in ([path, path.parent] if path is not None else []) + \
            ([Path(base)] if base is not None else []):
        try:
            if candidate is not None and candidate.is_dir():
                return candidate
        except OSError:
            continue
    return Path.home()
