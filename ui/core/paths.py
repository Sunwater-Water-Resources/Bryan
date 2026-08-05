"""Path handling for sims-list values, and the small file primitives.

The one thing to keep straight here: sims lists hold **Windows** relative paths
(``sims_mc\\results\\x.csv``). Two different jobs need them:

- the UI resolving a path to stat it, which must work on Linux too, so the
  separators are normalised first (``resolve_value``);
- the run copy the UI writes for Bryan, which must keep the original string
  **verbatim** - Bryan runs on Windows and the value is its own to interpret.

Never normalise on the way into a run copy. ``runwriter`` copies cell values
untouched for exactly this reason.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path, PurePosixPath, PureWindowsPath

import pandas as pd

# Run ids end up inside the ``simulation_list`` string, which lib/RunLog.py:29
# turns into a log path with a substring ``.replace('.xlsx', '_log.csv')``. A
# run id containing '.xlsx' would corrupt that, so the charset is constrained.
RUN_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def is_blank(value) -> bool:
    """True for the several ways a sims-list cell says 'nothing here'."""
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        pass  # arrays and the like are never blank
    return isinstance(value, str) and not value.strip()


def cell_text(value) -> str:
    """A sims-list cell as the string Bryan would see, or '' when blank."""
    return "" if is_blank(value) else str(value).strip()


def normalise_sep(value: str) -> str:
    """A sims-list path with its separators made local.

    Windows accepts forward slashes, so going to the local separator is safe
    in both directions. Used only for resolving and stat-ing - never for a
    value written back out for Bryan.
    """
    text = str(value).strip().strip('"')
    if not text:
        return ""
    if os.sep == "/":
        # A Windows path being resolved on POSIX: split on either separator.
        return PureWindowsPath(text.replace("/", "\\")).as_posix()
    return str(PureWindowsPath(PurePosixPath(text.replace("\\", "/"))))


def resolve_value(base: Path, value) -> Path | None:
    """Resolve a sims-list path cell against ``base``. None when blank.

    Absolute values are returned as-is, matching ``os.path.join``'s behaviour
    in Bryan (``Main.py:49``).
    """
    text = cell_text(value)
    if not text:
        return None
    local = Path(normalise_sep(text))
    if local.is_absolute():
        return local
    return (Path(base) / local).resolve() if base else local


def safe_run_id(candidate: str) -> str:
    """Validate a run id, raising rather than silently producing a bad path."""
    if not RUN_ID_RE.match(candidate):
        raise ValueError(
            f"run id {candidate!r} must match {RUN_ID_RE.pattern} - it becomes "
            f"part of the 'simulation_list' string that lib/RunLog.py turns "
            f"into a log filename by substring replacement"
        )
    if ".xlsx" in candidate.lower():
        raise ValueError(f"run id {candidate!r} must not contain '.xlsx'")
    return candidate


def atomic_write_json(path: Path, data) -> None:
    """Write JSON via a sibling temp file and ``os.replace``.

    Two UI servers on different ports do not share an in-process lock, so the
    replace is what stops a half-written run.json being read.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2, default=str)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def read_json(path: Path, default=None):
    """Best-effort JSON read - a half-written or absent file is not an error."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def stat_or_none(path: Path | None) -> os.stat_result | None:
    if path is None:
        return None
    try:
        return Path(path).stat()
    except OSError:
        return None
