"""Reading the end of a file that is still being written to.

Bryan's console output is redirected to a file rather than a pipe, so there is
no reader thread to keep fed and nothing deadlocks on a full pipe buffer. The
cost is that the UI has to go and look, which is what this does - seeking near
the end rather than reading a log that may be hundreds of megabytes.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

DEFAULT_TAIL_BYTES = 64 * 1024


@dataclass(frozen=True)
class Tail:
    text: str
    size: int
    mtime: float
    truncated: bool

    @property
    def exists(self) -> bool:
        return self.size > 0 or bool(self.text)

    def lines(self) -> list:
        return self.text.splitlines()


def tail(path, max_bytes: int = DEFAULT_TAIL_BYTES) -> Tail:
    """The last ``max_bytes`` of ``path``, decoded leniently.

    Bryan prints a lot and its output can carry anything a model wrote, so
    decoding never raises - a mangled byte is not worth losing the log over.
    """
    path = Path(path)
    try:
        stat = path.stat()
    except OSError:
        return Tail("", 0, 0.0, False)

    size = stat.st_size
    start = max(0, size - max_bytes)
    try:
        with path.open("rb") as stream:
            if start:
                stream.seek(start, os.SEEK_SET)
            data = stream.read()
    except OSError:
        return Tail("", size, stat.st_mtime, False)

    text = data.decode("utf-8", errors="replace")
    if start:
        # Drop the partial first line so nothing is shown half-decoded.
        newline = text.find("\n")
        text = text[newline + 1:] if newline >= 0 else text
    return Tail(text, size, stat.st_mtime, truncated=bool(start))


def head_and_tail(path, max_bytes: int = DEFAULT_TAIL_BYTES) -> str:
    """The start and end of a file, for a log the user wants to eyeball.

    The head matters for Bryan's per-simulation logs: reservoir routing opens
    each one with a SIMULATION INPUTS block recording every input path, size
    and mtime, which is exactly what someone checking a result wants to see.
    """
    path = Path(path)
    try:
        size = path.stat().st_size
    except OSError:
        return ""
    if size <= max_bytes * 2:
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""

    try:
        with path.open("rb") as stream:
            start = stream.read(max_bytes).decode("utf-8", errors="replace")
            stream.seek(-max_bytes, os.SEEK_END)
            end = stream.read().decode("utf-8", errors="replace")
    except OSError:
        return ""
    skipped = size - 2 * max_bytes
    return f"{start}\n\n... {skipped:,} bytes not shown ...\n\n{end}"
