"""The narrow seam onto Bryan itself.

The UI imports exactly two Bryan modules, and only through here:

- ``lib.RunLog``   - the run log format, so the UI reads what Bryan writes
- ``lib.LogFiles`` - ``resolve_duplicates``, so predicted log paths are the
  ones Bryan will actually use

Both import only os/platform/datetime and pandas, so the UI environment needs
no scipy and no matplotlib. Nothing else in ``lib/`` is importable that cheaply,
and nothing else should be imported at all - the allow-list is pinned by
``ui/tests/test_dependency_direction.py``.

Bryan has no packaging, so the repo root has to go on ``sys.path`` first.
"""

from __future__ import annotations

import sys
from pathlib import Path

BRYAN_ROOT = Path(__file__).resolve().parents[2]

ALLOWED_MODULES = ("lib.RunLog", "lib.LogFiles")


def ensure_importable() -> Path:
    """Put the Bryan repo root on sys.path. Idempotent."""
    root = str(BRYAN_ROOT)
    if root not in sys.path:
        sys.path.insert(0, root)
    return BRYAN_ROOT


def run_log():
    ensure_importable()
    from lib import RunLog
    return RunLog


def log_files():
    ensure_importable()
    from lib import LogFiles
    return LogFiles


def run_log_columns() -> list:
    """The run log's columns - the schema core/progress.py reads."""
    return list(run_log().COLUMNS)


def statuses() -> dict:
    module = run_log()
    return {
        "completed": module.COMPLETED,
        "failed": module.FAILED,
        "missing_input": module.FAILED_MISSING_INPUT,
    }


def resolve_duplicate_logs(frame, capture=False):
    """Bryan's own per-batch log dedupe, applied across a whole plan.

    ``lib/LogFiles.resolve_duplicates`` only dedupes within one call, so two
    chunks each holding a row that names ``sims_mc\\log\\TFD_run.log`` would
    both keep it and both open it with "w". Calling it once over every chunk's
    rows in plan order gives global uniqueness *and* exactly the names a single
    sequential run would have produced - so Bryan's own call, when it runs each
    chunk, becomes a no-op.

    The function announces its renames with ``print``. A UI has no console to
    print to, so ``capture=True`` returns ``(frame, text)`` instead of letting
    it escape to wherever the server's stdout happens to point.
    """
    if not capture:
        return log_files().resolve_duplicates(frame)

    import contextlib
    import io

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        resolved = log_files().resolve_duplicates(frame)
    return resolved, buffer.getvalue()
