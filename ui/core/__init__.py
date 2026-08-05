"""Bryan-facing logic for the run launcher UI.

Nothing in this package may import ``nicegui``. Everything here is a plain,
testable module so the later result viewer - and any CLI - can reuse it. The
rule is pinned by ``ui/tests/test_dependency_direction.py``.

The package touches Bryan itself through a deliberately narrow allow-list:
``lib.RunLog`` and ``lib.LogFiles``. Both import only os/platform/datetime and
pandas, so the UI environment needs no scipy and no matplotlib, and the log
plumbing is Bryan's own rather than a copy that drifts.
"""
