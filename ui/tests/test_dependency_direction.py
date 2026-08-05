"""Pin the three binding design rules for ui/.

Modelled on hydraulic_model/tests/ui/test_dependency_direction.py.

1. nothing under lib/ or the top-level scripts references ui;
2. nothing under ui/core/ imports nicegui;
3. ui/ touches Bryan only through the allow-list in core/bryan.py.

Rule 3 is what keeps the UI environment free of scipy and matplotlib, and
what stops the UI quietly coupling itself to Bryan's internals.
"""

from __future__ import annotations

import re
from pathlib import Path

from core.bryan import ALLOWED_MODULES, BRYAN_ROOT

UI_ROOT = BRYAN_ROOT / "ui"
CORE = UI_ROOT / "core"

IMPORT_RE = re.compile(
    r"^\s*(?:from\s+(lib[\w.]*)\s+import|import\s+(lib[\w.]*))", re.MULTILINE)


def python_files(root: Path):
    return [path for path in root.rglob("*.py")
            if "__pycache__" not in path.parts]


def test_bryan_never_references_the_ui():
    offenders = []
    for path in python_files(BRYAN_ROOT / "lib"):
        if re.search(r"\bui\.", path.read_text(encoding="utf-8")):
            offenders.append(str(path.relative_to(BRYAN_ROOT)))
    for name in ("Main.py", "RouteFlows.py"):
        script = BRYAN_ROOT / name
        if script.is_file() and "ui." in script.read_text(encoding="utf-8"):
            offenders.append(name)
    assert not offenders, f"Bryan must not depend on the UI: {offenders}"


def test_core_never_imports_nicegui():
    """core/ must stay usable by the tests, a CLI, and the result viewer."""
    offenders = [
        str(path.relative_to(UI_ROOT))
        for path in python_files(CORE)
        if re.search(r"^\s*(from|import)\s+nicegui", path.read_text(encoding="utf-8"),
                     re.MULTILINE)
    ]
    assert not offenders, f"core must not import nicegui: {offenders}"


def test_ui_imports_bryan_only_through_the_allow_list():
    allowed = set(ALLOWED_MODULES)
    offenders = []
    for path in python_files(UI_ROOT):
        if path.name in ("bryan.py", "fake_main.py"):
            continue          # the seam itself, and the deliberate stand-in
        text = path.read_text(encoding="utf-8")
        for match in IMPORT_RE.finditer(text):
            module = match.group(1) or match.group(2)
            if module not in allowed and module != "lib":
                offenders.append(f"{path.relative_to(UI_ROOT)}: {module}")
    assert not offenders, (
        f"import Bryan through core/bryan.py, not directly: {offenders}")


def test_the_allow_list_is_actually_importable():
    from core.bryan import log_files, run_log

    assert hasattr(run_log(), "COLUMNS")
    assert hasattr(log_files(), "resolve_duplicates")


def test_the_allow_list_stays_cheap():
    """RunLog and LogFiles must not drag scipy or matplotlib into the UI env."""
    import subprocess
    import sys

    code = (
        "import sys; sys.path.insert(0, %r);"
        "from lib import RunLog, LogFiles;"
        "bad = [m for m in ('scipy', 'matplotlib') if m in sys.modules];"
        "print(','.join(bad))" % str(BRYAN_ROOT)
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "", (
        f"the allow-listed modules pulled in {result.stdout.strip()}")
