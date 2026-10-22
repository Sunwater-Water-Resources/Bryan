"""No source file in Bryan compiles with an invalid escape sequence.

``sep="\\s+"`` written without the ``r`` worked, but Python 3.12 prints a
SyntaxWarning for it on every first import, and the downstream log opened with
two of them from ``lib/URBSmodel.py`` - read by a user as an error. A later
Python makes it one. The files are compiled, not imported, so this checks
``lib/`` without the UI environment needing anything ``lib/`` imports.
"""

import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SKIP = {".venv", "venv", "__pycache__", "site-packages", ".git"}


def _sources():
    for path in sorted(ROOT.rglob("*.py")):
        if not SKIP.intersection(path.relative_to(ROOT).parts):
            yield path


def test_no_source_has_an_invalid_escape_sequence():
    found = []
    for path in _sources():
        source = path.read_text(encoding="utf-8", errors="replace")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            compile(source, str(path), "exec")
        found += [f"{path.relative_to(ROOT)}:{w.lineno} {w.message}" for w in caught
                  if issubclass(w.category, SyntaxWarning)]
    assert not found, "write these as raw strings:\n" + "\n".join(found)
