"""The user manual, served by the launcher from this copy of Bryan.

Each page's help button opens its section of ``Manual/SubDocs/ui.md`` - the
manual of the Bryan the launcher belongs to, so it always matches it, and
without GitHub access or a network. The Markdown is shown as it is written;
here it only gets an anchor at each heading, in the manual's own style (lower
case, letters and digits: "Dam inputs" is ``#daminputs``), and its links to the
other docs are pointed at the launcher's copies.
"""

from __future__ import annotations

import re
from pathlib import Path

from .bryan import BRYAN_ROOT

MANUAL_DIR = BRYAN_ROOT / "Manual" / "SubDocs"
HOME = "ui"
_NAME = re.compile(r"^[A-Za-z0-9_-]+$")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
_DOC_LINK = re.compile(r"\]\(([A-Za-z0-9_-]+)\.md(#[^)]*)?\)")

# Each launcher page's section: (doc, heading). Pages are named by the title
# page_frame is given, as the menu highlights them.
SECTIONS = {
    "Study": (HOME, "The study"),
    "Simulations": (HOME, "Opening a sims list"),
    "Select": (HOME, "Choosing simulations"),
    "Run": (HOME, "What is checked before a run"),
    "Results": (HOME, "Viewing results"),
    "Ensemble": (HOME, "Viewing ensemble results"),
    "Events": (HOME, "Choosing representative events"),
    "Lake levels": (HOME, "Lake level frequency"),
    "Downstream": ("downstream_storms", ""),
    "Report": (HOME, "Report tables"),
    "PMF": (HOME, "The PMF"),
    "Figures": (HOME, "Report figures"),
    "Lake record": (HOME, "Lake record"),
}


def slug(text: str) -> str:
    """A heading's anchor, as the manual's own links write them."""
    return re.sub(r"[^a-z0-9]", "", text.lower())


def doc_path(name: str) -> Path | None:
    """A doc of the manual by name ('ui', 'downstream_storms'), or None."""
    if not name or not _NAME.match(name):
        return None
    path = MANUAL_DIR / f"{name}.md"
    return path if path.is_file() else None


def prepare(text: str) -> tuple[str, list]:
    """The Markdown with an anchor before each heading and its doc links made the
    launcher's; and the headings, as (level, text, anchor), for a contents list."""
    headings, lines, fenced = [], [], False
    for line in text.splitlines():
        if line.lstrip().startswith("```") and line.count("```") % 2 == 1:
            fenced = not fenced                     # never a heading inside code
        found = None if fenced else _HEADING.match(line)
        if found:
            level, title = len(found.group(1)), found.group(2)
            anchor = slug(title)
            headings.append((level, title, anchor))
            lines += [f'<a id="{anchor}"></a>', ""]
        lines.append(line)
    body = "\n".join(lines)
    body = _DOC_LINK.sub(lambda m: f"](/manual/{m.group(1)}{m.group(2) or ''})", body)
    return body, headings


def section_for(title: str) -> str | None:
    """The address of a page's section of the manual, or None when it has none."""
    for label, (doc, heading) in SECTIONS.items():
        if title.startswith(label):
            return f"/manual/{doc}" + (f"#{slug(heading)}" if heading else "")
    return None
