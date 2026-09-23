"""A report table as Word will take it off the clipboard.

The Report page puts two forms of one table on the clipboard at once: HTML,
which Word pastes as a real table with the house formatting, and tab-separated
text, which is what Excel, a plain-text editor or Word's "keep text only" paste
option take instead. Both come from one ``ReportTable`` so they cannot differ.

The formatting defaults are read off the Callide design hydrology report
(CLD_Design_Hydrology_Report_DRAFT_v1.docx): header rows in ``SunwaterTableHeading``
- white, bold, 10 pt, on the brand cyan ``#00B0CA`` - group rows bold on the
cyan tint ``#C1F7FF``, body text Rubik Light 10 pt, half-point cyan rules. A
study may override any of them under its ``"word"`` key.

These are colours for a Word document, not for the window, so they do not come
from ``palette.py`` and are not held to its contrast rules.

Nothing here imports nicegui.
"""

from __future__ import annotations

import html
from dataclasses import dataclass, field

DATA = "data"
SECTION = "section"

WORD_DEFAULTS = {
    "font": "Rubik Light",
    "size_pt": 10,
    "header_fill": "#00B0CA",
    "header_text": "#FFFFFF",
    "section_fill": "#C1F7FF",
    "rule": "#00B0CA",
}

# What a header or cell may carry beyond plain text: a line break, and the
# superscript in m3 so it pastes as m³ in Word and reads m3 as text.
SUPERSCRIPTS = {"m3/s": ("m<sup>3</sup>/s", "m3/s"),
                "m³/s": ("m<sup>3</sup>/s", "m3/s")}


@dataclass
class Row:
    cells: list
    kind: str = DATA          # DATA, or SECTION - one label across the table
    bold: bool = False
    shaded: bool = False      # the section tint on a data row (T1's DCF row)


@dataclass
class ReportTable:
    header: list
    rows: list = field(default_factory=list)
    align: list = field(default_factory=list)      # per column: left/right/center
    footnotes: list = field(default_factory=list)
    title: str = ""
    problems: list = field(default_factory=list)   # what could not be filled, and why

    @property
    def width(self) -> int:
        return len(self.header)

    def add(self, cells, **kwargs) -> Row:
        row = Row(cells=[str(cell) for cell in cells], **kwargs)
        self.rows.append(row)
        return row

    def section(self, label: str) -> Row:
        row = Row(cells=[str(label)], kind=SECTION, bold=True)
        self.rows.append(row)
        return row

    def alignment(self, column: int) -> str:
        if column < len(self.align):
            return self.align[column]
        return "left" if column == 0 else "right"


def word_style(overrides: dict | None = None) -> dict:
    style = dict(WORD_DEFAULTS)
    for key, value in (overrides or {}).items():
        if key in style and value not in (None, ""):
            style[key] = value
    return style


# -- text ----------------------------------------------------------------------

def _plain(text: str) -> str:
    out = str(text)
    for key, (_, plain) in SUPERSCRIPTS.items():
        out = out.replace(key, plain)
    return out.replace("\n", " ")


def to_text(table: ReportTable) -> str:
    """Tab-separated, one line per row; a section row is its label alone."""
    lines = ["\t".join(_plain(cell) for cell in table.header)]
    for row in table.rows:
        lines.append(_plain(row.cells[0]) if row.kind == SECTION
                     else "\t".join(_plain(cell) for cell in row.cells))
    lines.extend(_plain(note) for note in table.footnotes)
    return "\r\n".join(lines) + "\r\n"


# -- HTML ----------------------------------------------------------------------

def _rich(text: str) -> str:
    out = html.escape(str(text))
    for key, (rich, _) in SUPERSCRIPTS.items():
        out = out.replace(html.escape(key), rich)
    return out.replace("\n", "<br>")


def to_html(table: ReportTable, style: dict | None = None) -> str:
    """A table fragment Word pastes as a table, with the footnotes after it.

    Styles are inline and in points, because Word ignores a stylesheet on the
    clipboard and reads CSS pixels at 96 dpi.
    """
    s = word_style(style)
    font = f"font-family:'{s['font']}',sans-serif;font-size:{s['size_pt']}pt"
    rule = f"border:0.5pt solid {s['rule']}"
    pad = "padding:1.5pt 5pt"

    def cell(tag, text, column, extra=""):
        align = table.alignment(column)
        return (f"<{tag} style=\"{rule};{pad};{font};text-align:{align};"
                f"vertical-align:bottom;{extra}\">{_rich(text)}</{tag}>")

    parts = [f"<table style=\"border-collapse:collapse;{font}\">", "<thead><tr>"]
    header_style = (f"background:{s['header_fill']};color:{s['header_text']};"
                    f"font-weight:bold")
    parts += [cell("th", text, column, header_style)
              for column, text in enumerate(table.header)]
    parts.append("</tr></thead><tbody>")
    for row in table.rows:
        if row.kind == SECTION:
            parts.append(
                f"<tr><td colspan=\"{table.width}\" style=\"{rule};{pad};{font};"
                f"background:{s['section_fill']};font-weight:bold;text-align:left\">"
                f"{_rich(row.cells[0])}</td></tr>")
            continue
        extra = ""
        if row.shaded:
            extra += f"background:{s['section_fill']};"
        if row.bold:
            extra += "font-weight:bold;"
        parts.append("<tr>" + "".join(cell("td", text, column, extra)
                                      for column, text in enumerate(row.cells)) + "</tr>")
    parts.append("</tbody></table>")
    for note in table.footnotes:
        parts.append(f"<p style=\"{font};margin:2pt 0 0 0\">{_rich(note)}</p>")
    return "".join(parts)


def clipboard_document(fragment: str) -> str:
    """The HTML clipboard payload: a whole document, which Word expects."""
    return ("<html><head><meta charset=\"utf-8\"></head><body>"
            f"<!--StartFragment-->{fragment}<!--EndFragment--></body></html>")


def clipboard_script(html_document: str, text: str) -> str:
    """JavaScript that puts both forms on the clipboard and says which it managed.

    ``ClipboardItem`` is what carries HTML; a browser without it (or a webview
    that refuses it) gets the text, and the page says so rather than claiming a
    formatted copy it did not make. Needs a secure context, which localhost is.
    """
    import json
    return (
        "(async () => {"
        f"const html = {json.dumps(html_document)}; const text = {json.dumps(text)};"
        "try {"
        "  if (window.ClipboardItem && navigator.clipboard && navigator.clipboard.write) {"
        "    await navigator.clipboard.write([new ClipboardItem({"
        "      'text/html': new Blob([html], {type: 'text/html'}),"
        "      'text/plain': new Blob([text], {type: 'text/plain'})})]);"
        "    return 'html';"
        "  }"
        "  await navigator.clipboard.writeText(text); return 'text';"
        "} catch (error) {"
        "  try { await navigator.clipboard.writeText(text); return 'text'; }"
        "  catch (second) { return 'failed: ' + second; }"
        "}})()")
