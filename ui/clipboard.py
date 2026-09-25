"""Copy for Word, wherever a page shows a table.

The Report page's tables, the Results page's and the Ensemble page's critical
durations all go to the clipboard the same way (core/wordtable.py): as a
formatted table in the report's style, which Word pastes, and as tab-separated
text, which Excel takes. The style is the open study's ``"word"`` settings,
or the house defaults without a study.
"""

from __future__ import annotations

from nicegui import ui

from core import wordtable
from state import STATE


def table_from_rows(columns, rows, *, title="") -> wordtable.ReportTable:
    """A table as a page shows it in ``ui.table`` - its column labels and its rows
    of already formatted cells - as a report table."""
    table = wordtable.ReportTable(header=[str(column["label"]) for column in columns],
                                  title=title)
    for row in rows:
        table.add([row.get(column["field"], "") for column in columns])
    return table


def word_style() -> dict | None:
    return STATE.study.extra.get("word") if STATE.study is not None else None


async def copy_table(table: wordtable.ReportTable, *, rich: bool, warnings=()) -> None:
    """``warnings`` - what is out of date in the results the table is drawn from -
    are said again as it is copied, so a stale table never reaches Word unsaid."""
    text = wordtable.to_text(table)
    if not rich:
        ui.clipboard.write(text)
        _said("Copied as text", warnings)
        return
    fragment = wordtable.to_html(table, word_style())
    outcome = await ui.run_javascript(
        wordtable.clipboard_script(wordtable.clipboard_document(fragment), text),
        timeout=5.0)
    if outcome == "html":
        _said("Copied - paste into Word", warnings)
    elif outcome == "text":
        ui.notify("This browser would only take text; copied as text", type="warning")
    else:
        ui.notify(f"Could not copy: {outcome}", type="negative")


def _said(done: str, warnings) -> None:
    if warnings:
        ui.notify(f"{done} - but {warnings[0]}"
                  + (f" (and {len(warnings) - 1} more)" if len(warnings) > 1 else ""),
                  type="warning", multi_line=True, timeout=10000, close_button=True)
    else:
        ui.notify(done)


def copy_buttons(get_table, *, mark: str, get_warnings=None) -> None:
    """Copy for Word and Copy as text, for the table ``get_table()`` returns now
    (None when there is nothing to copy)."""
    async def copy(rich: bool) -> None:
        table = get_table()
        if table is None or not table.rows:
            ui.notify("Nothing to copy", type="warning")
            return
        await copy_table(table, rich=rich,
                         warnings=get_warnings() if get_warnings is not None else ())

    ui.button("Copy for Word", icon="content_copy", on_click=lambda: copy(True)) \
        .props("dense no-caps").mark(f"copy-word-{mark}")
    ui.button("Copy as text", icon="notes", on_click=lambda: copy(False)) \
        .props("flat dense no-caps").mark(f"copy-text-{mark}") \
        .tooltip("Tab-separated, for Excel")
