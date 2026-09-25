"""Manual: this Bryan's user manual, opened at a page's section by its help button."""

from __future__ import annotations

from nicegui import ui

from core import manual
from layout import page_frame


def manual_page(doc: str = manual.HOME) -> None:
    with page_frame("Manual"):
        path = manual.doc_path(doc)
        if path is None:
            ui.label(f"There is no manual page called {doc!r}.").classes("text-muted")
            ui.link("The launcher's manual", f"/manual/{manual.HOME}")
            return
        # A heading scrolled to sits below the menu bar, not under it.
        ui.add_css(".nicegui-markdown a[id] { display: block; scroll-margin-top: 72px; }")
        body, headings = manual.prepare(path.read_text(encoding="utf-8"))
        tops = [(title, anchor) for level, title, anchor in headings if level == 2]
        with ui.card().classes("w-full"):
            ui.label(f"From {path}").classes("mono text-xs text-muted")
            if tops:
                with ui.row().classes("gap-x-4 gap-y-1 flex-wrap").mark("manual-contents"):
                    for title, anchor in tops:
                        ui.link(title, f"#{anchor}").classes("text-sm")
            ui.markdown(body).classes("w-full").mark("manual-body")
        # The browser looks for the #section before the page is drawn, so look again.
        ui.timer(0.3, lambda: ui.run_javascript(
            "if (location.hash) { const target = document.getElementById("
            "decodeURIComponent(location.hash.slice(1))); "
            "if (target) target.scrollIntoView(); }"), once=True)
