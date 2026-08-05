"""Page registry.

Pages are plain functions bound to routes by ``register_all()`` rather than by
import-time ``@ui.page`` decorators, so tests can re-register them onto a fresh
NiceGUI app. Same arrangement as swe2d_ui/pages/__init__.py.
"""

from __future__ import annotations

from nicegui import ui


def register_all() -> None:
    from pages import edit, history, project, run, select

    for route, handler in [
        ("/", project.project_page),
        ("/select", select.select_page),
        ("/run", run.run_page),
        ("/edit", edit.edit_page),
        ("/history", history.history_page),
    ]:
        ui.page(route)(handler)
