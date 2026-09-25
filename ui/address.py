"""Addresses that carry the view: ``/results?run=E013%20RFSL&group=...&type=level``.

What a page is showing - the run, the group, the result type, the report table -
is kept in its address as it changes, so Back returns to it, and a view can be
bookmarked or sent. A link sent to a colleague names the run by the study's name
for it, so it opens that run rather than whichever one they last had open.

The address is read once, as the page is drawn (``param``); changes are written
back with ``keep``, which replaces the address rather than adding a history
entry per click, so Back goes to the previous page, not the previous tick box.
"""

from __future__ import annotations

from pathlib import Path
from urllib.parse import urlencode

from nicegui import ui

from state import STATE

_VIEW = "_bryan_view"          # the client's current query, as ``keep`` last wrote it


def _client():
    try:
        return ui.context.client
    except (AttributeError, RuntimeError):
        return None


def _view() -> dict:
    client = _client()
    if client is None:
        return {}
    view = getattr(client, _VIEW, None)
    if view is None:
        try:
            view = dict(client.request.query_params)
        except (AttributeError, RuntimeError, AssertionError):
            view = {}
        setattr(client, _VIEW, view)
    return view


def _path() -> str:
    client = _client()
    try:
        return client.request.url.path if client is not None else "/"
    except (AttributeError, RuntimeError, AssertionError):
        return "/"


def param(name: str, default: str = "") -> str:
    """One value from the page's address, as it was when the page was drawn."""
    return str(_view().get(name) or default)


def url(path: str, **params) -> str:
    kept = {key: str(value) for key, value in params.items() if value not in (None, "")}
    return f"{path}?{urlencode(kept)}" if kept else path


def keep(**params) -> None:
    """Put these in the address (a blank or None value takes one out)."""
    view = _view()
    for key, value in params.items():
        if value in (None, ""):
            view.pop(key, None)
        else:
            view[key] = str(value)
    _replace(url(_path(), **view))


def _replace(address: str) -> None:
    try:
        ui.navigate.history.replace(address)
    except (AttributeError, RuntimeError):
        pass


def follow_run() -> str:
    """Open the study run the address names, when another one is open.

    Returns why it could not be, or ''. Called by ``page_frame`` before a page
    draws, so every page shows the run its address names.
    """
    name = param("run")
    if not name or STATE.study is None:
        return ""
    path = STATE.study.run_config_path(name)
    if path is None:
        return f"The address names a run, {name}, that the study does not have."
    current = STATE.project.config.config_path if STATE.project is not None else None
    if current is not None and Path(current).resolve() == Path(path).resolve():
        return ""
    try:
        STATE.open_project(path)
    except Exception as exc:                      # noqa: BLE001 - said on the page
        return f"The address names {name}, which could not be opened: {exc}"
    return ""


def keep_run(name: str) -> None:
    """The open run's study name in the address (blank when it is not a study run)."""
    keep(run=name)


def switch_run_to(name: str) -> str:
    """The current page's address with another run in it, for the runs panel."""
    view = dict(_view())
    view["run"] = name
    return url(_path(), **view)
