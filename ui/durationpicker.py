"""Which storm durations a critical duration is taken over - the two pickers.

The Tables and Figures pages both read a group's design curve as the envelope
over its storm durations. The study keeps a default set, and each table and
figure may follow it, take every duration, or choose its own; the rule and its
storage are ``core/study.py``'s (``durations_for``). These are the controls.
"""

from __future__ import annotations

from nicegui import ui

from core import study as studies

FOLLOW, ALL, CHOOSE = "follow", "all", "choose"


def _options(study, chosen=()) -> dict:
    """Every duration the study's runs hold, with any already chosen, as '24 h'."""
    hours = sorted(set(studies.durations_in(study)) | set(studies.parse_hours(chosen)))
    return {hour: f"{hour:g} h" for hour in hours}


class Deferred:
    """A multi-select's value, applied once its list closes rather than per click.

    Applying the study default rebuilds every open table or figure, so applying
    it on each tick - and, on the Tables page, rebuilding the select under the
    cursor - made picking several durations a fight. While the list is open the
    picks are only held; they are applied when it closes, and only if they
    differ from what was last applied. A change with the list shut (a chip's
    cross, the clear button) is one action, so it applies at once.
    """

    def __init__(self, value, apply) -> None:
        self.applied = studies.parse_hours(value)
        self.value = list(self.applied)
        self.apply = apply
        self.is_open = False

    def opened(self) -> None:
        self.is_open = True

    def closed(self) -> None:
        self.is_open = False
        self._flush()

    def changed(self, value) -> None:
        self.value = studies.parse_hours(value)
        if not self.is_open:
            self._flush()

    def _flush(self) -> None:
        if self.value != self.applied:
            self.applied = list(self.value)
            self.apply(list(self.value))


def default_card(study, on_change) -> None:
    """The study's default: a card for the top of the Tables and Figures pages.

    ``on_change()`` runs after the default is stored and saved, so the page can
    rebuild what was drawn under the old one.
    """
    chosen = studies.default_durations(study)

    def store(hours) -> None:
        studies.set_default_durations(study, hours)
        try:
            study.save()
        except OSError as exc:
            ui.notify(f"Could not save {study.path.name}: {exc}", type="warning")
            return
        on_change()

    pending = Deferred(chosen, store)

    with ui.card().classes("w-full gap-1").mark("study-durations"):
        with ui.row().classes("w-full items-center gap-4 no-wrap"):
            with ui.column().classes("gap-0 w-80 shrink-0"):
                ui.label("Storm durations").classes("font-bold")
                ui.label("The critical duration of every table and figure is taken over "
                         "these, unless the table or figure chooses its own. None chosen "
                         "is every duration. Applied when the list closes."
                         ).classes("text-xs text-muted")
            ui.select(_options(study, chosen), multiple=True, value=chosen,
                      label="Study default (blank: all durations)",
                      on_change=lambda e: pending.changed(e.value)) \
                .on("popup-show", lambda _: pending.opened()) \
                .on("popup-hide", lambda _: pending.closed()) \
                .classes("grow").props("dense use-chips clearable") \
                .mark("study-durations-select")


def mode_of(spec: dict) -> str:
    value = spec.get(studies.DURATIONS_KEY)
    if isinstance(value, str) and value.strip().lower() == studies.ALL_DURATIONS:
        return ALL
    return FOLLOW if studies.follows_default(spec) else CHOOSE


def spec_picker(study, spec: dict, *, what: str = "table") -> None:
    """Follow the study, all durations, or chosen ones - for one table or figure.

    Writes ``spec["durations"]``: None, ``"all"`` or a list of hours. Choosing
    starts from what the table reads now, so switching to it changes nothing until
    a duration is taken off.
    """
    holder = ui.column().classes("w-full gap-1")

    def draw() -> None:
        holder.clear()
        mode = mode_of(spec)
        default = studies.default_durations(study)
        with holder:
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                ui.label("Critical duration over").classes("text-sm w-44 shrink-0")
                ui.toggle({FOLLOW: f"Study default ({studies.durations_text(default)})",
                           ALL: "All durations", CHOOSE: "Choose"},
                          value=mode, on_change=lambda e: set_mode(e.value)) \
                    .props("dense no-caps").mark(f"{what}-durations-mode")
            if mode == CHOOSE:
                chosen = studies.parse_hours(spec.get(studies.DURATIONS_KEY))
                ui.select(_options(study, chosen), multiple=True, value=chosen,
                          label=f"Durations this {what} takes the critical duration over",
                          on_change=lambda e: spec.update(
                              {studies.DURATIONS_KEY: studies.parse_hours(e.value)})) \
                    .classes("w-full").props("dense use-chips") \
                    .mark(f"{what}-durations-select")
                ui.label("None chosen follows the study default.") \
                    .classes("text-xs text-muted")

    def set_mode(mode) -> None:
        if mode == mode_of(spec):
            return
        if mode == FOLLOW:
            spec[studies.DURATIONS_KEY] = None
        elif mode == ALL:
            spec[studies.DURATIONS_KEY] = studies.ALL_DURATIONS
        else:
            current = studies.durations_for(study, spec) or studies.durations_in(study)
            spec[studies.DURATIONS_KEY] = list(current)
        draw()

    draw()


def describe(study, spec: dict) -> str:
    """'' when a table or figure reads every duration; otherwise which, and whose."""
    hours = studies.durations_for(study, spec)
    if not hours:
        return "" if mode_of(spec) != ALL or not studies.default_durations(study) \
            else "all durations"
    whose = "study default" if mode_of(spec) == FOLLOW else "own choice"
    return f"{studies.durations_text(hours)} ({whose})"
