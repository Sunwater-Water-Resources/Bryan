"""Study: the study file the report pages read, opened here.

A **study file** (``bryan_study.json``) sits at the top of a dam's design flood
study and holds what is about the dam and its report rather than one run: the
runs the report draws on, the report tables, the PMF, the figures and the lake
record. Report, PMF, Figures and Lake record all read it.

The study the launcher last had open is reopened when it starts
(``AppState.reopen_last_study``); this page says so when that could not be done,
opens or creates another, and shows what the open one holds. Its paths are kept
relative to its own folder, so a study copied to another disk still opens.
"""

from __future__ import annotations

from nicegui import ui

from core import lakerecord, study as studies
from core.paths import clean_path_text
from layout import page_frame, severity_banner
from state import STATE

# The pages that read the study, and what each keeps in it.
READERS = (("Report", "/report", "the runs and the report tables"),
           ("PMF", "/pmf", "the PMF groups and the adopted notional AEP"),
           ("Figures", "/figures", "the report figures"),
           ("Lake record", "/lake-record", "the lake level record and its four steps"))


def study_page() -> None:
    with page_frame("Study"):
        _StudyView().build()


def _plural(count, word) -> str:
    return f"{count} {word}{'' if count == 1 else 's'}"


def contents(study) -> dict:
    """What the study holds, counted, for the page and its tests."""
    record = lakerecord.settings(study)["homogenise"]
    return {
        "runs": [run["name"] for run in study.runs],
        "tables": len(study.tables),
        "pmf_groups": len(studies.pmf_entries(study)),
        "figures": len(study.extra.get("figures") or []),
        "gauges": len([g for g in record["gauges"] if str(g).strip()]),
    }


class _StudyView:
    def __init__(self) -> None:
        self.body = None

    def build(self) -> None:
        self.body = ui.column().classes("w-full gap-4")
        self.redraw()

    def redraw(self) -> None:
        self.body.clear()
        with self.body:
            if STATE.study_problem:
                severity_banner("warn", STATE.study_problem,
                                "Open it again below once it can be reached, or open "
                                "another study.")
            if STATE.study is None:
                self._open_card(first=True)
                return
            self._study_card()
            self._contents_card()
            with ui.expansion("Open or create another study", icon="folder_open") \
                    .classes("w-full").mark("another-study"):
                self._open_card(first=False)

    # -- opening -----------------------------------------------------------------

    def _open_card(self, *, first: bool) -> None:
        with ui.card().classes("w-full"):
            if first:
                ui.label("Open a study").classes("text-lg font-bold")
                ui.label("A study file names the runs a report is written from, and "
                         "keeps the report tables, the PMF, the figures and the lake "
                         "record. Keep it at the top of the study folder: its paths "
                         "are stored relative to it, so the study can be moved.") \
                    .classes("text-sm text-body")
            default = STATE.settings.last_study or ""
            if not default and STATE.project is not None:
                default = str(STATE.project.config.config_path.parent
                              / studies.DEFAULT_NAME)
            path = ui.input("Study file (or its folder)", value=default) \
                .classes("w-full").props("dense").mark("study-path")
            name = ui.input("Name, for a new study", value="").classes("w-full") \
                .props("dense").mark("study-name")
            with ui.row().classes("gap-2"):
                ui.button("Open", icon="folder_open",
                          on_click=lambda: self._open(path.value)).mark("open-study")
                ui.button("New study", icon="add",
                          on_click=lambda: self._open(path.value, create=True,
                                                      name=name.value)) \
                    .props("outline").mark("new-study")
            recent = [entry for entry in STATE.settings.recent_studies if entry]
            if recent:
                ui.label("Recent").classes("text-sm text-muted pt-2")
                for entry in recent:
                    ui.button(entry, on_click=lambda _, e=entry: self._open(e)) \
                        .props("flat dense no-caps align=left").classes("mono text-sm")

    def _open(self, path, *, create=False, name="") -> None:
        text = clean_path_text(path or "")
        if not text:
            ui.notify("Give the study file's path", type="warning")
            return
        try:
            STATE.open_study(text, create=create, name=name)
        except (studies.StudyError, OSError) as exc:
            ui.notify(str(exc), type="negative", multi_line=True)
            return
        self.redraw()

    def _close(self) -> None:
        STATE.close_study()
        self.redraw()

    # -- the open study ------------------------------------------------------------

    def _study_card(self) -> None:
        study = STATE.study
        with ui.card().classes("w-full").mark("study-card"):
            with ui.row().classes("w-full items-start justify-between no-wrap"):
                with ui.column().classes("gap-0 grow"):
                    name = ui.input("Study", value=study.name) \
                        .props("dense borderless").classes("text-lg font-bold w-full") \
                        .mark("study-title")
                    name.on("blur", lambda: self._rename(name.value))
                    ui.label(str(study.path)).classes("mono text-xs text-muted")
                    ui.label(f"Paths in it are kept relative to {study.folder}.") \
                        .classes("text-xs text-muted")
                ui.button("Close", icon="close", on_click=self._close) \
                    .props("flat").mark("close-study")

    def _rename(self, value) -> None:
        study = STATE.study
        if study is None or value == study.name:
            return
        study.name = value
        try:
            study.save()
        except OSError as exc:
            ui.notify(f"Could not save {study.path.name}: {exc}", type="warning")

    def _contents_card(self) -> None:
        found = contents(STATE.study)
        counts = {
            "/report": (_plural(len(found["runs"]), "run") + ", "
                        + _plural(found["tables"], "table")),
            "/pmf": _plural(found["pmf_groups"], "PMF group"),
            "/figures": _plural(found["figures"], "figure"),
            "/lake-record": _plural(found["gauges"], "gauge export"),
        }
        with ui.card().classes("w-full").mark("study-contents"):
            ui.label("What it holds").classes("text-lg font-bold")
            for label, target, what in READERS:
                with ui.row().classes("w-full items-center gap-3 no-wrap"):
                    ui.button(label, on_click=lambda _, t=target: ui.navigate.to(t)) \
                        .props("flat dense no-caps align=left").classes("w-32")
                    ui.label(counts[target]).classes("w-48 text-body")
                    ui.label(what).classes("text-sm text-muted")
            if found["runs"]:
                ui.label("Runs: " + ", ".join(found["runs"])) \
                    .classes("text-sm text-body pt-1").mark("study-runs")
