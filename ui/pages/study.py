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

from core import dam as dams, lakerecord, study as studies
from core.paths import clean_path_text
from layout import page_frame, severity_banner
from state import STATE
from widgets import EITHER, path_input, path_list

# The pages that read the study, and what each keeps in it.
READERS = (("Report", "/report", "the runs and the report tables"),
           ("PMF", "/pmf", "the PMF groups and the adopted notional AEP"),
           ("Figures", "/figures", "the report figures"),
           ("Lake record", "/lake-record", "the lake level record and its four steps"))


def study_page() -> None:
    with page_frame("Study"):
        _StudyView().build()


def _number(value, default=None):
    try:
        return float(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return default


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
            self._dam_card()
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
            path = path_input("Study file (or its folder)", default, expect=EITHER,
                              suffixes=(".json",), mark="study-path")
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

    # -- the dam inputs ------------------------------------------------------------

    def _dam_card(self) -> None:
        """What more than one analysis reads about the dam, entered once.

        Kept under ``"dam"`` in the study (``core/dam.py``). The Lake record page
        shows them and links back here rather than editing them itself.
        """
        study = STATE.study
        dam = dams.settings(study)
        with ui.card().classes("w-full").mark("dam-card"):
            ui.label(lakerecord.RECORD_CARD).classes("text-lg font-bold")
            ui.label("What more than one analysis reads about the dam, entered once. "
                     "The Lake record page's steps read these. Paths are kept relative "
                     "to the study.").classes("text-xs text-muted")

            ui.label("Lake level record").classes("text-sm font-bold pt-1")
            path_list("Gauge exports (WMIP / Hydstra), one per line, in the order the "
                      "gauges operated", dams.gauges(dam), base=study.folder,
                      base_name="the study folder", mark="dam-gauges",
                      on_commit=lambda lines: self._set_dam(dam, "gauges", [
                          studies.portable(study.folder, line) for line in lines]))
            overlay = dam.get("overlay")
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                ui.checkbox("Overlay gauge below a level", value=bool(overlay),
                            on_change=lambda e: self._toggle_overlay(dam, e.value)) \
                    .mark("dam-overlay")
                if overlay:
                    with ui.element("div").classes("grow"):
                        self._dam_path("Overlay gauge export", overlay, dam, "file")
                    ui.input("below (m)", value=f"{overlay.get('below') or ''}") \
                        .classes("w-28").props("dense") \
                        .on("blur", lambda e: self._set_dam(
                            overlay, "below", _number(e.sender.value), dam=dam))
                    ui.input("reconnect margin (m)",
                             value=f"{overlay.get('reconnect_margin', 0.10)}") \
                        .classes("w-40").props("dense") \
                        .on("blur", lambda e: self._set_dam(
                            overlay, "reconnect_margin", _number(e.sender.value, 0.10),
                            dam=dam))

            ui.label("Storage and release").classes("text-sm font-bold pt-1")
            with ui.row().classes("w-full gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._dam_path("Storage table (.els: EL, A, V)", dam, dam, "storage",
                                   mark="dam-storage", suffixes=(".els", ".csv"))
                with ui.element("div").classes("grow"):
                    self._dam_path("Ratings: a register (.xlsx), or one rating for the "
                                   "whole record (.rat, level,flow .csv or .sq)",
                                   dam, dam, "register", mark="dam-register", redraw=True,
                                   suffixes=dams.REGISTER_SUFFIXES + (".rat", ".csv", ".sq"))
                if dams.single_rating(dam):
                    ui.input("Its full supply level (m AHD; blank: the file's)",
                             value=f"{dam['register_fsl']:g}" if dam["register_fsl"] else "") \
                        .classes("w-64").props("dense").mark("dam-register-fsl") \
                        .on("blur", lambda e: self._set_dam(dam, "register_fsl",
                                                            _number(e.sender.value)))
            if dams.single_rating(dam):
                ui.label("One rating for the whole record, for a dam whose spillway has "
                         "not changed. Homogenising then routes the record through much the "
                         "rating it was made under, so the homogenised levels follow the "
                         "recorded ones; the antecedent storage and the inflow record are "
                         "what it is for. A .sq that does not state its full supply level "
                         "needs it given.").classes("text-xs text-muted") \
                    .mark("dam-single-rating")
            with ui.row().classes("w-full items-start gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._dam_path("Evaporation (SILO Data Drill)", dam, dam, "evaporation",
                                   mark="dam-evaporation")
                ui.input("Pan factors, Jan to Dec",
                         value=" ".join(f"{v:g}" for v in dam["pan_factors"])) \
                    .classes("w-[28rem]").props("dense").mark("dam-pan") \
                    .on("blur", lambda e: self._set_pan(dam, e.sender.value))

            ui.label("Catchment").classes("text-sm font-bold pt-1")
            with ui.row().classes("w-full items-start gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._dam_path("Catchment shapefile (.shp)", dam, dam, "shapefile",
                                   mark="dam-shapefile", suffixes=(".shp",))
                ui.input("Field (blank: every polygon)", value=dam["field"]) \
                    .classes("w-48").props("dense") \
                    .on("blur", lambda e: self._set_dam(dam, "field", e.sender.value.strip()))
                ui.input("equal to", value=dam["value"]).classes("w-36").props("dense") \
                    .on("blur", lambda e: self._set_dam(dam, "value", e.sender.value.strip()))
                ui.input("Catchment area (km2)",
                         value=f"{dam['catchment_km2']:g}" if dam["catchment_km2"] else "") \
                    .classes("w-40").props("dense").mark("dam-area") \
                    .on("blur", lambda e: self._set_dam(dam, "catchment_km2",
                                                        _number(e.sender.value)))

            with ui.row().classes("w-full items-end gap-3 no-wrap pt-1"):
                ui.select({month + 1: name for month, name in enumerate(lakerecord.MONTHS)},
                          value=int(dam["water_year_start"]), label="Water year starts",
                          on_change=lambda e: self._set_dam(dam, "water_year_start",
                                                            int(e.value))) \
                    .classes("w-44").props("dense").mark("dam-water-year")
                ui.label("The water year every analysis labels its annual maxima by.") \
                    .classes("text-xs text-muted")

    def _dam_path(self, label, holder, dam, key, *, mark="", redraw=False, suffixes=()):
        """A path input; ``redraw`` for one whose value changes what else is asked."""
        folder = STATE.study.folder
        return path_input(label, holder.get(key) or "", base=folder,
                          base_name="the study folder", suffixes=suffixes, mark=mark,
                          on_commit=lambda text: self._set_dam(
                              holder, key, studies.portable(folder, text)
                              if str(text).strip() else "", dam=dam, redraw=redraw))

    def _set_dam(self, holder, key, value, *, dam=None, redraw=False) -> None:
        """Change one dam input and save the study - ``dam`` is the whole section
        when ``holder`` is a part of it, such as the overlay gauge."""
        dam = holder if dam is None else dam
        if holder.get(key) == value:
            return
        holder[key] = value
        self._save_dam(dam)
        if redraw:
            self.redraw()

    def _save_dam(self, dam) -> None:
        study = STATE.study
        dams.store(study, dam)
        try:
            study.save()
        except OSError as exc:
            ui.notify(f"Could not save {study.path.name}: {exc}", type="warning")

    def _toggle_overlay(self, dam, on) -> None:
        dam["overlay"] = ({"file": "", "below": None, "reconnect_margin": 0.10}
                          if on else None)
        self._save_dam(dam)
        self.redraw()

    def _set_pan(self, dam, text) -> None:
        values = [_number(token) for token in str(text).replace(",", " ").split()]
        if len(values) != 12 or None in values:
            ui.notify("Give twelve pan factors, January to December", type="warning")
            return
        self._set_dam(dam, "pan_factors", values)

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
