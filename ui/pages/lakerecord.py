"""Lake record: catchment rainfall, homogenised levels, antecedent storage, inflow.

Three steps that feed each other and end in the ``lake_config.json`` files the
Monte Carlo runs sample their antecedent storage from, and a fourth off the same
record:

1. **Catchment rainfall** from a catchment shapefile and the folder of daily
   AWAP / AWRA-L grids, area-weighted - run here, in the launcher.
2. **Homogenisation**: the recorded lake levels re-routed through each target
   rating - ``util/HomogeniseLakeLevels.py`` under Bryan's interpreter.
3. **Antecedent storage**: the storm behind each year's maximum, the lake volume
   when it started, the S-curve and the lake configs -
   ``util/AntecedentStorage.py`` under Bryan's interpreter.
4. **Inflow record**: the inflow step 2 derives, as annual maximum peak inflow
   and burst volumes (for an inflow flood frequency analysis) and as event
   hydrographs (for calibration) - ``util/InflowRecord.py``.

What steps 2 to 4 all read - the gauge exports, the storage table, the rating
register, the evaporation, the longest step and the water year - is set once,
in "The lake level record" card at the top.

Everything is kept in the study file (open it on the Report page) and saved as
it changes; each run leaves its job file beside its outputs. See
``core/lakerecord.py``.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pandas as pd
from nicegui import app, run, ui

from core import awap, dam as dams, lakerecord, wordtable
from core.paths import clean_path_text
from layout import no_study, page_frame, severity_banner
from state import STATE
from theme import house_echart


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


def lake_record_page() -> None:
    with page_frame("Lake record"):
        _LakeRecordView().build()


def _float(value, default=None):
    try:
        return float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return default


class _LakeRecordView:
    def __init__(self) -> None:
        self.study = STATE.study
        self.section = lakerecord.settings(self.study) if self.study else None

    def build(self) -> None:
        if self.study is None:
            no_study("The lake record is kept in the study file.")
            return
        ui.label(f"Kept in {self.study.path}; paths are stored relative to it."
                 ).classes("text-xs text-muted")
        self._record_card()
        self._rainfall_card()
        self._homogenise_card()
        self._antecedent_card()
        self._inflow_card()

    def save(self) -> None:
        lakerecord.store(self.study, self.section)
        try:
            self.study.save()
        except OSError as exc:
            ui.notify(f"Could not save {self.study.path.name}: {exc}", type="warning")

    def _path_input(self, label, holder, key, *, mark=""):
        box = ui.input(label, value=holder.get(key) or "").classes("w-full").props("dense")
        box.on("blur", lambda e: self._set_path(holder, key, e.sender.value))
        if mark:
            box.mark(mark)
        return box

    def _set_path(self, holder, key, text) -> None:
        value = lakerecord.keep_path(self.study, text) if str(text).strip() else ""
        if holder.get(key) != value:
            holder[key] = value
            self.save()

    def _set(self, holder, key, value) -> None:
        if holder.get(key) != value:
            holder[key] = value
            self.save()

    def _blocked(self, step) -> bool:
        problems = lakerecord.problems_before_running(self.study, self.section, step)
        if problems:
            ui.notify("Not ready:\n" + "\n".join(problems), type="warning", multi_line=True)
        return bool(problems)

    # -- 1. rainfall ---------------------------------------------------------

    def _rainfall_card(self) -> None:
        r = self.section["rainfall"]
        with ui.card().classes("w-full"):
            ui.label("1. Catchment rainfall").classes("text-lg font-bold")
            ui.label("The daily catchment average of the AWAP / AWRA-L grids (one netCDF "
                     "file per year, as downloaded). A day D is the 24 hours to 9 am on D, "
                     "as the grids stamp it and as the antecedent search expects. The "
                     "series is kept in the study folder and ships with the model; the "
                     "grids stay on this computer, so anyone without them uses the series "
                     "as it is.").classes("text-xs text-muted")
            polygon = (f", where {r['field']} is {r['value']}" if r["field"] else
                       ", every polygon")
            ui.label(f"Catchment: {r['shapefile'] or 'not given'}{polygon} - one of the dam "
                     f"inputs above.").classes("text-sm text-body").mark("rain-catchment")
            with ui.row().classes("w-full gap-2 no-wrap"):
                with ui.column().classes("grow gap-1"):
                    ui.input("Folder of daily grids on this computer (yours, not the "
                             "study's)", value=STATE.settings.awap_folder) \
                        .classes("w-full").props("dense").mark("rain-grids") \
                        .on("blur", lambda e: self._set_grids(e.sender.value))
                    with ui.row().classes("w-full gap-2 no-wrap"):
                        ui.input("Files", value=r["pattern"]).classes("w-32").props("dense") \
                            .on("blur", lambda e: self._set(r, "pattern", e.sender.value.strip()
                                                            or "*.nc"))
                        ui.input("Variable (blank: the one daily grid)", value=r["variable"]) \
                            .classes("grow").props("dense") \
                            .on("blur", lambda e: self._set(r, "variable", e.sender.value.strip()))
                        ui.select({awap.AREA: "Area-weighted", awap.CENTRE: "Cell centres"},
                                  value=r["weighting"], label="Weighting",
                                  on_change=lambda e: self._set(r, "weighting", e.value)) \
                            .classes("w-40").props("dense")
            self._path_input("Write the series to", r, "output")
            with ui.row().classes("items-center gap-2"):
                ui.button("Make the rainfall series", icon="water_drop",
                          on_click=self._make_rainfall).mark("run-rainfall")
                self.rain_status = ui.label("").classes("text-sm text-body")
            self.rain_box = ui.column().classes("w-full")
            self._draw_rainfall()

    def _set_grids(self, text) -> None:
        text = clean_path_text(text or "")
        if text != STATE.settings.awap_folder:
            STATE.settings.awap_folder = text
            STATE.settings.save()

    async def _make_rainfall(self) -> None:
        grids = STATE.settings.awap_folder
        problems = lakerecord.problems_before_running(self.study, self.section, "rainfall",
                                                      grids)
        if problems:
            ui.notify("Not ready:\n" + "\n".join(problems), type="warning", multi_line=True)
            return
        r = self.section["rainfall"]
        self.rain_status.text = "reading the grids..."
        progress = {"text": ""}

        def work():
            catchment = awap.read_catchment(lakerecord.path_of(self.study, r["shapefile"]),
                                            r["field"], r["value"])
            files = awap.rainfall_files(grids, r["pattern"] or "*.nc")
            series = awap.catchment_rainfall(
                files, catchment, r["variable"], r["weighting"],
                progress=lambda n, total, name: progress.update(text=f"{n + 1} of {total}: {name}"))
            awap.write_series(series, lakerecord.path_of(self.study, r["output"]),
                              source=f"catchment {Path(r['shapefile']).name}"
                                     f"{' ' + r['field'] + '=' + r['value'] if r['field'] else ''}"
                                     f", {r['weighting']} weighting")
            return series

        timer = ui.timer(0.5, lambda: setattr(self.rain_status, "text", progress["text"]))
        try:
            series = await _off_thread(work)
        except awap.AwapError as exc:
            self.rain_status.text = ""
            ui.notify(str(exc), type="negative", multi_line=True)
            return
        finally:
            timer.cancel()
        self.rain_status.text = (f"{len(series.frame):,} days, {series.cells} cells - "
                                 f"written")
        for problem in series.problems[:5]:
            ui.notify(problem, type="warning")
        self._draw_rainfall()

    def _draw_rainfall(self) -> None:
        self.rain_box.clear()
        path = lakerecord.path_of(self.study, self.section["rainfall"]["output"])
        if not path or not path.is_file():
            return
        try:
            rain = awap.read_series(path)
        except Exception:                              # noqa: BLE001 - a view, never crash
            return
        totals = rain.groupby(rain.index.year).sum()
        with self.rain_box:
            ui.label(f"{path.name}: {rain.index.min():%d %b %Y} to {rain.index.max():%d %b %Y}"
                     f", mean {totals.mean():,.0f} mm a year").classes("text-sm")
            house_echart({
                "tooltip": {"trigger": "axis"},
                "grid": {"left": 60, "right": 20, "top": 20, "bottom": 30},
                "xAxis": {"type": "category", "data": [str(y) for y in totals.index]},
                "yAxis": {"type": "value", "name": "mm a year", "nameLocation": "middle",
                          "nameGap": 45},
                "series": [{"type": "bar", "name": "annual total",
                            "data": [round(float(v), 1) for v in totals]}],
            }).classes("w-full h-56").mark("rain-chart")

    # -- the record, shared by steps 2-4 -------------------------------------------

    def _record_card(self) -> None:
        """The dam inputs the steps read, shown - they are edited on the Study page.

        Each with whether it is there, so a missing file is seen before a step is
        run rather than after.
        """
        h = self.section["homogenise"]
        r, i = self.section["rainfall"], self.section["inflow"]
        with ui.card().classes("w-full").mark("record-card"):
            with ui.row().classes("w-full items-center justify-between no-wrap"):
                ui.label(lakerecord.RECORD_CARD).classes("text-lg font-bold")
                ui.button("Edit on the Study page", icon="edit",
                          on_click=lambda: ui.navigate.to("/study")) \
                    .props("flat dense no-caps").mark("edit-dam")
            ui.label("The recorded lake levels, what turns them into storage and release, "
                     "and the catchment. Every step below reads some of these.") \
                .classes("text-xs text-muted")
            gauges = [g for g in h["gauges"] if str(g).strip()]
            rows = [(f"Gauge export {n + 1}" if len(gauges) > 1 else "Gauge export", g)
                    for n, g in enumerate(gauges)] or [("Gauge exports", "")]
            overlay = h.get("overlay")
            if overlay:
                rows.append((f"Overlay gauge, below {overlay.get('below')} m",
                             overlay.get("file")))
            single = dams.single_rating(h)
            rows += [("Storage table", h["storage"]),
                     ("One rating" + (f", FSL {h['register_fsl']:g} m" if single and
                                      h.get("register_fsl") else "")
                      if single else "Rating register", h["register"]),
                     ("Evaporation (SILO)", h["evaporation"]),
                     ("Catchment shapefile", r["shapefile"])]
            with ui.element("div").classes("w-full grid gap-x-4 gap-y-0") \
                    .style("grid-template-columns: 15rem 1fr"):
                for label, value in rows:
                    ui.label(label).classes("text-sm text-muted")
                    self._file_state(value)
                ui.label("Catchment area").classes("text-sm text-muted")
                ui.label(f"{i['catchment_km2']:g} km2" if i["catchment_km2"] else "not given") \
                    .classes("text-sm text-body")
                ui.label("Water year starts").classes("text-sm text-muted")
                ui.label(f"{lakerecord.MONTHS[lakerecord.shared_water_year(self.section) - 1]}"
                         f", unless a step below says otherwise") \
                    .classes("text-sm text-body").mark("record-water-year")
            self.water_year_box = ui.column().classes("w-full")
            self._draw_water_year_note()
            with ui.row().classes("w-full items-end gap-2 no-wrap pt-2"):
                ui.input("Longest step", value=h["step"]).classes("w-28").props("dense") \
                    .on("blur", lambda e: self._set(h, "step", e.sender.value.strip() or "1h"))
                ui.label("Gaps in the level record longer than this are filled, for "
                         "steps 2 to 4.").classes("text-xs text-muted")

    # -- the water year of each step ---------------------------------------------------

    def _water_year_select(self, step) -> None:
        """The step's own water year, or the study's - blank means the study's."""
        shared = lakerecord.MONTHS[lakerecord.shared_water_year(self.section) - 1]
        options = {0: f"{shared} (the study's)",
                   **{month + 1: name for month, name in enumerate(lakerecord.MONTHS)}}
        ui.select(options, value=int(self.section[step].get("water_year") or 0),
                  label="Water year starts",
                  on_change=lambda e: self._set_water_year(step, e.value)) \
            .classes("w-48").props("dense").mark(f"water-year-{step}")

    def _set_water_year(self, step, value) -> None:
        self._set(self.section[step], "water_year", int(value) or None)
        self._draw_water_year_note()

    def _draw_water_year_note(self) -> None:
        self.water_year_box.clear()
        note = lakerecord.water_year_note(self.section)
        if note:
            with self.water_year_box:
                severity_banner("warn", note, "Each step's water year is chosen on its "
                                              "card; blank takes the study's.")

    def _file_state(self, value) -> None:
        """A dam input's path, and whether it is there."""
        with ui.row().classes("items-center gap-2 no-wrap"):
            if not str(value or "").strip():
                ui.label("not given").classes("text-sm text-muted")
                return
            path = lakerecord.path_of(self.study, value)
            found = path is not None and path.is_file()
            ui.label(str(value)).classes("mono text-sm text-body")
            ui.label("found" if found else "not found") \
                .classes(f"text-xs {'text-positive' if found else 'text-negative'}")

    # -- 2. homogenisation ---------------------------------------------------------

    def _homogenise_card(self) -> None:
        h = self.section["homogenise"]
        with ui.card().classes("w-full"):
            ui.label("2. Homogenisation").classes("text-lg font-bold")
            ui.label("The recorded levels, re-routed through each target rating: the net "
                     "inflow is derived against the rating in force at each step (the "
                     "register), then routed through the target. It reads the dam inputs at "
                     "the top of the page.").classes("text-xs text-muted")
            with ui.row().classes("w-full items-end gap-4"):
                ui.checkbox("Recession correction", value=bool(h["recession_correction"]),
                            on_change=lambda e: self._set(h, "recession_correction", e.value))
                self._water_year_select("homogenise")
            ui.label("Target ratings").classes("text-sm font-bold pt-1")
            self.targets_box = ui.column().classes("w-full gap-1")
            self._draw_targets()
            self._path_input("Write the results under", h, "out")
            with ui.row().classes("items-center gap-2"):
                ui.button("Homogenise", icon="timeline", on_click=self._homogenise) \
                    .mark("run-homogenise")
                self.h_status = ui.label("").classes("text-sm text-body")
            self.h_box = ui.column().classes("w-full")
            self._draw_homogenised(lakerecord.last_summary(self.study, self.section,
                                                           "homogenise"))

    def _draw_targets(self) -> None:
        targets = self.section["homogenise"]["targets"]
        self.targets_box.clear()
        with self.targets_box:
            for index, target in enumerate(targets):
                with ui.row().classes("w-full items-center gap-2 no-wrap"):
                    ui.input("Name", value=target.get("name", "")).classes("w-32") \
                        .props("dense").mark(f"target-name-{index}") \
                        .on("blur", lambda e, t=target: self._set(t, "name",
                                                                  e.sender.value.strip()))
                    with ui.element("div").classes("grow"):
                        self._path_input("Rating (.sq, or level,flow .csv)", target, "rating")
                    ui.input("FSL (m AHD)", value=f"{target.get('fsl') or ''}") \
                        .classes("w-28").props("dense") \
                        .on("blur", lambda e, t=target: self._set(t, "fsl", _float(e.sender.value)))
                    ui.button(icon="close", on_click=lambda _, i=index: self._drop_target(i)) \
                        .props("flat dense round")
            ui.button("Add a target rating", icon="add", on_click=self._add_target) \
                .props("flat dense no-caps").mark("add-target")

    def _add_target(self) -> None:
        self.section["homogenise"]["targets"].append({"name": "", "rating": "", "fsl": None})
        self.save()
        self._draw_targets()

    def _drop_target(self, index) -> None:
        self.section["homogenise"]["targets"].pop(index)
        self.save()
        self._draw_targets()

    async def _homogenise(self) -> None:
        if self._blocked("homogenise"):
            return
        self.h_status.text = "homogenising - about half a minute per 500,000 steps..."
        result = await _off_thread(lakerecord.homogenise, self.study,
                                   copy.deepcopy(self.section), STATE.settings.bryan_python)
        self.h_status.text = ""
        if not result.ok:
            self._failed(self.h_box, result)
            return
        self._draw_homogenised(result.summary)

    def _draw_homogenised(self, summary) -> None:
        self.h_box.clear()
        if not summary:
            return
        with self.h_box:
            record = summary.get("record", {})
            ui.label(f"{record.get('steps', 0):,} steps, {record.get('start')} to "
                     f"{record.get('end')}, from {', '.join(record.get('gauges', []))}"
                     ).classes("text-sm").mark("homogenised-summary")
            for note in summary.get("notes", []):
                severity_banner("info", note)
            series, years = [], None
            for name, ams in lakerecord.annual_maxima(summary).items():
                years = [str(y) for y in ams["WaterYear"]]
                if not series:
                    series.append({"type": "scatter", "name": "recorded",
                                   "data": [round(float(v), 3) for v in ams["Level_max"]]})
                series.append({"type": "line", "name": name, "symbolSize": 4,
                               "data": [round(float(v), 3) for v in ams["New_Level_max"]]})
            if series:
                house_echart({
                    "tooltip": {"trigger": "axis"}, "legend": {"top": 0},
                    "grid": {"left": 60, "right": 20, "top": 30, "bottom": 30},
                    "xAxis": {"type": "category", "data": years},
                    "yAxis": {"type": "value", "scale": True, "name": "Annual maximum (m AHD)",
                              "nameLocation": "middle", "nameGap": 45},
                    "series": series}).classes("w-full h-72").mark("ams-chart")

    # -- 3. antecedent storage --------------------------------------------------------

    def _antecedent_card(self) -> None:
        a = self.section["antecedent"]
        s = a["settings"]
        with ui.card().classes("w-full"):
            ui.label("3. Antecedent storage").classes("text-lg font-bold")
            ui.label("For each year's maximum, the rarest 1-5 day burst in the days before "
                     "it, the homogenised lake volume where the burst and its pre-burst "
                     "started, and the S-curve fitted to them with the ceiling at the full "
                     "supply volume - written as the lake configs the Monte Carlo runs "
                     "sample.").classes("text-xs text-muted")
            with ui.row().classes("w-full items-end gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._path_input("Rainfall series (blank: step 1's)", a, "rainfall")
                with ui.element("div").classes("grow"):
                    self._path_input("IFD (duration_h against '1 in X' columns, mm)", a, "ifd",
                                     mark="ifd")
                self._water_year_select("antecedent")
            with ui.row().classes("w-full items-end gap-2 no-wrap"):
                for key, label, width in (("window_days", "Search window (d)", "w-32"),
                                          ("threshold_fraction", "Significance (x 1 in 2)", "w-40"),
                                          ("preburst_mm", "Pre-burst edge (mm/d)", "w-40"),
                                          ("cunnane_a", "Cunnane a", "w-24"),
                                          ("round_ml", "Round floor to (ML)", "w-36")):
                    ui.input(label, value=f"{s[key]}").classes(width).props("dense") \
                        .on("blur", lambda e, k=key: self._set(s, k, _float(e.sender.value, s[k])))
                ui.input("Restriction factors by duration (d = factor)",
                         value=", ".join(f"{d} = {f}" for d, f in s["restriction"].items())) \
                    .classes("grow").props("dense") \
                    .on("blur", lambda e: self._set_restriction(e.sender.value))
            with ui.row().classes("items-center gap-4"):
                ui.label("Series").classes("text-sm")
                for basis, label in (("burst", "burst - storm without its pre-burst"),
                                     ("storm", "storm - with the pre-burst")):
                    ui.checkbox(label, value=basis in a["bases"],
                                on_change=lambda e, b=basis: self._set_basis(b, e.value))
            self._path_input("Write the results under", a, "out")
            with ui.row().classes("items-center gap-2"):
                ui.button("Antecedent storage", icon="show_chart", on_click=self._antecedent) \
                    .mark("run-antecedent")
                self.a_status = ui.label("").classes("text-sm text-body")
            self.a_box = ui.column().classes("w-full")
            self._draw_antecedent(lakerecord.last_summary(self.study, self.section,
                                                          "antecedent"))

    def _set_restriction(self, text) -> None:
        pairs = {}
        for part in str(text).split(","):
            if "=" in part:
                day, _, factor = part.partition("=")
                if _float(factor) is not None and day.strip().isdigit():
                    pairs[day.strip()] = _float(factor)
        if pairs:
            self._set(self.section["antecedent"]["settings"], "restriction", pairs)
            self._set(self.section["antecedent"]["settings"], "durations_d",
                      sorted(int(day) for day in pairs))

    def _set_basis(self, basis, on) -> None:
        bases = [b for b in self.section["antecedent"]["bases"] if b != basis]
        if on:
            bases.append(basis)
        self._set(self.section["antecedent"], "bases", sorted(bases))

    async def _antecedent(self) -> None:
        if self._blocked("antecedent"):
            return
        self.a_status.text = "homogenising and searching the rainfall..."
        result = await _off_thread(lakerecord.antecedent, self.study,
                                   copy.deepcopy(self.section), STATE.settings.bryan_python)
        self.a_status.text = ""
        if not result.ok:
            self._failed(self.a_box, result)
            return
        self._draw_antecedent(result.summary)

    def _draw_antecedent(self, summary) -> None:
        self.a_box.clear()
        if not summary:
            return
        with self.a_box:
            for target in summary.get("targets", []):
                ui.label(f"{target['name']}: {target['qualified']} of {target['years']} water "
                         f"years give a sample; full supply volume {target['fsv_ML']:,.0f} ML"
                         ).classes("text-sm font-bold").mark(f"antecedent-{target['name']}")
                if target.get("years_outside_rainfall"):
                    severity_banner("warn", f"{target['years_outside_rainfall']} annual "
                                    f"maxima fall outside the rainfall series and cannot give "
                                    f"a sample.")
                rows = [{"basis": basis, **{k: fit[k] for k in ("n", "k", "z0", "H")},
                         "Vf": f"{fit['Vf']:,.0f}", "Vc": f"{fit['Vc']:,.0f}",
                         "median": f"{fit['mean_z0_ML']:,.0f}",
                         "config": target["configs"][basis]}
                        for basis, fit in target["fits"].items()]
                columns = [{"name": key, "label": label, "field": key, "align": "left"}
                           for key, label in (("basis", "Series"), ("n", "n"), ("k", "k"),
                                              ("z0", "z0"), ("H", "H"), ("Vf", "Vf (ML)"),
                                              ("Vc", "Vc (ML)"), ("median", "At z = 0 (ML)"),
                                              ("config", "Lake config"))]
                ui.table(columns=columns, rows=rows, row_key="basis") \
                    .props("dense flat bordered").classes("w-full")
                series = []
                for basis in target["fits"]:
                    try:
                        samples, curve = lakerecord.scurve_points(target, basis)
                    except (OSError, KeyError, ValueError):
                        continue
                    series.append({"type": "scatter", "name": f"{basis} samples",
                                   "data": [[round(z, 4), round(v, 1)] for z, v in samples]})
                    series.append({"type": "line", "name": f"{basis} fit", "symbol": "none",
                                   "data": [[round(z, 4), round(v, 1)] for z, v in curve]})
                if series:
                    house_echart({
                        "tooltip": {"trigger": "item"}, "legend": {"top": 0},
                        "grid": {"left": 70, "right": 20, "top": 30, "bottom": 40},
                        "xAxis": {"type": "value", "name": "Standard normal variate",
                                  "nameLocation": "middle", "nameGap": 25},
                        "yAxis": {"type": "value", "name": "Antecedent storage (ML)",
                                  "nameLocation": "middle", "nameGap": 55},
                        "series": series}).classes("w-full h-72") \
                        .mark(f"scurve-{target['name']}")

    # -- 4. inflow record ---------------------------------------------------------------

    def _inflow_card(self) -> None:
        i = self.section["inflow"]
        with ui.card().classes("w-full"):
            ui.label("4. Inflow record").classes("text-lg font-bold")
            ui.label("Step 2's derived inflow - the change in storage plus the release "
                     "through the rating in force at each step - as each water year's peak "
                     "inflow and largest burst volumes, and as event hydrographs for "
                     "calibration. It reads the dam inputs at the top of the page, but none "
                     "of step 2's target ratings.") \
                .classes("text-xs text-muted")
            ui.label("On a recession above full supply the balance often goes negative: the "
                     "gates released more than the rating says. The recession correction "
                     "books that as release, so the corrected inflow sits at zero there; "
                     "each hydrograph also keeps the uncorrected inflow and flags the "
                     "intervals where the release is uncertain. The rising limb and the peak "
                     "are unaffected either way.").classes("text-xs text-muted")
            with ui.row().classes("w-full items-end gap-2 no-wrap"):
                ui.checkbox("Keep the lake evaporation in the inflow",
                            value=bool(i["evaporation"]),
                            on_change=lambda e: self._set(i, "evaporation", e.value)) \
                    .mark("inflow-evaporation")
                ui.checkbox("Recession correction", value=bool(i["recession_correction"]),
                            on_change=lambda e: self._set(i, "recession_correction", e.value))
                self._water_year_select("inflow")
                ui.input("Peak averaged over", value=i["smoothing"] or "") \
                    .classes("w-36").props("dense") \
                    .on("blur", lambda e: self._set(i, "smoothing", e.sender.value.strip()))
                ui.input("Burst durations (h)",
                         value=" ".join(f"{d:g}" for d in i["durations_h"])) \
                    .classes("w-40").props("dense") \
                    .on("blur", lambda e: self._set_durations(e.sender.value))
                ui.label(f"Runoff depths over {i['catchment_km2']:g} km2, the catchment area "
                         f"in the dam inputs." if i["catchment_km2"] else
                         "No catchment area in the dam inputs, so no runoff depths.") \
                    .classes("text-xs text-muted").mark("inflow-area")
            with ui.row().classes("w-full items-end gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._path_input("Rainfall series, beside each burst (blank: step 1's)",
                                     i, "rainfall")
                for key, label in (("top_events", "Hydrographs of the largest"),
                                   ("before_days", "from days before"),
                                   ("after_days", "to days after")):
                    ui.input(label, value=f"{i[key]:g}").classes("w-40").props("dense") \
                        .on("blur", lambda e, k=key: self._set(
                            i, k, _float(e.sender.value, i[k])))
            ui.textarea("More hydrographs: name, start, end - one a line "
                        "(e.g. Jan 2013, 2013-01-20, 2013-02-05 12:00)",
                        value=lakerecord.events_text(i["events"])) \
                .classes("w-full").props("dense autogrow").mark("inflow-events") \
                .on("blur", lambda e: self._set_events(e.sender.value))
            self._path_input("Write the results under", i, "out")
            with ui.row().classes("items-center gap-2"):
                ui.button("Inflow record", icon="waves", on_click=self._inflow) \
                    .mark("run-inflow")
                self.i_status = ui.label("").classes("text-sm text-body")
            self.i_box = ui.column().classes("w-full")
            self._draw_inflow(lakerecord.last_summary(self.study, self.section, "inflow"))

    def _set_durations(self, text) -> None:
        values = [_float(token) for token in str(text).replace(",", " ").split()]
        if not values or None in values or min(values) <= 0:
            ui.notify("Give the burst durations in hours, e.g. 24 36 48 72", type="warning")
            return
        self._set(self.section["inflow"], "durations_h", sorted({int(v) for v in values}))

    def _set_events(self, text) -> None:
        events, problems = lakerecord.parse_events(text)
        for problem in problems:
            ui.notify(f"Hydrographs, {problem}", type="warning")
        if not problems:
            self._set(self.section["inflow"], "events", events)

    async def _inflow(self) -> None:
        if self._blocked("inflow"):
            return
        self.i_status.text = "deriving the inflow..."
        result = await _off_thread(lakerecord.inflow, self.study,
                                   copy.deepcopy(self.section), STATE.settings.bryan_python)
        self.i_status.text = ""
        if not result.ok:
            self._failed(self.i_box, result)
            return
        self._draw_inflow(result.summary)

    def _draw_inflow(self, summary) -> None:
        self.i_box.clear()
        if not summary:
            return
        ams = lakerecord.inflow_maxima(summary)
        with self.i_box:
            record = summary.get("record", {})
            ui.label(f"{summary.get('years', 0)} water years "
                     f"({summary.get('complete_years', 0)} complete), {record.get('start')} to "
                     f"{record.get('end')}").classes("text-sm").mark("inflow-summary")
            for note in summary.get("notes", []):
                severity_banner("info", note)
            if ams is not None and len(ams):
                complete = ams["Complete"].astype(str).str.lower() == "true"
                house_echart({
                    "tooltip": {"trigger": "axis"}, "legend": {"top": 0},
                    "grid": {"left": 60, "right": 20, "top": 30, "bottom": 30},
                    "xAxis": {"type": "category", "data": [str(p) for p in ams["Period"]]},
                    "yAxis": {"type": "value", "name": "Annual maximum inflow (m3/s)",
                              "nameLocation": "middle", "nameGap": 45},
                    "series": [{"type": "bar", "name": name, "stack": "peak",
                                "data": [round(float(v), 1) if keep == want else None
                                         for v, keep in zip(ams["Peak_inflow_m3s"], complete)]}
                               for name, want in (("complete year", True),
                                                  ("part year", False))]}) \
                    .classes("w-full h-64").mark("inflow-ams-chart")
                self._ams_actions(summary, ams)
                self._depth_check(ams, summary.get("settings", {}))
            events = summary.get("hydrographs") or []
            if events:
                ui.label("Event hydrographs").classes("text-sm font-bold pt-2")
                names = {index: f"{item['name']} ({item['peak_m3s']:,.0f} m3/s)"
                         for index, item in enumerate(events)}
                ui.select(names, value=0, label="Hydrograph",
                          on_change=lambda e: self._draw_hydrograph(holder, events[e.value])) \
                    .classes("w-96").props("dense").mark("inflow-hydrograph-select")
                holder = ui.column().classes("w-full")
                self._draw_hydrograph(holder, events[0])
            if lakerecord.intervals_path(summary) and lakerecord.intervals_path(summary).is_file():
                self._record_window(summary, ams)

    def _ams_actions(self, summary, ams) -> None:
        path = Path(summary["ams"])
        durations = summary.get("settings", {}).get("durations_h") or []
        with ui.row().classes("items-center gap-2"):
            ui.button("Download AMS CSV", icon="download",
                      on_click=lambda: ui.download.file(path, path.name)) \
                .props("outline dense no-caps").mark("inflow-ams-download")
            ui.button("Copy AMS for Word", icon="content_copy",
                      on_click=lambda: self._copy_ams(ams, durations, rich=True)) \
                .props("outline dense no-caps").mark("inflow-ams-copy")
            ui.button("Copy as text", on_click=lambda: self._copy_ams(ams, durations, rich=False)) \
                .props("flat dense no-caps")
            ui.label(f"Written to {path}").classes("text-xs text-muted")

    async def _copy_ams(self, ams, durations, *, rich: bool) -> None:
        table = lakerecord.ams_table(ams, durations)
        text = wordtable.to_text(table)
        if not rich:
            ui.clipboard.write(text)
            ui.notify("Copied as text")
            return
        fragment = wordtable.to_html(table, self.study.extra.get("word"))
        outcome = await ui.run_javascript(
            wordtable.clipboard_script(wordtable.clipboard_document(fragment), text),
            timeout=5.0)
        if outcome == "html":
            ui.notify("Copied - paste into Word")
        elif outcome == "text":
            ui.notify("This browser would only take text; copied as text", type="warning")
        else:
            ui.notify(f"Could not copy: {outcome}", type="negative")

    # the record between two dates

    def _record_window(self, summary, ams) -> None:
        i = self.section["inflow"]
        window = i.setdefault("window", {"start": "", "end": "", "step": ""})
        if not window.get("start") and ams is not None and len(ams):
            peak = pd.Timestamp(ams.loc[ams["Peak_inflow_m3s"].idxmax(), "Peak_time"])
            window.update(start=f"{(peak - pd.Timedelta(days=3)).floor('D'):%Y-%m-%d}",
                          end=f"{(peak + pd.Timedelta(days=7)).ceil('D'):%Y-%m-%d}")
        record = summary.get("record", {})
        ui.label("The record between two dates").classes("text-sm font-bold pt-2")
        ui.label(f"Anywhere from {record.get('start')} to {record.get('end')}. Blank time step: "
                 "the record's own intervals, with the inflow averaged as the peaks are. A "
                 "time step (15min, 1h, 1D): the mean over each step, stamped at its end as a "
                 "model hydrograph is, so each step keeps its volume.") \
            .classes("text-xs text-muted")
        with ui.row().classes("items-end gap-2"):
            start = ui.input("Start", value=window["start"]).classes("w-44").props("dense") \
                .mark("window-start")
            end = ui.input("End", value=window["end"]).classes("w-44").props("dense") \
                .mark("window-end")
            step = ui.input("Time step (blank: as recorded)", value=window.get("step", "")) \
                .classes("w-52").props("dense").mark("window-step")
            ui.button("Plot", icon="show_chart",
                      on_click=lambda: self._window(summary, start.value, end.value, step.value,
                                                    holder, save=False)) \
                .props("dense").mark("window-plot")
            ui.button("Save CSV", icon="download",
                      on_click=lambda: self._window(summary, start.value, end.value, step.value,
                                                    holder, save=True)) \
                .props("outline dense").mark("window-save")
        holder = ui.column().classes("w-full")

    async def _window(self, summary, start, end, step, holder, *, save: bool) -> None:
        try:
            first, last, delta = lakerecord.parse_window(start, end, step)
        except ValueError as exc:
            ui.notify(str(exc), type="warning")
            return
        self._set(self.section["inflow"], "window",
                  {"start": str(start).strip(), "end": str(end).strip(),
                   "step": str(step or "").strip()})
        path = lakerecord.intervals_path(summary)

        def work():
            table = lakerecord.record_window(lakerecord.read_record(path), first, last, delta)
            written = None
            if save and not table.empty:
                written = lakerecord.write_window(
                    table, lakerecord.window_file(path.parent, first, last, delta))
            return table, written

        try:
            table, written = await _off_thread(work)
        except OSError as exc:
            ui.notify(f"Could not write the CSV: {exc}", type="negative")
            return
        holder.clear()
        if table.empty:
            ui.notify("The record has nothing between those dates", type="warning")
            return
        with holder:
            volume = table["Volume_ML"].sum()
            ui.label(f"{len(table):,} rows, {table.index[0]:%d %b %Y %H:%M} to "
                     f"{table.index[-1]:%d %b %Y %H:%M}; peak {table['Inflow_m3s'].max():,.0f} "
                     f"m3/s, inflow volume {volume:,.0f} ML"
                     + (f"; written to {written}" if written else "")) \
                .classes("text-xs text-muted").mark("window-note")
            self._flow_chart(table, "window-chart")
            if written:
                ui.download.file(written, written.name)
                ui.notify(f"Wrote {written}")

    def _depth_check(self, ams, settings) -> None:
        """Runoff depth beside the catchment rain: the one check the inflow did not make."""
        durations = settings.get("durations_h") or []
        if not durations:
            return
        hours = max(durations)
        depth, rain = f"Depth_{hours}h_mm", f"Rain_{hours}h_mm"
        if depth not in ams or rain not in ams:
            return
        over = ams[ams[depth] > ams[rain]]
        if len(over):
            severity_banner("warn", f"{len(over)} water years show more {hours} h runoff than "
                            f"rain fell ({', '.join(map(str, over['Period'][:6]))}) - "
                            f"check the level record there.")
        house_echart({
            "tooltip": {"trigger": "axis"}, "legend": {"top": 0},
            "grid": {"left": 60, "right": 20, "top": 30, "bottom": 30},
            "xAxis": {"type": "category", "data": [str(p) for p in ams["Period"]]},
            "yAxis": {"type": "value", "name": f"Largest {hours} h burst (mm)",
                      "nameLocation": "middle", "nameGap": 45},
            "series": [{"type": "bar", "name": "catchment rain (the burst's days and the day before)",
                        "data": [None if v != v else round(float(v), 1) for v in ams[rain]]},
                       {"type": "bar", "name": "runoff depth",
                        "data": [None if v != v else round(float(v), 1) for v in ams[depth]]}]}) \
            .classes("w-full h-56").mark("inflow-depth-chart")

    def _draw_hydrograph(self, holder, item) -> None:
        holder.clear()
        frame = lakerecord.hydrograph(item)
        if frame is None or frame.empty:
            return
        share = float((frame["Release_uncertain"].astype(str).str.lower() == "true").mean())
        if "Inflow_smoothed_m3s" in frame:
            frame = frame.assign(Inflow_m3s=frame["Inflow_smoothed_m3s"])
        with holder:
            ui.label(f"{Path(item['file']).name}: {item['start']} to {item['end']}; the release "
                     f"is uncertain over {share:.0%} of it"
                     ).classes("text-xs text-muted").mark("inflow-hydrograph-note")
            self._flow_chart(frame, "inflow-hydrograph")

    def _flow_chart(self, table, mark) -> None:
        """Inflow, release and level; the uncorrected inflow dashed where it differs."""
        points = lakerecord.chart_points(table)
        stamps = [f"{t:%Y-%m-%d %H:%M:%S}" for t in points.index]
        uncertain = (points["Uncertain"].astype(bool).tolist() if "Uncertain" in points
                     else [False] * len(points))

        def pairs(values, keep=None):
            keep = keep or [True] * len(stamps)
            return [[t, None if (v != v or not k) else round(float(v), 2)]
                    for t, v, k in zip(stamps, values, keep)]

        house_echart({
            "tooltip": {"trigger": "axis"}, "legend": {"top": 0},
            "grid": {"left": 60, "right": 60, "top": 30, "bottom": 60},
            "xAxis": {"type": "time"},
            "yAxis": [{"type": "value", "name": "m3/s", "nameLocation": "middle",
                       "nameGap": 45},
                      {"type": "value", "name": "Level (m AHD)", "scale": True,
                       "nameLocation": "middle", "nameGap": 45, "splitLine": {"show": False}}],
            "dataZoom": [{"type": "inside"}, {"type": "slider", "height": 18, "bottom": 8}],
            "series": [
                {"type": "line", "name": "inflow", "symbol": "none",
                 "data": pairs(points["Inflow_m3s"])},
                {"type": "line", "name": "inflow, uncorrected where the release is uncertain",
                 "symbol": "none", "lineStyle": {"type": "dashed"},
                 "data": pairs(points["Inflow_uncorrected_m3s"], uncertain)},
                {"type": "line", "name": "release", "symbol": "none",
                 "data": pairs(points["Release_m3s"])},
                {"type": "line", "name": "level", "symbol": "none", "yAxisIndex": 1,
                 "lineStyle": {"width": 1}, "data": pairs(points["Level_m"])},
            ]}).classes("w-full h-80").mark(mark)

    def _failed(self, box, result) -> None:
        box.clear()
        with box:
            errors = [line for line in result.output.splitlines()
                      if line.startswith("ERROR")]
            severity_banner("block", errors[0] if errors else "The run failed - the log is "
                                                              "below.")
            with ui.expansion("Log").classes("w-full"):
                ui.code(result.output[-20000:] or "(no output)").classes("w-full text-xs")
