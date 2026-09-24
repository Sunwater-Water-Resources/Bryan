"""Lake record: catchment rainfall, homogenised levels and antecedent storage.

Three steps that feed each other and end in the ``lake_config.json`` files the
Monte Carlo runs sample their antecedent storage from:

1. **Catchment rainfall** from a catchment shapefile and the folder of daily
   AWAP / AWRA-L grids, area-weighted - run here, in the launcher.
2. **Homogenisation**: the recorded lake levels re-routed through each target
   rating - ``util/HomogeniseLakeLevels.py`` under Bryan's interpreter.
3. **Antecedent storage**: the storm behind each year's maximum, the lake volume
   when it started, the S-curve and the lake configs -
   ``util/AntecedentStorage.py`` under Bryan's interpreter.

Everything is kept in the study file (open it on the Report page) and saved as
it changes; each run leaves its job file beside its outputs. See
``core/lakerecord.py``.
"""

from __future__ import annotations

import copy
from pathlib import Path

from nicegui import app, run, ui

from core import awap, lakerecord
from layout import page_frame, severity_banner
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
            with ui.card().classes("w-full items-center p-8"):
                ui.icon("description").classes("text-5xl text-muted")
                ui.label("Open a study first - the lake record is kept in the study file."
                         ).classes("text-muted")
                ui.button("Go to Report", on_click=lambda: ui.navigate.to("/report"))
            return
        ui.label(f"Kept in {self.study.path}; paths are stored relative to it."
                 ).classes("text-xs text-muted")
        self._rainfall_card()
        self._homogenise_card()
        self._antecedent_card()

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
            with ui.row().classes("w-full gap-2 no-wrap"):
                with ui.column().classes("grow gap-1"):
                    self._path_input("Catchment shapefile (.shp)", r, "shapefile",
                                     mark="rain-shapefile")
                    with ui.row().classes("w-full gap-2 no-wrap"):
                        ui.input("Field (blank: every polygon)", value=r["field"]) \
                            .classes("grow").props("dense") \
                            .on("blur", lambda e: self._set(r, "field", e.sender.value.strip()))
                        ui.input("equal to", value=r["value"]).classes("grow").props("dense") \
                            .on("blur", lambda e: self._set(r, "value", e.sender.value.strip()))
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
        text = str(text or "").strip().strip('"')
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

    # -- 2. homogenisation ----------------------------------------------------------

    def _homogenise_card(self) -> None:
        h = self.section["homogenise"]
        with ui.card().classes("w-full"):
            ui.label("2. Homogenisation").classes("text-lg font-bold")
            ui.label("The recorded levels, re-routed through each target rating: the net "
                     "inflow is derived against the rating in force at each step (the "
                     "register), then routed through the target.").classes("text-xs text-muted")
            ui.textarea("Gauge exports (WMIP / Hydstra), one per line, in the order the gauges "
                        "operated", value="\n".join(h["gauges"])) \
                .classes("w-full").props("dense autogrow").mark("gauges") \
                .on("blur", lambda e: self._set(h, "gauges", [
                    lakerecord.keep_path(self.study, line)
                    for line in str(e.sender.value).splitlines() if line.strip()]))
            overlay = h.get("overlay")
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                ui.checkbox("Overlay gauge below a level", value=bool(overlay),
                            on_change=lambda e: self._toggle_overlay(e.value))
                if overlay:
                    with ui.element("div").classes("grow"):
                        self._path_input("Overlay gauge export", overlay, "file")
                    ui.input("below (m)", value=f"{overlay.get('below', '')}") \
                        .classes("w-28").props("dense") \
                        .on("blur", lambda e: self._set(overlay, "below", _float(e.sender.value)))
                    ui.input("reconnect margin (m)",
                             value=f"{overlay.get('reconnect_margin', 0.10)}") \
                        .classes("w-40").props("dense") \
                        .on("blur", lambda e: self._set(overlay, "reconnect_margin",
                                                        _float(e.sender.value, 0.10)))
            with ui.row().classes("w-full gap-2 no-wrap"):
                for key, label in (("storage", "Storage table (.els: EL, A, V)"),
                                   ("register", "Rating register (.xlsx)"),
                                   ("evaporation", "Evaporation (SILO Data Drill)")):
                    with ui.element("div").classes("grow"):
                        self._path_input(label, h, key)
            with ui.row().classes("w-full items-end gap-2 no-wrap"):
                ui.input("Pan factors, Jan to Dec",
                         value=" ".join(f"{v:g}" for v in h["pan_factors"])) \
                    .classes("grow").props("dense") \
                    .on("blur", lambda e: self._set_pan(e.sender.value))
                ui.input("Longest step", value=h["step"]).classes("w-28").props("dense") \
                    .on("blur", lambda e: self._set(h, "step", e.sender.value.strip() or "1h"))
                ui.select({month + 1: name for month, name in enumerate(lakerecord.MONTHS)},
                          value=int(h["water_year_start"]), label="Water year starts",
                          on_change=lambda e: self._set(h, "water_year_start", e.value)) \
                    .classes("w-40").props("dense")
                ui.checkbox("Recession correction", value=bool(h["recession_correction"]),
                            on_change=lambda e: self._set(h, "recession_correction", e.value))
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

    def _toggle_overlay(self, on) -> None:
        self.section["homogenise"]["overlay"] = ({"file": "", "below": None,
                                                  "reconnect_margin": 0.10} if on else None)
        self.save()
        ui.navigate.reload()

    def _set_pan(self, text) -> None:
        values = [_float(token) for token in str(text).replace(",", " ").split()]
        if len(values) != 12 or None in values:
            ui.notify("Give twelve pan factors, January to December", type="warning")
            return
        self._set(self.section["homogenise"], "pan_factors", values)

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
            with ui.row().classes("w-full gap-2 no-wrap"):
                with ui.element("div").classes("grow"):
                    self._path_input("Rainfall series (blank: step 1's)", a, "rainfall")
                with ui.element("div").classes("grow"):
                    self._path_input("IFD (duration_h against '1 in X' columns, mm)", a, "ifd",
                                     mark="ifd")
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

    def _failed(self, box, result) -> None:
        box.clear()
        with box:
            errors = [line for line in result.output.splitlines()
                      if line.startswith("ERROR")]
            severity_banner("block", errors[0] if errors else "The run failed - the log is "
                                                              "below.")
            with ui.expansion("Log").classes("w-full"):
                ui.code(result.output[-20000:] or "(no output)").classes("w-full text-xs")
