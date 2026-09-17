"""Lake level frequency: the recorded annual maxima against the design floods.

Point it at the dam's headwater level record - one or more Hydstra or WMIP
exports, in the order the gauges operated - or at an annual maximum series
someone has already derived, and it shows the annual maxima on a frequency
axis. **Fit curves** runs util/LakeLevelFrequency.py under Bryan's interpreter
for the curves and their resampled bands, with the Monte Carlo design floods of
one group laid over them. The figure export is formatted for an A4 report page,
with or without the design floods; the CSV export is the annual maximum series
with what produced it written above it.

Everything on the page is saved as it changes to ``lake_frequency.json`` beside
the sims_config.json, so the analysis comes back the way it was left and can be
handed to someone else with the project. See ``core/lakefreq.py``.
"""

from __future__ import annotations

from pathlib import Path

from nicegui import app, run, ui

from core import events, lakefreq
from layout import page_frame, require_project, severity_banner
from state import STATE

MONTHS = {index + 1: name for index, name in enumerate(lakefreq.RECORD.MONTH_NAMES)}


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


def lake_levels_page() -> None:
    with page_frame("Lake levels"):
        project = require_project()
        if project is None:
            return
        _LakeLevelsView(project).build()


def _number(value):
    try:
        return None if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return None


def _reference_text(levels) -> str:
    return "\n".join(f"{item.get('label', '')} = {item.get('level')}" for item in levels)


def _parse_references(text) -> list:
    out = []
    for line in str(text or "").splitlines():
        if "=" not in line:
            continue
        label, _, value = line.rpartition("=")
        level = _number(value.strip())
        if level is not None:
            out.append({"label": label.strip(), "level": level})
    return out


class _LakeLevelsView:
    def __init__(self, project) -> None:
        self.project = project
        self.config_path = project.config.config_path
        self.settings = lakefreq.load_settings(self.config_path)
        self.groups = events.sources_by_group(project)
        if self.settings["design"]["group"] not in self.groups and self.groups:
            self.settings["design"]["group"] = next(iter(self.groups))
        self.view = None
        self.results = None
        self.show_design = True
        self.whole_record = False
        self._pending = None

    # -- settings ----------------------------------------------------------

    def save(self) -> None:
        try:
            lakefreq.save_settings(self.config_path, self.settings)
        except OSError as exc:
            ui.notify(f"Could not save {lakefreq.SETTINGS_NAME}: {exc}", type="warning")

    def set(self, section, key, value, *, reread=False) -> None:
        target = self.settings if section is None else self.settings[section]
        target[key] = value
        self.save()
        if reread:
            self.refresh_record_soon()
        else:
            self.redraw()

    def plan(self) -> lakefreq.JobPlan:
        return lakefreq.build_job(self.config_path, self.settings, self.groups)

    # -- layout ------------------------------------------------------------

    def build(self) -> None:
        ui.label(f"Settings are kept in {lakefreq.settings_path(self.config_path)}"
                 ).classes("text-xs text-gray-500")
        with ui.row().classes("w-full items-start gap-4 no-wrap"):
            with ui.column().classes("w-96 shrink-0 gap-2"):
                self._inputs()
            with ui.column().classes("grow gap-2 min-w-0"):
                with ui.row().classes("gap-2 items-center"):
                    self.fit_button = ui.button("Fit curves", icon="show_chart",
                                                on_click=self.fit).mark("fit-curves")
                    ui.button("Export figure", icon="image",
                              on_click=self.export_dialog).props("outline")
                    ui.button("Export AMS CSV", icon="table_view",
                              on_click=self.export_csv).props("outline")
                    self.spinner = ui.spinner(size="md")
                    self.spinner.set_visibility(False)
                    ui.switch("Show design floods", value=self.show_design,
                              on_change=lambda e: (setattr(self, "show_design", e.value),
                                                   self.redraw()))
                    ui.switch("Whole record", value=self.whole_record,
                              on_change=lambda e: (setattr(self, "whole_record", e.value),
                                                   self.redraw()))
                self.messages = ui.column().classes("w-full gap-1")
                self.chart = ui.echart({"series": []}).classes("w-full h-[34rem]") \
                    .mark("lake-chart")
                with ui.expansion("Annual maxima", icon="list").classes("w-full"):
                    self.table_area = ui.column().classes("w-full")
                with ui.expansion("Water year options", icon="date_range"
                                  ).classes("w-full") as self.water_year_panel:
                    self._water_year_panel()
        self.refresh_record_soon()

    def _inputs(self) -> None:
        settings = self.settings
        with ui.card().classes("w-full").props("flat bordered"):
            ui.label("Level record").classes("font-bold")
            ui.textarea(
                "Hydstra / WMIP exports, one per line, earliest gauge first",
                value="\n".join(settings["record"]["files"]),
                on_change=lambda e: self.set(
                    "record", "files",
                    [line.strip() for line in (e.value or "").splitlines() if line.strip()],
                    reread=True),
            ).classes("w-full").props("autogrow dense")
            ui.input("…or an annual maximum CSV", value=settings["record"]["ams_csv"],
                     on_change=lambda e: self.set("record", "ams_csv", e.value or "",
                                                  reread=True)
                     ).classes("w-full").props("dense")
            ui.label("Relative paths are read from the sims_config.json folder."
                     ).classes("text-xs text-gray-500")

        with ui.card().classes("w-full").props("flat bordered"):
            ui.label("Annual maxima").classes("font-bold")
            ui.select(MONTHS, value=int(settings["water_year_start"]),
                      label="Water year starts",
                      on_change=lambda e: self.set(None, "water_year_start", int(e.value),
                                                   reread=True)
                      ).classes("w-full").props("dense")
            with ui.row().classes("w-full no-wrap gap-2"):
                ui.number("Carry-over window (days)", value=settings["carryover_days"],
                          min=0, step=0.5, format="%g",
                          on_change=lambda e: self.set(None, "carryover_days",
                                                       _number(e.value) or 0.0, reread=True)
                          ).classes("grow").props("dense")
                ui.number("Minimum coverage", value=settings["min_coverage"], min=0,
                          max=1, step=0.05, format="%g",
                          on_change=lambda e: self.set(None, "min_coverage",
                                                       _number(e.value) or 0.0, reread=True)
                          ).classes("grow").props("dense")
            ui.checkbox("Include water years not fully covered",
                        value=bool(settings["include_incomplete"]),
                        on_change=lambda e: self.set(None, "include_incomplete", e.value,
                                                     reread=True)).props("dense")

        with ui.card().classes("w-full").props("flat bordered"):
            ui.label("Levels").classes("font-bold")
            with ui.row().classes("w-full no-wrap gap-2"):
                ui.number("Full supply level (m AHD)", value=settings["fsl"], format="%g",
                          on_change=lambda e: self.set(None, "fsl", _number(e.value))
                          ).classes("grow").props("dense")
                ui.input("Label", value=settings["fsl_label"],
                         on_change=lambda e: self.set(None, "fsl_label", e.value or "FSL")
                         ).classes("w-24").props("dense")
            ui.textarea("Other levels to mark, 'label = level' per line",
                        value=_reference_text(settings["reference_levels"]),
                        on_change=lambda e: self.set(None, "reference_levels",
                                                     _parse_references(e.value))
                        ).classes("w-full").props("autogrow dense")

        fit = settings["fit"]
        with ui.card().classes("w-full").props("flat bordered"):
            ui.label("Curve").classes("font-bold")
            ui.select(lakefreq.FORM_LABELS, value=fit["form"], label="Form",
                      on_change=lambda e: self.set("fit", "form", e.value)
                      ).classes("w-full").props("dense")
            with ui.row().classes("w-full no-wrap gap-2"):
                ui.number("Plateau tolerance (mm, 0 = none)",
                          value=fit["plateau_tolerance"] * 1000,
                          min=0, format="%g",
                          on_change=lambda e: self.set("fit", "plateau_tolerance",
                                                       (_number(e.value) or 0) / 1000)
                          ).classes("grow").props("dense")
                ui.number("Plateau gap (mm)", value=fit["plateau_gap"] * 1000, min=0,
                          format="%g",
                          on_change=lambda e: self.set("fit", "plateau_gap",
                                                       (_number(e.value) or 0) / 1000)
                          ).classes("grow").props("dense")
            with ui.row().classes("w-full no-wrap gap-2"):
                ui.number("Degree below FSL", value=fit["degree"], min=1, max=6, step=1,
                          format="%d",
                          on_change=lambda e: self.set("fit", "degree",
                                                       int(_number(e.value) or 4))
                          ).classes("grow").props("dense")
                ui.number("Degree above FSL", value=fit["upper_degree"], min=1, max=4,
                          step=1, format="%d",
                          on_change=lambda e: self.set("fit", "upper_degree",
                                                       int(_number(e.value) or 1))
                          ).classes("grow").props("dense")
            ui.select(lakefreq.UPPER_JOIN_LABELS, value=fit["upper_join"],
                      label="Curve above FSL starts",
                      on_change=lambda e: self.set("fit", "upper_join", e.value)
                      ).classes("w-full").props("dense")
            ui.label(lakefreq.CURVE_HELP).classes("text-xs text-gray-500")
            ui.number("Resamples", value=fit["draws"], min=50, step=50, format="%d",
                      on_change=lambda e: self.set("fit", "draws",
                                                   int(_number(e.value) or 400))
                      ).classes("w-full").props("dense")
            ui.checkbox("Fit the storm-driven maxima as well",
                        value=bool(fit["storm_driven"]),
                        on_change=lambda e: self.set("fit", "storm_driven", e.value)
                        ).props("dense")

        design = settings["design"]
        with ui.card().classes("w-full").props("flat bordered"):
            ui.label("Design floods").classes("font-bold")
            if not self.groups:
                ui.label("No Monte Carlo database found in this project's sims list."
                         ).classes("text-xs text-gray-500")
            else:
                ui.checkbox("Compare with the Monte Carlo results",
                            value=bool(design["include"]),
                            on_change=lambda e: self.set("design", "include", e.value)
                            ).props("dense")
                self.duration_select = None

                def durations_for(group):
                    return {choice.duration: choice.label
                            for choice in lakefreq.design_choices(self.groups, group)}

                def on_group(event) -> None:
                    self.settings["design"]["durations"] = None
                    self.set("design", "group", event.value)
                    options = durations_for(event.value)
                    self.duration_select.set_options(options, value=list(options))

                ui.select(list(self.groups), value=design["group"], label="Group",
                          on_change=on_group).classes("w-full").props("dense")
                options = durations_for(design["group"])
                chosen = list(options) if design["durations"] is None else [
                    value for value in options if value in set(design["durations"])]
                self.duration_select = ui.select(
                    options, value=chosen, multiple=True, label="Durations",
                    on_change=lambda e: self.set("design", "durations",
                                                 sorted(float(v) for v in e.value))
                ).classes("w-full").props("dense use-chips")
            ui.number("Axis reaches (1 in X)", value=settings["axes"]["rare_aep_1_in_x"],
                      min=10, format="%g",
                      on_change=lambda e: self.set("axes", "rare_aep_1_in_x",
                                                   _number(e.value) or 2000)
                      ).classes("w-full").props("dense")
            with ui.row().classes("w-full no-wrap gap-2"):
                ui.number("Level axis from", value=settings["axes"]["level_min"],
                          format="%g",
                          on_change=lambda e: self.set("axes", "level_min", _number(e.value))
                          ).classes("grow").props("dense clearable")
                ui.number("to", value=settings["axes"]["level_max"], format="%g",
                          on_change=lambda e: self.set("axes", "level_max", _number(e.value))
                          ).classes("grow").props("dense clearable")

    # -- the record --------------------------------------------------------

    def refresh_record_soon(self) -> None:
        if lakefreq.is_quick(self.plan().job):
            # Drawn with the page rather than after it, when nothing slow is needed.
            self.view = None
            plan = self.plan()
            if plan.can_read:
                try:
                    self.view = lakefreq.read_view(plan.job)
                except (OSError, ValueError) as exc:
                    ui.notify(f"Could not read the record: {exc}", type="negative")
            self.redraw()
            return
        # A path is typed a character at a time; read the export once typing stops.
        if self._pending is not None:
            self._pending.cancel()
        self._pending = ui.timer(0.6, self.refresh_record, once=True)

    async def refresh_record(self) -> None:
        self._pending = None
        plan = self.plan()
        self.view = None
        if not plan.can_read:
            self.redraw()
            return
        self.spinner.set_visibility(True)
        try:
            if lakefreq.is_quick(plan.job):
                self.view = lakefreq.read_view(plan.job)
            else:
                self.view = await _off_thread(lakefreq.read_view, plan.job)
        except (OSError, ValueError) as exc:
            self.view = None
            ui.notify(f"Could not read the record: {exc}", type="negative")
        finally:
            self.spinner.set_visibility(False)
        self.redraw()

    def redraw(self) -> None:
        plan = self.plan()
        self.messages.clear()
        with self.messages:
            for problem in plan.problems:
                severity_banner("warn", problem)
            if self.view is not None:
                for note in self.view.notes:
                    severity_banner("info", note)
        if self.view is None:
            self.chart.options.clear()
            self.chart.options.update({"series": []})
            self.chart.update()
            self.table_area.clear()
            return

        self.results = lakefreq.cached_results(self.config_path, plan.job) \
            if plan.can_fit else None
        with self.messages:
            for problem in plan.fit_problems:
                severity_banner("warn", problem)
            if plan.can_fit and self.results is None \
                    and plan.job["fit"]["form"] != "none":
                ui.label("The curves have not been fitted for these settings - "
                         "press Fit curves.").classes("text-sm text-gray-600")
            for name, block in ((self.results or {}).get("fits") or {}).items():
                for problem in (block.get("error"), block.get("warning")):
                    if problem:
                        severity_banner(
                            "warn", f"{'Storm-driven' if name == 'storm' else 'All'}"
                                    f" maxima: {problem}")
        options = lakefreq.chart_options(self.view.positions, self.results, self.settings,
                                         show_design=self.show_design,
                                         whole_record=self.whole_record)
        self.chart.options.clear()
        self.chart.options.update(options)
        self.chart.update()
        self._table()

    def _table(self) -> None:
        self.table_area.clear()
        frame = self.view.positions.sort_values("water_year")
        columns = [("period", "Water year"), ("level", "Level (m AHD)"),
                   ("level_at", "When"), ("days_into_year", "Day of year"),
                   ("carried_over", "Carried over"), ("coverage", "Coverage"),
                   ("reading_interval_min", "Reading interval (min)"),
                   ("aep_1_in_x", "1 in X")]
        rows = []
        for record in frame.to_dict("records"):
            rows.append({
                "period": record["period"], "level": f"{record['level']:.3f}",
                "level_at": "" if record["level_at"] != record["level_at"]
                else str(record["level_at"])[:16],
                "days_into_year": "" if record["days_into_year"] != record["days_into_year"]
                else f"{record['days_into_year']:.1f}",
                "carried_over": "yes" if record["carried_over"] else "",
                "coverage": f"{record['coverage']:.0%}",
                "reading_interval_min": "" if record["reading_interval_min"]
                != record["reading_interval_min"] else f"{record['reading_interval_min']:.0f}",
                "aep_1_in_x": f"{record['aep_1_in_x']:.1f}",
            })
        with self.table_area:
            ui.table(columns=[{"name": key, "label": label, "field": key, "align": "left"}
                              for key, label in columns],
                     rows=rows, row_key="period").classes("w-full").props("dense flat")

    # -- fitting -------------------------------------------------------------

    async def fit(self) -> None:
        plan = self.plan()
        problem = lakefreq.interpreter_problem(STATE.settings.bryan_python)
        blockers = plan.problems + plan.fit_problems + ([problem] if problem else [])
        if blockers:
            ui.notify(blockers[0], type="warning")
            return
        job_path = lakefreq.write_job(self.config_path, plan.job)
        results = lakefreq.results_path(self.config_path, plan.job)
        argv = lakefreq.command(STATE.settings.bryan_python, job_path, results=results)
        self.fit_button.set_enabled(False)
        self.spinner.set_visibility(True)
        ui.notify("Fitting - the resampling takes of the order of a minute")
        try:
            done = await _off_thread(lakefreq.run, argv, (results,))
        finally:
            self.fit_button.set_enabled(True)
            self.spinner.set_visibility(False)
        if not done.ok:
            with self.messages:
                severity_banner("block", f"The fit failed (exit {done.returncode}).",
                                done.output[-1500:])
            return
        self.redraw()

    # -- exports -------------------------------------------------------------

    def export_dialog(self) -> None:
        folder, name = lakefreq.default_export(self.config_path, self.settings)
        with ui.dialog() as dialog, ui.card().classes("min-w-[36rem]"):
            ui.label("Export the report figure").classes("text-lg font-bold")
            ui.label("Runs util/LakeLevelFrequency.py with Bryan's own interpreter. "
                     "6.3 x 4.0 in at 300 dpi, for an A4 page. Uses the fitted "
                     "curves when they are current, and fits them first when not."
                     ).classes("text-xs text-gray-500")
            folder_input = ui.input("Output folder", value=str(folder)
                                    ).classes("w-full").props("dense")
            name_input = ui.input("Base name", value=name).classes("w-full").props("dense")
            with_design = ui.checkbox("Include the design floods",
                                      value=bool(self.settings["design"]["include"]))
            preview = ui.column().classes("w-full gap-0")
            report = ui.column().classes("w-full gap-1").mark("export-report")

            def paths():
                return lakefreq.export_paths(folder_input.value.strip(),
                                             name_input.value, with_design=with_design.value)

            def refresh() -> None:
                preview.clear()
                built = paths()
                with preview:
                    for problem in built.problems:
                        severity_banner("block", problem)
                    if built.png:
                        exists = " - exists, will be overwritten" if built.png.is_file() else ""
                        ui.label(f"{built.png}{exists}").classes(
                            "text-xs " + ("text-orange-700" if exists else "text-gray-500"))

            async def export() -> None:
                built = paths()
                plan = self.plan()
                problem = lakefreq.interpreter_problem(STATE.settings.bryan_python)
                blockers = built.problems + plan.problems + plan.fit_problems + (
                    [problem] if problem else [])
                report.clear()
                if blockers:
                    with report:
                        severity_banner("block", blockers[0])
                    return
                self.settings["export"] = {"folder": folder_input.value.strip(),
                                           "name": name_input.value.strip()}
                self.save()
                job_path = lakefreq.write_job(self.config_path, plan.job)
                argv = lakefreq.command(
                    STATE.settings.bryan_python, job_path,
                    results=lakefreq.results_path(self.config_path, plan.job),
                    png=built.png, without_design=not with_design.value)
                run_button.set_enabled(False)
                with report:
                    ui.label("Exporting...").classes("text-xs text-gray-500")
                done = await _off_thread(lakefreq.run, argv, (built.png,))
                run_button.set_enabled(True)
                report.clear()
                with report:
                    if done.ok and done.written:
                        severity_banner("info", f"Wrote {built.png}")
                    else:
                        severity_banner("block", f"The export failed (exit {done.returncode}).",
                                        done.output[-1500:])
                self.redraw()

            with ui.row().classes("justify-end gap-2 w-full"):
                ui.button("Close", on_click=dialog.close).props("flat")
                run_button = ui.button("Export", on_click=export).mark("export-run")
            for element in (folder_input, name_input):
                element.on("blur", refresh)
            with_design.on_value_change(refresh)
            refresh()
        dialog.open()

    def export_csv(self) -> None:
        if self.view is None:
            ui.notify("Read a record first.", type="warning")
            return
        folder, name = lakefreq.default_export(self.config_path, self.settings)
        built = lakefreq.export_paths(folder, name, with_design=False)
        if built.problems:
            ui.notify(built.problems[0], type="warning")
            return
        try:
            path = lakefreq.write_ams_csv(self.view, self.plan().job, built.csv)
        except OSError as exc:
            ui.notify(f"Could not write the CSV: {exc}", type="negative")
            return
        ui.notify(f"Wrote {path}")

    # -- water year options --------------------------------------------------

    def _water_year_panel(self) -> None:
        ui.label(lakefreq.SCORE_HELP).classes("text-xs text-gray-500")
        with ui.row().classes("items-end gap-2"):
            fsl = self.settings["fsl"]
            self.high_level = ui.number(
                "Count carried-over maxima above (m AHD)",
                value=None if fsl is None else round(fsl - 1.5, 2), format="%g"
            ).classes("w-64").props("dense")
            self.rise = ui.number("A rise at the cut is more than (m)", value=0.10,
                                  min=0, step=0.05, format="%g").classes("w-56").props("dense")
            ui.button("Score the start months", on_click=self.score).props("outline") \
                .mark("score-water-years")
        self.score_area = ui.column().classes("w-full gap-2")

    async def score(self) -> None:
        if self.view is None or self.view.level is None:
            ui.notify("Scoring needs the level record itself, not an annual maximum CSV.",
                      type="warning")
            return
        job = self.plan().job
        record = lakefreq.RECORD

        def compute():
            scores = record.water_year_scores(
                self.view.level, carryover_days=float(job["carryover_days"]),
                min_coverage=float(job["min_coverage"]),
                high_level=_number(self.high_level.value), full_supply=job["fsl"],
                rise_m=_number(self.rise.value) or 0.10,
                include_incomplete=bool(job["include_incomplete"]))
            return (scores, record.monthly_levels(self.view.level),
                    record.maxima_by_month(self.view.ams))

        self.spinner.set_visibility(True)
        try:
            scores, monthly, counts = await _off_thread(compute)
        finally:
            self.spinner.set_visibility(False)
        rows = lakefreq.score_rows(scores, int(job["water_year_start"]))
        self.score_area.clear()
        with self.score_area:
            table = ui.table(
                columns=[{"name": key, "label": label, "field": key, "align": "left"}
                         for key, label in lakefreq.SCORE_LABELS.items()],
                rows=rows, row_key="start_month").classes("w-full").props("dense flat") \
                .mark("water-year-scores")
            table.add_slot("body", r'''
                <q-tr :props="props" :class="props.row.adopted ? 'bg-blue-1 text-weight-bold' : ''">
                  <q-td v-for="col in props.cols" :key="col.name" :props="props">
                    {{ col.value === null ? '-' : col.value }}
                  </q-td>
                </q-tr>''')
            ui.label("Highlighted: the start month in use. The bars count the storm-"
                     "driven maxima under that start, so they describe the choice "
                     "rather than justify it; the level by month does not depend on it."
                     ).classes("text-xs text-gray-500")
            ui.echart(lakefreq.seasonal_chart(monthly, counts)).classes("w-full h-72")
