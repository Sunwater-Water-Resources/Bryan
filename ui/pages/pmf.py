"""PMF: the ensemble results of the probable maximum flood, and a notional AEP.

Each **PMF group** pairs an ensemble group (the PMF run of one climate horizon)
with the Monte Carlo group of the same horizon. The page shows the spread of the
temporal patterns by duration, picks the PMF - the highest event, as the report
takes it, with the median-pattern answer beside it - and places the PMF level on
the upper tail of the Monte Carlo realisations to give it a notional AEP, with
the fit laid out and a grid of how the answer moves with the fit's settings.

The adopted notional AEP is the analyst's call, typed in at the top and kept in
the study file with a note of how it was reached. Everything on the page is
saved to the study file (``bryan_study.json``, under ``"pmf"``) as it changes;
the study is opened on the Report page. See ``core/ensemble.py``.
"""

from __future__ import annotations

import math

from nicegui import app, run, ui

from core import ensemble, pmfchart, reporttables, staleness, study as studies
from core.results import format_aep
from layout import confirm, no_study, page_frame, severity_banner
from state import STATE
from theme import house_echart

CAVEAT = ("The grid is the sensitivity of the fit only. How far the Monte Carlo "
          "realisations reach past the AEP of the PMP - and so where the PMF sits in "
          "them - is set by how the storm config extrapolates the rainfall (a GEV at "
          "Callide) and by the rating and storage curves, which no fit setting touches.")


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


def pmf_page() -> None:
    with page_frame("PMF"):
        _PmfView().build()


def _aep_text(aep) -> str:
    return f"1 in {aep:,.0f}" if aep and math.isfinite(aep) else "-"


def _grouped(value) -> str:
    """A number as the inputs show it, '500,000'; blank when there is none."""
    number = _number(value)
    return f"{number:,.0f}" if number else ""


def _number(value):
    try:
        return None if value in (None, "") else float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None


class _PmfView:
    def __init__(self) -> None:
        self.study = STATE.study
        self.section = ensemble.settings(self.study) if self.study else None
        self.selected = 0
        self.result = "level"
        self.body = None

    def build(self) -> None:
        if self.study is None:
            no_study("The PMF page keeps its settings in the study file.")
            return
        self._adopted_card()
        self.body = ui.column().classes("w-full gap-4")
        self.redraw()

    def save(self) -> None:
        ensemble.store(self.study, self.section)
        try:
            self.study.save()
        except OSError as exc:
            ui.notify(f"Could not save {self.study.path.name}: {exc}", type="warning")

    # -- the study-wide settings -------------------------------------------

    def _adopted_card(self) -> None:
        with ui.card().classes("w-full"):
            ui.label("Notional AEP of the PMF").classes("text-lg font-bold")
            with ui.row().classes("w-full items-start gap-4 no-wrap"):
                ui.input("Adopted (1 in x)", value=_grouped(self.section["adopted_aep"])) \
                    .classes("w-48").props("dense").mark("adopted-aep") \
                    .on("blur", lambda e: self._set_adopted(e.sender.value))
                ui.input("How it was reached", value=self.section["adopted_note"]) \
                    .classes("grow").props("dense") \
                    .on("blur", lambda e: self._set("adopted_note", e.sender.value))
                ui.input("AEP of the PMP (1 in x)", value=_grouped(self.section["pmp_aep"])) \
                    .classes("w-48").props("dense") \
                    .on("blur", lambda e: self._set("pmp_aep",
                                                    _number(e.sender.value) or
                                                    ensemble.DEFAULT_PMP_AEP, redraw=True))
            ui.textarea("Reference levels on the box plot, one 'label = level' per line",
                        value=reporttables.levels_text(self.section["reference_levels"])) \
                .classes("w-full").props("dense autogrow") \
                .on("blur", lambda e: self._set(
                    "reference_levels", reporttables.parse_levels(e.sender.value),
                    redraw=True))

    def _set_adopted(self, text) -> None:
        self._set("adopted_aep", _number(text))

    def _set(self, key, value, *, redraw=False) -> None:
        if self.section.get(key) == value:
            return
        self.section[key] = value
        self.save()
        if redraw:
            self.redraw()

    # -- the groups ----------------------------------------------------------

    def redraw(self) -> None:
        self.body.clear()
        with self.body:
            groups = self.section["groups"]
            with ui.row().classes("w-full items-center gap-2"):
                if groups:
                    options = {i: entry.get("label") or f"group {i + 1}"
                               for i, entry in enumerate(groups)}
                    self.selected = min(self.selected, len(groups) - 1)
                    ui.toggle(options, value=self.selected,
                              on_change=lambda e: self._select(e.value)) \
                        .props("no-caps").mark("pmf-groups")
                ui.button("Add PMF group", icon="add", on_click=self._add_dialog) \
                    .props("outline dense no-caps").mark("add-pmf-group")
                if groups:
                    ui.button("Summary of all", icon="table_view",
                              on_click=self._summary).props("flat dense no-caps") \
                        .mark("pmf-summary")
            if not groups:
                ui.label("Add a PMF group: the ensemble group of one climate horizon's "
                         "PMF run, and the Monte Carlo group of the same horizon.") \
                    .classes("text-sm text-muted")
                return
            self.summary_box = ui.column().classes("w-full")
            _GroupView(self, groups[self.selected]).build()

    def _select(self, index) -> None:
        self.selected = index
        self.redraw()

    def _add_dialog(self) -> None:
        runs = self.study.run_names()
        entry = {"label": "", "ensemble": {"run": runs[0] if runs else "", "group": ""},
                 "mc": {"run": runs[0] if runs else "", "group": ""}}
        with ui.dialog() as dialog, ui.card().classes("min-w-[44rem]"):
            ui.label("Add a PMF group").classes("text-lg font-bold")
            if not runs:
                ui.label("Add the runs on the Report page first - the PMF ensemble run "
                         "and the Monte Carlo run.").classes("text-muted")
            label = ui.input("Label, e.g. 'Near-term RFSL'").classes("w-full").props("dense")
            box = ui.column().classes("w-full gap-1")

            def draw() -> None:
                box.clear()
                with box:
                    ui.label("PMF ensemble group").classes("text-sm text-muted")
                    _source_picker(self.study, entry["ensemble"], draw, mark="ensemble")
                    ui.label("Monte Carlo group of the same horizon - its realisations "
                             "give the AEP").classes("text-sm text-muted")
                    _source_picker(self.study, entry["mc"], draw, mark="mc")

            draw()

            def add() -> None:
                if not entry["ensemble"]["group"]:
                    ui.notify("Choose the PMF ensemble group", type="warning")
                    return
                entry["label"] = label.value or entry["ensemble"]["group"]
                self.section["groups"].append(entry)
                self.selected = len(self.section["groups"]) - 1
                self.save()
                dialog.close()
                self.section = ensemble.settings(self.study)
                self.redraw()

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=dialog.close).props("flat")
                ui.button("Add", on_click=add).mark("confirm-add-pmf")
        dialog.open()

    def remove(self, entry) -> None:
        self.section["groups"] = [e for e in self.section["groups"] if e is not entry]
        self.selected = 0
        self.save()
        self.redraw()

    async def _summary(self) -> None:
        self.summary_box.clear()
        with self.summary_box:
            ui.spinner()
        rows = await _off_thread(_summarise_all, self.study, self.section)
        self.summary_box.clear()
        with self.summary_box:
            columns = [{"name": key, "label": label, "field": key,
                        "align": "left" if key == "label" else "right"}
                       for key, label in (("label", "PMF group"), ("level", "PMF level"),
                                          ("duration", "Duration (h)"),
                                          ("pattern", "Pattern"),
                                          ("median", "Median-pattern level"),
                                          ("aep", "Notional AEP (this fit)"),
                                          ("range", "Range over the grid"))]
            ui.table(columns=columns, rows=rows, row_key="label") \
                .props("dense flat bordered").classes("w-full").mark("pmf-summary-table")


def _source_picker(study, holder, redraw, *, mark) -> None:
    runs = study.run_names()

    def on_run(event) -> None:
        holder["run"] = event.value
        holder["group"] = ""
        redraw()

    groups = []
    if holder.get("run"):
        try:
            groups = study.open_run(holder["run"]).groups()
        except studies.StudyError:
            groups = []
    with ui.row().classes("w-full gap-2 no-wrap"):
        ui.select(runs, value=holder.get("run") if holder.get("run") in runs else None,
                  label="Run", on_change=on_run).classes("w-48").props("dense") \
            .mark(f"run-{mark}")
        ui.select(groups, value=holder.get("group") if holder.get("group") in groups else None,
                  label="Group", with_input=True,
                  on_change=lambda e: holder.update(group=e.value)) \
            .classes("grow").props("dense").mark(f"group-{mark}")


def _duration_for(study, entry, highest):
    """The Monte Carlo duration the fit reads: chosen, or the PMF event's own."""
    databases = ensemble.mc_databases(study, entry["mc"]["run"], entry["mc"]["group"])
    chosen = entry.get("duration") or ""
    if chosen and chosen in databases:
        return chosen, databases
    if highest is not None and highest.found:
        return ensemble.nearest_duration(list(databases), highest.duration), databases
    return (next(iter(databases)) if databases else None), databases


def _estimate(study, entry, result="level"):
    """Everything the AEP card draws, for one PMF group. Runs off the event loop."""
    frame = ensemble.load(study, entry["ensemble"]["run"], entry["ensemble"]["group"])
    highest = ensemble.pick(frame, ensemble.HIGHEST)
    duration, databases = _duration_for(study, entry, highest)
    if duration is None:
        raise studies.StudyError("the Monte Carlo group has no results database on disk")
    real = ensemble.read_realisations(databases[duration], result)
    target = getattr(highest, result)
    estimate = ensemble.fit(real, target, entry.get("lower_aep") or ensemble.DEFAULT_LOWER_AEP,
                            entry.get("upper_aep"), int(entry.get("degree") or 1))
    grid = ensemble.sensitivity(real, target,
                                entry.get("lower_aep") or ensemble.DEFAULT_LOWER_AEP,
                                entry.get("upper_aep"))
    return {"frame": frame, "highest": highest, "duration": duration,
            "durations": list(databases), "real": real, "estimate": estimate, "grid": grid}


def _summarise_all(study, section) -> list:
    rows = []
    for entry in section["groups"]:
        row = {"label": entry.get("label", "")}
        try:
            out = _estimate(study, entry)
        except studies.StudyError as exc:
            row.update(level=str(exc))
            rows.append(row)
            continue
        highest, estimate, grid = out["highest"], out["estimate"], out["grid"]
        median = ensemble.pick(out["frame"], ensemble.MEDIAN)
        values = [v for v in grid.to_numpy().ravel() if math.isfinite(v)]
        row.update(level=f"{highest.level:.2f}", duration=f"{highest.duration:g}",
                   pattern=highest.pattern, median=f"{median.level:.2f}",
                   aep=_aep_text(estimate.aep) if estimate.ok else estimate.refused,
                   range=(f"{min(values) / 1e6:.1f}M - {max(values) / 1e6:.1f}M"
                          if values else "-"))
        rows.append(row)
    return rows


class _GroupView:
    def __init__(self, page: _PmfView, entry: dict) -> None:
        self.page = page
        self.study = page.study
        self.entry = entry

    def build(self) -> None:
        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between"):
                ui.label(self.entry.get("label", "")).classes("text-lg font-bold")
                with ui.row().classes("gap-2 items-center"):
                    ui.toggle({"level": "Level", "inflow": "Inflow", "outflow": "Outflow"},
                              value=self.page.result, on_change=self._on_result) \
                        .props("dense no-caps").mark("pmf-result")
                    ui.button(icon="delete", on_click=self._ask_remove) \
                        .props("flat dense round").tooltip("Remove this PMF group") \
                        .mark("remove-pmf-group")
            ui.label(f"{self.entry['ensemble']['group']} in {self.entry['ensemble']['run']}; "
                     f"realisations from {self.entry['mc']['group'] or '(none)'} in "
                     f"{self.entry['mc']['run'] or '(none)'}").classes("text-xs text-muted")
            self.stale_box = ui.column().classes("w-full gap-2")
            self.box = ui.column().classes("w-full gap-2")
            with self.box:
                ui.spinner()
        self.fit_card = ui.card().classes("w-full")
        ui.timer(0.01, self.fill, once=True)

    def _ask_remove(self) -> None:
        confirm(f"Remove the PMF group '{self.entry.get('label') or '(unlabelled)'}'?",
                "Its runs and results are not touched; the adopted notional AEP is kept.",
                lambda: self.page.remove(self.entry))

    def _on_result(self, event) -> None:
        self.page.result = event.value
        self.page.redraw()

    async def fill(self) -> None:
        result = self.page.result
        try:
            out = await _off_thread(_estimate, self.study, self.entry, result)
        except studies.StudyError as exc:
            self.box.clear()
            with self.box:
                severity_banner("warn", str(exc))
            self.fit_card.clear()
            return
        if self.box.is_deleted:
            return
        self._draw_ensemble(out, result)
        self._draw_fit(out)
        stale = await _off_thread(staleness.for_sources, self.study,
                                  [self.entry["ensemble"], self.entry["mc"]]) or []
        if stale and not self.stale_box.is_deleted:
            with self.stale_box:
                severity_banner("warn", "\n".join(stale),
                                "Re-run the group before quoting this PMF.")

    def _draw_ensemble(self, out, result) -> None:
        frame, highest = out["frame"], out["highest"]
        median = ensemble.pick(frame, ensemble.MEDIAN)
        stats = ensemble.by_duration(frame, result)
        self.box.clear()
        with self.box:
            with ui.row().classes("w-full gap-8"):
                with ui.column().classes("gap-0"):
                    ui.label("PMF - highest event").classes("text-xs text-muted")
                    ui.label(f"{highest.level:.2f} m AHD, {highest.duration:g} h, "
                             f"{highest.pattern}").classes("text-lg font-bold") \
                        .mark("pmf-highest")
                    ui.label(f"inflow {highest.inflow:,.0f} m³/s, outflow "
                             f"{highest.outflow:,.0f} m³/s").classes("text-sm")
                with ui.column().classes("gap-0"):
                    ui.label("Median pattern (the convention for other ensemble results)"
                             ).classes("text-xs text-muted")
                    ui.label(f"{median.level:.2f} m AHD, {median.duration:g} h, "
                             f"{median.pattern}").classes("text-lg")
                    ui.label(f"{highest.level - median.level:+.2f} m to the highest event"
                             ).classes("text-sm text-muted")
            house_echart(pmfchart.box_chart(
                stats, ensemble.pattern_points(frame, result), result, highest=highest,
                reference_levels=self.page.section["reference_levels"])) \
                .classes("w-full h-80").mark("pmf-box")
            columns = [{"name": k, "label": l, "field": k, "align": "right"}
                       for k, l in (("duration", "Duration (h)"), ("n", "Patterns"),
                                    ("median", "Median"), ("median_pattern", "Median pattern"),
                                    ("max", "Highest"), ("max_pattern", "Highest pattern"))]
            decimals = 2 if result == "level" else 0
            rows = [{"duration": f"{hours:g}", "n": int(row["n"]),
                     "median": f"{row['median']:,.{decimals}f}",
                     "median_pattern": row["median_pattern"],
                     "max": f"{row['max']:,.{decimals}f}", "max_pattern": row["max_pattern"]}
                    for hours, row in stats.iterrows()]
            ui.table(columns=columns, rows=rows, row_key="duration") \
                .props("dense flat bordered").classes("w-full")

    def _draw_fit(self, out) -> None:
        entry, estimate, grid = self.entry, out["estimate"], out["grid"]
        self.fit_card.clear()
        with self.fit_card:
            ui.label("Notional AEP from the Monte Carlo realisations").classes(
                "text-lg font-bold")
            with ui.row().classes("w-full items-end gap-3 no-wrap"):
                ui.select({"": f"PMF's own ({out['duration']})",
                           **{label: label for label in out["durations"]}},
                          value=entry.get("duration") or "", label="Duration",
                          on_change=lambda e: self._set("duration", e.value)) \
                    .classes("w-48").props("dense").mark("pmf-duration")
                ui.input("Window from 1 in", value=_grouped(entry.get("lower_aep"))) \
                    .classes("w-40").props("dense").mark("pmf-lower") \
                    .on("blur", lambda e: self._set("lower_aep", _number(e.sender.value)
                                                    or ensemble.DEFAULT_LOWER_AEP))
                ui.input("to 1 in (blank: top of the sample)",
                         value=_grouped(entry.get("upper_aep"))) \
                    .classes("w-56").props("dense") \
                    .on("blur", lambda e: self._set("upper_aep", _number(e.sender.value)))
                ui.toggle({1: "Degree 1", 2: "2", 3: "3"},
                          value=int(entry.get("degree") or 1),
                          on_change=lambda e: self._set("degree", e.value)) \
                    .props("dense no-caps").mark("pmf-degree")
            real = out["real"]
            ui.label(f"{out['duration']} realisations reach 1 in {real.top_aep:,.0f} at "
                     f"{real.top_value:,.2f}; {estimate.n} in the window.") \
                .classes("text-sm text-body")
            if estimate.ok:
                ui.label(f"The PMF, {estimate.target:,.2f}, is at "
                         f"{_aep_text(estimate.aep)}").classes("text-xl font-bold") \
                    .mark("pmf-aep")
            else:
                severity_banner("block", f"No estimate: {estimate.refused}")
            for warning in estimate.warnings:
                severity_banner("warn", warning)
            house_echart(pmfchart.fit_chart(real, estimate,
                                            pmp_aep=self.page.section["pmp_aep"])) \
                .classes("w-full h-96").mark("pmf-fit")
            ui.label("How the answer moves with the fit (1 in x, millions)").classes(
                "text-sm font-bold")
            columns = [{"name": "window", "label": "Window from 1 in", "field": "window",
                        "align": "left"}] + [
                {"name": str(d), "label": f"Degree {d}", "field": str(d), "align": "right"}
                for d in grid.columns]
            rows = [{"window": format_aep(lower),
                     **{str(d): (f"{grid.loc[lower, d] / 1e6:.1f}"
                                 if math.isfinite(grid.loc[lower, d]) else "-")
                        for d in grid.columns}}
                    for lower in grid.index]
            ui.table(columns=columns, rows=rows, row_key="window") \
                .props("dense flat bordered").mark("pmf-grid")
            ui.label(CAVEAT).classes("text-xs text-muted")

    def _set(self, key, value) -> None:
        if self.entry.get(key) == value:
            return
        self.entry[key] = value
        self.page.save()
        self.page.redraw()
