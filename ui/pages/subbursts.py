"""The Results page's Sub-bursts tab: the neutrality check, and calibrated pattern weights.

``core/subburst.py`` holds the logic and says what the check is; this is the view. The
tab is drawn when first shown, because a chart built inside a hidden tab panel measures
zero and stays blank - the same reason the Groups tab waits (see pages/results.py) - and
because finding which groups recorded sub-bursts reads the header of every database,
which should not slow the Results page for someone who never opens this tab.
"""

from __future__ import annotations

from nicegui import app, run, ui

from core import subburst
from core.results import format_aep
from layout import severity_banner
from state import STATE
from theme import house_echart
from widgets import Running

NONE = "(none)"


def _durations_of(results, role) -> dict:
    return (((results or {}).get("groups") or {}).get(role) or {}).get("durations") or {}


async def _off_thread(function, *args):
    """``run.io_bound`` returns None on cancellation; run inline unless stopping."""
    result = await run.io_bound(function, *args)
    if result is None and not app.is_stopping:
        return function(*args)
    return result


class SubBurstView:
    def __init__(self, project) -> None:
        self.project = project
        self.config_path = project.config.config_path
        self.candidates = {}
        self.recorded = []
        self.tested = None
        self.compare = None
        self.low = self.high = None
        self.parent = None
        self.level = None
        self.results = None
        self.chart = None
        self._started = False

    # -- building --------------------------------------------------------------

    def build(self) -> None:
        self.root = ui.column().classes("w-full gap-4")

    def _build(self) -> None:
        self.candidates = subburst.candidates(self.project)
        self.recorded = [g for g, c in self.candidates.items() if c.recorded]
        self.tested = self.recorded[0] if self.recorded else None
        with self.root:
            self._build_controls()

    def _build_controls(self) -> None:
        if not self.recorded:
            with ui.card().classes("w-full items-center p-8").mark("subburst-empty"):
                ui.icon("grain").classes("text-5xl text-muted")
                ui.label("No group here has recorded sub-burst depths.").classes("text-muted")
                ui.label("Bryan records them in every Monte Carlo database it writes; a database "
                         "made by an older version has none. Re-running the storms records them - "
                         "'Run models: storms only' does it in minutes, without the hydrologic model."
                         ).classes("text-xs text-muted max-w-lg text-center")
            return
        with ui.card().classes("w-full"):
            ui.label("Are the embedded bursts in the design storms neutral? For each storm duration, "
                     "the wettest window of each shorter duration is taken through the total "
                     "probability theorem and divided by the IFD for that window at the same AEP. "
                     "Above 1, embedded bursts of that duration occur more often than the rainfall "
                     "statistics allow.").classes("text-sm text-body")
            names = subburst.display_names(self.candidates)
            recorded = {group: names.get(group, group) for group in self.recorded}
            with ui.row().classes("w-full items-end gap-4"):
                ui.select(recorded, value=self.tested, label="Group to test",
                          on_change=lambda e: self._choose(tested=e.value)) \
                    .classes("min-w-[18rem]").mark("subburst-tested")
                ui.select({NONE: NONE, **recorded}, value=NONE, label="Compare with",
                          on_change=lambda e: self._choose(compare=e.value)) \
                    .classes("min-w-[18rem]").mark("subburst-compare")
                self.check_button = ui.button("Check neutrality", icon="rule",
                                              on_click=lambda: self._run(False)) \
                    .props("no-caps").mark("subburst-check")
                self.calibrate_button = ui.button("Calibrate pattern weights", icon="tune",
                                                  on_click=lambda: self._run(True)) \
                    .props("no-caps outline").mark("subburst-calibrate")
            ui.label("Checking takes seconds per duration. Calibrating down-weights the patterns "
                     "whose embedded bursts occur too often until the group is neutral, and shows "
                     "what that does to the design level curve - a few minutes for a full group."
                     ).classes("text-xs text-muted")
            self.running_box = ui.row().classes("w-full")
        self.messages = ui.column().classes("w-full gap-2")
        self.body = ui.column().classes("w-full gap-4").mark("subburst-body")

    def activate(self) -> None:
        if not self._started:
            self._started = True
            self._build()
            if self.recorded:
                self._load()
                self.draw()
        if self.chart is not None:
            try:
                self.chart.run_chart_method("resize")
            except Exception:        # noqa: BLE001 - never worth failing a page over
                pass

    # -- state -----------------------------------------------------------------

    def _choose(self, tested=None, compare=None) -> None:
        if tested is not None:
            self.tested = tested
        if compare is not None:
            self.compare = None if compare == NONE else compare
        self._load()
        self.draw()

    def _load(self) -> None:
        tested = self.candidates.get(self.tested)
        compare = self.candidates.get(self.compare) if self.compare else None
        self.results = subburst.shown(self.project, tested, compare) if tested else None
        self.parent = None

    # -- running ---------------------------------------------------------------

    async def _run(self, calibrate: bool) -> None:
        tested = self.candidates.get(self.tested)
        compare = self.candidates.get(self.compare) if self.compare else None
        plan = subburst.build_job(self.project, tested, compare, calibrate=calibrate)
        problem = subburst.interpreter_problem(STATE.settings.bryan_python)
        blockers = plan.problems + ([problem] if problem else [])
        if blockers:
            ui.notify(blockers[0], type="warning")
            return
        job_path = subburst.write_job(self.config_path, plan)
        results = subburst.results_path(self.config_path, plan.fingerprint)
        argv = subburst.command(STATE.settings.bryan_python, job_path, results)
        what = ("Calibrating the pattern weights, duration by duration" if calibrate
                else "Checking the sub-burst neutrality")
        for button in (self.check_button, self.calibrate_button):
            button.set_enabled(False)
        self.messages.clear()
        running = Running(self.running_box, what, mark="running-subburst")
        try:
            done = await _off_thread(subburst.run, argv, running.cancel)
        finally:
            for button in (self.check_button, self.calibrate_button):
                button.set_enabled(True)
            running.done()
        if done.cancelled:
            ui.notify("Cancelled - the results are as they were")
            return
        saved = subburst.cached(self.config_path, plan)
        if done.returncode != 0 or saved is None:
            with self.messages:
                severity_banner("block", f"The check failed (exit {done.returncode}).",
                                done.output[-1500:])
            return
        self.results = saved
        self.draw()

    # -- drawing ---------------------------------------------------------------

    def draw(self) -> None:
        self.body.clear()
        self.chart = None
        with self.body:
            if self.results is None:
                ui.label("Not checked yet for this choice of groups - press Check neutrality."
                         ).classes("text-sm text-muted").mark("subburst-none")
                return
            problems = self.results.get("problems") or []
            if problems:
                severity_banner("warn", "\n".join(problems[:6])
                                + (f"\n... and {len(problems) - 6} more" if len(problems) > 6 else ""))
            aeps = subburst.window(self.results, self.low, self.high)
            severity, message, hint = subburst.verdict(self.results, aeps)
            severity_banner(severity, message, hint)
            self._window_controls()
            self._matrix(aeps)
            self._chart_card(aeps)
            if (self.results.get("calibration") or {}).get("durations"):
                self._calibration()
            self._files()

    def _window_controls(self) -> None:
        options = {a: f"1 in {format_aep(a)}" for a in subburst.standard_aeps(self.results)}
        current = subburst.window(self.results, self.low, self.high)
        with ui.row().classes("items-end gap-4"):
            ui.select(options, value=current[0] if current else None, label="Worst margin from",
                      on_change=lambda e: self._set_window(low=e.value)).classes("w-40")
            ui.select(options, value=current[-1] if current else None, label="to",
                      on_change=lambda e: self._set_window(high=e.value)).classes("w-40")

    def _set_window(self, low=None, high=None) -> None:
        if low is not None:
            self.low = low
            self.high = self.high if self.high is not None else subburst.window(self.results)[-1]
        if high is not None:
            self.high = high
            self.low = self.low if self.low is not None else subburst.window(self.results)[0]
        self.draw()

    def _matrix(self, aeps) -> None:
        subs, rows = subburst.matrix_rows(self.results, aeps)
        compared = bool(self.compare) and bool(self.results["groups"].get(subburst.COMPARE))
        with ui.card().classes("w-full").mark("subburst-matrix"):
            ui.label("Worst margin in the AEP range, by storm duration (rows) and window duration "
                     "(columns)" + (" - tested → compared" if compared else "")
                     ).classes("text-sm font-medium")
            with ui.element("div").classes("overflow-x-auto w-full"):
                with ui.grid(columns=len(subs) + 1).classes("gap-x-4 gap-y-1 items-center"):
                    ui.label("Storm").classes("text-xs text-muted")
                    for sub in subs:
                        ui.label(sub).classes("text-xs text-muted")
                    for row in rows:
                        ui.label(row["parent"]).classes("text-sm font-medium")
                        for sub in subs:
                            cell = ui.label(row.get(sub, "")).classes("text-sm tabular-nums")
                            if row.get(f"{sub}_breach"):
                                cell.classes("text-attention font-medium")
            ui.label(f"Above {subburst.BREACH:.2f} (in the attention colour) the patterns produce "
                     "embedded bursts of that duration more often than the IFD allows.").classes(
                "text-xs text-muted")

    def _chart_card(self, aeps) -> None:
        parents = subburst.parents(self.results)
        if not parents:
            return
        worst = subburst.worst_overall(self.results, subburst.TESTED, aeps)
        if self.parent not in parents:
            self.parent = worst[1] if worst else parents[-1]
        with ui.card().classes("w-full"):
            with ui.row().classes("items-end gap-4"):
                ui.select(parents, value=self.parent, label="Storm duration",
                          on_change=lambda e: self._set_parent(e.value)).classes("w-40") \
                    .mark("subburst-parent")
                compared = bool(_durations_of(self.results, subburst.COMPARE))
                ui.label(("Solid: the tested group. Dashed: the compared group. " if compared else "")
                         + "The line at 1 is the IFD.").classes("text-xs text-muted")
            self.chart = house_echart(subburst.margin_chart(self.results, self.parent)) \
                .classes("w-full h-96").mark("subburst-chart")

    def _set_parent(self, value) -> None:
        self.parent = value
        self.draw()

    def _calibration(self) -> None:
        rows = subburst.calibration_rows(self.results)
        with ui.card().classes("w-full").mark("subburst-calibration"):
            ui.label("Calibrated pattern weights").classes("text-sm font-medium")
            ui.label("A duration that did not converge could not be brought to neutral by weighting "
                     "alone; patterns at the floor are all but excluded. Either way, filtering is the "
                     "robust choice. A weights file can be given to a row's 'TP weights' column to "
                     "sample with those weights.").classes("text-xs text-muted")
            columns = [{"name": k, "label": label, "field": k, "align": "left"}
                       for k, label in (("duration", "Storm"), ("converged", "Neutral"),
                                        ("iterations", "Iterations"), ("margin", "Worst margin"),
                                        ("floored", "Patterns at the floor"),
                                        ("weights_file", "Weights file"))]
            ui.table(columns=columns, rows=rows, row_key="duration").classes("w-full") \
                .props("dense flat").mark("subburst-weights")
            if subburst.has_curves(self.results):
                with ui.row().classes("items-end gap-4"):
                    ui.number("AEP of a level (m AHD)", value=self.level, format="%.2f",
                              on_change=lambda e: self._set_level(e.value)) \
                        .classes("w-56").props("dense").mark("subburst-level")
                    ui.label("The dam crest, say: its AEP on each curve.").classes("text-xs text-muted")
                cases, level_rows = subburst.level_rows(self.results, self.level)
                columns = [{"name": "aep", "label": "AEP (1 in X)", "field": "aep", "align": "left"}]
                columns += [{"name": c, "label": subburst.CASE_LABELS[c], "field": c, "align": "right"}
                            for c in cases]
                ui.table(columns=columns, rows=level_rows, row_key="aep").classes("w-full") \
                    .props("dense flat").mark("subburst-levels")

    def _set_level(self, value) -> None:
        self.level = value if value not in ("", None) else None
        self.draw()

    def _files(self) -> None:
        with ui.expansion("Databases read").classes("w-full"):
            for role in (subburst.TESTED, subburst.COMPARE):
                entry = (self.results.get("groups") or {}).get(role)
                if not entry:
                    continue
                ui.label(entry.get("label", role)).classes("text-sm font-medium")
                for duration, data in entry.get("durations", {}).items():
                    ui.label(f"{duration}: {data.get('database', '')}").classes("text-xs text-muted")
