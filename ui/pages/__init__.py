"""Page registry.

Pages are plain functions bound to routes by ``register_all()`` rather than by
import-time ``@ui.page`` decorators, so tests can re-register them onto a fresh
NiceGUI app. Same arrangement as swe2d_ui/pages/__init__.py.
"""

from __future__ import annotations

from nicegui import ui


def register_all() -> None:
    from pages import (downstream, edit, ensembleresults, events, figures, history,
                       lakefreq, lakerecord, manual, pmf, report, results, run, select, simulations,
                       study)

    for route, handler in [
        ("/", simulations.simulations_page),
        ("/study", study.study_page),
        ("/select", select.select_page),
        ("/run", run.run_page),
        ("/results", results.results_page),
        ("/ensemble", ensembleresults.ensemble_page),
        ("/events", events.events_page),
        ("/lake-levels", lakefreq.lake_levels_page),
        ("/downstream", downstream.downstream_page),
        ("/report", report.report_page),
        ("/pmf", pmf.pmf_page),
        ("/figures", figures.figures_page),
        ("/lake-record", lakerecord.lake_record_page),
        ("/edit", edit.edit_page),
        ("/history", history.history_page),
        ("/manual", manual.manual_page),
        ("/manual/{doc}", manual.manual_page),
    ]:
        ui.page(route)(handler)
