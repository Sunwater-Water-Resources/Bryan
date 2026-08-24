"""Run the real util/CriticalDurationAnalysis.py, through the export planner.

The point of shelling out rather than recomputing is that the exported files
match what the study post-processing already produces - which is only true if
the script actually accepts what the UI builds for it. Skipped where Bryan's
own dependencies (scipy, matplotlib) are not importable, as the UI environment
deliberately does not have them.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from core import critexport, results
from test_critexport import sources
from test_real_bryan import BRYAN_PYTHON, needs_bryan
from test_results import STANDARD_AEPS, level_curve, mc_row, write_quantiles

pytestmark = needs_bryan


def perc_files(base: Path, values: dict) -> None:
    """The confidence percentile files MCScheme.smooth_percentiles writes."""
    for suffix in ("_perc", "_perc_smooth"):
        pd.DataFrame({"AEP": list(values),
                      "z": [results.normal_variate(aep) for aep in values],
                      "5%": [value * 0.9 for value in values.values()],
                      "95%": [value * 1.1 for value in values.values()],
                      }).to_csv(f"{base}{suffix}.csv", index=False)


def export(tmp_path, found, keys=("level",), **kwargs):
    plan = critexport.plan(found, list(keys), folder=tmp_path / "out",
                           base_name="CLD", python=BRYAN_PYTHON, **kwargs)
    assert plan.can_run, plan.problems
    return plan, [critexport.run(job) for job in plan.jobs]


@needs_bryan
def test_the_script_writes_the_table_and_the_plot(tmp_path):
    found = sources(tmp_path, [12, 24, 48, 96])
    plan, done = export(tmp_path, found)

    assert all(result.ok for result in done), done[0].output
    table = pd.read_csv(plan.jobs[0].csv, index_col=0)
    assert table.index.name == "aep (1 in x)"
    assert list(table.columns[:4]) == ["12h", "24h", "48h", "96h"]
    assert "max" in table.columns and "critical_duration" in table.columns
    assert plan.jobs[0].png.is_file()

    # whole durations stay whole in the labels, as the util scripts write them
    assert set(table["critical_duration"]) <= {"12h", "24h", "48h", "96h"}


@needs_bryan
def test_the_exported_critical_durations_match_what_the_page_showed(tmp_path):
    """The plot and the export must not disagree about which duration wins."""
    found = sources(tmp_path, [12, 24, 48, 96])
    plan, done = export(tmp_path, found)
    assert done[0].ok, done[0].output

    analysis = results.analyse(results.compare(found["level"]))
    exported = pd.read_csv(plan.jobs[0].csv, index_col=0)

    for aep, owner in analysis.critical.items():
        assert exported.loc[aep, "critical_duration"] == owner
        assert exported.loc[aep, "max"] == pytest.approx(analysis.envelope[aep])


@needs_bryan
def test_results_without_percentile_files_still_export(tmp_path):
    """Reservoir routing writes quantiles but no _perc_smooth files.

    UtilModule.analyse_percentiles used to pd.concat an empty list, so every
    re-routed result failed here instead of exporting a table with no
    confidence limits.
    """
    found = sources(tmp_path, [24, 48])
    plan, done = export(tmp_path, found)

    assert done[0].ok, done[0].output
    assert "no confidence limits" in done[0].output.lower()
    table = pd.read_csv(plan.jobs[0].csv, index_col=0)
    assert list(table.columns) == ["24h", "48h", "max", "critical_duration"]


@needs_bryan
def test_percentile_columns_come_through_when_the_files_are_there(tmp_path):
    found = sources(tmp_path, [24, 48])
    for duration in (24, 48):
        base = tmp_path / f"CLD_mc_{duration}h_E010_GWL0p3_level"
        perc_files(base, level_curve(duration))
    plan, done = export(tmp_path, found)

    assert done[0].ok, done[0].output
    table = pd.read_csv(plan.jobs[0].csv, index_col=0)
    assert "5%" in table.columns and "95%" in table.columns


@needs_bryan
def test_dropping_an_aep_that_is_not_there_does_not_fail(tmp_path):
    """The AEP set depends on the config and the AEP of the PMP."""
    found = sources(tmp_path, [24, 48])
    plan, done = export(tmp_path, found, drop_aeps=(2, 999999999))
    assert done[0].ok, done[0].output
    assert plan.jobs[0].png.is_file()


@needs_bryan
def test_a_missing_quantile_file_fails_loudly(tmp_path):
    """Rather than quietly analysing the durations that happen to be there."""
    found = sources(tmp_path, [24, 48])
    Path(found["level"][0].path).unlink()
    plan = critexport.plan(found, ["level"], folder=tmp_path / "out",
                           base_name="CLD", python=BRYAN_PYTHON)
    result = critexport.run(plan.jobs[0])

    assert not result.ok
    assert "not found" in result.output
    assert not plan.jobs[0].csv.is_file()


@needs_bryan
def test_an_aep_no_duration_reached_does_not_lose_the_analysis(tmp_path):
    """A dam that does not spill at frequent AEPs has no outflow quantile there.

    compute_std_quantiles interpolates in log space, so those rows come out
    empty for every duration - and idxmax raises on an all-NA row rather than
    returning NaN. It cost the whole outflow export on the real Callide E010
    results.
    """
    for duration in (24, 48):
        name = f"CLD_mc_{duration}h_E010_GWL0p3"
        write_quantiles(Path(f"{tmp_path / name}_outflow.csv"), "outflow",
                        {2: float("nan"), 5: float("nan"),
                         10: 100.0 * duration, 100: 500.0 * duration})
    frame = pd.DataFrame([mc_row(tmp_path, f"CLD_mc_{d}h_E010_GWL0p3", d)
                          for d in (24, 48)])
    found = results.sources_for_rows(frame, tmp_path)
    plan, done = export(tmp_path, found, keys=("outflow",))

    assert done[0].ok, done[0].output
    table = pd.read_csv(plan.jobs[0].csv, index_col=0)
    assert pd.isna(table.loc[2, "critical_duration"])     # nothing reached it
    assert table.loc[10, "critical_duration"] == "48h"


@needs_bryan
def test_a_volume_export_reads_the_column_inside_the_file(tmp_path):
    found = sources(tmp_path, [24, 48], key="inflowVol24h", column="Vol24h")
    plan, done = export(tmp_path, found, keys=("inflowVol24h",))

    assert done[0].ok, done[0].output
    assert plan.jobs[0].csv.name == "CLD_inflowVol24h_critical.csv"
    table = pd.read_csv(plan.jobs[0].csv, index_col=0)
    assert list(table.columns[:2]) == ["24h", "48h"]
