"""Reading a realisation's stored hydrographs, and saying what shape it is.

The point of the shape is the choice it supports: a representative event is
wanted simple - one rise, one peak, one recession - and the mcdf records the
peak but nothing about what happens either side of it. So the tests here are
mostly about the peak counting being usable: a real second flood counts, and
the steps and shoulders every routed outflow has do not.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest

from core import events, hydrographs
from state import Project
from test_events import build


@pytest.fixture(autouse=True)
def clear_cache():
    hydrographs.forget_cached()
    events.forget_cached()
    yield
    hydrographs.forget_cached()
    events.forget_cached()


def series(values, dt=1.0):
    index = [step * dt for step in range(len(values))]
    return pd.Series(values, index=index)


# -- the shape ---------------------------------------------------------------

def test_a_plain_flood_is_single_peaked():
    shape = hydrographs.shape_of(series([0, 5, 20, 60, 100, 70, 40, 20, 8, 2]))
    assert shape.single
    assert shape.peaks == 1
    assert shape.summary == "single peaked"
    assert shape.peak_time == pytest.approx(4.0)


def test_a_second_flood_counts():
    """Down to a fifth and back up to nine tenths - two events, not one."""
    shape = hydrographs.shape_of(
        series([0, 40, 100, 60, 20, 45, 90, 50, 20, 5]))
    assert shape.peaks == 2
    assert shape.summary == "2 peaks"


def test_a_step_on_the_recession_is_not_a_second_peak():
    """Every routed outflow has these - a gate, a fuse plug, a spillway lip.

    Counting bare local maxima would call almost every event multi-peaked and
    the column would say nothing at all.
    """
    shape = hydrographs.shape_of(
        series([0, 50, 100, 80, 78, 79, 60, 30, 10, 2]))
    assert shape.single


def test_a_small_early_rise_is_not_a_peak():
    """Below a fifth of the flood it is not what makes an event awkward."""
    shape = hydrographs.shape_of(
        series([0, 12, 4, 30, 100, 70, 30, 10, 2, 0]))
    assert shape.single


def test_a_shoulder_on_the_rising_limb_is_not_a_peak():
    shape = hydrographs.shape_of(
        series([0, 30, 55, 54, 70, 100, 60, 25, 8, 1]))
    assert shape.single


def test_an_empty_or_dry_hydrograph_says_so():
    assert hydrographs.shape_of(series([])).summary == "nothing stored"
    assert hydrographs.shape_of(series([0, 0, 0])).summary == "no flow"


def test_a_flat_top_is_one_peak():
    """A dam spilling at capacity holds its peak for hours."""
    shape = hydrographs.shape_of(
        series([0, 40, 100, 100, 100, 100, 60, 20, 4, 0]))
    assert shape.single


# -- reading the files -------------------------------------------------------

SIMS = ("sim_00000", "sim_00001", "sim_00002")


def stored(tmp_path, kinds=("inflows", "outflows", "levels")):
    """A monte carlo project whose runs stored their hydrographs."""
    return Project.open(write_hydrographs(build(tmp_path), tmp_path, kinds))


def write_hydrographs(config, tmp_path, kinds=("inflows", "outflows", "levels")):
    """The files a run with 'Store hydrographs' leaves in ../urbs_results."""
    folder = tmp_path / "sims_mc" / "urbs_results"
    folder.mkdir(parents=True, exist_ok=True)
    times = [step * 0.5 for step in range(41)]
    for duration in (12, 24, 48, 120):
        for kind in kinds:
            frame = pd.DataFrame({"Time": times})
            for position, name in enumerate(SIMS):
                peak = 8.0 + position
                base = 200.0 if kind == "levels" else 0.0
                frame[name] = [base + 100.0 / (1.0 + abs(hour - peak))
                               for hour in times]
            frame.set_index("Time").to_csv(
                folder / f"TFD_mc_{duration}h_GWL1p3_{kind}.csv")
    return config


def source_of(project, label="24h"):
    return next(source for source in events.sources_for_rows(project)
                if source.label == label)


def test_the_files_are_found_where_bryan_writes_them(tmp_path):
    project = stored(tmp_path)
    found = hydrographs.paths_for(project, source_of(project))
    assert set(found) == {"inflows", "outflows", "levels"}
    assert all(path.is_file() for path in found.values())


def test_a_run_that_stored_nothing_has_no_files(tmp_path):
    project = Project.open(build(tmp_path))
    assert hydrographs.paths_for(project, source_of(project)) == {}


def test_one_simulation_comes_back_with_its_three_series(tmp_path):
    project = stored(tmp_path)
    series, notes = hydrographs.series_for(project, source_of(project), 1)
    assert set(series) == {"inflows", "outflows", "levels"}
    assert notes == []
    assert series["outflows"].idxmax() == pytest.approx(9.0)


def test_a_run_without_stored_hydrographs_says_which_switch(tmp_path):
    project = Project.open(build(tmp_path))
    series, notes = hydrographs.series_for(project, source_of(project), 1)
    assert series == {}
    assert any("Store hydrographs" in note for note in notes)


def test_a_simulation_that_was_not_stored_is_a_note(tmp_path):
    project = stored(tmp_path)
    series, notes = hydrographs.series_for(project, source_of(project), 99)
    assert series == {}
    assert any("sim_00099" in note for note in notes)


def test_the_shapes_of_a_whole_table_come_from_one_read(tmp_path):
    project = stored(tmp_path)
    shapes, notes = hydrographs.shapes_for(project, source_of(project), [0, 1, 2])
    assert set(shapes) == {0, 1, 2}
    assert all(shape.single for shape in shapes.values())
    assert notes == []


def test_missing_candidates_are_counted_not_dropped_silently(tmp_path):
    project = stored(tmp_path)
    shapes, notes = hydrographs.shapes_for(project, source_of(project), [0, 99])
    assert set(shapes) == {0}
    assert any("1 of the candidates" in note for note in notes)


def test_a_re_run_is_picked_up(tmp_path):
    """The cache is keyed on mtime and size, as the mcdf cache is."""
    project = stored(tmp_path)
    source = source_of(project)
    first, _ = hydrographs.series_for(project, source, 0)

    path = hydrographs.paths_for(project, source)["outflows"]
    frame = pd.read_csv(path, index_col=0)
    frame["sim_00000"] = frame["sim_00000"] * 2
    frame.to_csv(path)

    again, _ = hydrographs.series_for(project, source, 0)
    assert again["outflows"].max() == pytest.approx(first["outflows"].max() * 2)


def test_plotting_thins_a_long_series_without_losing_its_ends():
    long = series([float(step % 50) for step in range(5000)], dt=0.1)
    times, values = hydrographs.for_plot(long, limit=200)
    assert len(times) <= 250
    assert times[0] == pytest.approx(0.0)
    assert not any(math.isnan(value) for value in values)
