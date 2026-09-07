"""The representative event analysis, against a synthetic Monte Carlo database.

The fixture carries the column names ``lib/MCScheme.py:26-31`` and
``lib/Simulator.py`` actually write, including the mixed units - ``rain_aep`` in
'1 in X' and the TPT columns as probabilities - because confusing those two is
the failure this module is most exposed to and it would not look wrong on a
plot.
"""

from __future__ import annotations

import json
import math

import pandas as pd
import pytest

from lib import RepresentativeEvents as events

TARGET = 1000.0


def realisation(sim_id, rain_aep, level_aep, **overrides):
    """One mcdf row, with every column a real database has."""
    row = {
        "m": sim_id // 10, "n": sim_id % 10,
        "rain_z": events.normal_variate(rain_aep),
        "rain_aep": rain_aep,
        "mean_rain_mm": 250.0,
        "tp": 3, "storm_method": "ARR point", "tp_frequency": "rare",
        "il_p": 0.5, "il_scaling": 1.0,
        "preburst_p": 0.5, "preburst_proportion": 0.15, "preburst_mm": 37.5,
        "initial_loss": 30.0,
        "cl_p": 0.5, "cl_scaling": 1.0, "continuing_loss": 2.5,
        "residual_depth": 7.5,
        "lake_z": 0.0,
        "embedded_bursts": events.NO_EMBEDDED_BURSTS,
        "ADV": 210_000.0,
        "subburst_2h": 60.0, "ifd_2h": 100.0,
        "subburst_6h": 120.0, "ifd_6h": 160.0,
        "inflow": 3000.0, "level": 220.0, "outflow": 2500.0,
        "level_aep": 1.0 / level_aep,          # probability, as the TPT writes it
        "inflow_aep": 1.0 / level_aep,
        "outflow_aep": 1.0 / level_aep,
    }
    row.update(overrides)
    return row


@pytest.fixture
def mcdf():
    """Seven realisations against a 1 in 1,000 lake level loading.

    0 is the answer: rainfall and flood both 1 in 1,000, nothing flagged.
    """
    rows = {
        0: realisation(0, 1000, 1000),
        1: realisation(1, 100, 1000, level=220.5),              # coincidence
        2: realisation(2, 1000, 100, level=214.0),              # rain right, flood not
        3: realisation(3, 900, 1100, level=220.1,               # embedded burst
                       embedded_bursts="Unfiltered embedded bursts: 2h burst exceeds by 20.0%",
                       subburst_2h=120.0),
        4: realisation(4, 900, 1100, level=220.1, preburst_p=0.88,
                       preburst_proportion=0.6, preburst_mm=150.0),
        5: realisation(5, 900, 1100, level=220.1, lake_z=2.0, ADV=260_000.0),
        6: realisation(6, 5_000_000, 1000, level=221.0),        # beyond the PMP
    }
    frame = pd.DataFrame.from_dict(rows, orient="index")
    return frame


@pytest.fixture
def scored(mcdf):
    return events.score(events.prepare(mcdf, "level"), TARGET)


LEVEL_CURVE = {10: 214.0, 100: 216.0, 1000: 220.0, 10_000: 223.0}


# -- units -------------------------------------------------------------------

def test_the_tpt_column_is_a_probability_and_rain_aep_is_not(mcdf):
    """1/p for the result, as read for the rainfall. Mixing them is invisible."""
    prepared = events.prepare(mcdf, "level")
    assert prepared.loc[0, "result_aep"] == pytest.approx(1000.0)
    assert prepared.loc[0, "rain_aep"] == pytest.approx(1000.0)
    # Both describe the same rarity, so their variates agree.
    assert prepared.loc[0, "z_result"] == pytest.approx(prepared.loc[0, "z_rain"])


def test_a_database_without_the_analysis_says_so(mcdf):
    raw = mcdf.drop(columns=["level_aep"])
    with pytest.raises(ValueError, match="Analyse results"):
        events.prepare(raw, "level")


# -- ranking -----------------------------------------------------------------

def test_the_aep_neutral_event_wins(scored):
    ranking = events.rank(scored, count=3)
    assert ranking.candidates.index[0] == 0


def test_a_coincidence_ranks_below_a_neutral_event(scored):
    """Row 1 reaches the right level off 1 in 100 rainfall."""
    order = list(events.rank(scored, count=7).candidates.index)
    assert order.index(0) < order.index(1)


def test_distance_is_measured_in_z(scored):
    """Not in '1 in X', where the rare end swamps everything.

    Row 1 and row 2 are each one decade of AEP away from the target, one in
    rainfall and one in flood. In z they are comparable; in '1 in X' row 2 (out
    by 900) would look ten times worse than row 1 (out by 900 the other way).
    """
    assert scored.loc[1, "delta_z"] == pytest.approx(scored.loc[2, "delta_z"], rel=0.35)


def test_the_pmp_cap_drops_rainfall_beyond_it(scored):
    ranking = events.rank(scored, events.Filters(aep_of_pmp=1e6), count=7)
    assert 6 not in ranking.candidates.index
    assert sum(ranking.excluded.values()) == 1


def test_without_the_cap_nothing_is_dropped_for_rarity(scored):
    ranking = events.rank(scored, count=7)
    assert 6 in ranking.candidates.index
    assert not ranking.excluded


def test_a_distance_limit_reports_what_it_removed(scored):
    ranking = events.rank(scored, events.Filters(max_delta_z=0.5), count=7)
    assert 1 not in ranking.candidates.index
    assert any("in z" in reason for reason in ranking.excluded)


# -- the flags ---------------------------------------------------------------

def test_an_embedded_burst_is_flagged_by_comment_and_by_ratio(scored):
    flags = events.flags_for(scored.loc[3])
    assert any("Unfiltered embedded bursts" in flag for flag in flags)
    assert any("sub-burst 1.20x" in flag for flag in flags)


def test_a_high_preburst_is_flagged(scored):
    assert any("pre-burst percentile 0.88" in flag
               for flag in events.flags_for(scored.loc[4]))


def test_an_unusual_antecedent_storage_is_flagged(scored):
    flags = events.flags_for(scored.loc[5])
    assert any("high antecedent storage" in flag for flag in flags)


def test_a_clean_event_carries_no_flags(scored):
    assert events.flags_for(scored.loc[0]) == ()


def test_flags_are_reported_by_default_not_excluded(scored):
    """The whole point: a flagged event is still offered, with its reason."""
    ranking = events.rank(scored, count=7)
    assert {3, 4, 5} <= set(ranking.candidates.index)
    assert ranking.candidates.loc[4, "flags"]


def test_excluding_embedded_bursts_removes_only_those(scored):
    ranking = events.rank(scored, events.Filters(exclude_embedded=True), count=7)
    assert 3 not in ranking.candidates.index
    assert {4, 5} <= set(ranking.candidates.index)


def test_a_high_preburst_is_not_an_embedded_burst(scored):
    """'pre-burst' contains 'burst'. Ask the data, never the flag text."""
    assert events.has_embedded_burst(scored.loc[3])
    assert not events.has_embedded_burst(scored.loc[4])
    assert not events.has_embedded_burst(scored.loc[0])


def test_excluding_everything_flagged_leaves_the_clean_events(scored):
    ranking = events.rank(scored, events.Filters(exclude_flagged=True), count=7)
    assert not ({3, 4, 5} & set(ranking.candidates.index))
    assert 0 in ranking.candidates.index


def test_an_older_database_without_subburst_columns_still_works(mcdf):
    """The ratio is unknown; the comment still says whether there was one."""
    older = mcdf.drop(columns=["subburst_2h", "ifd_2h", "subburst_6h", "ifd_6h"])
    scored = events.score(events.prepare(older, "level"), TARGET)
    assert scored["subburst_ratio"].isna().all()
    flags = events.flags_for(scored.loc[3])
    assert any("Unfiltered embedded bursts" in flag for flag in flags)
    assert not any("sub-burst" in flag for flag in flags)


def test_the_vectorised_masks_agree_with_the_row_by_row_rules(scored):
    """rank() filters a whole database at once; flags_for writes the words.

    Two implementations of one rule, for speed - so they are held together
    here rather than left to drift.
    """
    assert list(events.embedded_mask(scored)) == [
        events.has_embedded_burst(row) for _, row in scored.iterrows()]
    assert list(events.flag_mask(scored)) == [
        bool(events.flags_for(row)) for _, row in scored.iterrows()]


def test_the_masks_still_agree_with_the_bands_moved(scored):
    filters = events.Filters(preburst_band=0.1, lake_band=0.5, loss_band=0.05,
                             subburst_limit=0.5)
    assert list(events.flag_mask(scored, filters)) == [
        bool(events.flags_for(row, filters)) for _, row in scored.iterrows()]


def test_ranking_a_full_size_database_is_quick(mcdf):
    """The page re-ranks on every control change, and an mcdf is m x n rows."""
    import time

    big = pd.concat([mcdf] * 1500, ignore_index=True)      # 10,500 realisations
    prepared = events.prepare(big, "level")

    start = time.perf_counter()
    for _ in range(5):
        ranking = events.rank(events.score(prepared, TARGET), count=10)
    elapsed = (time.perf_counter() - start) / 5

    assert len(ranking.candidates) == 10
    assert elapsed < 0.25, f"{elapsed:.3f}s per re-rank of {len(big):,} rows"


# -- a lake level as the target ----------------------------------------------

def test_a_level_on_the_curve_becomes_an_aep():
    lookup = events.aep_for_level(LEVEL_CURVE, 220.0)
    assert lookup.found
    assert lookup.aep == pytest.approx(1000.0, rel=0.01)


def test_a_level_between_points_interpolates():
    lookup = events.aep_for_level(LEVEL_CURVE, 221.0)
    assert 1000 < lookup.aep < 10_000


def test_a_level_above_the_curve_is_reported_not_invented():
    lookup = events.aep_for_level(LEVEL_CURVE, 230.0)
    assert not lookup.found
    assert lookup.above_curve
    assert "above the top of the curve" in lookup.note


def test_a_level_above_the_curve_ranks_by_the_highest_events(scored):
    ranking = events.rank(scored, count=3, above_curve=True)
    assert ranking.candidates.index[0] == 6          # level 221.0, the highest
    assert ranking.candidates["level"].is_monotonic_decreasing


def test_an_aep_beyond_the_reach_of_a_variate_does_not_raise():
    """aep_of_variate goes infinite once the upper tail underflows.

    A level curve that flattens off puts the interpolation there, and an
    exception out of the middle of a page redraw is not the answer.
    """
    assert events.normal_variate(float("inf")) != events.normal_variate(float("inf"))
    flat = {10: 214.0, 100: 216.0, 1000: 220.0, 10_000: 220.000000001}
    lookup = events.aep_for_level(flat, 220.0000000005)
    assert lookup.above_curve or (lookup.aep and lookup.aep < float("inf"))


def test_a_level_below_the_curve_says_so():
    lookup = events.aep_for_level(LEVEL_CURVE, 100.0)
    assert not lookup.found
    assert not lookup.above_curve


# -- targets on disk ---------------------------------------------------------

def test_a_target_round_trips_through_the_selection_file(tmp_path):
    targets = [events.Target(kind="level", value=220.5, result_type="level",
                             source="TFD_mc_48h", count=5, picked=17),
               events.Target(kind="aep", value=2000, rain_aep=1000)]
    path = tmp_path / events.SELECTION_FILE
    path.write_text(json.dumps(events.selection_payload(targets, {"group": "GWL1.3"})))

    loaded, settings = events.read_selection(path)
    assert [t.value for t in loaded] == [220.5, 2000]
    assert loaded[0].picked == 17
    assert loaded[1].rain_aep == 1000
    assert settings["group"] == "GWL1.3"


def test_a_missing_selection_file_is_empty_not_an_error(tmp_path):
    assert events.read_selection(tmp_path / "nothing.json") == ([], {})


def test_target_labels_read_as_loadings():
    assert events.Target(kind="aep", value=2000).label == "1 in 2,000"
    assert events.Target(kind="level", value=220.5).label == "220.5 m AHD"


# -- the AEP of the PMP ------------------------------------------------------

def test_the_pmp_aep_is_read_from_the_ifd_files_config(tmp_path):
    """Where a Monte Carlo run gets it - not from the method config file."""
    (tmp_path / "ifd_files.json").write_text(json.dumps({"AEP_of_PMP": 2_000_000}))
    storm = tmp_path / "storm_config.json"
    storm.write_text(json.dumps({"file_paths": {"rare_ifds": "ifd_files.json"}}))
    assert events.pmp_aep_from_storm_config(storm) == 2_000_000


def test_a_windows_separator_in_the_chain_still_resolves(tmp_path):
    folder = tmp_path / "ifds"
    folder.mkdir()
    (folder / "ifd_files.json").write_text(json.dumps({"AEP_of_PMP": 1e7}))
    storm = tmp_path / "storm_config.json"
    storm.write_text(json.dumps({"file_paths": {"rare_ifds": r"ifds\ifd_files.json"}}))
    assert events.pmp_aep_from_storm_config(storm) == 1e7


def test_a_broken_chain_returns_none_rather_than_failing(tmp_path):
    storm = tmp_path / "storm_config.json"
    storm.write_text(json.dumps({"file_paths": {"rare_ifds": "missing.json"}}))
    assert events.pmp_aep_from_storm_config(storm) is None
    assert events.pmp_aep_from_storm_config(None) is None


# -- the tie back to the hydrographs -----------------------------------------

def test_the_sim_label_matches_the_stored_hydrograph_columns():
    """URBSmodel.py:653 zero-pads to five."""
    assert events.sim_label(0) == "sim_00000"
    assert events.sim_label(42) == "sim_00042"
    assert events.sim_label(12345) == "sim_12345"


# -- reading the stored series -----------------------------------------------

def _hydrograph_frame():
    return pd.DataFrame({"sim_00000": [1.0, 2.0, 3.0],
                         "sim_00001": [4.0, 5.0, 6.0]},
                        index=pd.Index([0.0, 0.5, 1.0], name="time"))


def test_hydrographs_read_from_a_csv(tmp_path):
    path = tmp_path / "flows_inflows.csv"
    _hydrograph_frame().to_csv(path)
    frame = events.read_hydrographs(path)
    assert list(frame.columns) == ["sim_00000", "sim_00001"]
    assert frame["sim_00001"].tolist() == [4.0, 5.0, 6.0]


def test_hydrographs_read_from_a_parquet(tmp_path):
    """A routed row's inflows come from the sims list, and are often parquet.

    Reading one with read_csv raises UnicodeDecodeError on the parquet magic,
    which used to end the whole plotting run rather than one panel.
    """
    pytest.importorskip("pyarrow")
    path = tmp_path / "flows_inflows.parquet"
    # index=False is how they are written: the time axis arrives as an
    # ordinary first column and has to be promoted back.
    _hydrograph_frame().reset_index().to_parquet(path, index=False)

    frame = events.read_hydrographs(path)
    assert list(frame.columns) == ["sim_00000", "sim_00001"]
    assert frame.index.tolist() == [0.0, 0.5, 1.0]
    assert frame["sim_00000"].tolist() == [1.0, 2.0, 3.0]


def test_a_parquet_read_as_a_csv_is_the_failure_being_fixed(tmp_path):
    pytest.importorskip("pyarrow")
    path = tmp_path / "flows_inflows.parquet"
    _hydrograph_frame().reset_index().to_parquet(path, index=False)
    # The exact exception is pandas', not ours - only that it cannot be done.
    with pytest.raises((UnicodeDecodeError, ValueError, pd.errors.ParserError)):
        pd.read_csv(path, index_col=0)


# -- ranking on the result rather than on neutrality -------------------------
#
# Hitting the loading is often what the event is for - a gate operation or a
# dambreak run needs the lake at a level - and AEP neutrality is the thing
# traded against it. So the two are offered as orders, and whichever is not
# ranked on is still in the table.

def test_the_design_value_is_read_back_off_the_curve():
    """value_for_aep inverts aep_for_level on the same interpolation."""
    assert events.value_for_aep(LEVEL_CURVE, 1000) == pytest.approx(220.0)
    level = events.value_for_aep(LEVEL_CURVE, 3000)
    assert 220.0 < level < 223.0
    assert events.aep_for_level(LEVEL_CURVE, level).aep == pytest.approx(3000, rel=0.01)


def test_an_aep_off_the_curve_has_no_design_value():
    """Reported as NaN, not extrapolated - the caller says so and falls back."""
    assert math.isnan(events.value_for_aep(LEVEL_CURVE, 1_000_000))
    assert math.isnan(events.value_for_aep(LEVEL_CURVE, 2))


def test_ranking_on_the_result_takes_the_event_that_reaches_the_level(mcdf, scored):
    """Row 1 reaches 220.5 m off 1 in 100 rainfall; row 0 is the neutral one.

    Ranked on neutrality row 0 wins, which is the default and right when the
    AEP is the loading. Ranked on the level, the event that actually gets there
    wins and its rainfall is left to the flags to argue about.
    """
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.5)
    assert events.rank(on_level, count=3, order=events.RESULT).candidates.index[0] == 1
    assert events.rank(on_level, count=3).candidates.index[0] == 0


def test_the_result_ranking_breaks_its_ties_on_neutrality(mcdf):
    """Four events reach the same level; the least neutral one comes last."""
    mcdf.loc[1, "level"] = 220.1                  # row 1's rainfall is 1 in 100
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.1)
    order = list(events.rank(on_level, count=7, order=events.RESULT).candidates.index)
    tied = [sim for sim in order if sim in (1, 3, 4, 5)]
    assert tied == [3, 4, 5, 1], "the 1 in 100 rainfall should lose the tie"


def test_without_a_design_value_the_result_ranking_says_what_it_used(scored):
    """No level to measure against - the result axis in z is the same order."""
    ranking = events.rank(scored, count=7, order=events.RESULT)
    assert any("no design value" in note for note in ranking.notes)
    distances = ranking.candidates["d_z_result"].abs()
    assert distances.is_monotonic_increasing


def test_the_delta_to_the_loading_is_reported_either_way(mcdf, scored):
    """The column is there whether or not it is ranked on - and NaN, not zero,
    where there is no design value, so nothing reads as a perfect match."""
    assert scored["delta_value"].isna().all()
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.5)
    assert on_level.loc[1, "delta_value"] == pytest.approx(0.0)
    assert on_level.loc[2, "delta_value"] == pytest.approx(6.5)


def test_a_band_makes_near_enough_levels_equal_and_neutrality_decide(mcdf):
    """Row 0 is 90 mm from the loading and neutral; rows 3-5 are 10 mm off it.

    Without a band the 10 mm wins on a difference the rating curve cannot
    resolve. Inside a 100 mm band all four count as reaching the loading, and
    the AEP neutral one comes first.
    """
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.09)
    close = events.rank(on_level, count=7, order=events.RESULT)
    assert close.candidates.index[0] in (3, 4, 5)

    banded = events.rank(on_level, count=7, order=events.RESULT, band=0.1)
    assert banded.candidates.index[0] == 0


def test_the_band_only_groups_what_is_inside_it(mcdf):
    """A 20 mm band leaves the 90 mm event where it was - behind the closer ones."""
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.09)
    ranking = events.rank(on_level, count=7, order=events.RESULT, band=0.02)
    order = list(ranking.candidates.index)
    assert order[0] in (3, 4, 5)
    assert order.index(0) > order.index(3)


def test_no_band_is_the_plain_distance(mcdf):
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.09)
    distances = events.banded(on_level["delta_value"], 0, 0)
    assert distances.equals(on_level["delta_value"])


def test_the_band_survives_a_missing_distance(mcdf):
    """NaN stays NaN - a band is not a reason to call an unknown result a match."""
    distances = events.banded(pd.Series([float("nan"), 0.015, 0.045]), 0.02)
    assert math.isnan(distances.iloc[0])
    assert distances.iloc[1] == 0
    assert distances.iloc[2] > 0


# -- the band and the rounding are two different questions --------------------

def test_the_band_is_what_counts_as_reaching_the_loading():
    """Everything inside it is one group, whatever the rounding says."""
    distances = pd.Series([0.000, 0.012, 0.020, 0.021])
    key = events.banded(distances, band=0.02, rounding=0.01)
    assert list(key[:3]) == [0, 0, 0]
    assert key.iloc[3] > 0


def test_the_rounding_groups_what_is_outside_the_band():
    """32 mm and 34 mm are the same distance away to any defensible precision."""
    distances = pd.Series([0.032, 0.034, 0.055])
    key = events.banded(distances, band=0.02, rounding=0.01)
    assert key.iloc[0] == key.iloc[1]
    assert key.iloc[2] != key.iloc[0]


def test_a_coarse_rounding_never_swallows_the_band():
    """Outside the band is outside it, however coarse the grid."""
    key = events.banded(pd.Series([0.010, 0.030]), band=0.02, rounding=0.5)
    assert key.iloc[0] == 0
    assert key.iloc[1] > 0


def test_rounding_alone_still_ties_the_near_enough_events():
    """No band: the closest event leads, but 32 mm and 34 mm still tie."""
    key = events.banded(pd.Series([0.032, 0.034, 0.004]), band=0, rounding=0.01)
    assert key.iloc[0] == key.iloc[1]
    assert key.iloc[2] < key.iloc[0]


def test_the_rounding_decides_between_equally_distant_events(mcdf):
    """Rows 3-5 are 40 mm out and row 0 is 60 mm; on a 50 mm grid they tie.

    Row 0 is the AEP neutral one, so it comes first once they do - which it
    does not on the raw distance, where 20 mm of lake level decides it.
    """
    on_level = events.score(events.prepare(mcdf, "level"), TARGET,
                            target_value=220.06)
    raw = events.rank(on_level, count=7, order=events.RESULT)
    assert raw.candidates.index[0] in (3, 4, 5)

    rounded = events.rank(on_level, count=7, order=events.RESULT, rounding=0.05)
    assert rounded.candidates.index[0] == 0


# -- the design curve and the realisations are two different answers ----------

def _run(levels, variates):
    return pd.DataFrame({"result_value": levels, "z_result": variates})


def test_a_level_is_placed_where_this_run_reached_it():
    """Not where the design curve puts it - the two are not the same number."""
    run = _run([214.0, 220.0, 226.0], [2.0, 3.0, 4.0])
    assert events.variate_at_value(run, 220.0) == pytest.approx(3.0)
    assert events.variate_at_value(run, 223.0) == pytest.approx(3.5)


def test_a_level_outside_the_run_has_no_place_in_it():
    """NaN, not the nearest end: nothing in the run reached it."""
    run = _run([214.0, 220.0, 226.0], [2.0, 3.0, 4.0])
    assert math.isnan(events.variate_at_value(run, 230.0))
    assert math.isnan(events.variate_at_value(run, 210.0))
    assert math.isnan(events.variate_at_value(run, None))


def test_a_run_with_nothing_analysed_places_nothing():
    assert math.isnan(events.variate_at_value(_run([], []), 220.0))
    assert math.isnan(events.variate_at_value(
        _run([220.0], [float("nan")]), 220.0))


def test_repeated_levels_are_placed_at_the_first_one_to_reach_it():
    """Exceedance, not equality: everything rarer reaches the level too.

    So a flat spot in the curve is entered at its frequent end - the AEP of
    'this run reaches 220 m' is the first realisation that does.
    """
    run = _run([214.0, 220.0, 220.0, 226.0], [2.0, 3.0, 3.4, 4.0])
    assert events.variate_at_value(run, 220.0) == pytest.approx(3.0)
