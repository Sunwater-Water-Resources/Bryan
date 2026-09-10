"""The very rare - extreme spatial transition follows the AEP, not the storm method.

Under ``interpolate_weights`` the spatial pattern is interpolated from the gridded
ARR pattern at the lower smoothing bound towards the GSDM/GTSMR pattern at the
upper one, as a function of the standard normal variate alone.

In the AEP changeover zone the storm method is sampled: ARR half the time, GSDM or
GTSMR the other half (``sample_storm_method``). That coin toss is about *temporal*
patterns, and it was deciding the spatial pattern as well - ``get_depth_z`` returned
the unsmoothed gridded depths whenever the sampled method was ARR, so half the Monte
Carlo realisations in the zone skipped the transition entirely, and an ensemble run
with ``interim_for_ensemble = "arr"`` skipped all of them. It also left the mean
spatial pattern with a step at 1 in 2,000, which is the discontinuity this smoothing
method exists to remove.

The extreme pattern the weights head for is now sampled in its own right
(``StormBurst.sample_spatial_method``), on the same duration rule an extreme storm
uses, and passed to ``get_depth_z`` as ``spatial_method``. A caller that supplies
none falls back to the storm method - so every database written before
10 September 2026 rebuilds as the run that wrote it.

Run with Bryan's own interpreter, ``python -m pytest tests -q``: ``lib/RainfallScheme``
imports scipy.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
pd = pytest.importorskip("pandas")

from scipy.special import ndtri                                     # noqa: E402

from lib.RainfallScheme import ifdCurves                            # noqa: E402
from lib.StormGenerator import StormBurst                           # noqa: E402

DURATION = 24
SUBCATCHMENTS = ['sub_A', 'sub_B']
AREAS = pd.Series([50.0, 50.0], index=SUBCATCHMENTS)
# Area-weighted mean of each pattern is 1.0, so the catchment average is preserved.
LOWER_BOUND = pd.Series([1.2, 0.8], index=SUBCATCHMENTS)
PMP_GTSMR = pd.Series([0.6, 1.4], index=SUBCATCHMENTS)
GRIDDED = pd.Series([120.0, 80.0], index=SUBCATCHMENTS)             # catchment average 100 mm


class FakeDurationCurve:
    """Stands in for a DurationCurve: the ARF-adjusted depths by subcatchment."""

    def __init__(self, depths):
        self.depths = depths

    def get_depth_z(self, z):
        return self.depths


def curves_with_weights(tmp_path, bounds=(100, 2000)):
    """An ``ifdCurves`` on the weights method, with the depths and patterns injected.

    Everything the rainfall import would have built is faked - this is a test of the
    branch in ``get_depth_z``, not of the IFD file reading.
    """
    config = tmp_path / 'ifd_files.json'
    config.write_text(json.dumps({
        'AEP_of_PMP': 1e7,
        'extreme_rainfall_interpolation_method': 'GEV',
        'extreme_spatial_smoothing_method': {'interpolate_weights': list(bounds)},
    }))
    curves = ifdCurves(str(config))
    curves.arr_ifd_curves[DURATION] = FakeDurationCurve(GRIDDED)
    curves.catchment_ave_curves[DURATION] = FakeDurationCurve(pd.Series([300.0]))
    pmp_pattern = pd.DataFrame({'GTSMR': PMP_GTSMR, f'GSDM_{DURATION}': PMP_GTSMR})
    curves.smoothed_weights.set_pattern(pd.DataFrame({DURATION: LOWER_BOUND}), pmp_pattern, AREAS)
    return curves


def z_of(aep):
    return ndtri(1 - 1 / aep)


def catchment_average(depths):
    return float((depths * AREAS).sum() / AREAS.sum())


# ---------------------------------------------------------------- get_depth_z

def test_an_arr_realisation_gets_the_interpolated_pattern(tmp_path):
    curves = curves_with_weights(tmp_path)
    z = z_of(2000)                                  # top of the smoothing zone: the PMP pattern
    depths = curves.get_depth_z(z=z, duration=DURATION, storm_method='ARR areal',
                                spatial_method='GTSMR', print_msg=False)
    assert list(depths) == pytest.approx([60.0, 140.0])


def test_the_transition_is_the_same_whichever_temporal_method_was_drawn(tmp_path):
    curves = curves_with_weights(tmp_path)
    z = z_of(1000)
    arr = curves.get_depth_z(z=z, duration=DURATION, storm_method='ARR areal',
                             spatial_method='GTSMR', print_msg=False)
    extreme = curves.get_depth_z(z=z, duration=DURATION, storm_method='GTSMR',
                                 spatial_method='GTSMR', print_msg=False)
    assert list(arr) == pytest.approx(list(extreme))


def test_the_catchment_average_depth_is_unchanged_by_the_transition(tmp_path):
    # The weights are normalised on area, so smoothing moves rain between subcatchments
    # without changing what the catchment receives - mean_rain_mm, the losses and the
    # embedded burst filter all see the same number.
    curves = curves_with_weights(tmp_path)
    z = z_of(1500)
    smoothed = curves.get_depth_z(z=z, duration=DURATION, storm_method='ARR areal',
                                  spatial_method='GTSMR', print_msg=False)
    assert catchment_average(smoothed) == pytest.approx(catchment_average(GRIDDED))


def test_without_a_spatial_method_an_arr_realisation_keeps_the_gridded_pattern(tmp_path):
    # The pre-10 September 2026 behaviour, kept for callers that sample no spatial method
    # so that an old database rebuilds as the run that wrote it.
    curves = curves_with_weights(tmp_path)
    depths = curves.get_depth_z(z=z_of(1000), duration=DURATION, storm_method='ARR areal',
                                print_msg=False)
    assert list(depths) == pytest.approx(list(GRIDDED))


def test_below_the_lower_bound_there_is_no_transition(tmp_path):
    curves = curves_with_weights(tmp_path, bounds=(1000, 2000))
    depths = curves.get_depth_z(z=z_of(500), duration=DURATION, storm_method='ARR areal',
                                spatial_method='GTSMR', print_msg=False)
    assert list(depths) == pytest.approx(list(GRIDDED))


def test_the_filtering_lookup_still_gets_the_gridded_depths(tmp_path):
    # The embedded burst and pre-burst filters only want catchment average depths.
    curves = curves_with_weights(tmp_path)
    depths = curves.get_depth_z(z=z_of(2000), duration=DURATION, storm_method='ARR areal',
                                spatial_method='GTSMR', print_msg=False, for_filtering=True)
    assert list(depths) == pytest.approx(list(GRIDDED))


def test_an_arr_method_above_1_in_2000_no_longer_falls_over(tmp_path):
    # Reachable when aep_changeover_to_extreme runs past 1 in 2,000: the extreme branch has
    # no gridded fallback, and asking WeightInterpolatedSpatialPattern for an ARR pattern
    # raised. With the spatial method sampled in its own right there is always one to use.
    curves = curves_with_weights(tmp_path)
    depths = curves.get_depth_z(z=z_of(5000), duration=DURATION, storm_method='ARR areal',
                                spatial_method='GTSMR', print_msg=False)
    assert list(depths) == pytest.approx([180.0, 420.0])            # 300 mm on the PMP pattern


# -------------------------------------------------- StormBurst.sample_spatial_method

def storm_burst(tmp_path, smoothing_method, changeover=(9, 18)):
    config = tmp_path / 'storm_config.json'
    config.write_text(json.dumps({
        'file_paths': {},
        'storm_method_config': {'aep_changeover_to_extreme': [100, 2000],
                                'gsdm_gtsmr_changover_duration': list(changeover),
                                'extreme_pattern_for_ensemble': 'GTSMR'},
    }))
    storm = StormBurst(str(config), generate_storms=False)
    storm.rainfall = curves_with_weights(tmp_path)
    storm.rainfall.extreme_spatial_method = smoothing_method
    return storm


def test_the_spatial_method_follows_the_duration_for_an_arr_realisation(tmp_path):
    storm = storm_burst(tmp_path, 'interpolate_weights')
    assert storm.sample_spatial_method('ARR point', duration=6) == 'GSDM'
    assert storm.sample_spatial_method('ARR areal', duration=24) == 'GTSMR'


def test_the_spatial_method_is_sampled_in_the_duration_changeover_band(tmp_path):
    storm = storm_burst(tmp_path, 'interpolate_weights')
    drawn = {storm.sample_spatial_method('ARR areal', duration=12) for _ in range(50)}
    assert drawn == {'GSDM', 'GTSMR'}


def test_the_ensemble_takes_its_overlap_pattern_from_the_config(tmp_path):
    storm = storm_burst(tmp_path, 'interpolate_weights')
    for _ in range(10):
        assert storm.sample_spatial_method('ARR areal', duration=12,
                                           sim_method='ensemble') == 'GTSMR'


def test_an_extreme_realisation_heads_for_its_own_method(tmp_path):
    storm = storm_burst(tmp_path, 'interpolate_weights')
    assert storm.sample_spatial_method('GSDM', duration=12) == 'GSDM'
    assert storm.sample_spatial_method('GTSMR', duration=12) == 'GTSMR'


def test_the_legacy_smoothing_method_is_left_alone(tmp_path):
    # interpolate_depths applies no spatial pattern below 1 in 2,000, so there is nothing
    # to head for and the storm method is returned unchanged.
    storm = storm_burst(tmp_path, 'interpolate_depths')
    assert storm.sample_spatial_method('ARR areal', duration=24) == 'ARR areal'


def test_a_run_that_never_reaches_the_extreme_rainfall_samples_nothing(tmp_path):
    storm = storm_burst(tmp_path, 'interpolate_weights')
    storm.rainfall.smoothed_weights.lower_bound_pattern = None       # set_up_extreme_rain not called
    assert storm.sample_spatial_method('ARR areal', duration=24) == 'ARR areal'


def test_a_blank_spatial_method_is_not_a_method(tmp_path):
    # A Replicate file from a run that recorded none leaves NaN in the column.
    curves = curves_with_weights(tmp_path)
    depths = curves.get_depth_z(z=z_of(1000), duration=DURATION, storm_method='ARR areal',
                                spatial_method=float('nan'), print_msg=False)
    assert list(depths) == pytest.approx(list(GRIDDED))
