"""The pre-burst extension must terminate.

``prepend_excess_preburst`` walks a storm backwards, timestep by timestep, until
the enveloping-burst curve reaches the depth the filtering needs. That loop had no
bound, so a curve that never gets there hung the run outright - no error, no
traceback, and a log whose last line is mid-sentence somewhere else entirely. It
took a 30 minute ``timeout`` to notice.

The curve that produced it came from one line above. ``filter_temppat`` picks the
duration with the shallowest slope *from* the main burst, but sliced the filter
curve with ``.loc[main_burst_duration:]``, which includes the main burst duration
itself - and the slope of a point to itself divides by ``log(d) - log(d) == 0``.
That is ``-inf`` whenever the point sits below 100, so ``idxmin`` chose it every
time, the two-point interpolation became one point repeated, and extrapolating it
gives ``exp(-inf) == 0`` at every duration. Zero is never greater than the target.
"""
import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from lib.StormGenerator import StormBurst


def burst(timesteps=1.0):
    obj = StormBurst.__new__(StormBurst)
    obj.timesteps = timesteps
    return obj


def curve_from(points):
    index = np.log(list(points))
    values = np.log(list(points.values()))
    return scipy.interpolate.interp1d(index, values, fill_value='extrapolate')


def test_a_rising_curve_still_terminates_and_prepends():
    """The normal case: the enveloping burst grows with duration and is reached."""
    pattern = pd.Series({-3.0: 2.0, -2.0: 3.0, -1.0: 4.0})
    curve = curve_from({48: 100.0, 96: 130.0, 240: 200.0})
    extended, extension = burst().prepend_excess_preburst(pattern, 48.0, 12.0, curve)
    assert extension > 0
    assert len(extended) > len(pattern)
    assert extended.tail(len(pattern)).equals(pattern)


def test_a_degenerate_curve_is_reported_rather_than_spun_on():
    """The hang: a curve that returns 0 at every duration never reaches the target."""
    zero_everywhere = lambda _: np.array(-np.inf)
    pattern = pd.Series({-3.0: 2.0, -2.0: 3.0, -1.0: 4.0})
    with pytest.raises(Exception, match='not converging'):
        burst().prepend_excess_preburst(pattern, 48.0, 12.0, zero_everywhere)


def test_a_flat_curve_is_reported_rather_than_spun_on():
    pattern = pd.Series({-3.0: 2.0, -2.0: 3.0, -1.0: 4.0})
    flat = curve_from({48: 100.0, 4800: 100.0})
    with pytest.raises(Exception, match='not converging'):
        burst().prepend_excess_preburst(pattern, 48.0, 50.0, flat)


def test_the_slope_search_cannot_select_the_main_burst_duration():
    """What made the degenerate curve. The slope of the main burst duration to
    itself is -inf when it sits below 100, which wins every idxmin."""
    main = 48.0
    tabulated = pd.Series({48.0: 96.0, 72.0: 103.0, 96.0: 108.0})
    whole = tabulated.loc[main:]
    with np.errstate(divide='ignore', invalid='ignore'):
        including = ((np.log(whole) - np.log(100))
                     / (np.log(whole.index.to_numpy()) - np.log(main)))
    assert including.idxmin() == main, 'the -inf that used to be selected'

    longer = whole.loc[whole.index > main]
    excluding = ((np.log(longer) - np.log(100))
                 / (np.log(longer.index.to_numpy()) - np.log(main)))
    assert excluding.idxmin() == 72.0
    assert np.isfinite(excluding).all()


def test_the_source_slices_strictly_longer_than_the_main_burst():
    import pathlib
    source = (pathlib.Path(__file__).resolve().parents[1] / 'lib' / 'StormGenerator.py').read_text()
    assert 'longer = curve.loc[curve.index > main_burst_duration]' in source, \
        'the slope search is back to including the main burst duration itself'
