"""TotalProbTheorem.assign_aep_all must agree with assign_aep, value for value.

The vectorised path is what the Monte Carlo analysis runs; the scalar one is
kept for util/CalibrateTpWeights.py, which calls it with a single IFD depth and
swaps the mcdf underneath it. They are two implementations of one piece of
arithmetic, so they need pinning together.
"""
import numpy as np
import pandas as pd
import pytest

from lib.MCScheme import TotalProbTheorem


def divisions(m, lower_aep=2, upper_aep=100000):
    from scipy.special import ndtri
    return np.linspace(ndtri(1 - 1 / lower_aep), ndtri(1 - 1 / upper_aep), m + 1)


def sample(m=6, n=200, seed=0, nan_frac=0.0, weighted=False):
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame({
        'm': np.repeat(np.arange(m), n),
        'inflow': rng.lognormal(mean=6.0, sigma=0.8, size=m * n),
    })
    # Make the peaks grow with the division, as a real sample does.
    frame['inflow'] *= 1.0 + 0.35 * frame['m']
    if nan_frac:
        spoil = rng.random(len(frame)) < nan_frac
        frame.loc[spoil, 'inflow'] = np.nan
    if weighted:
        frame['tp_w'] = rng.choice([0.05, 0.1, 0.3], size=len(frame))
    return frame


def compare(frame, m=6, n=200):
    tpt = TotalProbTheorem(m, n, divisions(m), frame)
    values = frame['inflow'].to_numpy(dtype=float)
    fast = tpt.assign_aep_all(values, 'inflow')
    slow = np.array([tpt.assign_aep(v, 'inflow') for v in values])
    return fast, slow


def test_unweighted_matches_the_scalar_path():
    fast, slow = compare(sample())
    np.testing.assert_allclose(fast, slow, rtol=1e-12, atol=0)


def test_weighted_matches_the_scalar_path():
    """The 'tp_w' branch: weights applied in the analysis, not the sampling."""
    fast, slow = compare(sample(weighted=True))
    np.testing.assert_allclose(fast, slow, rtol=1e-12, atol=0)


@pytest.mark.parametrize('weighted', [False, True])
def test_a_failed_realisation_never_counts_as_an_exceedance(weighted):
    """NaN peaks must not exceed anything.

    np.sort puts NaN last, so a searchsorted over the raw sorted division counts
    every NaN as above every value. That is the bug this pins: it inflates pH,
    and the effect is largest in the top division, which is where the rare tail
    of the frequency curve comes from.
    """
    frame = sample(nan_frac=0.15, weighted=weighted)
    assert frame['inflow'].isna().any()
    fast, slow = compare(frame)
    np.testing.assert_allclose(fast, slow, rtol=1e-12, atol=0)


def test_a_nan_value_gets_the_same_answer_as_the_scalar_path():
    frame = sample()
    tpt = TotalProbTheorem(6, 200, divisions(6), frame)
    fast = tpt.assign_aep_all(np.array([np.nan]), 'inflow')
    slow = tpt.assign_aep(np.nan, 'inflow')
    np.testing.assert_allclose(fast, [slow], rtol=1e-12, atol=0)


def test_ties_are_counted_as_not_exceeding():
    """assign_aep counts strictly greater, so a value equal to a peak is not an
    exceedance. searchsorted has to be side='right' to agree."""
    frame = pd.DataFrame({'m': [0, 0, 1, 1], 'inflow': [10.0, 10.0, 10.0, 20.0]})
    tpt = TotalProbTheorem(2, 2, divisions(2), frame)
    fast = tpt.assign_aep_all(np.array([10.0]), 'inflow')
    slow = tpt.assign_aep(10.0, 'inflow')
    np.testing.assert_allclose(fast, [slow], rtol=1e-12, atol=0)


def test_an_empty_division_is_reported_rather_than_dropped():
    """The one place the vectorised path deliberately does not copy the scalar one.

    A division with no realisations makes the scalar path's groupby miss that row
    of compute_df, leaving NaN - which pandas' skipna then drops straight back out
    of the sum, so the division's share of the probability space silently goes
    missing and every AEP comes out low but entirely plausible. A stratified
    sample cannot have an empty division, so this only ever means the mcdf is not
    the sample the scheme config describes.
    """
    frame = pd.DataFrame({'m': [0, 0, 2, 2], 'inflow': [1.0, 2.0, 3.0, 4.0]})
    tpt = TotalProbTheorem(3, 2, divisions(3), frame)
    with pytest.raises(Exception, match='no realisations'):
        tpt.assign_aep_all(np.array([1.5]), 'inflow')
    # What it does today: the middle division's 0.42 of the sample space vanishes.
    assert tpt.assign_aep(1.5, 'inflow') == pytest.approx(0.3902260193782334)
