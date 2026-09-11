"""get_preburst_pattern draws from the patterns it actually found.

The selection takes the three pre-burst patterns whose total depth is nearest the
sampled proportion, then rescales the winner to that proportion. It is not always
three: ``rank()`` averages ties, so a set of patterns with equal totals - or totals
equal to floating point - can put every rank above 3 and return one pattern, or
none. The draw used to be ``randint(3)`` over that list, which raised IndexError.

This is not hypothetical. The example project's synthetic pre-burst patterns were
normalised to a common total, and 12 and 48 hour Monte Carlo rows with pre-burst
enabled died on it while 24 hour rows happened to survive.
"""
import numpy as np
import pandas as pd
import pytest

from lib.TemporalPatterns import PreburstPatterns


def patterns_with(totals, length=6):
    """A PreburstPatterns holding one duration whose patterns have the given totals."""
    index = list(range(-length, 0))
    frame = pd.DataFrame(
        {i: np.full(length, total / length) for i, total in enumerate(totals)},
        index=index)
    obj = PreburstPatterns.__new__(PreburstPatterns)
    obj.preburst_patterns = {24: frame}
    return obj


def draw(totals, proportion=0.2, timestep=1.0):
    return patterns_with(totals).get_preburst_pattern(24, proportion, timestep)


def test_a_normal_spread_picks_one_of_the_three_nearest():
    totals = [0.03, 0.10, 0.18, 0.22, 0.30, 0.40]
    chosen = {draw(totals)[1] for _ in range(60)}
    assert chosen <= {2, 3, 4}, 'drew a pattern outside the three nearest to 0.2'
    assert len(chosen) == 3, 'the draw is not reaching all three'


def test_every_total_equal_does_not_raise():
    """The degenerate case: nothing to choose between the patterns."""
    pattern, sample_int = draw([0.2] * 6)
    assert sample_int in range(6)
    assert not pattern.empty


def test_two_candidates_does_not_raise():
    """A tie structure that leaves fewer than three under rank() <= 3."""
    pattern, sample_int = draw([0.2, 0.2, 0.2, 0.2, 0.9, 0.95])
    assert sample_int in range(6)
    assert not pattern.empty


def test_the_winner_is_rescaled_to_the_sampled_proportion():
    """Selection is on the nearest total; the depth comes from the proportion. The
    pattern is a percentage of the burst depth, so it sums to proportion * 100 -
    within the 3 dp each increment is rounded to on the way out.
    """
    for proportion in [0.05, 0.2, 0.33]:
        pattern, _ = draw([0.03, 0.10, 0.18, 0.22, 0.30, 0.40], proportion=proportion)
        assert pattern.sum() == pytest.approx(proportion * 100, abs=0.5e-3 * len(pattern))


def test_the_usual_three_take_one_number_from_the_random_stream():
    """The fix must not move the random stream where it used to work: with three
    candidates the draw is the same draw it always was."""
    totals = [0.03, 0.10, 0.18, 0.22, 0.30, 0.40]
    np.random.seed(7)
    mine = [draw(totals)[1] for _ in range(20)]
    np.random.seed(7)
    theirs = [[2, 3, 4][np.random.randint(3)] for _ in range(20)]
    assert mine == theirs
