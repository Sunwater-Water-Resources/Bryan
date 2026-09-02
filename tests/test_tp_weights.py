"""Temporal pattern weights: applied in the sampling, or in the analysis, never both.

There are two ways to give some temporal patterns more probability than others,
and they are alternatives:

- the sims-list ``TP weights`` column, which composes the weights into
  ``get_temporal_pattern_sample`` so the realisations are *drawn* in proportion
  to them. The TPT then counts them plainly, because a weighted sample already
  is the weighted probability. The run records what it used as ``tp_weight``.
- ``util/CalibrateTpWeights.py``, which attaches ``tp_w`` to a **uniformly**
  sampled mcdf and weights at analysis time, so a weight set can be tried
  without re-running the simulation.

Doing both weights the patterns twice, and nothing about the output would look
wrong - the curves stay smooth and plausible. Hence the guards, and hence these
tests, which also pin that the ordinary single-mechanism paths are untouched.

This is a **Bryan** test, not a UI one - it imports lib/MCScheme.py directly,
which the launcher is forbidden to do (it pulls in scipy and matplotlib). Run
it with Bryan's own interpreter:

    python -m pytest tests -q

It skips where those dependencies are missing rather than failing.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
pytest.importorskip("matplotlib")

from lib.MCScheme import SampleScheme                                # noqa: E402


def scheme():
    return SampleScheme(number_of_main_divisions=2, number_of_sub_divisions=3,
                        number_of_temporal_patterns=10, output_folder=".",
                        upper_aep=100, lower_aep=2)


def mcdf(**columns):
    frame = pd.DataFrame({"m": [0, 0, 0, 1, 1, 1],
                          "inflow": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})
    for name, values in columns.items():
        frame[name] = values
    return frame


def quantiles_of(frame):
    run = scheme()
    run.df = frame
    run.compute_std_quantiles("inflow")
    return run.quantiles["inflow"]["inflow"].tolist()


def calibrate_module():
    path = BRYAN_ROOT / "util" / "CalibrateTpWeights.py"
    spec = importlib.util.spec_from_file_location("CalibrateTpWeights", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# -- the guard ---------------------------------------------------------------

def test_both_columns_is_refused_rather_than_weighted_twice():
    with pytest.raises(Exception, match="already sampled|weight them again"):
        quantiles_of(mcdf(tp_w=1.0, tp_weight=0.5))


def test_the_message_names_both_mechanisms():
    """It has to say which one to drop, or it is just an obstacle."""
    with pytest.raises(Exception) as error:
        quantiles_of(mcdf(tp_w=1.0, tp_weight=0.5))
    message = str(error.value)
    assert "TP weights" in message
    assert "CalibrateTpWeights" in message


def test_calibrating_on_a_weighted_sample_is_refused():
    module = calibrate_module()
    frame = mcdf(tp_weight=0.5, tp=[0, 1, 2, 3, 4, 5], group="ARR point|rare")
    with pytest.raises(Exception, match="already sampled"):
        module.attach_weights(frame, {"ARR point|rare": [1.0] * 10})


# -- neither mechanism is disturbed -------------------------------------------

def test_a_weighted_sampling_run_analyses_exactly_as_an_unweighted_one():
    """'tp_weight' is provenance. The sample already carries the weighting."""
    plain = quantiles_of(mcdf())
    weighted_sample = quantiles_of(mcdf(tp_weight=0.5))
    assert weighted_sample == plain


def test_analysis_time_weights_still_change_the_answer():
    """The calibration path: 'tp_w' on a uniform sample must still be applied."""
    even = quantiles_of(mcdf(tp_w=1.0))
    lopsided = quantiles_of(mcdf(tp_w=[1.0, 1.0, 9.0, 1.0, 1.0, 9.0]))
    assert even != lopsided
