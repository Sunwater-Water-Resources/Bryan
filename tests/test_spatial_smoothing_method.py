"""The extreme spatial smoothing method is selected, not assumed.

``ifdCurves`` offers two methods for smoothing the very rare - extreme spatial
patterns, chosen by the ``extreme_spatial_smoothing_method`` key in the IFD
files config:

- ``interpolate_depths``, the legacy method, which applies the GSDM/GTSMR
  spatial patterns to the PMP only and interpolates between 1 in 2,000 and PMP;
- ``{"interpolate_weights": [lower, upper]}``, which interpolates the spatial
  scaling across the changeover zone instead.

Omitting the key falls back to ``interpolate_depths`` with a printed CHECK.
**That fallback was the only working way to get the legacy method**: the
explicit branch compared ``self.extreme_spatial_method`` to the string it meant
to assign, and since nothing else in ``__init__`` sets the attribute, naming the
method you were already getting by default raised ``AttributeError``. It looked
like a typo that could not matter, because the default and the explicit setting
select the same method - the config that says out loud what it does was the one
that crashed.

Selecting the method matters beyond the crash: the Callide E010 IFD config
carries no key at all and so ran on ``interpolate_depths``, while the E009 run
supplying its inflows used ``interpolate_weights``. Nothing in the outputs says
which was used, so these tests pin all three ways of asking and the error for a
value that is neither.

This is a **Bryan** test, not a UI one - ``lib/RainfallScheme.py`` imports scipy,
which the launcher is forbidden to pull in. Run it with Bryan's own interpreter:

    python -m pytest tests -q

It skips where those dependencies are missing rather than failing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

pytest.importorskip("scipy")
pytest.importorskip("pandas")

from lib.RainfallScheme import ifdCurves                             # noqa: E402


def ifd_config(tmp_path, **keys):
    """A minimal IFD files config. Only the keys under test are read here."""
    path = tmp_path / 'ifd_files.json'
    path.write_text(json.dumps(keys))
    return str(path)


def test_missing_key_defaults_to_the_legacy_method(tmp_path):
    curves = ifdCurves(ifd_config(tmp_path))
    assert curves.extreme_spatial_method == 'interpolate_depths'


def test_the_legacy_method_can_be_named_explicitly(tmp_path):
    # The regression: this raised AttributeError, so the only way to run the
    # legacy method was to leave the config silent about it.
    curves = ifdCurves(ifd_config(
        tmp_path, extreme_spatial_smoothing_method='interpolate_depths'))
    assert curves.extreme_spatial_method == 'interpolate_depths'


def test_weights_method_carries_its_bounds(tmp_path):
    curves = ifdCurves(ifd_config(
        tmp_path,
        extreme_spatial_smoothing_method={'interpolate_weights': [100, 2000]},
        AEP_of_PMP=1900000))
    assert curves.extreme_spatial_method == 'interpolate_weights'
    assert curves.smoothed_weights.bounds_aep == [100, 2000]
    # The bounds are held as standard normal variates, ascending.
    assert curves.smoothed_weights.z_lower < curves.smoothed_weights.z_upper


@pytest.mark.parametrize('value', ['nonsense', 'interpolate weights', 42, None])
def test_an_unrecognised_method_says_so(tmp_path, value):
    # A bare string other than 'interpolate_depths' used to reach .keys() and
    # raise AttributeError, which reads as a bug in Bryan rather than as a
    # config to fix.
    with pytest.raises(Exception, match='extreme_spatial_smoothing_method is invalid'):
        ifdCurves(ifd_config(tmp_path, extreme_spatial_smoothing_method=value))
