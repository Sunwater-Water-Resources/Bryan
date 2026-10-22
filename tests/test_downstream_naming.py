"""The launcher and the generator must name a storm file identically.

They cannot share the code. ``DownstreamStorms.py`` imports the depth assembly
and so pulls in scipy; ``ui/core/downstream.py`` runs in an environment that
deliberately has none, and Bryan may not import from ``ui`` in the other
direction. So the naming is written twice, and this is what stops the two
drifting: the launcher shows the user a filename before the run, and the
generator writes one after it. If they disagree, the page reports a file that
never appears and the run leaves one nobody was told about.

Bryan's own interpreter has both sides importable, which is why this test lives
here rather than in ``ui/tests``.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))
sys.path.insert(0, str(BRYAN_ROOT / 'ui'))

pytest.importorskip('scipy')
pytest.importorskip('pandas')

import DownstreamStorms as generator                                 # noqa: E402
from core import downstream as launcher                              # noqa: E402


@pytest.mark.parametrize('duration, gwl', [
    (24, 1.3), (6, 0.0), (120, 2.7), (4.5, 0.3), (96, 1.0),
])
@pytest.mark.parametrize('method', ['monte carlo', 'ensemble'])
@pytest.mark.parametrize('suffix', ['', '_RFSL'])
def test_the_two_sides_name_a_storm_file_the_same(duration, gwl, method, suffix):
    target = SimpleNamespace(picked=42)
    assert generator.storm_name(target, duration, gwl, suffix, method) ==         launcher.storm_filename(42, duration, gwl, suffix, method)


def test_an_ensemble_event_can_never_take_a_realisations_name():
    target = SimpleNamespace(picked=42)
    assert generator.storm_name(target, 9, 1.3, '_RFSL') == '000042_9h_GWL1p3_RFSL.9'
    assert generator.storm_name(target, 9, 1.3, '_RFSL', 'ensemble') == 'e00042_9h_GWL1p3_RFSL.9'


def test_the_launcher_reads_a_database_name_the_way_the_generator_does():
    name = 'CLD_mc_4p5h_E009_ebf_no-pbp_L105-1.7_GWL0p3_RGN__mcdf.parquet'
    assert generator.from_name(generator.DURATION, name, None, 'duration') == \
        launcher.read_from_name(launcher.DURATION, name)
    assert generator.from_name(generator.GWL, name, None, 'warming level') == \
        launcher.read_from_name(launcher.GWL, name)


def test_a_name_without_a_duration_stops_the_generator_and_flags_the_launcher():
    name = 'results__mcdf.parquet'
    assert launcher.read_from_name(launcher.DURATION, name) is None
    with pytest.raises(SystemExit):
        generator.from_name(generator.DURATION, name, None, 'duration')


def test_a_part_hour_storm_has_one_dot_so_its_run_name_is_whole():
    # '.4.5' left the run named 'e00022_4p5h_GWL2p7_RFSL.4', and URBS wrote no results for it
    import os
    target = SimpleNamespace(picked=22)
    name = generator.storm_name(target, 4.5, 2.7, '_RFSL', 'ensemble')
    assert name == 'e00022_4p5h_GWL2p7_RFSL.4p5'
    assert os.path.splitext(name)[0] == 'e00022_4p5h_GWL2p7_RFSL'
    # whole hours keep the extension every run so far was named with
    assert generator.storm_name(target, 9, 1.3, '_RFSL', 'ensemble') == 'e00022_9h_GWL1p3_RFSL.9'
