"""The 'stm' replicate key covers the spatial method, and says so when it cannot.

'stm' replicates the ARR/extreme storm method sampling. Since 10 September 2026 the
extreme *spatial* pattern the interpolate_weights transition heads for is sampled
separately and recorded as 'spatial_method', and 'stm' replicates that too - but only
where the source database carries the column. An older database falls through to a
fresh draw, which inside the gsdm_gtsmr_changover_duration band is a coin toss, so
the run reproduces neither the source nor another replicate of itself.

That fallback is the quiet kind of wrong, so opening such a file has to say so.
"""
import re
import pathlib

import pandas as pd
import pytest

SIMULATOR = pathlib.Path(__file__).resolve().parents[1] / 'lib' / 'Simulator.py'
SOURCE = SIMULATOR.read_text()


def test_stm_replicates_the_spatial_method_when_the_column_is_there():
    assert re.search(
        r"if self\.replicates\['storm_method'\] and 'spatial_method' in "
        r"self\.replicates_df\.columns:\s*\n\s*spatial_method = "
        r"self\.replicates_df\.loc\[sim_id, 'spatial_method'\]",
        SOURCE), 'the stm key no longer reads spatial_method out of the replication file'


def test_the_fallback_still_resamples():
    """Deliberate: an old database records no spatial method, and reproducing the
    run that wrote it would mean reproducing the bug it was written with."""
    assert 'spatial_method = storm.sample_spatial_method(storm_method, duration)' in SOURCE


def test_a_replication_file_without_the_column_is_warned_about_by_name():
    block = SOURCE.split('self.replicates_df = pd.read_csv(replicate_file, index_col=0)')[1]
    block = block.split('# Set up the lake')[0]
    assert "self.replicates['storm_method'] and 'spatial_method' not in self.replicates_df.columns" in block, \
        'nothing checks for the missing column when the replication file is opened'
    assert 'WARNING' in block, 'the missing column does not warn'
    assert 'replicate_file' in block, 'the warning does not name the replication file'
    assert 'RE-SAMPLED' in block, 'the warning does not say what happens instead'


def test_the_warning_is_raised_once_not_per_realisation():
    """It belongs where the file is opened. In the simulation loop it would print
    once per realisation and be scrolled away by the run it is warning about."""
    setup, _, loop = SOURCE.partition('def run_models')
    assert 'no "spatial_method" column' not in loop
    assert 'no "spatial_method" column' in setup
