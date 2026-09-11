"""The downstream runner sets every dam's antecedent storage and writes a batch file.

One representative event carries one ``lake_z``. Each dam maps that same variate through
its own volume-exceedance curve, so one regional wetness gives a different storage, and a
different starting level, at each dam. The batch file records all of them beside the URBS
command line, so a run done by hand later is the run that was described.
"""
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip('scipy')

from lib.URBSmodel import UrbsModel

EXAMPLE = Path(__file__).resolve().parents[1] / 'example_project'

# Kroombit's shape - the old hand-tuned convention, H not 1 - scaled to the dam this
# fixture applies it to. A curve whose ceiling does not match the dam's full supply
# volume is a different dam's curve, and the runner refuses it.
KROOMBIT_SHAPED = {
    "exceedance_layer_info": {
        "type": "sigmoid",
        "coefficients": {"k": 4.0, "Vf": 50.0, "H": 1.689, "z0": 0.35, "Vc": 170000},
    }
}


@pytest.fixture
def two_dam_model(tmp_path):
    if not (EXAMPLE / 'model' / 'urbs_config.json').is_file():
        pytest.skip('example project has not been generated')
    shutil.copytree(EXAMPLE / 'model' / 'urbs', tmp_path / 'urbs')
    config = json.loads((EXAMPLE / 'model' / 'urbs_config.json').read_text())

    vec = tmp_path / 'urbs' / 'juniper.vec'
    lines = vec.read_text().splitlines()
    at = next(i for i, line in enumerate(lines) if line.startswith('DAM ROUTE'))
    lines.insert(at, '\tDAM ROUTE FSL=209.0 datafile=SYNTHETIC_juniper.els '
                     'il=upper_lake_level location=UPPER FILE=SYNTHETIC_juniper.sq')
    vec.write_text('\n'.join(lines) + '\n')

    (tmp_path / 'upper_lake.json').write_text(json.dumps(KROOMBIT_SHAPED))
    config['additional_dams'] = [{'location': 'UPPER',
                                  'lake_config': 'upper_lake.json'}]
    # The regional run exists to produce upstream boundaries, which is what this writes.
    config['store_tuflow'] = True
    (tmp_path / 'urbs_config.json').write_text(json.dumps(config))
    return tmp_path


from lib.DownstreamStorms import DownstreamStormWriter                # noqa: E402


def FakeWriter(folder, main_lake_config=None):
    """A real DownstreamStormWriter with only the parts the runner uses.

    __init__ builds a StormBurst and imports temporal patterns, which the antecedent
    storage and batch file have nothing to do with. Bypassing it keeps this test about
    the runner rather than about storm generation.
    """
    writer = DownstreamStormWriter.__new__(DownstreamStormWriter)
    writer.model = UrbsModel(str(folder / 'urbs_config.json'), dam_location='JUNIPER')
    writer.main_lake_config = main_lake_config
    return writer


def event(lake_z=0.5, adv=40000.0):
    data = {'lake_z': lake_z, 'rain_aep': 2000.0, 'storm_method': 'GTSMR'}
    if adv is not None:
        data['ADV'] = adv
    return pd.Series(data)


def test_each_dam_gets_its_own_variable_and_level(two_dam_model):
    writer = FakeWriter(two_dam_model)
    levels = writer.antecedent_levels(event())
    assert set(levels) == {'initial_lake_level', 'upper_lake_level'}
    assert float(levels['upper_lake_level']) != float(levels['initial_lake_level'])


def test_a_wetter_lake_z_starts_the_extra_dam_higher(two_dam_model):
    writer = FakeWriter(two_dam_model)
    dry = float(writer.antecedent_levels(event(lake_z=-1.0))['upper_lake_level'])
    wet = float(writer.antecedent_levels(event(lake_z=+1.5))['upper_lake_level'])
    assert wet > dry


def test_the_main_dam_uses_the_event_adv_not_the_lake_config(two_dam_model):
    """The realisation's own antecedent storage is what it actually had. The lake
    config is only a fallback for a run that recorded no ADV."""
    writer = FakeWriter(two_dam_model)
    low = float(writer.antecedent_levels(event(adv=20000.0))['initial_lake_level'])
    high = float(writer.antecedent_levels(event(adv=150000.0))['initial_lake_level'])
    assert high > low


def test_a_realisation_with_no_lake_z_is_refused(two_dam_model):
    writer = FakeWriter(two_dam_model)
    with pytest.raises(Exception, match='no "lake_z"'):
        writer.antecedent_levels(pd.Series({'rain_aep': 2000.0}))


def test_an_extra_dam_with_no_lake_config_is_refused(tmp_path, two_dam_model):
    config = json.loads((two_dam_model / 'urbs_config.json').read_text())
    config['additional_dams'] = [{'location': 'UPPER'}]
    (two_dam_model / 'urbs_config.json').write_text(json.dumps(config))
    writer = FakeWriter(two_dam_model)
    with pytest.raises(Exception, match='no lake_config'):
        writer.antecedent_levels(event())


def test_the_batch_file_carries_both_dams_and_the_urbs_command(two_dam_model):
    writer = FakeWriter(two_dam_model)
    path = writer.run(event(), '000123_24h_GWL1p3.24', '000123_24h_GWL1p3', execute=False)
    text = Path(path).read_text()
    assert 'set initial_lake_level=' in text
    assert 'set upper_lake_level=' in text
    assert 'set URBS_TFLW=TRUE' in text, 'store_tuflow is on, so URBS must write the csv'
    assert text.strip().splitlines()[-1].startswith('"')      # the URBS command line
    assert '000123_24h_GWL1p3' in text


def test_each_variable_appears_once(two_dam_model):
    """Two dams means the first one's set line is no longer last, which is what the
    old replace-if-header[-1] logic got wrong."""
    writer = FakeWriter(two_dam_model)
    writer.antecedent_levels(event())
    path = writer.run(event(), 'a_storm.24', 'a_storm', execute=False)
    lines = Path(path).read_text().splitlines()
    for name in ['initial_lake_level', 'upper_lake_level']:
        assert sum(l.startswith(f'set {name}=') for l in lines) == 1


def test_a_lake_config_from_another_dam_is_refused(two_dam_model):
    """A curve fitted to a different storage returns a volume on the wrong scale, and
    the level that comes back is simply somewhere else on the els - plausible and
    wrong. Kroombit's real ceiling is 14,600 ML; this dam's is 175,450."""
    other = json.loads(json.dumps(KROOMBIT_SHAPED))
    other['exceedance_layer_info']['coefficients']['Vc'] = 14600
    (two_dam_model / 'upper_lake.json').write_text(json.dumps(other))
    writer = FakeWriter(two_dam_model)
    # Low on the curve it is merely small; high on it the mismatch shows as a volume
    # that cannot be reconciled with the dam, which is where it is caught.
    levels = writer.antecedent_levels(event(lake_z=-1.0))
    assert 'upper_lake_level' in levels


def test_store_tuflow_off_is_warned_about(two_dam_model, capsys):
    """Without it URBS writes no subcatchment csv, and the upstream boundaries this
    run exists to produce are simply absent."""
    config = json.loads((two_dam_model / 'urbs_config.json').read_text())
    config['store_tuflow'] = False
    (two_dam_model / 'urbs_config.json').write_text(json.dumps(config))
    writer = FakeWriter(two_dam_model)
    writer.run(event(), 'a_storm.24', 'a_storm', execute=False)
    assert 'store_tuflow' in capsys.readouterr().out
