"""A URBS model with more than one dam sets each one's starting level separately.

Every dam's DAM ROUTE line carries its full supply level, its storage curve and the
name of the environment variable URBS reads the starting level from. That last one is
per dam, and the regional Callide model shipped with ``il=initial_lake_level`` on both
lines - so one ``set`` statement drove both, and Kroombit started every event at
Callide's level. Nothing in the output says so: the model runs, and the second dam is
simply wrong.
"""
import json
import shutil
from pathlib import Path

import pytest

pytest.importorskip('scipy')

from lib.URBSmodel import UrbsModel

EXAMPLE = Path(__file__).resolve().parents[1] / 'example_project'


@pytest.fixture
def two_dams(tmp_path):
    """The example model with a second dam spliced into the vec."""
    if not (EXAMPLE / 'model' / 'urbs_config.json').is_file():
        pytest.skip('example project has not been generated')
    shutil.copytree(EXAMPLE / 'model' / 'urbs', tmp_path / 'urbs')
    config = json.loads((EXAMPLE / 'model' / 'urbs_config.json').read_text())

    vec_path = tmp_path / 'urbs' / 'juniper.vec'
    lines = vec_path.read_text().splitlines()
    at = next(i for i, line in enumerate(lines) if line.startswith('DAM ROUTE'))
    # An upstream dam on its own storage curve, with its own il= variable.
    lines.insert(at, 'DAM ROUTE FSL=209.0 datafile=SYNTHETIC_juniper.els '
                     'il=upper_lake_level location=UPPER FILE=SYNTHETIC_juniper.sq')
    vec_path.write_text('\n'.join(lines) + '\n')

    config['additional_dams'] = [{'location': 'UPPER'}]
    (tmp_path / 'urbs_config.json').write_text(json.dumps(config))
    return tmp_path


def model_for(folder, **kwargs):
    return UrbsModel(str(folder / 'urbs_config.json'), dam_location='JUNIPER', **kwargs)


def test_both_dams_are_found_with_their_own_variables(two_dams):
    model = model_for(two_dams)
    assert model.lake_level_variable == 'initial_lake_level'
    assert [d['variable'] for d in model.extra_dams] == ['upper_lake_level']
    assert model.extra_dams[0]['full_supply_level'] == 209.0


def test_each_dam_writes_its_own_set_statement(two_dams):
    model = model_for(two_dams)
    model.set_volume_below_fsl(5000.0)
    model.set_additional_dam_volume('UPPER', 2000.0)
    sets = [line for line in model.header if line.startswith('set ')]
    assert any(line.startswith('set initial_lake_level=') for line in sets)
    assert any(line.startswith('set upper_lake_level=') for line in sets)


def test_resetting_one_dam_replaces_its_line_rather_than_appending(two_dams):
    """The old code only inspected header[-1]. With two dams the first dam's line is no
    longer last, so re-setting it appended a duplicate and cmd.exe took the later one."""
    model = model_for(two_dams)
    model.set_volume_below_fsl(5000.0)
    model.set_additional_dam_volume('UPPER', 2000.0)
    model.set_volume_below_fsl(20000.0)          # the main dam again, now not last
    main = [line for line in model.header if line.startswith('set initial_lake_level=')]
    assert len(main) == 1, 'the main dam has two conflicting set statements'
    assert len([l for l in model.header if l.startswith('set upper_lake_level=')]) == 1


def test_a_lower_starting_volume_gives_a_lower_level(two_dams):
    model = model_for(two_dams)
    full = model.set_additional_dam_volume('UPPER', 0.0)
    drawn = model.set_additional_dam_volume('UPPER', 20000.0)
    assert drawn < full


def test_two_dams_sharing_one_variable_is_refused(two_dams):
    """The real fault. Both Callide dams shipped with il=initial_lake_level."""
    vec_path = two_dams / 'urbs' / 'juniper.vec'
    vec_path.write_text(vec_path.read_text().replace('il=upper_lake_level',
                                                     'il=initial_lake_level'))
    with pytest.raises(Exception, match='same variable'):
        model_for(two_dams)


def test_a_dam_with_no_variable_is_refused(two_dams):
    vec_path = two_dams / 'urbs' / 'juniper.vec'
    vec_path.write_text(vec_path.read_text().replace('il=upper_lake_level', 'il=208.5'))
    with pytest.raises(Exception, match='no il=<variable>'):
        model_for(two_dams)


def test_a_missing_dam_is_named(two_dams):
    config = json.loads((two_dams / 'urbs_config.json').read_text())
    config['additional_dams'] = [{'location': 'NOT_A_DAM'}]
    (two_dams / 'urbs_config.json').write_text(json.dumps(config))
    with pytest.raises(Exception, match='NOT_A_DAM'):
        model_for(two_dams)


def test_an_indented_dam_route_line_is_found(two_dams):
    """A dam nested inside a STORE./GET. block is indented, and a startswith on the raw
    line cannot see it. The regional Callide model routes Kroombit that way, so its
    DAM ROUTE line was invisible and the dam reported as missing from the vec."""
    vec_path = two_dams / 'urbs' / 'juniper.vec'
    lines = vec_path.read_text().splitlines()
    at = next(i for i, l in enumerate(lines) if 'location=UPPER' in l)
    lines[at] = '\t' + lines[at]
    vec_path.write_text('\n'.join(lines) + '\n')
    model = model_for(two_dams)
    assert [d['location'] for d in model.extra_dams] == ['UPPER']
    assert model.extra_dams[0]['variable'] == 'upper_lake_level'


def test_a_model_with_one_dam_is_unchanged(tmp_path):
    """No additional_dams key: everything behaves as it always has."""
    if not (EXAMPLE / 'model' / 'urbs_config.json').is_file():
        pytest.skip('example project has not been generated')
    shutil.copytree(EXAMPLE / 'model' / 'urbs', tmp_path / 'urbs')
    shutil.copy(EXAMPLE / 'model' / 'urbs_config.json', tmp_path / 'urbs_config.json')
    model = model_for(tmp_path)
    assert model.extra_dams == []
    assert model.lake_level_variable == 'initial_lake_level'
    model.set_volume_below_fsl(5000.0)
    assert model.header[-1].startswith('set initial_lake_level=')
