"""The dam's column of a regional inflow file, replaced by a dam failure model's outflow.

``DownstreamBoundaries.py`` is the last step of the downstream handover: URBS has written
``<run>_i.csv`` for each chosen event, and the dam's column in it is to come from a dam
failure model instead. What can go wrong quietly is the matching - the wrong event's
outflow under the right name - and a dam run that stops before the regional one does, so
those are what these pin, on a project small enough to read.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

BRYAN_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRYAN_ROOT))

import DownstreamBoundaries as boundaries                            # noqa: E402

DATABASE = 'CLD_mc_24h_E013_GWL1p3_RFSL__mcdf.csv'


def project(tmp_path, dam_hours=48.0, picked=42, scenarios=('',)):
    """A regional model with one event run, and a dam failure list holding that event."""
    model = tmp_path / 'regional' / 'model'
    model.mkdir(parents=True)
    config = tmp_path / 'regional' / 'urbs_config.json'
    config.write_text(json.dumps({'model_folder': 'model',
                                  'max_keys': {'outflow': 'DAM_OUT', 'level': 'DAM'}}))
    times = [0.0, 0.001] + [0.25 * i for i in range(1, 193)]            # to 48 h
    pd.DataFrame({'Time_hrs': times, 'TRIB': 5.0, 'DAM_OUT': 100.0}).to_csv(
        model / f'{picked:06d}_24h_GWL1p3_RFSL_i.csv', index=False)

    selection = tmp_path / 'RFSL_Design_GWL1p3_representative_events.json'
    selection.write_text(json.dumps({'settings': {}, 'targets': [dict(
        kind='aep', value=1000.0, result_type='level', rain_aep=None, source='',
        output_file='CLD_mc_24h', database=DATABASE, count=10, picked=picked, comment='')]}))

    dam = tmp_path / 'dam'
    (dam / 'results').mkdir(parents=True)
    rows = []
    for scenario in scenarios:
        for group, flow in (('RFSL_Design_GWL1p3', 250.0), ('FSL_Design_GWL1p3', 999.0)):
            name = f'RUN_{group}_{scenario or "gatesOK"}'
            steps = [0.25 * i for i in range(int(dam_hours / 0.25) + 1)]
            pd.DataFrame({'Time': steps, 'Total Outflow': [flow * t for t in steps]}).to_csv(
                dam / 'results' / f'{name}.csv', index=False)
            rows.append({'Simulation_name': name, 'Inflow_id': f'1_in_1_000_sim{picked:05d}',
                         'Output_name': f'results/{name}.csv', 'Gate_scenario': scenario,
                         'Bryan_group': group})
    pd.DataFrame(rows).to_excel(dam / 'Simulation_List.xlsx', index=False)
    return selection, config, dam / 'Simulation_List.xlsx', tmp_path / 'out'


def run(selection, config, sims, out, *extra):
    return boundaries.main(['--selection', str(selection), '--model', str(config),
                            '--dam-outflow', str(sims), '--out', str(out), '--suffix', '_RFSL',
                            *extra])


def test_the_dam_column_is_the_dam_models_and_the_rest_is_urbs(tmp_path):
    selection, config, sims, out = project(tmp_path)
    assert run(selection, config, sims, out) == 0
    written = pd.read_csv(out / '000042_24h_GWL1p3_RFSL_i.csv')
    assert list(written.columns) == ['Time_hrs', 'TRIB', 'DAM_OUT']
    assert (written['TRIB'] == 5.0).all()
    # the RFSL group's run, 250 m3/s an hour - not the FSL group's, which shares the realisation
    assert written['DAM_OUT'].iloc[-1] == pytest.approx(250.0 * 48.0)
    assert written.loc[written['Time_hrs'] == 12.0, 'DAM_OUT'].iloc[0] == pytest.approx(3000.0)
    assert len(written) == 194


def test_urbs_own_file_is_left_alone(tmp_path):
    selection, config, sims, out = project(tmp_path)
    run(selection, config, sims, out)
    original = pd.read_csv(config.parent / 'model' / '000042_24h_GWL1p3_RFSL_i.csv')
    assert (original['DAM_OUT'] == 100.0).all()


def test_a_dam_run_that_stops_early_is_refused(tmp_path):
    selection, config, sims, out = project(tmp_path, dam_hours=24.0)
    with pytest.raises(SystemExit, match='Extend the dam run'):
        run(selection, config, sims, out)


def test_the_gate_scenario_picks_the_run_and_names_the_file(tmp_path):
    selection, config, sims, out = project(tmp_path, scenarios=('', 'gate_fail_1'))
    run(selection, config, sims, out, '--gate-scenario', 'gate_fail_1')
    assert (out / '000042_24h_GWL1p3_RFSL_gate_fail_1_i.csv').is_file()
    index = pd.read_csv(out / 'RFSL_Design_GWL1p3_gate_fail_1_boundaries.csv')
    assert index['dam_simulation'].tolist() == ['RUN_RFSL_Design_GWL1p3_gate_fail_1']


def test_no_matching_dam_run_says_what_it_looked_for(tmp_path):
    selection, config, sims, out = project(tmp_path)
    with pytest.raises(SystemExit, match='0 simulations'):
        run(selection, config, sims, out, '--gate-scenario', 'not_a_scenario')


def test_a_location_not_marked_for_urbs_is_explained(tmp_path):
    selection, config, sims, out = project(tmp_path)
    config.write_text(json.dumps({'model_folder': 'model', 'max_keys': {'outflow': 'ELSEWHERE'}}))
    with pytest.raises(SystemExit, match='asterisk'):
        run(selection, config, sims, out)


def test_the_regional_run_has_to_exist(tmp_path):
    selection, config, sims, out = project(tmp_path)
    with pytest.raises(SystemExit, match='Run the regional model first'):
        boundaries.main(['--selection', str(selection), '--model', str(config),
                         '--dam-outflow', str(sims), '--out', str(out)])      # no suffix
