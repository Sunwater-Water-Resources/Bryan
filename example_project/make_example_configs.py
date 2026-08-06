"""
Build the config files, simulation lists and launchers for the Juniper Creek Dam example.

Kept separate from make_example_data.py so the two concerns stay apart: that script invents
rainfall and dam data, this one wires it into the config structure Bryan expects. Run it
after the data generator.

    python make_example_data.py
    python make_example_configs.py
"""
import json
import os
import stat

import pandas as pd

from make_example_data import (AEP_OF_PMP, AREA_BIN, AREAL_DURATIONS, DAM, DURATIONS,
                               FSL, HERE, REGION, SCHEME)


def write_json(relative_path, data):
    path = os.path.join(HERE, relative_path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    return path


def ifd_config():
    return {
        'col_suffix': '_AEP',
        'durations': [{'duration': d, 'filename': f'ifd/SYNTHETIC_ifd_{d:g}h.csv'}
                      for d in DURATIONS],
        'PMP_depths': 'ifd/SYNTHETIC_pmp_depths.csv',
        'PMP_scaling': 'ifd/SYNTHETIC_pmp_scaling.csv',
        'AEP_of_PMP': AEP_OF_PMP,
        'extreme_rainfall_interpolation_method': 'GEV',
        'extreme_spatial_smoothing_method': {'interpolate_weights': [1000, 2000]},
    }


def storm_config():
    return {
        'file_paths': {
            'ARR_datahub_file': 'SYNTHETIC_arr_datahub.txt',
            'rare_ifds': 'ifd_config.json',
            'point_patterns': 'patterns/SYNTHETIC_point_increments.csv',
            'areal_patterns': 'patterns/SYNTHETIC_areal_increments.csv',
            'gsdm_patterns': 'patterns/SYNTHETIC_gsdm.pat',
            'gtsmr_patterns': 'patterns/SYNTHETIC_gtsmr_~AREA~.pat',
            'preburst_patterns': 'patterns/SYNTHETIC_preburst_patterns.json',
        },
        'storm_method_config': {
            'aep_changeover_to_extreme': [100, 2000],
            'gsdm_gtsmr_changover_duration': [9, 18],
        },
        'CL_limit': {'apply': True, 'limit': 1.0},
    }


def mc_config():
    return {
        'scheme_config': {
            **SCHEME,
            'number_of_temporal_patterns': 10,
            'sample_method': 'normally distributed',
        },
        'aep_of_pmp': AEP_OF_PMP,
        'tpt_quantile_analysis': {
            'inflow': {'lower': 500, 'upper': 20000, 'step': 500},
            'level': {'lower': 209, 'upper': 218, 'step': 0.25},
            'outflow': {'lower': 500, 'upper': 20000, 'step': 500},
        },
    }


def enb_config():
    return {
        'storm_durations': AREAL_DURATIONS,
        'aep_list': [100, 1000, 2000, 10000, 100000, AEP_OF_PMP],
        'storm_method_config': {
            'interim_for_ensemble': 'both',
            'extreme_pattern_for_ensemble': 'GTSMR',
        },
    }


def lake_config():
    # Standard logistic (H = 1), so Vc is the exact ceiling - see LakeConfig.md
    return {
        'exceedance_layer_info': {
            'type': 'sigmoid',
            'coefficients': {'k': 2.0, 'Vf': 45000.0, 'H': 1, 'z0': 0.9, 'Vc': 190000.0},
        },
        'volume_cap': 'fsv',
    }


def climate_config():
    return {
        'NRM cluster': REGION,
        'KG classification': 'Cfa',
        'loss rates file': 'SYNTHETIC_rainfall_loss_rates_of_change.csv',
        'temporal pattern scaling': 'SYNTHETIC_temporal_pattern_scaling_factors.csv',
        'weighting files': {
            'ARR areal': 'Temporal_Pattern_Weighting/SYNTHETIC_areal.csv',
            'ARR point': 'Temporal_Pattern_Weighting/SYNTHETIC_point.csv',
            'GSDM': 'Temporal_Pattern_Weighting/SYNTHETIC_gsdm.csv',
            'GTSMR': 'Temporal_Pattern_Weighting/SYNTHETIC_gtsmr.csv',
        },
    }


def urbs_config():
    """Template only - needs a URBS licence, a .vec file and a ratings folder."""
    return {
        'model_type': 'URBS',
        'model_exe': 'C:\\URBS\\urbs.exe',
        'model_folder': 'urbs',
        'storms_folder': 'urbs/storms',
        'ratings_folder': 'urbs/ratings',
        'vec_file': 'juniper.vec',
        'result_prefix': 'JC',
        'store_tuflow': False,
        'time_increment': 0.25,
        'time_increment_override': {'2': 0.05, '3': 0.05, '6': 0.25, '9': 0.25},
        'alpha': 0.30,
        'beta': 1.5,
        'm_exponent': 0.8,
        'full_supply_volume': 190000.0,
        'simulation_periods': {str(d): max(2 * d, d + 48) for d in DURATIONS},
        'max_keys': {'inflow': 'JuniperDam', 'level': 'JuniperDam', 'outflow': 'JuniperDam'},
        'baseflow': {'apply': False},
    }


def rorb_config():
    """Template only - needs a RORB licence, a .catg and a .par file.

    RORB is the model type the storms-only demo uses: RorbModel does not check that the
    executable exists at construction, so the storm generation can be exercised on a machine
    without RORB installed. The models themselves are not run in that mode.
    """
    return {
        'model_type': 'RORB',
        'model_exe': 'C:\\RORB\\RORBWin.exe',
        'model_folder': 'rorb',
        'storms_folder': 'rorb/storms',
        'results_folder': 'rorb/results',
        'cat_file': 'juniper.catg',
        'par_file': 'juniper.par',
        'ADV_name': 'JuniperDam',
        'full_supply_volume': 190000.0,
        'simulation_periods': {str(d): max(2 * d, d + 48) for d in DURATIONS},
        'max_keys': {'inflow': 'JuniperDam', 'level': 'JuniperDam', 'outflow': 'JuniperDam'},
    }


STORMS_COLUMNS = [
    'Include', 'Method', 'Duration', 'Run models', 'Analyse results', 'Analyse sub-bursts',
    'Store hydrographs', 'Mop up files', 'ADV', 'Lake config', 'Focal subcatchments',
    'Config file', 'IL', 'CL', 'Replicates', 'Replicate file', 'Exclusions', 'GWL',
    'Log file', 'Output file', 'Comment',
]

ROUTING_COLUMNS = [
    'Include', 'Method', 'Run models', 'Analyse results', 'Store hydrographs',
    'Output suffix', 'Input MCDF', 'Inflow', 'ELS file', 'SQ file', 'FSL', 'ADV',
    'ADV source', 'Lake config', 'Config file', 'Log file', 'Output file',
    'Hydrographs folder', 'Results folder', 'Comment',
]


def storms_sims_list():
    common = {
        'Include': 'yes', 'Method': 'monte carlo', 'Run models': 'storms only',
        'Analyse results': 'no', 'Store hydrographs': 'no', 'Mop up files': 'no',
        'ADV': 'fsv', 'Lake config': 'sim_options/lake_conditions/lake_config.json',
        'Focal subcatchments':
            'sim_options/focal_locations/SYNTHETIC_dam_subcatchments.csv',
        'Config file': 'sim_options/mc_config.json',
        'IL': 30.0, 'CL': 2.0, 'Replicates': '', 'Replicate file': '',
    }
    rows = [
        {**common, 'Duration': 24, 'Analyse sub-bursts': 'yes', 'Exclusions': 'pb',
         'GWL': 0.0, 'Log file': 'results/storms_24h_gwl0_log.txt',
         'Output file': 'results/storms_24h_gwl0',
         'Comment': '24 h storms at present climate, no pre-burst'},
        {**common, 'Duration': 24, 'Analyse sub-bursts': 'yes', 'Exclusions': 'pb,ebf',
         'GWL': 0.0, 'Log file': 'results/storms_24h_gwl0_noebf_log.txt',
         'Output file': 'results/storms_24h_gwl0_noebf',
         'Comment': 'Same, with the embedded burst filter switched off for comparison'},
        {**common, 'Duration': 24, 'Analyse sub-bursts': 'no', 'Exclusions': 'pb',
         'GWL': 2.7, 'Log file': 'results/storms_24h_gwl27_log.txt',
         'Output file': 'results/storms_24h_gwl27',
         'Comment': '24 h storms at GWL 2.7 degC'},
        {**common, 'Include': 'no', 'Duration': 6, 'Analyse sub-bursts': 'yes',
         'Exclusions': 'pb', 'GWL': 0.0,
         'Log file': 'results/storms_6h_gwl0_log.txt',
         'Output file': 'results/storms_6h_gwl0',
         'Comment': 'A short duration, so GSDM patterns are sampled. Set Include to yes to run'},
    ]
    return pd.DataFrame(rows)[STORMS_COLUMNS]


def routing_sims_list():
    common = {
        'Include': 'yes', 'Method': 'reservoir routing', 'Run models': 'yes',
        'Analyse results': 'yes', 'Store hydrographs': 'no',
        'Input MCDF': 'inflows/SYNTHETIC_mc_24h__mcdf.csv',
        'Inflow': 'inflows/SYNTHETIC_mc_24h_inflows.csv',
        'ELS file': 'model/SYNTHETIC_juniper.els',
        'SQ file': 'model/SYNTHETIC_juniper.sq',
        'ADV': 'fsv', 'Lake config': 'sim_options/lake_conditions/lake_config.json',
        'Config file': 'sim_options/mc_config.json', 'Log file': '',
        'Hydrographs folder': 'results/routed_hydrographs',
        'Results folder': 'results',
    }
    rows = [
        {**common, 'Output suffix': '_existing', 'FSL': FSL, 'ADV source': 'mcdf',
         'Output file': 'routed_juniper',
         'Comment': 'Route the stored inflows through the existing dam'},
        {**common, 'Output suffix': '_lowered_fsl', 'FSL': FSL - 1.0, 'ADV source': 'mcdf',
         'Output file': 'routed_juniper',
         'Comment': 'The same inflows with the full supply level lowered by 1 m'},
        {**common, 'Output suffix': '_varying_adv', 'FSL': FSL, 'ADV source': 'lake_z',
         'Output file': 'routed_juniper',
         'Comment': 'Existing dam, antecedent volume resampled from lake_z via the lake config'},
    ]
    return pd.DataFrame(rows)[ROUTING_COLUMNS]


def write_launchers():
    """A .bat for Windows and a .sh companion, matching how the studies are driven."""
    for name, config in [('run_storms_only', 'sims_config_storms.json'),
                         ('run_routing', 'sims_config_routing.json')]:
        bat = os.path.join(HERE, f'{name}.bat')
        with open(bat, 'w', newline='\r\n') as f:
            f.write(f'''TITLE {DAM} example - {name}
cd /D "%~dp0"

:: Point these at your Bryan clone and its Python environment
set "BRYAN_ROOT=..\\.."
set "VENV_PY=%BRYAN_ROOT%\\env\\python.exe"

"%VENV_PY%" "%BRYAN_ROOT%\\Main.py" "{config}"

pause
''')
        sh = os.path.join(HERE, f'{name}.sh')
        with open(sh, 'w') as f:
            f.write(f'''#!/usr/bin/env bash
# {DAM} example - {name}. POSIX companion to {name}.bat.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

BRYAN_ROOT="${{BRYAN_ROOT:-..}}"
VENV_PY="${{VENV_PY:-$(command -v python3 || command -v python)}}"

exec "$VENV_PY" "$BRYAN_ROOT/Main.py" "{config}"
''')
        os.chmod(sh, os.stat(sh).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def main():
    write_json('storm_data/ifd_config.json', ifd_config())
    write_json('storm_data/storm_config.json', storm_config())
    write_json('sim_options/mc_config.json', mc_config())
    write_json('sim_options/enb_config.json', enb_config())
    write_json('sim_options/lake_conditions/lake_config.json', lake_config())
    write_json('climate_change/climate_config.json', climate_config())
    write_json('model/urbs_config.json', urbs_config())
    write_json('model/rorb_config.json', rorb_config())
    print('Config files written')

    write_json('sims_config_storms.json', {
        'simulation_list': 'SimsList_storms.xlsx',
        'test_runs': 0,
        'filepaths': {
            'model_config': 'model/rorb_config.json',
            'storm_config': 'storm_data/storm_config.json',
            'climate_config': 'climate_change/climate_config.json',
        },
    })
    write_json('sims_config_routing.json', {
        'simulation_list': 'SimsList_routing.xlsx',
        'test_runs': 0,
        'filepaths': {
            'model_config': 'model/rorb_config.json',
            'storm_config': 'storm_data/storm_config.json',
            'climate_config': 'climate_change/climate_config.json',
        },
    })

    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    storms_sims_list().to_excel(os.path.join(HERE, 'SimsList_storms.xlsx'), index=False)
    routing_sims_list().to_excel(os.path.join(HERE, 'SimsList_routing.xlsx'), index=False)
    print('Simulation lists written')

    write_launchers()
    print('Launchers written')


if __name__ == '__main__':
    main()
