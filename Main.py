"""
Script used to run URBS Simulations either as Monte Carlo or Ensemble events
Generally... units are: '1 in X' for AEP, hours for storm duration, and km² for area.
"""

import sys
import os
import json
import pandas as pd
import numpy as np
from lib.Simulator import MonteCarloSimulator, EnsembleSimulator
from datetime import datetime


# Get the config information
simulation_file = sys.argv[1]
print('Opening simulation configuration file: ', simulation_file)
f = open(simulation_file)
sim_config = json.load(f)
f.close()
print(sim_config)

# Get the simulation information
project_folder = os.path.dirname(simulation_file)
if 'project_folder' in sim_config:
    if not sim_config['project_folder'].lower() == 'default':
        project_folder = sim_config['project_folder']

print('\nProject folder is: ', os.path.abspath(project_folder))
sim_list_file = os.path.join(project_folder, sim_config['simulation_list'])
print('\nOpening simulation list: ', sim_list_file)
sim_df = pd.read_excel(sim_list_file, sheet_name=0)
print(sim_df[sim_df['Include'] == 'yes'])

# Set up a dataframe for logging the simulations


# Get config filepaths and convert paths from relative
filepaths = sim_config['filepaths']
folder = os.path.dirname(simulation_file)
if not os.path.isabs(folder):
    folder = os.path.join(os.getcwd(), folder)
for key, relpath in filepaths.items():
    path = os.path.join(folder, relpath)
    filepaths[key] = os.path.normpath(path)

print('\nClimate change config file: ', filepaths['climate_config'])

# Check if running in test mode
if "test_runs" in sim_config.keys():
    # this will limit the number of runs done - used for testing the code
    test_runs = sim_config['test_runs']
    print(f'Running in test mode with {test_runs} number of runs.')
else:
    test_runs = 0

# Run through the simulations
log = {}
for index, sim in sim_df.iterrows():
    if sim['Include'] == 'yes':
        log_data = {
            'Simulation': [sim['Output file']],
            'Start time': [datetime.now().strftime("%Y-%m-%d %H:%M")],
            'End time': [''],
            'Computer': [os.environ['COMPUTERNAME']],
            'Method': [sim['Method']],
            'Run models': [sim['Run models']],
            'Analyse results': [sim['Analyse results']],
            'Executable': [sys.argv[0]],
            'Config file': [sys.argv[1]],
            'Model config': [filepaths['model_config']],
            'Storm config': [filepaths['storm_config']],
            'Climate config': [filepaths['climate_config']],
            'Method config': [sim['Config file']]
        }

        print(sim['Config file'])
        model_method = sim['Method']  # Method is either: monte carlo | ensemble
        if model_method == 'monte carlo':
            sim = MonteCarloSimulator(sim, filepaths, test_runs)
        elif model_method == 'ensemble':
            sim = EnsembleSimulator(sim, filepaths, test_runs)
        else:
            raise Exception(f'Modelling method {model_method} not recognised. Use either: monte carlo | ensemble')

        # Output the logged info
        log_data['End time'] = [datetime.now().strftime("%Y-%m-%d %H:%M")]
        df = pd.DataFrame(log_data)
        log_filepath = str(sim_config['simulation_list']).replace('.xlsx', '_log.csv')
        if os.path.isfile(log_filepath):
            print('Appending to log file:', log_filepath)
            df.to_csv(log_filepath, mode='a', index=False, header=False)
        else:
            print('Creating log file:', log_filepath)
            df.to_csv(log_filepath, index=False)

