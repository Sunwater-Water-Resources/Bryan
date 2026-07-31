"""
Script used to run URBS Simulations either as Monte Carlo or Ensemble events
Generally... units are: '1 in X' for AEP, hours for storm duration, and km² for area.
"""

import sys
import os
import json
import traceback
import pandas as pd
import numpy as np
from lib.Simulator import MonteCarloSimulator, EnsembleSimulator
from lib.ReservoirRouting import ReservoirRoutingSimulator
from lib import RunLog
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
log_filepath = RunLog.log_filepath(sim_config)
included = sim_df[sim_df['Include'] == 'yes']
failures = []
for index, sim in included.iterrows():
    log_entry = RunLog.start_entry(sim, filepaths, executable=sys.argv[0], config_file=sys.argv[1])

    print(sim['Config file'])
    model_method = sim['Method']  # Method is: monte carlo | ensemble | reservoir routing
    try:
        if model_method == 'monte carlo':
            simulator = MonteCarloSimulator(sim, filepaths, test_runs)
        elif model_method == 'ensemble':
            simulator = EnsembleSimulator(sim, filepaths, test_runs)
        elif model_method == 'reservoir routing':
            simulator = ReservoirRoutingSimulator(sim, filepaths, test_runs)
        else:
            raise Exception(f'Modelling method {model_method} not recognised. Use one of: monte carlo | ensemble | reservoir routing')
        status, error = RunLog.COMPLETED, ''
    except Exception as exception:
        # Keep going with the rest of the list: one simulation with a missing input
        # file should not throw away the simulations queued behind it. The failure is
        # flagged in the run log and summarised at the end.
        status = RunLog.status_for(exception)
        error = f'{type(exception).__name__}: {exception}'
        failures.append((sim['Output file'], status, error))
        print(f'\n{"!" * 60}')
        print(f'{status}: {sim["Output file"]}')
        traceback.print_exc()
        print(f'Moving on to the next simulation in the list.')
        print(f'{"!" * 60}')

    # Output the logged info. The simulators redirect stdout to their own log file, so
    # put it back first, else this lands in the log of the simulation that just ran.
    sys.stdout = sys.__stdout__
    RunLog.write_entry(RunLog.close_entry(log_entry, status, error), log_filepath)

RunLog.print_summary(failures, len(included))
if failures:
    sys.exit(1)

