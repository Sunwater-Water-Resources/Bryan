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
from lib import RunLog
import multiprocessing


def run_simulation(sim_df, filepaths, test_runs, log_filepath, lock, executable, config_file):
    # Run through the simulations
    failures = []
    for index, sim in sim_df.iterrows():
        log_entry = RunLog.start_entry(sim, filepaths, executable=executable, config_file=config_file)

        print(sim['Config file'])
        model_method = sim['Method']  # Method is either: monte carlo | ensemble
        try:
            if model_method == 'monte carlo':
                simulator = MonteCarloSimulator(sim, filepaths, test_runs)
            elif model_method == 'ensemble':
                simulator = EnsembleSimulator(sim, filepaths, test_runs)
            else:
                raise Exception(f'Modelling method {model_method} not recognised. Use either: monte carlo | ensemble')
            status, error = RunLog.COMPLETED, ''
        except Exception as exception:
            # Carry on with the rest of this process's share of the list - see Main.py
            status = RunLog.status_for(exception)
            error = f'{type(exception).__name__}: {exception}'
            failures.append((sim['Output file'], status, error))
            print(f'\n{"!" * 60}')
            print(f'{status}: {sim["Output file"]}')
            traceback.print_exc()
            print(f'Moving on to the next simulation in the list.')
            print(f'{"!" * 60}')

        # The simulators redirect stdout to their own log file - put it back first
        sys.stdout = sys.__stdout__
        RunLog.write_entry(RunLog.close_entry(log_entry, status, error), log_filepath, lock)

    RunLog.print_summary(failures, len(sim_df))
    if failures:
        # Non-zero exit code so the parent process can report the failures
        sys.exit(1)


def main():
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
    sim_df = sim_df[sim_df['Include'] == 'yes']
    print(sim_df)

    number_of_concurrent_runs = 1
    if 'multiprocessing' in sim_config:
        number_of_concurrent_runs = sim_config['multiprocessing']
        print(f'Found multiprocessing key, using {number_of_concurrent_runs} processes')
    else:
        print('Multiprocessing not defined, using a single process')

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

    # Split up the simulations and run them
    print(f'\nSplitting the simulation list into {number_of_concurrent_runs} parts...')
    # Split on the row positions rather than handing the dataframe itself to
    # np.array_split, which returns plain arrays (not dataframes) on newer numpy.
    sim_parts = [sim_df.iloc[positions] for positions
                 in np.array_split(np.arange(len(sim_df)), number_of_concurrent_runs)]
    sim_parts = [part for part in sim_parts if len(part)]
    print(sim_parts)
    # One run log for the whole list, so the processes take turns to write to it
    log_filepath = RunLog.log_filepath(sim_config)
    lock = multiprocessing.Lock()
    processes = []
    for sim_part in sim_parts:
        p = multiprocessing.Process(target=run_simulation,
                                    args=(sim_part, filepaths, test_runs, log_filepath, lock,
                                          sys.argv[0], simulation_file))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    if any(p.exitcode for p in processes):
        print(f'\n{"!" * 60}')
        print('Some simulations did not complete - see the "Status" and "Error" columns in')
        print(f'the run log: {log_filepath}')
        print(f'{"!" * 60}')
        sys.exit(1)
    print(f'\nAll {len(sim_df)} simulations completed.')


if __name__ == "__main__":
    main()
