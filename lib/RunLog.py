"""
lib/RunLog.py

The run log written beside the simulation list: one row per simulation, recording
when it ran, what it ran, and how it went.

Used by Main.py. A simulation that fails - most often because one of its input
files is not where the sims list says it is - is written to the log with a FAILED
status and the error, and the batch carries on with the next simulation rather
than stopping part way through the list.
"""
import os
import platform
from datetime import datetime

import pandas as pd

COLUMNS = ['Simulation', 'Start time', 'End time', 'Computer', 'Method', 'Run models',
           'Analyse results', 'Executable', 'Config file', 'Model config', 'Storm config',
           'Climate config', 'Method config', 'Status', 'Error']

COMPLETED = 'completed'
FAILED = 'FAILED'
FAILED_MISSING_INPUT = 'FAILED - missing input'


def log_filepath(sim_config):
    """The run log sits beside the simulation list, with the same name."""
    return str(sim_config['simulation_list']).replace('.xlsx', '_log.csv')


def start_entry(sim, filepaths, executable, config_file):
    """Open a log entry for a simulation that is about to run."""
    return {
        'Simulation': sim['Output file'],
        'Start time': datetime.now().strftime("%Y-%m-%d %H:%M"),
        'End time': '',
        # COMPUTERNAME is a Windows variable - Bryan's home ground, but not the only place
        # the code gets run.
        'Computer': os.environ.get('COMPUTERNAME', platform.node()),
        'Method': sim['Method'],
        'Run models': sim['Run models'],
        'Analyse results': sim['Analyse results'],
        'Executable': executable,
        'Config file': config_file,
        'Model config': filepaths.get('model_config', ''),
        'Storm config': filepaths.get('storm_config', ''),
        'Climate config': filepaths.get('climate_config', ''),
        'Method config': sim['Config file'],
        'Status': '',
        'Error': '',
    }


def close_entry(entry, status, error=''):
    """Finish a log entry with how the simulation went."""
    entry['End time'] = datetime.now().strftime("%Y-%m-%d %H:%M")
    entry['Status'] = status
    entry['Error'] = error
    return entry


def status_for(error):
    """Flag the missing input files separately - they are the usual reason for a
    simulation not running, and the usual thing to go and fix."""
    if isinstance(error, (FileNotFoundError, NotADirectoryError)):
        return FAILED_MISSING_INPUT
    return FAILED


def write_entry(entry, filepath):
    """Append one entry to the run log, creating it if needed."""
    df = pd.DataFrame([entry], columns=COLUMNS)
    if not os.path.isfile(filepath):
        print('Creating log file:', filepath)
        df.to_csv(filepath, index=False)
        return

    existing_columns = pd.read_csv(filepath, nrows=0).columns.tolist()
    if existing_columns == COLUMNS:
        print('Appending to log file:', filepath)
        df.to_csv(filepath, mode='a', header=False, index=False)
        return

    # A log file from an older version of Bryan (before the Status and Error columns).
    # Appending to it as-is would leave ragged rows, so bring the whole file up to the
    # current columns first. Anything the user has added of their own is kept.
    print(f'Updating the columns in the existing log file: {filepath}')
    columns = existing_columns + [c for c in COLUMNS if c not in existing_columns]
    old = pd.read_csv(filepath).reindex(columns=columns)
    pd.concat([old, df.reindex(columns=columns)]).to_csv(filepath, index=False)


def print_summary(failures, total):
    """Report the simulations that did not run at the end of a batch."""
    if not failures:
        print(f'\nAll {total} simulations completed.')
        return
    print(f'\n{"!" * 60}')
    print(f'{len(failures)} of {total} simulations did not complete:')
    for name, status, error in failures:
        print(f'  {name}: {status}')
        print(f'    {error}')
    print(f'{"!" * 60}')
