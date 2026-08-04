"""
lib/LogFiles.py

Keeps every simulation in a batch on its own log file.

The 'Log file' column of a sims list is often filled in per *group* rather than
per row - the Callide list, for instance, gives 92 simulations 7 log files, one
per GWL. Each simulation opens its log with "w", so the runs sharing a path
overwrite one another and only one survives. Worse, whether a run's output lands
in its own file or a later one's depends on when the previous log's buffer is
flushed, which is why a second or third run's log could come out holding an
earlier run's output.

Rather than silently renaming everything, only the repeats are disambiguated:
a path used by one simulation is left exactly as the sims list asked for it.
"""
import os

import pandas as pd


def _disambiguate(path, output_file, used):
    """A log path close to `path` that no other simulation in the batch has."""
    stem, extension = os.path.splitext(path)
    tag = str(output_file).strip().replace(os.sep, '_').replace('/', '_')
    candidate = f'{stem}_{tag}{extension}' if tag else path
    suffix = 2
    while candidate.lower() in used:
        candidate = f'{stem}_{tag}_{suffix}{extension}' if tag else f'{stem}_{suffix}{extension}'
        suffix += 1
    return candidate


def resolve_duplicates(included, column='Log file', name_column='Output file'):
    """Return `included` with a unique 'Log file' for every row.

    Rows whose log path is already unique keep it untouched. Repeats get the
    simulation's output name folded into the filename, and what changed is
    printed so the new names are not a surprise when someone goes looking.
    """
    if column not in included.columns:
        return included

    included = included.copy()
    used = set()
    changes = []
    for index, row in included.iterrows():
        path = row[column]
        if pd.isna(path) or not str(path).strip():
            continue
        path = str(path)
        if path.lower() in used:
            new_path = _disambiguate(path, row.get(name_column, ''), used)
            included.at[index, column] = new_path
            changes.append((path, new_path))
            path = new_path
        used.add(path.lower())

    if changes:
        print(f'\n{len(changes)} simulation(s) in this batch share a "{column}" with an '
              f'earlier one.\nEach needs its own file, or they overwrite each other - '
              f'renaming the repeats:')
        for old, new in changes:
            print(f'  {old}\n    -> {new}')
    return included
