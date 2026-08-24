"""Critical duration analysis over a set of storm durations.

For one result type, reads the standard-AEP quantile file of each duration,
takes the maximum at every AEP, records which duration produced it, and writes
the table and a plot of the duration curves:

    <output folder>/<output name>.csv
    <output folder>/<output name>_durations.png

The analysis itself lives in ``util/UtilModule.py`` and is unchanged - this is
the entry point onto it. ``CrticalDurationAnalysis.py`` beside this file is the
same thing driven by a hard-coded block of paths, for a study you are sitting
in front of; this one takes them as arguments, which is what lets the run
launcher's Results page export exactly what it is showing.

    python util/CriticalDurationAnalysis.py --result-type level \\
        --output-folder .../results --output-name CLD_mc_E010_level_critical \\
        --sim 24 .../CLD_mc_24h_E010_level.csv \\
        --sim 36 .../CLD_mc_36h_E010_level.csv

Run it with Bryan's own interpreter: UtilModule needs scipy and matplotlib.

The confidence percentile columns come from the ``_perc_smooth.csv`` files
written beside each quantile file. The reservoir routing method does not write
those, so results routed by it produce a table without those columns rather
than a failure.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from UtilModule import MonteCarloSimulation, MonteCarloSimulationGroup


def run_group(simulations, result_type, output_folder, output_name,
              drop_aeps=None, storage_filepath=None, plot=True):
    """Analyse one result type across durations.

    ``simulations`` is a sequence of (duration, quantile file path). Returns the
    group, whose ``critical_durations`` holds the table that was written.
    """
    drop_aeps = list(drop_aeps or [])
    group = MonteCarloSimulationGroup(output_folder, output_name=output_name,
                                      drop_aeps=drop_aeps)
    if storage_filepath:
        group.set_volume_curve(storage_filepath)

    for duration, filepath in simulations:
        group.add_simulation(
            MonteCarloSimulation(filepath, result_type, duration), duration)

    group.compute_critical_durations(plot=plot)
    return group


def _aep(text):
    """An AEP argument, kept whole when it is whole.

    The quantile files index the standard AEPs as integers, so a float 2.0 asks
    pandas to match across dtypes when the plot drops it.
    """
    value = float(text)
    return int(value) if value.is_integer() else value


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='CriticalDurationAnalysis.py',
        description='Critical duration analysis across storm durations.')
    parser.add_argument('--sim', nargs=2, action='append', required=True,
                        metavar=('DURATION', 'QUANTILE_FILE'),
                        help='one storm duration and its quantile file; repeat '
                             'for each duration')
    parser.add_argument('--result-type', required=True,
                        help="the column to read: inflow, level, outflow, or a "
                             "volume column such as Vol24h")
    parser.add_argument('--output-folder', required=True)
    parser.add_argument('--output-name', required=True,
                        help='without an extension - .csv and _durations.png '
                             'are added')
    parser.add_argument('--drop-aep', type=_aep, action='append', default=[],
                        metavar='AEP',
                        help='an AEP (1 in X) to leave off the PLOT; the table '
                             'still holds every AEP. Repeatable')
    parser.add_argument('--storage', default=None, metavar='ELS_FILE',
                        help='a storage curve, to derive volume results from a '
                             'level analysis')
    parser.add_argument('--no-plot', action='store_true')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    simulations = []
    for duration, filepath in args.sim:
        try:
            value = float(duration)
        except ValueError:
            print(f'ERROR: --sim duration must be a number, not {duration!r}')
            return 2
        # A whole duration stays whole: the labels go into the table's
        # critical_duration column, and the convention there is '24h', not
        # '24.0h'. Callide's 4.5 hour storm keeps its fraction.
        simulations.append((int(value) if value.is_integer() else value, filepath))
    simulations.sort()

    missing = [path for _, path in simulations if not os.path.isfile(path)]
    if missing:
        # Say so here rather than leaving UtilModule to warn per file and then
        # produce a table quietly missing those durations.
        print('ERROR: quantile files not found:')
        for path in missing:
            print(' ', path)
        return 2

    os.makedirs(args.output_folder, exist_ok=True)
    print('Critical duration analysis')
    print('  result type :', args.result_type)
    print('  durations   :', ', '.join(f'{d:g}h' for d, _ in simulations))
    print('  output      :', os.path.join(args.output_folder,
                                          f'{args.output_name}.csv'))

    group = run_group(simulations, args.result_type, args.output_folder,
                      args.output_name, drop_aeps=args.drop_aep,
                      storage_filepath=args.storage, plot=not args.no_plot)

    if group.critical_durations is None:
        print('ERROR: no quantiles could be read - nothing was written.')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
