"""
lib/EnbAnalysis.py

The ensemble result analysis: for each AEP, the median temporal pattern of every
storm duration, the box plots behind it, and the critical duration taken as the
duration whose median is the largest.

Extracted verbatim from EnsembleSimulator.analyse_results so that the
`reservoir routing` method can re-analyse a re-routed ensemble without going
through a Simulator. The body only ever depended on `self.outputfile`, so it
takes that plus the results dataframe.

Outputs, all relative to the folder holding `outputfile`:
    plots/<outputfile>_<aep>_bp.png   box plots of the three result types
    csv/<outputfile>_<aep>_med.csv    the median of each duration
    csv/<outputfile>_critical.csv     the critical duration for each AEP
"""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


RESULT_TYPES = ['inflow', 'level', 'outflow']
RESULT_LABELS = {'inflow': 'Flow (m³/s)',
                 'level': 'Elevation (m AHD)',
                 'outflow': 'Flow (m³/s)'}
RESULT_TITLES = {'inflow': 'Dam inflow',
                 'level': 'Lake level',
                 'outflow': 'Dam outflow'}


def analyse_ensemble(df, outputfile, result_types=None,
                     result_labels=None, result_titles=None):
    """Analyse an ensemble results dataframe.

    df           : one row per simulation, with rain_aep, duration, tp,
                   storm_method and the columns being analysed.
    outputfile   : the output basename, as used by the simulators - the plots
                   and csv sub-folders are created beside it.
    result_types : the columns to analyse, defaulting to the three peaks. The
                   inflow volume analysis of lib/ReservoirRouting.py passes its
                   volume columns instead, so it cannot drift from this one.
    result_labels / result_titles : the y-axis label and box plot title of each,
                   keyed by column name. Needed with result_types.
    """
    result_types = list(result_types) if result_types is not None else RESULT_TYPES
    result_labels = result_labels if result_labels is not None else RESULT_LABELS
    result_titles = result_titles if result_titles is not None else RESULT_TITLES
    # A database read from parquet brings storm_method back as a Categorical
    # (dictionary-encoded), which will not concatenate with a string when the
    # composite pattern label is built. Reading the same table from csv gives
    # object dtype, which is what this was written against.
    df = df.copy()
    df['storm_method'] = df['storm_method'].astype(str)
    aep_lst = df['rain_aep'].unique()
    duration_lst = df['duration'].unique()
    positions = range(len(duration_lst))
    pattern_ids = df['tp'].unique()
    print(df)

    # Set up the critical results dataframe
    headers = []
    for result_type in result_types:
        headers.append(result_type)
        headers.append(f'{result_type}_duration')
        headers.append(f'{result_type}_tp')
    criticals = pd.DataFrame(columns=headers, index=aep_lst)

    # Create sub-folders for outputs
    output_path = os.path.join(os.path.dirname(outputfile), 'plots')
    os.makedirs(output_path, exist_ok=True)
    output_path = os.path.join(os.path.dirname(outputfile), 'csv')
    os.makedirs(output_path, exist_ok=True)

    # Analyse each AEP
    for aep in aep_lst:
        print(f'Working on 1 in {aep} AEP events')
        aep_df = df[df['rain_aep'] == aep]

        # Set up the box plot
        # squeeze=False so a single result type still indexes as ax[ind]
        fig, axes = plt.subplots(nrows=1, ncols=len(result_types),
                                 figsize=(3 * len(result_types), 4), squeeze=False)
        ax = axes[0]
        fig.suptitle(f'1 in {aep} AEP')

        # Set up the median dataframe
        medians = pd.DataFrame(columns=headers, index=duration_lst)

        # Collate the data for each type: inflows, outflows, and levels
        for ind, result_type in enumerate(result_types):
            # extract the data for each duration into the result dataframe
            storm_types = aep_df['storm_method'].unique()
            number_of_storm_types = storm_types.shape[0]
            print(f'Found {number_of_storm_types} storm types:', storm_types)
            result_index = []
            for storm_type in storm_types:
                for pattern_id in pattern_ids:
                    index_label = f'{storm_type}: {pattern_id}'
                    result_index.append(index_label)
            result = pd.DataFrame(index=result_index)
            for duration_ind, duration in enumerate(duration_lst):
                print(f'Working on {duration} hour storm duration')
                dur_df = aep_df[aep_df['duration'] == duration].copy()
                dur_df['composite_index'] = dur_df['storm_method'] + ': ' + dur_df['tp'].astype(str)
                dur_df.set_index('composite_index', inplace=True)
                print(dur_df)
                result[duration] = dur_df[result_type]
                print(result)

                local = dur_df[[result_type]].sort_values(by=result_type, ascending=True)
                median_id = int(np.around(local.shape[0] / 2, 0))
                print('\nGetting median from position:', median_id)
                median = local[result_type].iloc[median_id]
                median_tp = local.index[median_id]
                medians.loc[duration, result_type] = median
                medians.loc[duration, f'{result_type}_tp'] = median_tp
                medians.loc[duration, f'{result_type}_duration'] = duration
                print(f'\n{result_type} results for 1 in {aep} AEP with {duration} hour duration:')
                print('Median value is {} of {} for pattern {}'.format(result_labels[result_type],
                                                                       median,
                                                                       median_tp))
                # Plot the local data points
                print(f'\nDuration index: {duration_ind}')
                x = np.full(result.shape[0], duration_ind)
                print(x)
                ax[ind].plot(x, result[duration], '.', mfc='b', markeredgewidth=0.0)

            # Plot the box data
            print(result)
            plot_data = []
            for res in result.columns:
                local_result = result[res].dropna()
                plot_data.append(local_result)
            ax[ind].boxplot(plot_data, positions=positions, usermedians=medians[result_type])
            # Format the plot
            ax[ind].set_xticks(positions)
            ax[ind].set_xticklabels(duration_lst)
            ax[ind].set_xlabel("Storm duration (hours)")
            ax[ind].set_ylabel(result_labels[result_type])
            ax[ind].set_title(result_titles[result_type])
            ax[ind].grid()

        # Output the plot
        output_file = f'{outputfile}_{aep}_bp.png'
        output_path = os.path.join(os.path.dirname(output_file),
                                   'plots',
                                   os.path.basename(output_file))
        print('\nWriting figure to:', output_path)
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close(fig)

        # Output the medians
        output_file = f'{outputfile}_{aep}_med.csv'
        output_path = os.path.join(os.path.dirname(output_file),
                                   'csv',
                                   os.path.basename(output_file))
        medians.index.name = 'aep (1 in x)'
        medians.to_csv(output_path)

        # Get the critical events
        for result_type in result_types:
            local_df = medians[result_type].astype(float)
            crit_index = local_df.squeeze().argmax()
            criticals.loc[aep, result_type] = medians[result_type].iloc[crit_index]
            criticals.loc[aep, f'{result_type}_tp'] = medians[f'{result_type}_tp'].iloc[crit_index]
            criticals.loc[aep, f'{result_type}_duration'] = medians[f'{result_type}_duration'].iloc[crit_index]

    print('\nFinal results of critical events from median patterns:')
    print(criticals)
    output_file = f'{outputfile}_critical.csv'
    output_path = os.path.join(os.path.dirname(output_file),
                               'csv',
                               os.path.basename(output_file))
    criticals.to_csv(output_path)
    return criticals
