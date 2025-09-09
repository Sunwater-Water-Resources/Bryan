import numpy as np
import pandas as pd
import os
from scipy.special import ndtri, ndtr
import matplotlib.pyplot as plt
from scipy import interpolate

folder = r'directory containing your spreadsheet that lists the representative event queries'  # change this to your folder
filename = 'spreadsheet filename'  # change this to your file
filepath = os.path.join(folder, filename)


def main():
    # Extract representative events for the list of AEPs
    print('Opening file:', filepath)
    df = pd.read_excel(filepath, sheet_name='AEP_list', dtype=pd.StringDtype())
    df = df.loc[df['Include'] == 'yes']
    if df.shape[0] == 0:
        print('\nIt looks like no events have been provided.')
    else:
        pass
        find_events(df, find_aep=False)

    # Extract representative events for the list of loads
    print('Opening file:', filepath)
    df = pd.read_excel(filepath, sheet_name='Load_list', dtype=pd.StringDtype())
    df = df.loc[df['Include'] == 'yes']
    if df.shape[0] == 0:
        print('It looks like no flood loads have been provided.')
    else:
        df['Target AEP'] = np.nan
        find_events(df, find_aep=True)


def find_events(df, find_aep=False):
    # loop through the list of results provided
    for ind, event in df.iterrows():
        print(event)
        filepath = os.path.join(event['Input folder'], event['Input filename'])
        print('\nOpening file:', filepath)
        results = pd.read_csv(filepath, index_col=0)

        # find the AEP if needed
        if find_aep:
            level_df = results[['level', 'level_aep']].astype(float)
            level_df.dropna(inplace=True)
            snvs = ndtri(1 - level_df['level_aep'])
            log_levels = np.log10(level_df['level'])
            f = interpolate.interp1d(log_levels, snvs)
            level = float(event['Level'])
            print(f'Searching for the AEP for a level of {level} m AHD')
            log_level = np.log10(level)
            try:
                do_get_the_max = False
                snv = f(log_level)
                aep = int(np.around(1 / (1 - ndtr(snv)), 0))
                print(f'Assigning AEP of 1 in {aep}')
            except:
                do_get_the_max = True
                print('Could not interpolate the AEP!')
                aep = np.nan
            event['Target AEP'] = aep
            # Add this step to find the AEPS

        # Set the AEP for the axes
        target_aep = float(event['Target AEP'])
        target_z = ndtri(1 - 1 / target_aep)
        rain_aep = target_aep
        if 'Rain AEP' in df.columns:
            if pd.notna(event['Rain AEP']):
                rain_aep = float(event['Rain AEP'])
        rain_z = ndtri(1 - 1 / rain_aep)

        # Compute the distance of the AEPs from the target AEP
        result_key = '{}_aep'.format(event['Result type'])
        results['result_aep'] = 1 / results[result_key]
        print(f'\nWorking on 1 in {rain_aep} AEP for {result_key}')
        results['z'] = ndtri(1 - results[result_key])
        # results['d_rain'] = results['rain_z'] - target_z
        results['d_rain'] = results['rain_aep'] - rain_aep
        # results['d_results'] = results['z'] - target_z
        results['d_results'] = 1 / results[result_key] - target_aep
        results['delta'] = np.sqrt(results['d_rain']**2 + results['d_results']**2)

        # Get the best matches
        results = results.loc[results['rain_aep'] < float(event['AEP of PMP']) * 1.1]
        filter_number = int(event['Filter'])
        results.sort_values(by=['delta', 'd_results'], inplace=True)
        if do_get_the_max:
            results.sort_values(by=['level', 'rain_aep'], ascending=False, inplace=True)
        results = results.iloc[: filter_number]

        # Output the results
        print()
        print(results)
        # filepath = os.path.join(event['Output folder'], event['Output filename'])
        # print('\nOpening file:', filepath)
        # results.to_csv(filepath)
        # more logic like lake_z > 0.0; preburst_p close to 0.5; il_p close to 0.5

        # Plot the results
        fig, ax = plt.subplots(figsize=(7, 6), dpi=120)
        aep_factor = 0.1
        ax.axvspan(target_aep * (1 - aep_factor), target_aep * (1 + aep_factor), color='gray', alpha=0.2, lw=0)
        ax.plot(target_aep, rain_aep, marker='o', markersize=15, color='r', alpha=0.5)
        ax.grid(which='both')
        ax.axvline(target_aep, color='k', alpha=0.5)
        ax.axhline(rain_aep, color='k', alpha=0.5)
        # aep_deviation = 1 / (1 - ndtr(results['delta']))
        # aep_deviation = results['delta'] / target_aep * 100
        aep_deviation = np.sqrt(((results['d_rain'] - rain_aep) / rain_aep)**2 + ((results['d_results'] - target_aep) / target_aep)**2)
        im = ax.scatter(results['result_aep'], results['rain_aep'], c=aep_deviation, zorder=1)
        for i in range(results.shape[0]):
            indx = results.index[i]
            ax.annotate(indx, (results.loc[indx, 'result_aep'], results.loc[indx, 'rain_aep']))
        # max_aep = min([results['result_aep'].max(), results['rain_aep'].max()])
        # ax.set_ylim(top=max_aep)
        ax.set_xlabel('AEP for {} (1 in X)'.format(event['Result type']))
        ax.set_ylabel('AEP for the rainfall (1 in X)')
        plt.xticks(rotation=90)
        plt.colorbar(im, label="Deviation from 1 in {} AEP (%)".format(event['Target AEP']), orientation="vertical")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(event['Output filename'])
        plt.tight_layout
        # plt.show()
        filepath = os.path.join(event['Output folder'], event['Output filename'])
        plt.savefig(f'{filepath}.png')
        plt.close()

        # Set up an output folder
        folder_name = event['Output filename']
        try:
            folder = os.path.join(event['Output folder'], folder_name)
            os.mkdir(folder)
            print('\nFolder created:', folder)
        except:
            print('\nFolder already exists:', folder)

        # Get the URBS results
        main_time_period = float(event['Main burst period'])
        urbs_folder = os.path.join(event['Output folder'], event['URBS results'])
        result_types = ['inflows', 'levels', 'outflows']
        sims = ['sim_{}'.format(str(n).zfill(5)) for n in results.index]
        print('\nGetting UBS results for the simulations:')
        print(sims)
        all_results = {}
        for result_type in result_types:
            if '_mcdf' in event['Input filename']:
                urbs_result_file = event['Input filename'].replace('_mcdf', result_type)
            else:
                urbs_result_file = event['Input filename'].replace('.csv', f'_{result_type}.csv')
            urbs_path = os.path.join(urbs_folder, urbs_result_file)
            print('\nOpening URBS results:', urbs_path)
            urbs_df = pd.read_csv(urbs_path)
            all_sims = []
            for sim in sims:
                print(f'\nGetting {result_type} for sim {sim}')
                sim_df = urbs_df[['Time', sim]]
                sim_df.dropna(inplace=True)
                end_time = sim_df['Time'].iloc[-1]
                time_shift = end_time - main_time_period
                print(f'Shifting time by {time_shift} hours')
                sim_df['Time'] = sim_df['Time'] - time_shift
                sim_df.set_index('Time', inplace=True)
                all_sims.append(sim_df)
                print(sim_df)
            # urbs_df = urbs_df[sims]
            urbs_df = pd.concat(all_sims, axis=1).sort_index()
            print(urbs_df)

            # print(urbs_df)
            # urbs_df.to_csv(os.path.join(folder, urbs_result_file))
            all_results[result_type] = urbs_df

        for sim_ind, sim in enumerate(sims):
            sim_id = int(sim.strip('sim_'))
            rain_aep = np.around(results.loc[sim_id, 'rain_aep'], 1)
            level_aep = np.around(1 / results.loc[sim_id, 'level_aep'], 1)
            outflow_aep = np.around(1 / results.loc[sim_id, 'outflow_aep'], 1)
            plot_name = f'{sim} | Rain AEP: {rain_aep} | Level AEP: {level_aep} | Outflow AEP: {outflow_aep}'

            fig, axs = plt.subplots(nrows=2)
            fig.suptitle(plot_name)

            time = all_results['inflows'].index
            inflow = all_results['inflows'][sim]
            level = all_results['levels'][sim]
            outflow = all_results['outflows'][sim]
            axs[0].plot(time, inflow, label='Inflow', color='b')
            axs[0].plot(time, outflow, label='Outflow', color='r')
            # ax[1] = ax.twinx()
            axs[1].plot(time, level)
            axs[0].set_xlabel('Time (hours)')
            axs[1].set_xlabel('Time (hours)')
            axs[0].set_ylabel('Flow m³/s)')
            axs[1].set_ylabel('Lake level (m AHD)')
            # axs[1].text(0.25, 0.25, f'rain: {rain_aep}',
            #             horizontalalignment='left',
            #             verticalalignment='top')
            axs[0].legend()
            plt.tight_layout()
            plot_name = f'{str(sim_ind).zfill(2)}_{folder_name}_{sim}'
            plot_path = os.path.join(folder, f'{plot_name}.png')
            plt.savefig(plot_path)
            plt.close()

        # Store outcomes in an excel spreadsheet
        filepath = f'{filepath}.xlsx'
        with pd.ExcelWriter(filepath) as writer:
            print('\nOpening file:', filepath)
            results.to_excel(writer, sheet_name="mcdf")
            for result_type in result_types:
                all_results[result_type].to_excel(writer, sheet_name=result_type)


if __name__ == "__main__":
    main()
