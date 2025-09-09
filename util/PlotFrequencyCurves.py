"""
Script used to plot frequency curves of either/or/both the FFA outcomes from RMC Bestfit
and URS model results (ensemble and Monte Carlo methods).

WARNNING: This script was cobbled together over time and is not as tidy as other scripts
"""

import pandas as pd
import os
from scipy.special import ndtri
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle                          
from pathlib import Path

# Insert the folder and plot list filenames below
folder = r'folder containing your spreadsheets that list the required plots'
plot_list_files = ['filename of spreadsheet']
min_aep = 10  # can be specified in the spreadsheet or here
max_aep = 5000000  # can be specified in the spreadsheet or here
warnings = []


# The rest should not need to be edited
def main():
    
    for plot_list in plot_list_files:
        # Import all the data
        filepath = os.path.join(folder, plot_list)
        print('Opening plot list:', filepath)
        
        plot_info = pd.read_excel(filepath, 
                                  sheet_name='Plots',
                                  dtype=pd.StringDtype())
        plot_info.index = plot_info.index.map(str)
        plot_info = plot_info.loc[plot_info['Include'] == 'yes']
        
        urbs_info = pd.read_excel(filepath, 
                                  sheet_name='URBS',
                                  index_col='Key',
                                  dtype=pd.StringDtype())
        urbs_info.index = urbs_info.index.map(str)
        
        ffa_info = pd.read_excel(filepath, 
                                 sheet_name='FFA',
                                 index_col='Key',
                                 dtype=pd.StringDtype())
        ffa_info.index = ffa_info.index.map(str)
        
        # Separate out the data for each plot
        for ffa_plot in plot_info.iterrows():
            # Get the plot info for this plot
            print('\nCreating plot', ffa_plot[0])
            this_plot_info = ffa_plot[1]
            print(this_plot_info)
            if 'Min AEP' in this_plot_info.index:
                global min_aep
                min_aep = int(this_plot_info['Min AEP'])
                print(f'Using minimum AEP of 1 in {min_aep}')
            if 'Max AEP' in this_plot_info.index:
                global max_aep
                max_aep = int(this_plot_info['Max AEP'])
                print(f'Using maximum AEP of 1 in {max_aep}')
 
    
            # Get the URBS results keys for this plot               
            urbs_lst = string_to_list(this_plot_info['URBS'])
            print('URBS keys:', urbs_lst)
            
            # Get the FFA keys for this plot
            ffa_lst = string_to_list(this_plot_info['FFA'])               
            print('FFA keys:', ffa_lst)
            if ffa_lst:
                this_ffa_info = ffa_info.loc[ffa_lst]
            else:
                this_ffa_info = None
            
            # Create the plot
            create_plot(this_plot_info,
                        this_ffa_info,
                        urbs_info.loc[urbs_lst])
            plt.clf()
        ext = '.xlsx'
        if ext in filepath:
            warning_file = filepath.replace('.xlsx', 'errors.txt')
            with open(warning_file, 'w') as f:
                for warning in warnings:
                    f.write(warning+'\n')


def get_ffa_results(filepath):
    try:
        df = pd.read_csv(filepath)
    except Exception:
        warning = f'WARNING: failed to open: {filepath}'
        warnings.append(warning)
        print(warning)
        return None
    ffa = df[['Posterior Mode_x', 
              '90% Credible Intervals_y',
              'Posterior Mode_y',
              '90% Credible Intervals_y2',
              'Posterior Predictive_y']].dropna(how='all')
    ffa['AEP'] = 1 / ffa['Posterior Mode_x']
    ffa = ffa.loc[ffa['AEP'] >= min_aep]
    ffa['z'] = ndtri(1 - 1 / ffa['AEP'])
    ffa = ffa.loc[ffa['AEP'] >= min_aep]
    ffa = ffa.loc[ffa['AEP'] <= max_aep]                                    
    ffa.set_index('AEP', inplace=True)
    ffa.drop('Posterior Mode_x', axis=1, inplace=True)
    ffa.rename(columns={'90% Credible Intervals_y': '90% upper CI',
                        '90% Credible Intervals_y2': '90% lower CI',
                        'Posterior Mode_y': 'Mode',
                        'Posterior Predictive_y': 'Predictive'},
               inplace=True)
    print(ffa)
    
    print('\nGetting AMS data...')
    ams = df[['Systematic Data_x', 
              'Systematic Data_y']]
    ams['AEP'] = 1 / ams['Systematic Data_x']
    ams = ams.loc[ams['AEP'] >= min_aep]
    ams = ams.loc[ams['AEP'] <= max_aep]
    ams['z'] = ndtri(1 - 1 / ams['AEP'])
    ams.drop('Systematic Data_x', axis=1, inplace=True)
    ams.rename(columns={'Systematic Data_y': 'AMS'}, inplace=True)
    ams.set_index('AEP', inplace=True)
    print(ams)
    
    return [ffa, ams]


def string_to_list(input_str):
    if pd.isna(input_str):
        return None
    
    if ',' in str(input_str):
        input_lst = [n.strip() for n in input_str.split(',')]
    else:
        input_lst = [input_str]
    return input_lst


def get_urbs_results(filepathtemplate, durations, result_type, do_output=True):
    print('URBS result type:', result_type)
    all_durations = []
    print('\nFound durations:', durations)
    for duration in durations:
        filepath = filepathtemplate.replace('~DD~', duration)
        print('Opening URBS result file:', filepath)
        try:
            df = pd.read_csv(filepath)
            if df.iloc[0, 0] == 0:
                print('It looks like this is an older format from the TPT analysis')
                print('Remvoing the redundant index column')
                df = df.iloc[:, 1:]
            # df['z'] = ndtri(1 - 1 / df['aep (1 in x)'])
            # Below is a fix for missing label on the AEP (first) column                                                
            df.rename(columns={df.columns[0]: 'aep (1 in x)'}, inplace=True)
            df.set_index('aep (1 in x)', inplace=True)
            df.rename(columns={result_type: f'{duration}h'}, inplace=True)
            all_durations.append(df[[f'{duration}h']])
        except Exception:
            warning = f'WARNING: failed to open: {filepath}'
            warnings.append(warning)
            print(warning)
    if not all_durations:
        warning = f'WARNING: No URBS results found!'
        warnings.append(warning)
        print(warning)
        return None
    df = pd.concat(all_durations, axis=1)

    print('\nRecording the critical durations:')
    df['max'] = df.max(axis=1)
    df['critical_duration'] = np.nan

    print('\nGetting the upper and lower percentiles')
    lower_title = None
    upper_title = None
    percentile_method = 'discrete'
    if len(durations) > 1:
        df['critical_duration'] = df.idxmax(axis=1)
        if percentile_method == 'envelope':
            crit_durations = np.unique(df['critical_duration'].to_numpy())
            crit_durations = [dur.replace('h', '') for dur in crit_durations]
            print(crit_durations)
            lower_ls = []
            upper_ls = []
            for crit_duration in crit_durations:
                filepath = filepathtemplate.replace('~DD~', crit_duration).replace('.csv', '_perc_smooth.csv')
                print('\nOpening the percentiles file:', filepath)
                try:
                    perc_df = pd.read_csv(filepath, index_col=0)
                    print(perc_df)
                    lower_title = perc_df.columns[1]
                    upper_title = perc_df.columns[2]
                    lower_df = perc_df[[lower_title]]
                    lower_df.rename(columns={lower_title: f'{crit_duration}h'}, inplace=True)
                    lower_ls.append(lower_df)
                    upper_df = perc_df[[upper_title]]
                    upper_df.rename(columns={upper_title: f'{crit_duration}h'}, inplace=True)
                    upper_ls.append(upper_df)
                except Exception:
                    print('Unable to open the percentiles file')
            if len(lower_ls) > 0:
                lower_df = pd.concat(lower_ls, axis=1)
                lower_df[lower_title] = lower_df.min(axis=1)
                upper_df = pd.concat(upper_ls, axis=1)
                upper_df[upper_title] = upper_df.max(axis=1)
                # perc_df = pd.concat([lower_df[[lower_title]], upper_df[[upper_title]]], axis=1)
                df[lower_title] = lower_df[lower_title]
                df[upper_title] = upper_df[upper_title]

        elif percentile_method == 'discrete':
            filepath = filepathtemplate.replace('~DD~', '00')
            if '_mc_' in filepath:
                print('Getting the confidence limits')
                for aep in df.index:
                    crit_duration = df.loc[aep, 'critical_duration'].replace('h', '')
                    filepath = filepathtemplate.replace('~DD~', crit_duration).replace('.csv', '_perc_smooth.csv')
                    try:
                        perc_df = pd.read_csv(filepath, index_col=0)
                        lower_title = perc_df.columns[1]
                        upper_title = perc_df.columns[2]
                        df.loc[aep, lower_title] = perc_df.loc[aep, lower_title]
                        df.loc[aep, upper_title] = perc_df.loc[aep, upper_title]
                        # print(perc_df)
                    except Exception:
                        print('Unable to open the percentiles file')
                if upper_title is not None:
                    df = smooth_limits(df, lower_title, upper_title)
        print(df)
        if do_output:
            filepath = filepathtemplate.replace('~DD~', '00')
            print('Writing the critical durations:', filepath)
            df.to_csv(filepath)

    df['z'] = ndtri(1 - 1 / df.index)
    df = df.loc[df.index <= max_aep]
    df = df.loc[df.index >= min_aep]  
    df['max']=df['max'].replace(0,1)
                           
    return [df, lower_title, upper_title]


def smooth_limits(df_in, lower_title, upper_title):
    # split the data and remove zero/NA flows
    df = df_in[[lower_title, upper_title]]
    df['z'] = ndtri(1 - 1 / df.index)
    # df.reset_index(inplace=True)
    # df.set_index('z', inplace=True)
    upper = df[['z', upper_title]]
    lower = df[['z', lower_title]]
    upper[upper_title] = np.log10(upper[upper_title])
    lower[lower_title] = np.log10(lower[lower_title])
    upper[upper_title].replace([np.inf, -np.inf], np.nan, inplace=True)
    lower[lower_title].replace([np.inf, -np.inf], np.nan, inplace=True)
    upper.dropna(inplace=True)
    lower.dropna(inplace=True)
    lower = lower.loc[lower[lower_title] > 0]

    # Make monotonic
    lower = lower.iloc[::-1]
    x1 = upper['z']
    y1 = np.maximum.accumulate(upper[upper_title])
    x2 = lower['z']
    y2 = np.minimum.accumulate(lower[lower_title])

    # Build models of the data
    degrees = 5
    f1 = np.poly1d(np.polyfit(x1, y1, degrees))
    f2 = np.poly1d(np.polyfit(x2, y2, degrees))
    m1 = 10 ** (f1(x1))
    m2 = 10 ** (f2(x2))

    # build the smoothed dataframe
    new_lower_title = f'{lower_title}_smooth'
    lower[new_lower_title] = m2
    new_upper_title = f'{upper_title}_smooth'
    upper[new_upper_title] = m1
    return pd.concat([df_in, lower[[new_lower_title]], upper[[new_upper_title]]], axis=1)


def create_plot(plot_info, ffa_info, urbs_info):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
    ax.grid(which='both')
    
    # Plot the FFA
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    plot_cl = True
    print('\nFFA info:')
    print(ffa_info)
    if ffa_info is not None:
        for this_ffa_info in ffa_info.iterrows():
            this_ffa = this_ffa_info[1]
            print('\nFFA info:', this_ffa_info[0])
            print(this_ffa)
            filepath = os.path.join(this_ffa['Folder'], this_ffa['Filename'])
            ffa, ams = get_ffa_results(filepath)
            ffa_key = 'Mode'
            if pd.notna(this_ffa['Label']):
                amsl='AMS '+this_ffa['Label']
                ffal='FFA '+this_ffa['Label']
            else:
                amsl='AMS'
                ffal='FFA'
            if 'Central' in this_ffa.index:
                if pd.notna(this_ffa['Central']):
                    ffa_key = this_ffa['Central']
                    print(f'Using the posterior {ffa_key} as the central description')
                ax.plot(ffa['z'], ffa[ffa_key], next(linecycler), color='k', label=ffal)
            if plot_cl:
                ax.plot(ffa['z'], ffa['90% upper CI'], '--', color='k', alpha=0.2)
                ax.plot(ffa['z'], ffa['90% lower CI'], '--', color='k', alpha=0.2)
                plot_cl = False
            ax.plot(ams['z'], ams['AMS'], 'o', markeredgewidth=0, markersize=4,
                    label=amsl)
    
    print('\nPloting the URBS results...')
    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)
    min_urbs_flow = 100000

    for this_urbs_info in urbs_info.iterrows():
        this_urbs = this_urbs_info[1]
        print('\nURBS simulation info:', this_urbs_info[0])
        print(this_urbs)
        filepath = os.path.join(this_urbs['Folder'], this_urbs['Filename'])
        durations = string_to_list(this_urbs['Durations'])
        try:
            if this_urbs['Output durations'] == 'yes':
                do_output = True
        except:
            do_output = False
        urbs_results = get_urbs_results(filepath, durations, plot_info['Type'], do_output)
        if urbs_results is not None:
            urbs_df = urbs_results[0].fillna(1)
            lower_title = urbs_results[1]
            upper_title = urbs_results[2]
            ax.plot(urbs_df['z'], urbs_df['max'], next(linecycler),
                    label=this_urbs['Label'])
            if 'Uncertainty' in this_urbs.index:
                if pd.notna(this_urbs['Uncertainty']):
                    lower_label = '{} {}'.format(this_urbs['Label'], lower_title)
                    upper_label = '{} {}'.format(this_urbs['Label'], upper_title)
                    if this_urbs['Uncertainty'] == 'smooth':
                        lower_title = f'{lower_title}_smooth'
                        upper_title = f'{upper_title}_smooth'
                    ax.fill_between(urbs_df['z'],
                                    urbs_df[lower_title],
                                    urbs_df[upper_title],
                                    alpha=0.35, linewidth=1, color='k')
                    ax.plot(urbs_df['z'], urbs_df[lower_title], '--', color='gray',
                            label=lower_label)
                    ax.plot(urbs_df['z'], urbs_df[upper_title], '--', color='gray',
                            label=upper_label)
            min_urbs_flow_check = urbs_df['max'].min()
            print('Minimum flow in URBS results:', min_urbs_flow_check)
            if min_urbs_flow_check < min_urbs_flow:
                min_urbs_flow = min_urbs_flow_check
                print('Reducing minimum URBS flow:', min_urbs_flow)
    
    # Plotting AEP of the PMP
    key = 'AEP of PMP'
    if key in plot_info.index:
        aep_of_pmp = plot_info[key]
        if pd.notna(aep_of_pmp):
            print(f'AEP of the PMP is 1 in {aep_of_pmp}')
            z_of_pmp = ndtri(1 - 1 / float(aep_of_pmp))
            ax.axvline(z_of_pmp, color='k')

    # Plotting levels
    if plot_info['Type'] == 'level':
        x_aep = np.array([min_aep, max_aep])
        z_aep = ndtri(1 - 1 / x_aep)
        key = 'FSL'
        if key in plot_info.index:
            fsl = plot_info[key]
            if pd.notna(fsl):
                fsl = float(fsl)
                print(f'Full supply level is {aep_of_pmp} m AHD')
                # ax.axhline(fsl, linestyle='-', color='gold', alpha=0.5)
                ax.plot(z_aep, [fsl, fsl], '-', label='FSL', color='gold', alpha=0.5)
        key = 'DCL'
        if key in plot_info.index:
            dcl = plot_info[key]
            if pd.notna(dcl):
                dcl = float(dcl)
                print(f'Dam crest level is {aep_of_pmp} m AHD')
                # ax.axhline(dcl, linestyle='-', color='red', alpha=0.5)
                ax.plot(z_aep, [dcl, dcl], '-', label='DCL', color='red', alpha=0.5)
        if 'Level values' in plot_info.index:
            level_labels = string_to_list(plot_info['Level labels'])
            level_values = string_to_list(plot_info['Level values'])
            level_values = [float(level) for level in level_values]
            level_colours = string_to_list(plot_info['Level colours'])
            for level_ind, level_label in enumerate(level_labels):
                level_value = level_values[level_ind]
                level_colour = level_colours[level_ind]
                ax.plot(z_aep, [level_value, level_value],
                        '-', label=level_label,
                        color=level_colour, alpha=0.5)
    # Plot formatting
    if ffa_info is not None:
        min_ffa_flow = ffa['90% lower CI'].min()
        print('\nMinimum flow in FFA:', min_ffa_flow)
    else:
        min_ffa_flow = 100000

    min_flow = min([min_urbs_flow, min_ffa_flow])
    print('Overall minimum:', min_flow)
    print('\nThe minimum flow ')
    if not plot_info['Type'] == 'level':
        ax.set_yscale('log')
        if min_flow > 10000:
            bottom = 10000
        elif min_flow > 1000:
            bottom = 1000
        elif min_flow > 100:
            bottom = 100
        elif min_flow > 10:
            bottom = 10
        else:
            bottom = 1
        ax.set_ylim(bottom=bottom)
    # std_aeps = np.array([2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000])
    std_aeps = get_standard_aeps(min_aep, max_aep)                                         
    std_z = ndtri(1 - 1 / std_aeps)
    ax.set_xticks(std_z)
    ax.set_xticklabels(std_aeps, rotation=90)
    ax.set_xlabel("AEP (1 in X)")
    result_type = plot_info['Type']
    plot_labels = {'inflow': 'Flow (m³/s)',
                   'outflow': 'Flow (m³/s)',
                   'level': 'Lake level (m AHD)'}
    if result_type[0: 3] == 'Vol':
        duration = result_type[3: 6]
        duration = duration.replace('h', '')
        plot_labels = {result_type: f'{duration} hr Volume (m³)'}
    ax.set_ylabel(plot_labels[result_type])
    ax.legend(ncol=1, fontsize='x-small')
    if 'Title' in plot_info.index:
        if pd.notna(plot_info['Title']):
            plt.title(plot_info['Title'])
    Path(plot_info['Folder']).mkdir(parents=True, exist_ok=True)
    filepath = os.path.join(plot_info['Folder'], plot_info['Filename'])
    filepath = f'{filepath}.png'
    print('Writing plot to file:', filepath)
    plt.tight_layout()
    plt.savefig(filepath)


def get_standard_aeps(lower_aep, upper_aep):
    # Create a list of standard AEPS, e.g.: 2, 5, 10, 20, ...
    aep = lower_aep
    std_aeps = []
    # while aep < self.upper_aep/4:
    while aep <= upper_aep:
        # if aep >= 4 * self.lower_aep:
        if aep >= lower_aep:
            std_aeps.append(aep)

        if np.log10(aep / 2) % 1 == 0.0:  # is multiple of 20
            aep = aep * 5 // 2
        else:
            aep = aep * 2
    return np.array(std_aeps)


if __name__ == "__main__":
    main()