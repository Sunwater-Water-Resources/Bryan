# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:39:52 2024

@author: PanosotG
"""

import os
import pandas as pd
from scipy.special import ndtri
from RMCreader import rmc

TFD_SimsList = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\TFD_SimsList_validation.xlsx'
results_folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\sims_mc\results'

rmc_LPIII = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_inflows.csv'
rmc_genpareto = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_inflows_GenPareto.csv'
rmc_vol24 = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_volume_24hr_results.csv'
rmc_vol36 = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_volume_36hr_results.csv'
rmc_vol48 = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_volume_48hr_results.csv'
rmc_vol72 = r'C:\PythonProjects\TFD_2024\03_Design\runs\_Validation\FFA\TFD_FFA_volume_72hr_results.csv'


# Generate nested dictionaries grouping runs together. Customise for the way each dam is set up in the SimsList
def group_runs_list(site, result_type):
    if site in ['tinaroo_inflow', 'tinaroo_outflow']:
        revisions = {
            # '001': [
            #     # 'L50-2',
            #     'L20-2',
            #     'L20-1', 
            #     # 'L5-1',
            #     'L0-0'
            #     ],
            # '002': ['L20-1_ebf'],
            # '003': ['L20-1'],
            # 'C024': ['L20-1',
            #           'L20-3'],
            # 'C026': ['L20-1'],
            'C028': [
                # 'L20-3',
                     'L20-3_ADV_npb',
                      'L20-3_ADV_pbf'
                     ]
            }
        durations = [18, 24, 36, 48, 72, 96]
        # durations = [24, 36, 48, 72] #, 96]
        
        run_set = {}
        for rev, versions in revisions.items():
            for ver in versions:
                duration_group = {}
                for dur in durations:
                    csv = f'TFD_mc_{rev}_{dur}h_{ver}_{result_type}.csv'
                    duration_group[dur] = csv
                run_set[f'{rev} {ver}'] = duration_group
    elif site == 'picnic_xing':
        revisions = {'001': ['L20-1'],
                             # 'L30-1'],
                     'C028': ['L20-3']}
        durations = [18, 24, 36, 48, 72] #, 96]
        
        run_set = {}
        for rev, versions in revisions.items():
            for ver in versions:
                duration_group = {}
                for dur in durations:
                    csv = f'PNX_mc_{rev}_{dur}h_{ver}_{result_type}.csv'
                    duration_group[dur] = csv
                run_set[f'{rev} {ver}'] = duration_group
    else:
        raise Exception('Site not coded')
    
    
    return run_set

def max_quantiles(site, results_folder, result_type='inflow', run_set = None):
    if run_set is None:
        run_set = group_runs_list(site, result_type)

    quantiles = {}
    max_env = []
    for key, duration_group in run_set.items():
        data = []
        for dur, csv in duration_group.items():
            file = os.path.join(results_folder, csv)
            csv = pd.read_csv(file, 
                              usecols=['aep (1 in x)', result_type])
            # csv['aep (1 in x)'] = csv['aep (1 in x)'].astype(int)         # Older outputs may be type float
            csv.set_index('aep (1 in x)', inplace=True)
            col = csv[result_type]
            col.name = dur
            data.append(col)
        
        durations_df = pd.concat(data, axis=1)
        durations_df['Q_crit'] = durations_df.max(axis=1, skipna=False)
        durations_df['Dur_crit'] = durations_df.idxmax(axis = 1).astype(float)
        
        quantiles[key] = durations_df
        
        maxcol = durations_df['Q_crit']
        maxcol.name = key
        max_env.append(maxcol)
        
    max_df = pd.concat(max_env, axis=1).sort_index()
    all_df = pd.concat(quantiles, axis = 1, names = ['set', 'dur']).sort_index()
    return max_df, all_df

def plot_validation(max_df, rmc_file, title=None):
    max_df.reset_index(inplace = True)
    max_df['Zstd'] = ndtri(1 - 1/max_df['aep (1 in x)'])
    max_df.set_index('Zstd', inplace = True)
    max_df.sort_index(inplace=True)
    
    plot_df = max_df.loc[:3.540].copy()  # Excludes quantiles > 1 in 2000
    ax = plot_df.drop('aep (1 in x)', axis =1).plot(logy=True)
    
    ffa = rmc(rmc_file)         # flood frequency analysis output as CSV from RMC program
    ffa.annual_max_series.loc[0.8:].plot(color='black', linestyle='none', marker='o', 
                               label='AMS', legend=True, ax=ax)
    ffa.ffa.loc[0.8:].plot(color='black', label='FFA', legend=True, ax=ax)
    ffa.confidence_interval.loc[0.8:].plot(color='grey', linestyle='dashed', ax=ax, legend=False)
    
    if title:
        ax.set_title(title)
    
    ax.set_xticks(plot_df.index.to_numpy())
    ax.set_xticklabels(plot_df['aep (1 in x)'].to_numpy())
    ax.set_xlabel('aep (1 in x)')
    
    ax.grid(True, which='minor', linestyle='dotted')
    ax.grid(True, which='major')

max_df, quant_df = max_quantiles('tinaroo_inflow', results_folder)
plot_validation(max_df, rmc_genpareto, title='Tinaroo inflows (Generalised Pareto)')
plot_validation(max_df, rmc_LPIII, title='Tinaroo inflows (LP-III)')

Vol24h_run_set = group_runs_list('tinaroo_inflow', 'inflowVol24h')
Vol36h_run_set = group_runs_list('tinaroo_inflow', 'inflowVol36h')
Vol48h_run_set = group_runs_list('tinaroo_inflow', 'inflowVol48h')
Vol72h_run_set = group_runs_list('tinaroo_inflow', 'inflowVol72h')

maxVol24_df, quantVol24_df = max_quantiles('tinaroo_inflow', results_folder, result_type='Vol24h', run_set=Vol24h_run_set)
maxVol36_df, quantVol36_df = max_quantiles('tinaroo_inflow', results_folder, result_type='Vol36h', run_set=Vol36h_run_set)
maxVol48_df, quantVol48_df = max_quantiles('tinaroo_inflow', results_folder, result_type='Vol48h', run_set=Vol48h_run_set)
maxVol72_df, quantVol72_df = max_quantiles('tinaroo_inflow', results_folder, result_type='Vol72h', run_set=Vol72h_run_set)

plot_validation(maxVol24_df, rmc_vol24, title='Tinaroo dam 24h inflow volume')
plot_validation(maxVol36_df, rmc_vol36, title='Tinaroo dam 36h inflow volume')
plot_validation(maxVol48_df, rmc_vol48, title='Tinaroo dam 48h inflow volume')
plot_validation(maxVol72_df, rmc_vol72, title='Tinaroo dam 72h inflow volume')