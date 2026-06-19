# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:18:02 2024

@author: SharpeR
"""
import pandas as pd
from scipy.special import ndtri, ndtr
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import os

gwl = 1.3

method_options = {1: 'fitted polynomial',
                  2: 'Interpolation'}
degrees = 6
method = 2

# ----------------------------------------------------------------------
# The block below was added to interpolate up to PMF useing E012 results
# This is for the saddle dam aep interpolation which is above PMPF
add_pmf = True
pmf_all = {
       0: {'AEP': 16891955, 'level': 677.58},
       1.3: {'AEP': 25271837, 'level': 677.9},
       1.7: {'AEP': 24377139, 'level': 678.0},
       2.7: {'AEP': 32320611, 'level': 678.25}
       }
pmf = pmf_all[gwl]
# ----------------------------------------------------------------------


folder = r"C:\PythonProjects\TFD_2024\03_Design\runs\E012\US_E012\sims_mc\results"
# folder = r"C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\results"
filename = 'TFD_mc_C030_ebf_00h_L20-4_GWL-TEMP-_level.csv'
filepath = os.path.join(folder, filename)

gwl_str = str(gwl).replace('.', 'p')
filepath = filepath.replace('-TEMP-', gwl_str)

flood_levels = pd.DataFrame(index=['Flood 2', 'Flood 3', 'Flood 4'], 
                            data=[674.3, 676.31, 677.0],
                            columns=['Level'])

flood_levels = pd.DataFrame(index=['Main dam', 'Lamb St', 'Saddle Dam', 'Saddle top'], 
                            data=[674.31, 676.54, 677.79, 678.13],
                            columns=['Level'])

# flood_levels = pd.DataFrame(index=['PMF'], data=[677.9], columns=['Level'])
print('Opening:', filepath)
df = pd.read_csv(filepath, index_col=0)
df = df.loc[df.index > 10]
df = df.iloc[:-1]
# df.drop(['6h', '9h', '12h', '18h'], axis=1, inplace=True)
sub = [1000000, 1876173]
df.loc[sub, 'max'] = df.loc[sub, '24h']
if 'z' not in df.columns:
    df['z'] = ndtri(1 - 1 / df.index)
if add_pmf:
    df.loc[pmf['AEP'], 'max'] = pmf['level']
    df.loc[pmf['AEP'], 'z'] = ndtri(1 - 1 / pmf['AEP'])
print(df)

fig, ax = plt.subplots()
ax.plot(df['z'], df['5th percentile'],'o', color='gray')
ax.plot(df['z'], df['max'],'-', linewidth=2.0, color='k')
ax.plot(df['z'], df['95th percentile'],'o', color='gray')

# plt.yscale("log")
plt.ylim(672, 679)
plt.xlim(2, 6)

df.set_index('z', inplace=True)
max_log = np.log10(df[['max']])
max_log.replace([np.inf, -np.inf], np.nan, inplace=True)
max_log.dropna(inplace=True)
max_log = max_log.loc[max_log['max'] > 0]

x1 = max_log.index
y1 = np.maximum.accumulate(max_log['max'])

if method == 1:
    f1 = np.poly1d(np.polyfit(x1, y1, degrees))
    f2 = np.poly1d(np.polyfit(y1, x1, degrees))
    m1 = 10**(f1(x1)) 
    ax.plot(x1, m1, linestyle='dotted')
else:
    f1 = interpolate.interp1d(x1, y1)
    f2 = interpolate.interp1d(y1, x1)

flood_levels['z'] = f2(np.log10(flood_levels['Level']))
flood_levels['AEP'] = np.around(1 / (1 - ndtr(flood_levels['z'])), 0)

ax.plot(flood_levels['z'], flood_levels['Level'],'o', color='red')

plt.xlabel('Standard normal variate')
plt.xlabel('Flow / Level')

print('\nFlood levels and AEPs:')
print(flood_levels)

