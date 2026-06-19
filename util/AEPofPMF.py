# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 08:06:47 2025

@author: SharpeR
"""
import pandas as pd
from scipy.special import ndtri, ndtr
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import os

gwl = 2.7

folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\results'
folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E012\US_E012\sims_mc\results'
filename = 'TFD_mc_C030_ebf_24h_L20-4_GWL{}__mcdf.csv'.format(str(gwl).replace('.', 'p'))
filepath = os.path.join(folder, filename)

pmf_levels = {
    # 0: 677.47,
    # 1.3: 677.77,
    # 1.7: 677.86,
    # 2.7: 678.09
    0: 677.58,
    1.3: 677.90,
    1.7: 678.00,
    2.7: 678.25
    }

pmf_level = pmf_levels[gwl]
# pmf_level = 677.79
# pmf_level = 674.30
min_aep = 10
pmp_aep = 1876173
max_aep = 3500000
min_z = ndtri(1 - 1 / min_aep)
pmp_z = ndtri(1 - 1 / pmp_aep)
max_z = ndtri(1 - 1 / max_aep)

df = pd.read_csv(filepath, index_col=0)
print(df)

fig, ax = plt.subplots()
aeps = 1 / df['level_aep']
df['z'] = ndtri(1 - 1 / aeps)
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.set_index('z', inplace=True)
df = df.loc[df.index > min_z]
df = df.loc[df.index < max_z]
df.sort_index(inplace=True)
ax.plot(df.index, df['level'],'-', linewidth=2.0, color='k')
ax.axvline(pmp_z, color='r')

# plt.yscale("log")
# plt.ylim(672, 679)
# plt.xlim(2, 5)`

max_log = np.log10(df['level'])
max_log.replace([np.inf, -np.inf], np.nan, inplace=True)
max_log.dropna(inplace=True)

x1 = max_log.index
y1 = max_log

degrees = 10
f1 = np.poly1d(np.polyfit(x1, y1, degrees))
m1 = 10**(f1(x1))

ax.plot(x1, m1)

f2 = np.poly1d(np.polyfit(y1, x1, degrees))
pmf_z = f2(np.log10(pmf_level))
pmf_aep = int(np.around(1 / (1 - ndtr(pmf_z)), 0))
ax.plot(pmf_z, pmf_level,'o', color='gray')

plt.text(pmf_z - 1.2, pmf_level, f'1 in {pmf_aep:,}') #, fontsize = 22)
plt.title(filename)
plt.xlabel('Standard normal variate')
plt.ylabel('Level')
