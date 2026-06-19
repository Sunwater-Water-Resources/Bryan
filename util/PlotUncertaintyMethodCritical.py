# -*- coding: utf-8 -*-
"""
Created on Tue May 20 17:17:44 2025

@author: SharpeR
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np


folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E012\US_E012\sims_mc\results'

filename = 'TFD_mc_C030_ebf_00h_L20-4_GWL1p3_level.csv'

df = pd.read_csv(os.path.join(folder, filename)) 
df = df.iloc[:-1]

# make monotonic
# lower = df.loc[df['5th percentile_smooth'] > 0]
mon_inc = np.maximum.accumulate(df['5th percentile_smooth'])

fig, ax = plt.subplots()
ax.plot(df['z'], df['5th percentile'], 'o', markeredgewidth=0, 
        markersize=4, label='Initial', color='r', alpha=0.5)
ax.plot(df['z'], df['95th percentile'], 'o', markeredgewidth=0, 
        markersize=4, color='r', alpha=0.5)

ax.plot(df['z'], mon_inc, label='Final', color='k')
ax.plot(df['z'], df['95th percentile_smooth'], color='k')

ax.set_xticks(df['z'])
ax.set_xticklabels(df['aep (1 in x)'], rotation=90)
ax.set_xlabel("AEP (1 in X)")
ax.set_ylabel('Lake level (m AHD)')
ax.legend()
ax.set_title('Critical duration envelope')