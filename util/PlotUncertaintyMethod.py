# -*- coding: utf-8 -*-
"""
Created on Tue May 20 16:52:18 2025

@author: SharpeR
"""

import pandas as pd
import os
import matplotlib.pyplot as plt


folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E012\US_E012\sims_mc\results'

filename = 'TFD_mc_C030_ebf_24h_L20-4_GWL1p3_level_perc'

raw_path = os.path.join(folder, f'{filename}.csv') 
smooth_path = os.path.join(folder, f'{filename}_smooth.csv') 

raw = pd.read_csv(raw_path)
smooth = pd.read_csv(smooth_path)
smooth = smooth.iloc[:-1]


fig, ax = plt.subplots()
ax.plot(raw['z'], raw['5th percentile'], 'o', markeredgewidth=0, 
        markersize=4, label='Raw', color='r', alpha=0.5)
ax.plot(raw['z'], raw['95th percentile'], 'o', markeredgewidth=0, 
        markersize=4, color='r', alpha=0.5)

ax.plot(smooth['z'], smooth['5th percentile'], label='Smooth', color='k')
ax.plot(smooth['z'], smooth['95th percentile'], color='k')

ax.set_xticks(smooth['z'])
ax.set_xticklabels(smooth['AEP'], rotation=90)
ax.set_xlabel("AEP (1 in X)")
ax.set_ylabel('Lake level (m AHD)')
ax.legend()
ax.set_title('24-hour storm duration')
