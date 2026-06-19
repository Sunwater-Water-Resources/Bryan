# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:18:02 2024

@author: SharpeR
"""
import pandas as pd
from scipy.special import ndtri
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np


filepath = r"C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\results\TFD_mc_C030_ebf_00h_L20-4_GWL1p3_outflow.csv"
df = pd.read_csv(filepath, index_col=0)
df['z'] = ndtri(1 - 1 / df.index)

fig, ax = plt.subplots()
ax.plot(df['z'], df['5th percentile'],'o', color='gray')
ax.plot(df['z'], df['max'],'-', linewidth=2.0, color='k')
ax.plot(df['z'], df['95th percentile'],'o', color='gray')

# plt.yscale("log")
plt.ylim(100, 10000)

df.set_index('z', inplace=True)
upper = np.log10(df[['95th percentile']])
lower = np.log10(df[['5th percentile']])
upper.replace([np.inf, -np.inf], np.nan, inplace=True)
lower.replace([np.inf, -np.inf], np.nan, inplace=True)
upper.dropna(inplace=True)
lower.dropna(inplace=True)
lower = lower.loc[lower['5th percentile'] > 0]
# upper = upper.iloc[::-1]
lower = lower.iloc[::-1]

x1 = upper.index
y1 = np.maximum.accumulate(upper['95th percentile'])
x2 = lower.index
y2 = np.minimum.accumulate(lower['5th percentile'])

degrees = 5
f1 = np.poly1d(np.polyfit(x1, y1, degrees))
f2 = np.poly1d(np.polyfit(x2, y2, degrees))
m1 = 10**(f1(x1))
m2 = 10**(f2(x2))

ax.plot(x1, m1)
ax.plot(x2, m2)

plt.xlabel('Standard normal variate')
plt.xlabel('Flow / Level')

