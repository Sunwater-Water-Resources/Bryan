# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 12:23:43 2024
Read and convert RMC-BestFit outputs

@author: PanosotG
"""
# 
import pandas as pd
from scipy.special import ndtri

class rmc:
    def __init__(self, file):
        csv = pd.read_csv(file, index_col=0)
        ams_df = csv[['Systematic Data_y']].copy()
        ams_df['Zstd'] = ndtri(1 - csv['Systematic Data_x'])
        ams_df.set_index('Zstd', inplace=True)
        ams_df.rename(columns={'Systematic Data_y': 'Annual max series'}, inplace=True)
        ams_df.dropna(inplace=True)
        ams_df.sort_index(inplace=True)
        self.annual_max_series = ams_df['Annual max series']
        
        ci_df = csv[['90% Credible Intervals_y', '90% Credible Intervals_y2']].copy()
        ci_df['Zstd'] =  ndtri(1 - csv['90% Credible Intervals_x'])
        ci_df.set_index('Zstd', inplace=True)
        ci_df.dropna(inplace=True)
        ci_df.sort_index(inplace=True)
        self.confidence_interval = ci_df
        
        ffa = csv[['Posterior Mode_y']].copy()
        ffa['Zstd'] =  ndtri(1 - csv['Posterior Mode_x'])
        ffa.set_index('Zstd', inplace=True)
        ffa.rename(columns={'Posterior Mode_y': 'FFA'}, inplace=True)
        ffa.dropna(inplace=True)
        ffa.sort_index(inplace=True)
        self.ffa = ffa['FFA']
        
        
        
        