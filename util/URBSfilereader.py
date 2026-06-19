# -*- coding: utf-8 -*-
"""
Extract data from URBS output files and URBS formatted files

@author: PanosotG
"""
import os
import numpy as np
import pandas as pd

# Extract data from text file in URBS Rainfall Definition File, 
# including reading temporal patterns and spatial patterns into Pandas DataFrames
class urbs_storm():
    def __init__(self, filename):
        with open(filename, 'r') as f:
            self.lines = f.readlines()
        
        # self.stm_descript = self.all_lines[0]
        lines = self.lines[1:]
        for i, line in enumerate(lines):
            if line.title().startswith('Time Increment'):
                self.time_inc = float(line.split(':')[1].strip())
                continue
            elif line.title().startswith('Run Duration'):
                self.run_duration = float(line.split(':')[1].strip())
                continue
            elif line.title().startswith('Storm Duration'):
                self.storm_duration = float(line.split(':')[1].strip())
                continue
            elif line.title().startswith('Pluviograph.'):
                next_line = lines[i+1]
                self.data_interval = float(next_line.split(':')[1].strip())
                
                data =[]
                j = i+2
                while True:
                    next_line = lines[j]
                    if next_line[0].isalpha():
                        break
                    elif next_line.startswith('*'):
                        continue
                    else:
                        data.extend(next_line.split())
                    j+=1
                
                self.temppat = pd.Series(data).astype(float)
                continue
            
            elif line.title().startswith('Rain On Subareas'):
                data =[]
                j = i+1
                while True:
                    next_line = lines[j]
                    if next_line[0].isalpha():
                        break
                    elif next_line.startswith('*'):
                        continue
                    else:
                        data.extend(next_line.split())
                    j+=1
                
                self.subarea_depth = pd.Series(data).astype(float)
                continue
            
            elif line.upper().startswith('IL'):
                self.IL = float(line.split(':')[1].strip())
            elif line.upper().startswith('CL'):
                self.CL = float(line.split(':')[1].strip())

# Extract timeseries data from text file in URBS Pluviograph file format,
# such as basename.a and basename.e URBS output files
class urbs_pluvi():
    def __init__(self, filename):
        with open(filename, 'r') as f:
            self.lines = f.readlines()
        
        t_start, t_inc_seconds, p_num = self.lines[4].split()
        t_start = float(t_start)
        t_inc = float(t_inc_seconds) / 60 / 60         # Convert to hours
        p_num = int(p_num)
        
        t_index = pd.Index(np.linspace(t_start+t_inc, t_inc*p_num, p_num), dtype=float, name='Time step')
        
        sim_name = os.path.splitext(os.path.basename(filename))[0]
        self.pluvi = pd.read_csv(filename, skiprows=5, dtype=float, header=None, names=[sim_name])
        self.pluvi.index = t_index
        self.t_inc = t_inc
        self.p_num = p_num


# Extract timeseries data from URBS csv output file
# including rainfall, rain excess, river levels and flow rates
class urbs_csv():
    def __init__(self, filename):
        with open(filename, 'r') as f:
            all_lines = f.readlines()
        
        self.parameter = {}
        for i, line in enumerate(all_lines):
            if line[:4].isspace():          # First four characters are spaces
                continue
            elif line.title().startswith(r'"Gross Rain'):
                rain_0 = i
                continue
            elif line.title().startswith(r'"Effect. Rain'):
                rain_N = i-1
                excess_0 = i
                continue
            elif line.title().startswith(r'"River Levels'):
                excess_N = i-1
                level_0 = i
                continue
            elif line.title().startswith(r'"Flow Rates'):
                level_N = i-1
                flow_0 = i
                continue
            elif line.title().startswith(r'"Parameter Data'):
                flow_N = i-1
                self.parameter['storm_name'] = all_lines[i+1]
                self.parameter['run_info'] = all_lines[i+2]
                self.parameter['model_parameters'] = all_lines[i+3]
                break
        
        self.rain_df = pd.read_csv(filename, index_col=0, skiprows=rain_0, nrows=(rain_N-rain_0))
        self.excess_df = pd.read_csv(filename, index_col=0, skiprows=excess_0, nrows=(excess_N-excess_0))
        self.level_df = pd.read_csv(filename, index_col=0, skiprows=level_0, nrows=(level_N-level_0))
        self.flow_df = pd.read_csv(filename, index_col=0, skiprows=flow_0, nrows=(flow_N-flow_0))
        
        for df in [self.rain_df, self.excess_df, self.level_df, self.flow_df]:
            new = []
            for col in df.columns:
                new_column_name = col.split('(')[0].strip()
                new.append(new_column_name)
            df.columns = new
        return