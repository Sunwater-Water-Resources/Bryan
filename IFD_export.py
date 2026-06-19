# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 17:12:27 2024

@author: PanosotG
"""

import os
# import json
import pandas as pd
# import numpy as np
from scipy.special import ndtri
from lib.StormGenerator import StormBurst
from lib.RainfallScheme import ExtremeDurationCurve, PMP
from lib.ClimateChange import ClimateAdjustment


# filepath = r'C:\PythonProjects\TFD_2024\03_Design\storm_data\storm_config_02.json'
# filepath = r'C:\PythonProjects\FHD_2024\03_Design\storm_data\storm_config.json'
# focal_subcatchments = r'C:\PythonProjects\TFD_2024\03_Design\runs\E001\TFD_dam_subcatchment_areas.csv'
standard_durations = [6, 9, 12, 18, 24, 36, 48, 72, 96, 120, 144, 168]

def separate_ifd_sets(storm_config, subcatchment_data, 
                      storm_durations = standard_durations, 
                      aep_minmax = None, 
                      interpolation_method = None, 
                      pmp_file = None, 
                      climate_config = None, gwl = None):
    out_folder = os.path.dirname(storm_config)
    
    subarea_df = pd.read_csv(subcatchment_data, index_col=0)
    
    if pmp_file:
        areas_col = subarea_df.loc[:,'Area'].to_dict()
        reg_PMP = regional_PMP(pmp_file, areas_col)
    else:
        reg_PMP = [None] * len(subarea_df.index)
    
    subarea_df['Zero'] = 0.0
    zeroes = subarea_df['Zero']
    for i in subarea_df.index:
        series = zeroes.copy()
        series[i] = subarea_df.loc[i, 'Area']
        name = subarea_df.loc[i, 'Name']
        
        ifd_df = compile_ifd(storm_config, subcatchment_data = series,
                             storm_durations=storm_durations, aep_minmax=aep_minmax, 
                             interpolation_method=interpolation_method,
                             climate_config=climate_config, gwl=gwl, 
                              pmp_obj= reg_PMP[i])
        
        out_file = os.path.join(out_folder, f'Catchment_IFD_({name}).csv')
        ifd_df.to_csv(out_file)
        
        ax = ifd_df.T.plot(logx = True, grid=True, xlabel = 'AEP 1 in Y', ylabel = 'Rain depth (mm)', title = f'Catchment {name} IFD curves')
        ax.figure.savefig(f'{out_folder}/Catchment_IFD_({name}).png', dpi = 120)
    return

def compile_ifd(storm_config, subcatchment_data, 
                storm_durations = standard_durations, 
                aep_minmax = None, 
                interpolation_method = None, 
                pmp_obj = None, 
                climate_config = None,
                gwl = None):
    
    if climate_config:
        climate = ClimateAdjustment(config_file=climate_config, method='gwl', gwl=gwl)
    
    storm = StormBurst(storm_config)
    if type(subcatchment_data) is str:
        storm.load_subcatchment_areas(subcatchment_data)
    elif type(subcatchment_data) is pd.Series:
        storm.subcatch_area = subcatchment_data
        storm.area = subcatchment_data.sum()
    
    # storm.skip_extreme_methods([24, 24])
    storm.skip_method = 'gsdm'          # Hard coded.
    storm_method = 'GTSMR'              # Hard coded.
    
    print(storm.area)
    
    # Set up the rainfall depths
    storm.import_rare_rainfall()  # importing IFD depths more frequent than 1 in 2,000 AEP
    storm.apply_areal_reduction()  # applying the areal reduction factors - BEFORE extreme rainfall!
    
    if interpolation_method:
        storm.rainfall.interpolation_method = interpolation_method  # Override options: ['GEV', 'SiriwardenaWeinmann1998', 'HillAndOthers2000']
    
    if aep_minmax is None:
        if pmp_obj:
            storm.rainfall.pmp = pmp_obj
            for duration in storm_durations:
                gtsmr_df = storm.rainfall.arr_ifd_curves[duration].df.loc[[1000, 2000]].copy()
                gtsmr_df.loc[pmp_obj.aep_of_pmp] = storm.rainfall.pmp.df_gtsmr.loc[duration]
                gtsmr_duration_obj = ExtremeDurationCurve(gtsmr_df,
                                                          duration=duration,
                                                          interpolation_method=storm.rainfall.interpolation_method)
                storm.rainfall.gtsmr_ifd_curves[duration] = gtsmr_duration_obj
                print(gtsmr_df)
        else:
            storm.set_up_extreme_rainfall()  # importing the PMP depths and setting up extreme depths rarer than 1 in 2,000
        aep_of_pmp = storm.rainfall.pmp.aep_of_pmp
        aep_lower = storm.rainfall.aep_lower
        aep_list = list_std_aeps(aep_lower, aep_of_pmp)
    else:
        aep_list = list_std_aeps(aep_minmax[0], aep_minmax[1])
        if aep_minmax[1] > 2000:
            if pmp_obj:
                storm.rainfall.pmp = pmp_obj
                for duration in storm_durations:
                    gtsmr_df = storm.rainfall.arr_ifd_curves[duration].df.loc[[1000, 2000]].copy()
                    gtsmr_df.loc[pmp_obj.aep_of_pmp] = storm.rainfall.pmp.df_gtsmr.loc[duration]
                    gtsmr_duration_obj = ExtremeDurationCurve(gtsmr_df,
                                                              duration=duration,
                                                              interpolation_method=storm.rainfall.interpolation_method)
                    storm.rainfall.gtsmr_ifd_curves[duration] = gtsmr_duration_obj
                    print(gtsmr_df)
            else:
                storm.set_up_extreme_rainfall()  # importing the PMP depths and setting up extreme depths rarer than 1 in 2,000
    
    data = {}
    for aep in aep_list:
        rain_sample_z = get_rain_z(aep)
        curve = {}
        for duration in storm_durations:
            rain_depths = storm.rainfall.get_depth_z(z=rain_sample_z,
                                                     duration=duration,
                                                     storm_method=storm_method)
            ave_rain = storm.get_average_rain(rain_depths)
            
            if climate_config:
                rainfall_climate_adjustment = climate.get_rainfall_uplift_factor(duration=duration)
            else:
                rainfall_climate_adjustment = 1.0
                
            curve[duration] = ave_rain * rainfall_climate_adjustment
        # column = pd.Series(curve, name=aep)
        # data.append(column)
        data[aep] = curve
    df = pd.DataFrame(data)
    return df

def get_rain_z(aep):
    f = 1 - 1 / aep
    z = ndtri(f)
    return z

def list_std_aeps(min_aep, max_aep):
    if min_aep > max_aep:
        raise Exception('max_aep must be greater than min_aep')
    elif min_aep == max_aep:
        return [min_aep]
    
    template = [2, 5, 10]
    ls = []
    
    while template[-1] < min_aep:
        template = [a*10 for a in template]
    
    if template[1] < min_aep:
        ls.append(template[2])
    elif template[0] < min_aep:
        ls.extend(template[1:])
    else:
        ls.extend(template)
    
    template = [a*10 for a in template]
    
    while template[-1] < max_aep:
        ls.extend(template)
        template = [a*10 for a in template]
    
    if template[1] < max_aep:
        ls.extend(template[:2])
    elif template[0] < max_aep:
        ls.append(template[0])
        
    ls.append(max_aep)
    
    return ls
    
def regional_PMP(datafile, areas_col):
    reg_PMP = {}
    pmp_depths = pd.read_csv(datafile, index_col=0)
    pmp_depths.columns = pmp_depths.columns.astype(int)         # Convert columns to int
    for subcatchment in pmp_depths.columns:
        pmp_obj = PMP()
        pmp_obj.df_gtsmr = pmp_depths
        area = round(areas_col[subcatchment], 0)
        pmp_obj.set_aep_of_pmp(area)
        reg_PMP[subcatchment] = pmp_obj
    return reg_PMP
    
# class regional_PMP:
#     def __init__(self, datafile, areas):
#         self.pmp_depths = pd.read_csv(datafile, index_col=0)
#         self.aep_of_pmp = self.set_aep_of_pmp(areas)
        
#     def set_aep_of_pmp(self, areas):
#         aeps = {}
#         for name, area in areas.items():
#             if area < 100:  # 100 km²
#                 aep_of_pmp = 1e-7
#             elif area > 1e5:  # 100,000 km²
#                 aep_of_pmp = 1e-4
#             else:
#                 aep_of_pmp = 10 ** (np.log10(area) - 9)
    
#             aep_of_pmp = 1 / aep_of_pmp  # Express as 1 in Y
#             aeps[name] = aep_of_pmp
#         return aeps
        
    