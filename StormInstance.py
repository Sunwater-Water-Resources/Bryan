# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:13:47 2024

@author: PanosotG
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.special import ndtr
from lib.StormGenerator import StormBurst
from lib.ClimateChange import ClimateAdjustment
from lib.URBSmodel import UrbsModel
from lib.RORBmodel import RorbModel
from lib.Lake import LakeConditions


# project = 'Tinaroo_E005'
# project = 'Tinaroo_test'
# project = 'Teemburra'
project = 'FHD_E011'


# gwl = None; preburst_tp = None
# if project == 'Teemburra':
#     storm_config = r'C:\PythonProjects\teemburra\TEEM_2024\03_Design\storm_data\storm_config_H006A.json'
#     focal_subcatchments = r'C:\PythonProjects\teemburra\TEEM_2024\03_Design\run\H003\rorb\model\TEEM_dam_catchment_subareas_H003.csv'
#     # rain_sample_z = 3.1638723435887925; tp_sample = 9; storm_method = 'ARR point'; duration = 6
#     rain_sample_z = 2.5093784607472096; tp_sample = 9; storm_method = 'ARR point'; duration = 6
#     preburst_tp = 5
#     preburst_proportion = 1.24
# elif project == 'Tinaroo':
#     # storm_config = r'C:\PythonProjects\TFD_2024\03_Design\storm_data\storm_config_02.json'
#     climate_config = r'C:\PythonProjects\TFD_2024\03_Design\climate_change\climate_config.json'
#     storm_config = r'C:\PythonProjects\TFD_2024\03_Design\storm_data\storm_config_TFD_E002_03.json'
#     focal_subcatchments = r'C:\PythonProjects\TFD_2024\03_Design\sim_options\focal_locations\TFD_dam_subcatchment_areas.csv'
#     # rain_sample_z = 4.82563426095645; tp_sample = 2; storm_method = 'GSDM'; duration = 9
#     rain_sample_z = 1.9678216467457141; tp_sample = 9; storm_method = 'ARR point'; duration = 6; preburst_tp = 12; preburst_proportion = 1.67796
#     # rain_sample_z = 4.8205; tp_sample = 6; storm_method = 'GSDM'; duration = 6; preburst_tp = 3; preburst_proportion = 0.039
#     # rain_sample_z = 3.30536; tp_sample = 5; storm_method = 'GSDM'; duration = 6; preburst_tp = None; preburst_proportion = 0.51376
#     rain_sample_z = 2.0501871854063722; tp_sample = 9; storm_method = 'ARR point'; duration = 6; preburst_tp = 9; preburst_proportion = 0.61109; gwl = 1.7
# elif project == 'FHD':
#     # storm_config = r'C:\PythonProjects\FHD_2024\03_Design\storm_data\storm_config_07.json'
#     storm_config = r'C:\PythonProjects\FHD_2024\03_Design\storm_data\storm_config_06.json'
#     focal_subcatchments = r'C:\PythonProjects\FHD_2024\03_Design\sims_options\focal_locations\dam_focal_location.csv'
#     # rain_sample_z = 1.83810912019202; tp_sample = 4; storm_method = 'ARR areal'; duration = 96
#     # rain_sample_z = 1.90268; tp_sample = 4; storm_method = 'ARR point'; duration = 9; preburst_tp = 16
#     rain_sample_z = ndtri(1 - 1 / 10.0); tp_sample = 9; storm_method = 'ARR point'; duration = 9;

if project == 'Tinaroo_E005':
    sim_config = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\TFD_mc_sims_config_02.json'
    focal_subcatchments = r'C:\PythonProjects\TFD_2024\03_Design\sim_options\focal_locations\TFD_dam_subcatchment_areas.csv'
    events_list = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\representative_events\Downstream_storm_files\_event_list.xlsx'
    storm_output_folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\representative_events\Upstream_storm_files'
    do_filtering = True
    do_preburst_patterns = True
    do_baseflow = True
elif project == 'Tinaroo_test':
    sim_config = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\TFD_mc_sims_config_02.json'
    focal_subcatchments = r'C:\PythonProjects\TFD_2024\03_Design\sim_options\focal_locations\TFD_dam_subcatchment_areas.csv'
    events_list = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\representative_events\Upstream_storm_files\test_list.xlsx'
    storm_output_folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E005\US_E005\sims_mc\representative_events\Upstream_test'
    do_filtering = True
    do_preburst_patterns = True
    do_baseflow = True
elif project == 'FHD_E011':
    sim_config = r'C:\PythonProjects\FHD_2024\03_Design\run\E011\sims_config_DAM.json'
    focal_subcatchments = r'C:\PythonProjects\FHD_2024\03_Design\sims_options\focal_locations\dam_focal_location.csv'


def main():
    
    with open(sim_config) as f:
        config_data = json.load(f)
    
    filepaths = config_data['filepaths']
    folder = os.path.dirname(sim_config)
    if not os.path.isabs(folder):
        folder = os.path.join(os.getcwd(), folder)
    for key, relpath in filepaths.items():
        path = os.path.join(folder, relpath)
        filepaths[key] = os.path.normpath(path)
    
    # Get the event info
    event_df = pd.read_excel(events_list, sheet_name='event_list')
    event_df = event_df.loc[event_df['Include'] == 'yes']
    # print(event_df)
    
    model_config = filepaths['model_config']
    model_info = initialise_model(model_config, storm_output_folder)
    model = model_info['model']
    model.apply_baseflow = do_baseflow
    
    checks = []
    for index, row in event_df.iterrows():
        mcdf = pd.read_csv(os.path.join(row['Result folder'], row['Result filename']), index_col=0)
        sim_id = row['Event']
        sim = mcdf.loc[sim_id]
        check = generate_storm(row, sim, filepaths, model_info)
        checks.append(check)
    
    return pd.concat(checks, axis = 1)
    

def generate_storm(row, sim, filepaths, model_info):
    ##
    duration = row['Duration']
    gwl = row['GWL']
    method = row['Method'].lower()
    ensemble = method == 'ensemble'
    
    ##
    rain_sample_z = sim['rain_z']
    rain_sample_aep = 1 / (1 - ndtr(rain_sample_z))  # 1 in X AEP
    tp_sample = sim['tp']
    storm_method = sim['storm_method']
    
    ##
    preburst_proportion = sim['preburst_proportion']
    if ensemble:
        preburst_tp = 'median'
    else:
        preburst_tp = int(sim['preburst_tp'])
    
    ## 
    initial_loss = sim['initial_loss']
    continuing_loss = sim['continuing_loss']
    
    ####
    climate_config = filepaths['climate_config']
    storm_config = filepaths['storm_config']
    
    
    storm = StormBurst(storm_config)
    storm.load_subcatchment_areas(focal_subcatchments)
    storm_durations=[duration]
    
    storm.import_rare_rainfall()  # importing IFD depths more frequent than 1 in 2,000 AEP
    storm.apply_areal_reduction() 
    storm.skip_extreme_methods(storm_durations, do_preburst=True)
    storm.set_up_extreme_rainfall()
    
    if gwl:
        climate = ClimateAdjustment(config_file=climate_config, method='gwl', gwl=gwl)
        rainfall_climate_adjustment = climate.get_rainfall_uplift_factor(duration)
        rain_climate_adj = {}
        for dur in storm.rainfall.durations:
            rain_climate_adj[dur] = climate.get_rainfall_uplift_factor(duration=dur)
    else:
        rainfall_climate_adjustment = 1.0
        rain_climate_adj = None
    
    
    storm.import_arr_areal_patterns(storm_durations)
    storm.import_arr_point_patterns(storm_durations)
    storm.import_gtsmr_temporal_patterns(storm_durations)
    storm.import_gsdm_temporal_patterns(storm_durations)
    temporal_pattern = storm.get_temporal_pattern(storm_method=storm_method,
                                                  duration=duration,
                                                  tp_sample=tp_sample,
                                                  rain_sample_z=rain_sample_z)
    
    rain_depths = storm.rainfall.get_depth_z(z=rain_sample_z,
                                              duration=duration,
                                              storm_method=storm_method)
    
    ave_rain = storm.get_average_rain(rain_depths)
    
    rain_depths = rain_depths * rainfall_climate_adjustment
    ave_rain = ave_rain * rainfall_climate_adjustment
    
    if do_filtering:
        temporal_pattern, embedded_burst_comment = storm.filter_temppat(temporal_pattern, rain_sample_z,
                                                                        duration, ave_rain, storm_method,
                                                                        buffer = 1.1,  # hard coded buffer to reduce below burst depth
                                                                        climate_adjustment = rain_climate_adj)
    if ensemble:
        storm.import_preburst_patterns(median_only=True)
    else:
        storm.import_preburst_patterns()
    
    if preburst_proportion > 0:
        if preburst_tp == 'median':
            preburst_pattern, _ = storm.preburst_patterns.get_preburst_pattern(duration, preburst_proportion, storm.timesteps, sample_int = 'median')
        else:
            preburst_pattern, _ = storm.preburst_patterns.get_preburst_pattern(duration, preburst_proportion, storm.timesteps, sample_int = preburst_tp)
        
        preburst_pattern1, filter_comment = storm.filter_embedded_bursts_in_preburst(preburst_pattern, temporal_pattern, 
                                                                                      duration, rain_sample_z, ave_rain, 
                                                                                      storm_method, messages = True, 
                                                                                      climate_adjustment = rain_climate_adj) #, buffer=1.0, cap_duration=True)
        preburst_duration = preburst_pattern1.index.min() * -1.0
        temporal_pattern = pd.concat([preburst_pattern1, temporal_pattern], axis=0)      # prepend preburst pattern
    else:
        preburst_duration = 0.0
    
    model = model_info['model']
    model_type = model_info['model_type']
    max_keys = model_info['max_keys']
    # self.bfvf10 = 0
    
    lake_vol = sim['ADV']
    volume_below_fsl = model.full_supply_volume - lake_vol
    model.set_volume_below_fsl(volume_below_fsl)
    
    # create storm file and run the URBS model
    simulation_period = 240             # TODO: Hard coded
    if model_type == 'URBS':
        storm_duration = duration
        if model.apply_baseflow:
            model.insert_baseflow_into_vec(storm.get_bfvf10(), rain_sample_z)
        if do_preburst_patterns:
            rain_depths = rain_depths * (1+preburst_proportion)             # sub-area rain depths need to be scaled up for URBS but not RORB
            storm_duration = duration + preburst_duration                   # storm duration needs to be adjusted for preburst duration for URBS storm files (but not RORB)
            # simulation_period += preburst_duration

        duration_str = model.duration_string(duration)
        storm_filename = f'sim_{str(sim.name).zfill(5)}.{duration_str}h'
        print(f'DARN: {storm_filename}')
        model.create_storm_file(rainfall_depths=rain_depths,
                                temporal_pattern=temporal_pattern,
                                data_interval=storm.timesteps,
                                filename=storm_filename,
                                run_duration=simulation_period,
                                storm_duration=storm_duration,
                                storm_duration_excl_pb=duration,
                                initial_loss=initial_loss,
                                continuing_loss=continuing_loss,
                                ari=np.round(rain_sample_aep, 0),
                                ensemble=tp_sample)

        result_filename = f'sim_{str(sim.name).zfill(5)}_{duration_str}h'
        model.run_storm(storm_name=storm_filename,
                        result_name=result_filename)
        
        max_data = model.get_max_results(result_filename, max_keys)
        max_data['ave_rain'] = ave_rain
        print(max_data)
        print(sim[['mean_rain_mm', 'inflow', 'level', 'outflow']])
        check = pd.concat([sim[['mean_rain_mm', 'inflow', 'level', 'outflow']], pd.Series(max_data)])
        check.name = sim.name
        return check

def initialise_model(config_file, storm_output_folder):
    # Get the config data
    print(f'\nSetting up hydrological model from config file: {config_file}')
    f = open(config_file)
    config_data = json.load(f)
    f.close()
    
    # relpath = os.path.relpath(storm_output_folder, start=os.path.dirname(config_file))

    # Set up the URBS model
    model_type = config_data['model_type'].upper()
    max_keys = config_data['max_keys']
    if model_type == 'URBS':
        model = UrbsModel(config_file, dam_location=max_keys['level'])
    elif model_type == 'RORB':
        model = RorbModel(config_file, dam_location=config_data['ADV_name'])
    else:
        raise Exception(f'Model type must be "URBS" or "RORB". ({model_type} not recognised!)')
    
    model.storms_folder = os.path.join(storm_output_folder, 'storm')
    model.output_folder = os.path.join(storm_output_folder, 'output')
    os.makedirs(model.storms_folder, exist_ok=True)
    os.makedirs(model.output_folder, exist_ok=True)
    model.copy_catchment_data_file()
    
    if model_type == 'URBS':
        model.header = []
        model.header.append(f'cd {model.output_folder}')
        model.header.append(f'del {model.output_folder}\\urbserr.log')
        model.header.append(f'del {model.output_folder}\\urbsout.log')
        model.header.append('set URBS_LOGF=TRUE')
        model.header.append(f'set URBS_LOGD={model.output_folder}')
        model.header.append(f'set URBS_RATS={model.ratings_folder}')
    
    return {'model': model, 'model_type': model_type, 'max_keys': max_keys}