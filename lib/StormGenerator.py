"""
This file contains classes used to create the design storms: temporal patterns, areal reduction factors,
pre bursts, initial losses
"""

import os
import pandas as pd
import numpy as np
import json
from scipy.special import ndtri, ndtr
import scipy.interpolate
from lib.RainfallScheme import ifdCurves
from lib.TemporalPatterns import ArealPatterns, PointPatterns, GsdmPatterns, GtsmrPatterns, PreburstPatterns


class StormBurst:
    def __init__(self, json_file, generate_storms=True):
        # open the config file and get contents
        print('Opening the storm config file:', json_file)
        f = open(json_file)
        self.config_data = json.load(f)
        f.close()
        # self.filepaths = self.config_data['file_paths']
        
        if 'do_embedded_burst_filtering' in self.config_data.keys():
            print('Warning: embedded burst filter setting in config file ignored. Use setting in Simulation list to exclude/include embedded burst filter')
        
        self.set_filepaths(json_file)
        self.subcatch_area = None
        self.rainfall = None
        self.areal_reduction = None
        self.tp_df = pd.DataFrame()
        self.all_durations = []
        self.areal_temporal_patterns = None
        self.point_temporal_patterns = None
        self.gsdm_temporal_patterns = None
        self.gtsmr_temporal_patterns = None
        self.preburst_patterns = None
        self.pattern_type = ''
        self.timesteps = 0
        # Import the preburst data and set up initial loss object
        self.storm_initial_loss = None
        self.continuing_loss = None
        self.original_temporal_patterns = []
        self.filtered_temporal_patterns = []
        self.is_stored_hyetographs = False
        
        if generate_storms:
            if 'preburst_proportions' in self.filepaths.keys():
                self.pre_burst_obj = PreBursts(self.filepaths['preburst_proportions'], source='csv')
            else:
                self.pre_burst_obj = PreBursts(self.filepaths["ARR_datahub_file"])
            
        self.area = None
        if 'storm_method_config' in self.config_data.keys():
            self.storm_method_config = self.config_data['storm_method_config']
        else:
            print('WARNING: Zones for changing over from ARR to GTSMR/GSDM and from GSDM to GTSMR have  not been provided')
            print('         Reverting to defaults')
            self.storm_method_config = {"aep_changeover_to_extreme": [100, 2000],
                                        "gsdm_gtsmr_changover_duration": [9, 18]}
        self.do_embedded_burst_filtering = True
        # self.embedded_burst_comment = True
        self.areal_arr_pattern_area = 0.0
        self.gtsmr_pattern_area = 0.0
        self.duration_range = []  # this is the range that will be simulated [shortest, longest]
        self.preburst_envelope_durations = []       # range of storm depth-durations to load for computing pre-burst envelope
        self.skip_method = 'none'
        self.current_method = 'arr'
        if 'CL_limit' in self.config_data.keys():
            cal_cap_info = self.config_data['CL_limit']
            self.apply_cl_cap = cal_cap_info['apply']
            self.cl_cap = cal_cap_info['limit']
        else:
            self.apply_cl_cap = False
            self.cl_cap = None

    def get_aep_of_pmp(self):
        rainfall = ifdCurves(self.filepaths['rare_ifds'])
        return rainfall.get_aep_of_pmp()

    def do_the_cl_capping(self, rain_z):
        if self.continuing_loss < self.cl_cap:
            return self.continuing_loss
        cl_100 = self.continuing_loss
        aep_of_pmp = self.rainfall.pmp.aep_of_pmp
        pmp_z = ndtri(1 - 1 / aep_of_pmp)
        z_100 = ndtri(1 - 1 / 100)
        if rain_z > z_100:
            log_cl = (np.log10(self.cl_cap) - np.log10(cl_100)) / (pmp_z - z_100)
            log_cl = np.log10(cl_100) + (rain_z - z_100) * log_cl
            cl = 10 ** log_cl
            print(f'Capped continuing loss to {cl} mm/hr')
        else:
            cl = self.continuing_loss
        return cl

    def get_bfvf10(self):
        print('Getting the baseflow volume factor for 10% AEP (BFVF10) from ARR file...')
        arr_datahub_file = self.filepaths["ARR_datahub_file"]
        print('Opening the ARR file:', arr_datahub_file)
        baseflow_lines = []
        bfvf10 = 0
        with open(arr_datahub_file) as infile:
            copy = False
            for line in infile:
                if line.strip() == "[BASEFLOW]":
                    copy = True
                    continue
                elif line.strip() == "[END_BASEFLOW]":
                    copy = False
                    continue
                elif copy:
                    baseflow_lines.append(line)
        if not baseflow_lines:
            raise Exception('The baseflow information was not found in the ARR file - uses keyword "[BASEFLOW]"')
        for line in baseflow_lines:
            elements = line.split(',')
            first_element = elements[0]
            if first_element.strip() == 'Volume Factor':
                second_element = elements[1]
                bfvf10 = float(second_element.strip())
        if bfvf10 == 0:
            raise Exception('It looks like the BFVF10 was not found in the ARR file - check if the "Volume Factor" term is in the ARR file!')
        return bfvf10

    def skip_extreme_methods(self, duration_range, do_preburst=False):
        if do_preburst:
            longest_duration_to_load = max(self.rainfall.durations)
            
            if len(duration_range) > 1:
                shortest_duration_to_load = min(duration_range)
            else:
                if duration_range[0] < max(self.rainfall.durations):
                    shortest_duration_to_load = duration_range[0]
                else:
                    shortest_duration_to_load = sorted(self.rainfall.durations)[-2]    # 2nd largest duration
            
            for dur in self.rainfall.durations:
                if dur >= shortest_duration_to_load:
                    self.preburst_envelope_durations.append(dur)
                    
        else:
            shortest_duration_to_load = min(duration_range)
            longest_duration_to_load = max(duration_range)
        
        if shortest_duration_to_load >= 24:
            print('Skipping GSDM patterns... shortest duration is {} hours'.format(duration_range[0]))
            self.skip_method = 'gsdm'
            
        if longest_duration_to_load <= 6:
            print('Skipping GTSMR patterns... longest duration is {} hours'.format(duration_range[0]))
            self.skip_method = 'gtsmr'

    def set_filepaths(self, json_file):
        filepaths = self.config_data['file_paths']
        folder = os.path.dirname(json_file)
        for key, relpath in filepaths.items():
            path = os.path.join(folder, relpath)
            filepaths[key] = os.path.normpath(path)
        self.filepaths = filepaths
        
    def load_subcatchment_areas(self, filepath):
        # if filepath:
            # filename = self.filepaths["subcatchment_areas"]
        print('\nSubcatchment area data is used for weighting the average catchment rainfall')
        print('for the embedded burst filter and computing the preburst depth')
        # print('To exclude subcatchments (e.g. downstream of the dam) from the weighting, set subcatchment areas to zero in this file')
        print('Opening subcatchment area file', filepath)
        areas = pd.read_csv(filepath, index_col=0, usecols = ['Name', 'Area'])
        subcatch_area = areas['Area']
        # else:
        #     subcatch_area = None
        #     print('\nSubcatchment area data not given, mean rainfall will be calculated as a simple average')
        #     print('To weight average rainfall by subcatchment area, specify subcatchment area file in storm_config file')
        self.subcatch_area = subcatch_area
        self.area = subcatch_area.sum()
        print(f'\nThe catchment area to the focal point was computed as {self.area} km²')
        
    def get_average_rain(self, rain_depth, print_msg=True):
        if self.subcatch_area is None:
            average = rain_depth.mean()
            if print_msg:
                print(f'Simple average rainfall: {average:.3f} mm')
        else:
            areas = self.subcatch_area
            total_area = areas.sum()
            product = rain_depth * areas
            # print(product)
            average = product.sum() / total_area
            if print_msg:
                print(f'Weighted average rainfall: {average:.3f} mm')
        return average

    def point_tp_frequency_bins(self, aep):
        # Convert AEP from 1 in X to percent
        aep = 1 / aep * 100
        if aep < 3.2:
            frequency = 'rare'
        elif aep > 14.4:
            frequency = 'frequent'
        else:
            frequency = 'intermediate'
        return frequency

    def import_rare_rainfall(self):
        print('\nImporting rare rainfall')
        self.rainfall = ifdCurves(self.filepaths['rare_ifds'])
        self.rainfall.get_rare_rainfall_depths(area=self.area)

    def apply_areal_reduction(self):
        print('\nSetting up areal reduction factors')
        self.areal_reduction = ArealReduction(arr_datahub_file=self.filepaths["ARR_datahub_file"],
                                              area=self.area)
        self.rainfall.apply_areal_reduction(self.areal_reduction)

    def set_up_extreme_rainfall(self):
        self.rainfall.skip_method = self.skip_method
        print('\nSetting up the extreme rainfall')
        self.rainfall.set_up_extreme_rain()

    def import_gsdm_temporal_patterns(self, storm_durations):
        if not self.skip_method == 'gsdm':
            filepath = self.filepaths['gsdm_patterns']
            self.gsdm_temporal_patterns = GsdmPatterns(filepath, storm_durations)
            self.gsdm_temporal_patterns.import_patterns()

    def import_gtsmr_temporal_patterns(self, storm_durations):
        if not self.skip_method == 'gtsmr':
            filepath = self.filepaths['gtsmr_patterns']
            self.gtsmr_temporal_patterns = GtsmrPatterns(filepath, storm_durations, self.area)
            self.gtsmr_temporal_patterns.import_patterns()
            self.gtsmr_pattern_area = self.gtsmr_temporal_patterns.pattern_area

    def import_arr_areal_patterns(self, storm_durations):
        filepath = self.filepaths['areal_patterns']
        self.areal_temporal_patterns = ArealPatterns(filepath, storm_durations, self.area)
        self.areal_temporal_patterns.import_patterns()
        self.areal_arr_pattern_area = self.areal_temporal_patterns.pattern_area

    def import_arr_point_patterns(self, storm_durations='all'):
        filepath = self.filepaths['point_patterns']
        self.point_temporal_patterns = PointPatterns(filepath, storm_durations)
        self.point_temporal_patterns.import_patterns()
        
    def import_preburst_patterns(self, median_only=False):
        filepath = self.filepaths['preburst_patterns']
        self.preburst_patterns = PreburstPatterns(filepath)
        self.preburst_patterns.import_patterns(median_only)

    def sample_storm_method(self, rain_sample_z=None, rain_sample_aep=None, duration=None, sim_method='monte carlo', df_method=None):
        # The adopted storm method (ARR/GSDM/GTSMR) will depend on the AEP
        # Then using GSDM or GTSMR will depend on the duration
        # print(self.storm_method_config)
        aep_changeover = self.storm_method_config['aep_changeover_to_extreme']
        if rain_sample_z:
            aep = 1 / (1 - ndtr(rain_sample_z))
        elif rain_sample_aep:
            aep = rain_sample_aep
        else:
            raise Exception('rain_sample_z or rain_sample_aep must be provided')
        print(f'\nStorm AEP is 1 in {np.around(aep, 0)}')
        
        # if left of the overlapping range between ARR and GSDM/GTSMR temporal patterns
        if aep <= aep_changeover[0]:
            storm_method = self.apply_rare_method(duration)
            print(f'Storm AEP is less than threshold - using {storm_method}')
        # else if right of the overlapping range between ARR and GSDM/GTSMR temporal patterns
        elif aep > aep_changeover[1]:
            storm_method = self.apply_extreme_method(duration, sim_method)
            print(f'Storm AEP is greater than threshold - using {storm_method}')
        # else if in the overlapping range between ARR and GSDM/GTSMR temporal patterns
        else:
            if sim_method == 'ensemble':
                # For ensemble method just use ARR in the interim AEP range
                interim_method = self.storm_method_config['interim_for_ensemble']
                extreme_pattern = self.storm_method_config['extreme_pattern_for_ensemble']
                if interim_method == 'arr':
                    print('Using ARR pattern in AEP changeover zone')
                    storm_method = self.apply_rare_method(duration)
                elif interim_method == 'extreme':
                    print('Using extreme pattern in AEP changeover zone')
                    storm_method = self.apply_extreme_method(duration, sim_method, extreme_pattern)
                elif interim_method == 'both':
                    if df_method is None:
                        print(f'Using ARR pattern in AEP changeover zone')
                        storm_method = self.apply_rare_method(duration)
                    elif df_method == 'extreme':
                        print(f'Using extreme pattern in AEP changeover zone')
                        storm_method = self.apply_extreme_method(duration, sim_method, extreme_pattern)
                else:
                    print('Reading the "interim_for_ensemble" key in the storm config.')
                    raise Exception(f'The storm method "{interim_method}" is not recognised, use "arr", "extreme", or "both".')
            else:
                pattern_sample = np.random.randint(0, 2)
                if pattern_sample == 0:
                    storm_method = self.apply_rare_method(duration)
                else:
                    storm_method = self.apply_extreme_method(duration, sim_method)
            print(f'Storm AEP is in mixed range - using {storm_method}')
        return storm_method

    def apply_rare_method(self, duration):
        if self.area < 75:
            storm_method = 'ARR point'
        else:
            if duration < 12:
                storm_method = 'ARR point'
            else:
                storm_method = 'ARR areal'
        return storm_method

    def apply_extreme_method(self, duration, sim_method, overlap_method='GTSMR'):
        duration_changeover = self.storm_method_config['gsdm_gtsmr_changover_duration']
        if duration <= duration_changeover[0]:
            storm_method = 'GSDM'
            print(f'Storm duration is short ({duration} hours) - using {storm_method}')
        elif duration >= duration_changeover[1]:
            storm_method = 'GTSMR'
            print(f'Storm duration is long ({duration} hours) - using {storm_method}')
        else:
            if sim_method == 'ensemble':
                if overlap_method == 'GTSMR' or overlap_method == 'GSDM':
                    storm_method = overlap_method
                else:
                    raise Exception(f'{overlap_method} storm method for extreme temporal patterns not recognised!')
            else:
                pattern_sample = np.random.randint(0, 2)
                if pattern_sample == 0:
                    storm_method = 'GSDM'
                else:
                    storm_method = 'GTSMR'
            print(f'Storm duration is in mixed range ({duration} hours) - using {storm_method}')
        return storm_method

    # apply_extreme_method_old is superseded because now sampling randomly btw 6 and 24 hr duration.
    def apply_extreme_method_old(self, duration):
        if duration <= 9:
            storm_method = 'GSDM'
            print(f'Storm duration is short ({duration} hours) - using {storm_method}')
        else:
            storm_method = 'GTSMR'
            print(f'Storm duration is long ({duration} hours) - using {storm_method}')
        return storm_method

    def get_temporal_pattern(self, storm_method, duration, tp_sample, rain_sample_z):
        # The adopted temporal pattern (ARR/GSDM/GTSMR) will depend on the AEP
        # Then using GSDM or GTSMR will depend on the duration
        aep = 1 / (1 - ndtr(rain_sample_z))
        if storm_method == 'ARR areal':
            temporal_pattern = self.get_arr_areal_pattern(duration, tp_sample)
            self.timesteps = self.areal_temporal_patterns.timesteps[duration]
            self.pattern_type = 'ARR areal'
        elif storm_method == 'ARR point':
            frequency = self.point_tp_frequency_bins(aep)
            temporal_pattern = self.get_arr_point_pattern(duration, frequency, tp_sample)
            self.timesteps = self.point_temporal_patterns.timesteps[duration]
            self.pattern_type = 'ARR areal'
        elif storm_method == 'GSDM':
            temporal_pattern = self.get_gsdm_pattern(duration, tp_sample)
            self.timesteps = self.gsdm_temporal_patterns.timesteps[duration]
            self.pattern_type = 'GSDM'
        elif storm_method == 'GTSMR':
            # if duration < 24:
            #     print(f'Shortest duration available for GTSMR is 24 hours')
            #     print(f'Therefore, increasing duration of {duration} hours to 24 hours for temporal pattern extraction')
            #     duration = 24
            temporal_pattern = self.get_gtsmr_pattern(duration, tp_sample)
            self.timesteps = self.gtsmr_temporal_patterns.timesteps[duration]
            self.pattern_type = 'GTSMR'
        else:
            raise Exception(f'The {storm_method} storm method is not recognised!')
        return temporal_pattern

    def get_arr_areal_pattern(self, duration, pattern_number):
        # if duration < 12:
        #     print('increasing temporal pattern to 12 for ARR areal pattern')
        #     duration = 12
        tp_all = self.areal_temporal_patterns.temporal_patterns[duration]
        return tp_all[pattern_number]

    def get_arr_point_pattern(self, duration, frequency, pattern_number):
        # Also depend on the aep frequency
        # print(self.point_temporal_patterns.temporal_patterns)
        tp_all_all = self.point_temporal_patterns.temporal_patterns[duration]
        tp_all = tp_all_all[frequency]
        return tp_all[pattern_number]

    def get_gsdm_pattern(self, duration, pattern_number):
        tp_all = self.gsdm_temporal_patterns.temporal_patterns[duration]
        return tp_all[pattern_number]

    def get_gtsmr_pattern(self, duration, pattern_number):
        # if duration < 24:
        #     print('increasing temporal pattern to 24 for GTSMR pattern')
        #     duration = 24
        tp_all = self.gtsmr_temporal_patterns.temporal_patterns[duration]
        return tp_all[pattern_number]

    def get_preburst_proportion(self, rain_z, preburst_p, duration):
        preburst_proportion = self.pre_burst_obj.get_preburst_proportion(rain_z, preburst_p, duration)
        return preburst_proportion

    def get_areal_pattern_area(self, area):
        if area < 75.0:
            return 0
        elif area < 150.0:
            return 100
        elif area < 350.0:
            return 200
        elif area < 750.0:
            return 500
        elif area < 1750.0:
            return 1000
        elif area < 3750.0:
            return 2500
        elif area < 7500.0:
            return 5000
        elif area < 15000.0:
            return 10000
        elif area < 30000.0:
            return 20000
        else:
            return 40000

    def filter_temppat(self, temporal_pattern, 
                       z, main_burst_duration, ave_rain, storm_method, 
                       buffer=1.0, climate_adjustment = None):
        
        original_temporal_pattern = temporal_pattern.copy()
        
        # Construct temppat filter curve
        # For given Z, get depths for all durations < storm duration, return as a pd.Series
        data = {}
        for dur in self.rainfall.durations:
            if dur >= main_burst_duration:
                continue
            embedded_depth = self.rainfall.get_depth_z(z=z, duration=dur, storm_method=storm_method, print_msg=False)
            embedded_ave = self.get_average_rain(embedded_depth, print_msg=False)                       # Catchment average rain
            
            if climate_adjustment is not None:
                embedded_ave = embedded_ave * climate_adjustment[dur]
            
            normalised = embedded_ave / ave_rain * 100                              # Normalise as percentage of main burst depth
            if normalised > 100:
                print(f'{main_burst_duration}h rain depth: {ave_rain}mm')
                print(f'{dur}h rain depth: {embedded_ave}mm')
                error_msg = f'{dur}h embedded burst depth exceeds {main_burst_duration}h burst depth for z: {z}'
                return original_temporal_pattern, error_msg
            data[dur] = normalised        
        
        temppat_filter_curve = pd.Series(data)
        temppat_filter_curve.sort_index(ascending=True, inplace=True)
        
        embedded_burst_comment = "No embedded bursts"
        # storm_duration = temporal_pattern.index[-1]  # get storm duration 
        
        i = 0
        for dur in temppat_filter_curve.index:
            window = dur / self.timesteps
            if window.is_integer():
                window = int(window)
            else: 
                continue
                # error_msg = 'filter curve timestep indivisible'
                # return original_temporal_pattern, error_msg

            filter_depth = temppat_filter_curve[dur]
            stm_percentage = dur / main_burst_duration * 100
            if filter_depth < stm_percentage:
                error_msg = f'Error: IFD inconsistency: {dur}h intensity < {main_burst_duration}h intensity'
                return original_temporal_pattern, error_msg
            else:
                if buffer >= 1.0:
                    filter_depth = filter_depth * buffer
                else:
                    if filter_depth * buffer >= stm_percentage:
                        filter_depth = filter_depth * buffer
                    else:
                        filter_depth = stm_percentage

            while True:
                embedded_burst = temporal_pattern.rolling(window).sum().max()
                if embedded_burst <= filter_depth:
                    break
                else:
                    embedded_burst_comment = 'Embedded burst filtered'
                    i += 1
                    if i > 200:
                        print('filter iteration exceeds 200')
                        embedded_burst_comment = 'Embedded burst filter nonconvergent'
                        return original_temporal_pattern, embedded_burst_comment

                    # find index of embedded burst
                    time_index = temporal_pattern.rolling(window).sum().idxmax()
                    iloc_end = int(time_index / self.timesteps)
                    iloc_start = iloc_end - window

                    # Scale embedded burst down
                    scale_down = filter_depth / embedded_burst
                    temporal_pattern.iloc[iloc_start:iloc_end] = temporal_pattern.iloc[iloc_start:iloc_end] * scale_down

                    # Scale rest of design storm up
                    # if embedded_burst > 100:
                    #     self.embedded_burst_comment = 'Blowout'
                    #     return original_temporal_pattern
                    # scale_up = (100 - filter_depth) / (100 - embedded_burst)
                    # temporal_pattern.iloc[:iloc_start] = temporal_pattern.iloc[:iloc_start] * scale_up
                    # temporal_pattern.iloc[iloc_end:] = temporal_pattern.iloc[iloc_end:] * scale_up
                    
                    ## Rev 2 code: scale rest of design storm up
                    adjust = embedded_burst - temporal_pattern.iloc[iloc_start:iloc_end].sum()
                    pre_post_ = ~temporal_pattern.index.isin(temporal_pattern.index[iloc_start:iloc_end])  # index excl. embedded burst
                    max_val = temporal_pattern[pre_post_].max()
                    gap0 = (max_val - temporal_pattern[pre_post_]).sum()
                    if gap0 == 0:
                        temporal_pattern[pre_post_] = temporal_pattern[pre_post_] + adjust/len(temporal_pattern[pre_post_])
                    else:
                        scale_fact = (gap0 - adjust) / gap0
                        temporal_pattern[pre_post_] = (1-scale_fact) * max_val + scale_fact * temporal_pattern[pre_post_]
                    
                    # if temporal_pattern.isnull().values.any():
                    #     print(original_temporal_pattern)
                    #     print(temporal_pattern)
                    

        # Round and check still adds up to 100
        temporal_pattern = np.around(temporal_pattern, 2)
        adjust = 100.0 - temporal_pattern.sum()
        if abs(adjust) > 0.1:
            print (temporal_pattern, z, main_burst_duration, ave_rain, storm_method)
            error_msg = 'Error: sum of filtered temporal patterns is not in range of 99.9%–100.1%'
            return original_temporal_pattern, error_msg
        elif abs(adjust) > 0.0:
            temporal_pattern.iloc[-1] += adjust

        return temporal_pattern, embedded_burst_comment
    
    def check_embedded_burst(self, temporal_pattern,                        
                             z, main_burst_duration, ave_rain, storm_method, 
                             buffer=1.0, climate_adjustment = None):
        
        # Construct temppat filter curve
        # For given Z, get depths for all durations < storm duration, return as a pd.Series
        data = {}
        for dur in self.rainfall.durations:
            if dur >= main_burst_duration:
                continue
            embedded_depth = self.rainfall.get_depth_z(z=z, duration=dur, storm_method=storm_method, print_msg=False)
            embedded_ave = self.get_average_rain(embedded_depth, print_msg=False)                       # Catchment average rain
            
            if climate_adjustment is not None:
                embedded_ave = embedded_ave * climate_adjustment[dur]
            
            normalised = embedded_ave / ave_rain * 100                              # Normalise as percentage of main burst depth
            if normalised > 100:
                print(f'{main_burst_duration}h rain depth: {ave_rain}mm')
                print(f'{dur}h rain depth: {embedded_ave}mm')
                return f'Error: {dur}h embedded burst depth exceeds {main_burst_duration}h burst depth for z: {z}'
            data[dur] = normalised        
        
        temppat_filter_curve = pd.Series(data)
        temppat_filter_curve.sort_index(ascending=True, inplace=True)
        
        max_burst = {}
        for dur in temppat_filter_curve.index:
            window = dur / self.timesteps
            if window.is_integer(): 
                window = int(window)
            else: 
                continue
                # return f'Error: duration {dur} is not a multiple of timestep {self.timesteps}'
            
            max_burst[dur] = temporal_pattern.rolling(window).sum().max()
        
        curve_max = pd.Series(max_burst)
        condition = curve_max / temppat_filter_curve
        if condition.max() > 1:
            text = f'{condition.idxmax()}h burst exceeds by {condition.max()-1:.1%}'
            embedded_burst_comment = f'Unfiltered embedded bursts: {text}'
        else:
            embedded_burst_comment = 'No embedded bursts'
        return embedded_burst_comment
    
    # Not used
    def filter_preburst_pattern(self, temporal_pattern, initial_loss_normalised):
        
        preburst_pattern = temporal_pattern.loc[:0]
        mainburst_pattern = temporal_pattern.loc[0.001:]        # Index > 0
        
        if preburst_pattern.sum() <= initial_loss_normalised:
            return 'No embedded burst in pre-burst'
        else:           # Remove initial loss from pre-burst
            excess_pattern = temporal_pattern.copy()
            if preburst_pattern.iloc[0] <= initial_loss_normalised:
                idx0 = preburst_pattern[preburst_pattern.cumsum() <= initial_loss_normalised].index[-1]
                excess_pattern.loc[:idx0] = 0.0
                
            idx1 = preburst_pattern[preburst_pattern.cumsum() > initial_loss_normalised].index[0]
            excess_pattern.loc[idx1] = preburst_pattern.loc[:idx1].sum() - initial_loss_normalised
        
        for i in range(1, len(mainburst_pattern) + 1):
            timeofpeak = excess_pattern.rolling(i).sum().idxmax()
            if timeofpeak < self.timesteps * i:
                return 'Unfiltered embedded burst in pre-burst'
        
        return 'No embedded burst in pre-burst'
    
    # Test version only
    def uniform_preburst0(self, temporal_pattern):
        main_burst = temporal_pattern[temporal_pattern.index > 0]
        preburst_pattern = temporal_pattern.loc[:0]
        preburst_depth = preburst_pattern.sum()         # Normalised depth (i.e. as a percentage of main burst depth)
        
        intensity = []                                  # Identify minimum average intensity (to end of main burst) as preburst intensity limit
        for i in range(1, len(main_burst) + 1):
            limit = main_burst.iloc[-i:].mean()
            intensity.append(limit)
        
        min_timesteps = int(np.ceil(preburst_depth/ min(intensity)))     # Minimum timesteps required to limit preburst intensity
        num_timesteps = max(min_timesteps, len(preburst_pattern))   # Adopt original preburst duration if longer than the minimum
        
        pb_intensity = np.around(preburst_depth/num_timesteps, 2)
        uniform_preburst_pattern = pd.Series([pb_intensity] * num_timesteps)
        uniform_preburst_pattern.iloc[0] = preburst_depth - uniform_preburst_pattern.iloc[1:].sum()     # Correct for rounding errors
        
        index = pd.Index(np.linspace(-1 * num_timesteps * self.timesteps, -1 * self.timesteps, num_timesteps))
        uniform_preburst_pattern.index = index
        temporal_pattern = pd.concat([uniform_preburst_pattern, main_burst], axis=0)

        return temporal_pattern
    
    def uniform_preburst(self, duration, preburst_proportion, timesteps):
        # Arbitrary pre-burst duration
        preburst_duration_scale = {1.0: 4.0,
                                   6.0: 3.0,
                                   9.0: 2.75,
                                   12.0: 2.5,
                                   18.0: 2.25,
                                   24.0: 2.0,
                                   36.0: 1.75,
                                   48.0: 1.5,
                                   72.0: 1.25,
                                   96.0: 1.125,
                                   120.0: 1.0}
        
        preburst_depth = preburst_proportion * 100      # as percentage
        
        preburst_duration = preburst_duration_scale[duration] * duration
        num_timesteps = int(np.ceil(preburst_duration/timesteps))
        pb_intensity = np.around(preburst_depth/num_timesteps, 2)
        
        uniform_preburst_pattern = pd.Series([pb_intensity] * num_timesteps)
        uniform_preburst_pattern.iloc[0] = preburst_depth - uniform_preburst_pattern.iloc[1:].sum()     # Correct for rounding errors
        index = pd.Index(np.linspace(-1 * num_timesteps * timesteps, -1 * timesteps, num_timesteps))
        uniform_preburst_pattern.index = index
        
        return uniform_preburst_pattern
    
    def filter_embedded_bursts_in_preburst(self, preburst_pattern0, main_burst, main_burst_duration,
                                           z, ave_rain, storm_method, buffer=0.9, messages=False, 
                                           climate_adjustment = None, cap_duration = True):
        preburst_pattern = preburst_pattern0.copy()
        z_2000 = ndtri(1 - 1 / 2000.0)
        
        data = {}
        for dur in self.preburst_envelope_durations:
            if z <= z_2000:
                storm_method_1 = storm_method
            else:
                duration_changeover = self.storm_method_config['gsdm_gtsmr_changover_duration']
                if dur <= duration_changeover[0]:
                    storm_method_1 = 'GSDM'
                elif dur >= duration_changeover[1]:
                    storm_method_1 = 'GTSMR'
                else:
                    storm_method_1 = storm_method
            
            embedded_depth = self.rainfall.get_depth_z(z=z, duration=dur, storm_method=storm_method_1, print_msg=False)
            embedded_ave = self.get_average_rain(embedded_depth, print_msg=False)
            
            if climate_adjustment is not None:
                embedded_ave = embedded_ave * climate_adjustment[dur]
            
            normalised = embedded_ave / ave_rain * 100
            data[dur] = normalised
            
        temppat_filter_curve = pd.Series(data)
        temppat_filter_curve.sort_index(ascending=True, inplace=True)
        nonstd_durations = scipy.interpolate.interp1d(np.log(temppat_filter_curve.index), 
                                                      np.log(temppat_filter_curve.to_numpy()), fill_value='extrapolate')
        
        if cap_duration:
            if temppat_filter_curve.index.max() > main_burst_duration:
                # check if the longest storm duration rainfall depth (as proportion of main burst)
                # is less than the preburst rainfall depth... then cajole some extrapolation if so
                # because we need to check that the preburst depth is not greater than a burst of
                # the same duration.
                normalised_full_storm_burst = preburst_pattern.sum() + 100
                if temppat_filter_curve.iloc[-1] < normalised_full_storm_burst:
                    # get a slice of the depths (relative to main burst) with durations longer than the main burst duration
                    curve = temppat_filter_curve.loc[main_burst_duration:]              # Part of curve exceeding main_burst_duration
                    slope = (np.log(curve) - np.log(100)) / (np.log(curve.index.to_numpy()) - np.log(main_burst_duration))
                    min_slope = slope.idxmin()
                    
                    new_curve = scipy.interpolate.interp1d(np.log([main_burst_duration, min_slope]),
                                                           np.log([100, temppat_filter_curve[min_slope]]), fill_value='extrapolate')
                    nonstd_durations = new_curve
                    test = np.around(np.exp(new_curve(np.log(curve.index.to_series()))),4) > np.around(curve,4)
                    if test.all():
                        print('WARNING: Enveloping bursts may exceed design burst AEP')
        
        overall_tally = 0
        counter = 0
        while True:
            preburst_pattern, tally, unchanged = self.filter_max_intensitiy_preburst(preburst_pattern, 
                                                                                     main_burst, main_burst_duration, nonstd_durations, 
                                                                                     redistribute=True, buffer=buffer, messages=messages)
            if unchanged:
                if counter == 0:
                    preburst_filter_comment = 'No embedded bursts;'
                else: 
                    preburst_filter_comment = f'Embedded bursts filtered ({counter} passes);'
                break
            overall_tally += tally
            counter += 1
            if counter > 10: 
                while True:
                    preburst_pattern, tally, unchanged = self.filter_max_intensitiy_preburst(preburst_pattern, 
                                                                                             main_burst, main_burst_duration, nonstd_durations, 
                                                                                             redistribute=False, buffer=buffer, messages=messages)
                    if unchanged:
                        break
                    overall_tally += tally
                    counter += 1
                    if counter > 10: 
                        break
                        # print(f'rain_z: {z}; tp: {main_burst.name}; pb_pat: {preburst_pattern.name}; storm: {storm_method}')
                        # raise Exception('preburst filter error')
                preburst_filter_comment = f'Embedded bursts filtered ({counter} passes);'
                break
            
        preburst_pattern, tally, unchanged = self.filter_enveloping_burst(preburst_pattern, main_burst_duration, nonstd_durations)
        overall_tally += tally
        
        counter = 0
        if unchanged:
            preburst_filter_comment += ' ; ;'
        else:
            preburst_filter_comment += ' Enveloping burst filtered;'
            while True:
                preburst_pattern, tally, unchanged = self.filter_max_intensitiy_preburst(preburst_pattern, 
                                                                                         main_burst, main_burst_duration, nonstd_durations, 
                                                                                         redistribute=False, buffer=buffer, messages=messages)
                if unchanged:
                    if counter > 0:
                        preburst_filter_comment += f'Re-filtered embedded bursts ({counter} passes);'
                    else:
                        preburst_filter_comment += ' ;'
                    break
                overall_tally += tally
                counter += 1
                if counter > 10: 
                    while True:
                        preburst_pattern, tally, unchanged = self.filter_max_intensitiy_preburst(preburst_pattern, 
                                                                                                 main_burst, main_burst_duration, nonstd_durations, 
                                                                                                 redistribute=False, buffer=buffer, messages=messages)
                        if unchanged:
                            break
                        overall_tally += tally
                        counter += 1
                        if counter > 10: 
                            break
                            # print(f'rain_z: {z}; tp: {main_burst.name}; pb_pat: {preburst_pattern.name}; storm: {storm_method}')
                            # raise Exception('preburst filter error')
                    preburst_filter_comment = f'Re-filtered embedded bursts ({counter} passes);'
                    break
        
        if overall_tally > 0:
            preburst_pattern, extension = self.prepend_excess_preburst(preburst_pattern, main_burst_duration, overall_tally, nonstd_durations)
            preburst_filter_comment += f'Storm extended {extension}h.'
        
        if np.abs(preburst_pattern.sum() - preburst_pattern0.sum()) > 2:    # 2%
            change = preburst_pattern.sum() - preburst_pattern0.sum()
            print('Pre-burst percentage (before filtering)', preburst_pattern0.sum())
            print('Pre-burst percentage (after filtering)', preburst_pattern.sum())
            preburst_filter_comment = f'WARNING: pre-burst proportion changed by {change}% ;' + preburst_filter_comment
        
        return preburst_pattern, preburst_filter_comment
    
    def filter_max_intensitiy_preburst(self, preburst_pattern0, main_burst, main_burst_duration, 
                                       nonstd_durations, redistribute = True, buffer = 0.9, messages=False):
        # main_burst = temporal_pattern[temporal_pattern.index > 0]
        # preburst_pattern = temporal_pattern.loc[:0].copy()
        preburst_pattern = preburst_pattern0.copy()
        unchanged = True
        tally = 0
        # Scan pre-burst peak intensities for embedded burst
        # max_window = min(len(main_burst), len(preburst_pattern))
        max_window = len(preburst_pattern)
        for i in range(1, max_window + 1):
            if i <= len(main_burst):
                limit = main_burst.rolling(i).sum().max() * buffer
            else:
                window_dur = i * self.timesteps
                limit = np.exp(nonstd_durations(np.log(window_dur)).item())
                
            if preburst_pattern.rolling(i).sum().max() > limit:
                unchanged = False
                
                if messages:
                    print(f'Found embedded {i * self.timesteps}h burst')
                
                if i == 1:
                    idx = preburst_pattern.rolling(i).sum().idxmax()    # get index 
                    pos_n = preburst_pattern.index.get_loc(idx)         # get index location
                    pos_0 = pos_n
                    
                    excess = preburst_pattern.max() - limit
                    preburst_pattern.loc[idx] = limit
                else:
                    idx = preburst_pattern.rolling(i).sum().idxmax()    # get index 
                    pos_n = preburst_pattern.index.get_loc(idx)         # location of index at end of window
                    pos_0 = pos_n -i +1                                 # location of index at start of window
                    
                    excess = preburst_pattern.rolling(i).sum().max() - limit
                    if preburst_pattern.iloc[pos_0] < preburst_pattern.iloc[pos_n]:
                        preburst_pattern.iloc[pos_0] -= excess
                        if np.abs(preburst_pattern.iloc[pos_0]) < 0.001:
                            preburst_pattern.iloc[pos_0] = 0
                        if np.around(preburst_pattern.iloc[pos_0], 3) < 0: 
                        # if preburst_pattern.iloc[pos_0] < 0: 
                            print('ERROR filtering pre-burst pattern')
                            print(f'preburst_proportion: {preburst_pattern0.sum()/100.0}')
                            print(preburst_pattern0)
                            print(main_burst)
                            raise Exception()
                    else:
                        preburst_pattern.iloc[pos_n] -= excess
                        if np.abs(preburst_pattern.iloc[pos_n]) < 0.001:
                            preburst_pattern.iloc[pos_n] = 0
                        if np.around(preburst_pattern.iloc[pos_n], 3) < 0: 
                        # if preburst_pattern.iloc[pos_n] < 0: 
                            print('ERROR filtering pre-burst pattern')
                            print(f'preburst_proportion: {preburst_pattern0.sum()/100.0}')
                            print(preburst_pattern0)
                            print(main_burst)
                            raise Exception()
                
                if messages:
                    print(excess)
                
                tally = excess
                prev_capacity = None
                for j in range(i+1, max_window + 1):
                    # expand window by one until max_window, or tally exhausted
                    
                    if pos_0 == 0: 
                        pos_n += 1
                        min_pos = pos_n
                    elif pos_n == max_window -1: 
                        pos_0 -= 1
                        min_pos = pos_0
                    else:
                        if preburst_pattern.iloc[pos_0 -1] >= preburst_pattern.iloc[pos_n +1]:
                            pos_0 -= 1
                            min_pos = pos_0
                        else:
                            pos_n += 1
                            min_pos = pos_n
                        
                    # limit = main_burst.rolling(j).sum().max() * buffer
                    if j <= len(main_burst):
                        limit = main_burst.rolling(j).sum().max() * buffer
                    else:
                        window_dur = j * self.timesteps
                        limit = np.exp(nonstd_durations(np.log(window_dur)).item())
                    
                    excess = preburst_pattern.iloc[pos_0:pos_n+1].sum() - limit
                    if messages:
                        print(f'Expanding window at position {preburst_pattern.index[min_pos]}')
                        print('Tally: ', tally)
                        print('Excess: ', excess)
                    if excess == 0:
                        continue
                    elif excess > 0:
                        tally += excess
                        preburst_pattern.iloc[min_pos] -= excess
                        if np.abs(preburst_pattern.iloc[min_pos]) < 0.001:
                            preburst_pattern.iloc[min_pos] = 0
                        if np.around(preburst_pattern.iloc[min_pos], 3) < 0: 
                        # if preburst_pattern.iloc[min_pos] < 0: 
                            print('ERROR filtering pre-burst pattern')
                            print(f'preburst_proportion: {preburst_pattern0.sum()/100.0}')
                            print(preburst_pattern0)
                            print(main_burst)
                            print(f'Value: {preburst_pattern.iloc[min_pos]}')
                            raise Exception()
                        prev_capacity = None
                    elif redistribute:
                        capacity = -1 * excess
                        # tail = pos_0 - max_window      # negative number
                        # preburst_tail = preburst_pattern.iloc[pos_0:].sum()
                        # mainburst_tail = main_burst.iloc[tail:].sum() * buffer
                        # tail_capacity = mainburst_tail - preburst_tail
                        # capacity = min(capacity, tail_capacity)
                        
                        window_dur = (len(preburst_pattern) - pos_0) * self.timesteps + main_burst_duration
                        limit = np.exp(nonstd_durations(np.log(window_dur)).item()) - 100.0     # Subtract normalised main burst (= 100)
                        preburst_tail = preburst_pattern.iloc[pos_0:].sum()
                        tail_capacity = limit - preburst_tail
                        
                        if prev_capacity is None:
                            capacity = min(capacity, tail_capacity)
                        else:
                            capacity = min(capacity, tail_capacity, prev_capacity)
                        prev_capacity = capacity
                        
                        # if prev_capacity:
                        #     capacity = min(capacity, prev_capacity)
                        # prev_capacity = capacity
                        # capacity = min(capacity, tail_capacity)
                        
                        if messages:
                            print('Capacity: ', capacity)
                        
                        if capacity > 0:
                            if capacity < tally:
                                preburst_pattern.iloc[min_pos] += capacity
                                tally -= capacity
                            else:
                                preburst_pattern.iloc[min_pos] += tally
                                tally = 0
                                break
                break
        return preburst_pattern, tally, unchanged
    
    def filter_enveloping_burst(self, preburst_pattern0, main_burst_duration, nonstd_durations): #, buffer=1.0): 
        preburst_pattern = preburst_pattern0.copy()
        buffer = 1.0        # TODO: Code needs to be changed if using buffer < 1 to avoid negative depths
        unchanged = True
        # preburst_pattern = temporal_pattern.loc[:0].copy()
        tally = 0
        for i in range(1, len(preburst_pattern) + 1):
            env = preburst_pattern.iloc[-i:].sum() + 100          # Normalised main_burst depth = 100;
            env_dur = main_burst_duration + i * self.timesteps
            limit = np.exp(nonstd_durations(np.log(env_dur)).item()) * buffer
            if env > limit:
                unchanged = False
                excess = env - limit
                tally += excess
                preburst_pattern.iloc[-i] -= excess
                if preburst_pattern.iloc[-i] < 0:
                    print(f'Filtering preburst: Excess {excess} mm | Available (mm)', preburst_pattern.iloc[-i])
                    print('Trying to remove more rainfall than is available.')
                    print('WARNING: The filtering may not have worked as intended')
                    tally += preburst_pattern.iloc[-i]
                    preburst_pattern.iloc[-i] = 0.0
                    # raise Exception()
            else:
                if tally > 0:
                    capacity = limit - env
                    if tally > capacity:
                        tally -= capacity
                        preburst_pattern.iloc[-i] += capacity
                    else:
                        tally = 0
                        preburst_pattern.iloc[-i] += tally
        return preburst_pattern, tally, unchanged
    
    def prepend_excess_preburst(self, preburst_pattern0, main_burst_duration, tally, nonstd_durations): #, buffer=1.0):
        preburst_pattern = preburst_pattern0.copy()
        buffer = 1.0        # TODO
        # preburst_pattern = temporal_pattern.loc[:0]
        normalised = preburst_pattern.sum() + 100               # TODO: What if mainburst = 99.99 or 100.01?
        target = normalised + tally
        
        idx = preburst_pattern.index.min()
        storm_dur = len(preburst_pattern) * self.timesteps + main_burst_duration

        prepend = {idx: normalised}
        while True:
            idx -= self.timesteps
            storm_dur += self.timesteps
            normalised = np.exp(nonstd_durations(np.log(storm_dur)).item()) * buffer
            
            if normalised < target:
                prepend[idx] = normalised
            else:
                prepend[idx] = target
                break
        
        prepend = pd.Series(prepend)            # Convert to series
        prepend = prepend.diff().iloc[1:]       # Convert to incremental and drop first value (NaN)
        prepend.sort_index(inplace=True)        # Sort 
        extension = preburst_pattern.index.min() - idx 
        preburst_pattern1 = pd.concat([prepend, preburst_pattern], axis=0)
        
        return preburst_pattern1, extension

    def store_hyetographs(self, sim_id, original_temporal_pattern, filtered_temporal_pattern):
        print('Patterns filtered... storing temporal patterns')
        original_temporal_pattern.name = sim_id
        filtered_temporal_pattern.name = sim_id
        self.is_stored_hyetographs = True
        self.original_temporal_patterns.append(original_temporal_pattern.copy())
        self.filtered_temporal_patterns.append(filtered_temporal_pattern.copy())

    def output_hyetographs(self, output_file):
        if self.is_stored_hyetographs:
            output_path = f'{output_file}_hyeto_original.csv'
            print('\nOutputting original hyetographs', output_path)
            pd.concat(self.original_temporal_patterns, axis=1).to_csv(output_path)
            output_path = f'{output_file}_hyeto_filtered.csv'
            print('Outputting filtered hyetographs', output_path)
            pd.concat(self.filtered_temporal_patterns, axis=1).to_csv(output_path)


class PreBursts:
    def __init__(self, input_file, source='ARR_DataHub'):
        self.percentiles = [10, 25, 50, 75, 90]
        self.prebursts_by_percentile = {}
        self.prebursts_by_duration = {}

        # Get the preburst proportions
        print(f'\nReading pre-burst from {input_file}')
        for percentile in self.percentiles:
            print(f'\nPre-burst proportions for {percentile} percentile')
            if source=='ARR_DataHub':
                df = self.read_preburst_proportions(input_file, percentile)
            elif source=='csv':
                df = self.read_preburst_proportions_from_csv(input_file, percentile)
            print(df)
            self.prebursts_by_percentile[percentile] = df
        self.rearrange_preburst_data()

    def get_preburst_proportion(self, rain_z, preburst_p, duration):
        preburst_p = preburst_p * 100
        # print(f'Interpolating pre-burst percentile {np.around(preburst_p, 1)}')
        # preburst_p = np.around(ndtr(preburst_z) * 100, 1)
        # print(f'Getting preburst for rain AEP of 1 in {rain_aep} and preburst percentile of {preburst_p}%')
        if duration in self.prebursts_by_duration.keys():
            all_df = self.prebursts_by_duration[duration]
        elif duration > max(self.prebursts_by_duration.keys()):
            max_dur = max(self.prebursts_by_duration.keys())
            all_df = self.prebursts_by_duration[max_dur]
            print(f'{duration}h pre-burst ratios not provided')
            print('Applying pre-burst ratios from duration:', max_dur, 'h')
        elif duration < 1:
            all_df = self.prebursts_by_duration[1]
        elif duration == 9:
            all_df = (self.prebursts_by_duration[6] + self.prebursts_by_duration[12]) / 2
        elif duration == 4.5:
            all_df = (self.prebursts_by_duration[3] + self.prebursts_by_duration[6]) / 2
        else:
            durations = np.array(list(self.prebursts_by_duration.keys()))
            print('Pre-burst ratios available for durations:', durations)
            idx = (np.abs(durations - duration)).argmin()
            apply_duration = durations[idx]
            all_df = self.prebursts_by_duration[apply_duration]
            print('Applying pre-burst ratios from duration:', apply_duration)
            # raise Exception('Duration not in the preburst database!')
        # print(all_df)

        all_preburst_p = all_df.columns
        all_rain_z = all_df.index

        # Interpolate pre-burst percentile
        # print(f'\nInterpolating burst percentile {preburst_z}')
        all_df = all_df.transpose()
        # min_preburst_z = all_preburst_z[0]
        # max_preburst_z = all_preburst_z[-1]
        #if preburst_z < min_preburst_z:
        #     preburst_z = min_preburst_z
        # elif preburst_z > max_preburst_z:
        #     preburst_z = max_preburst_z
        if preburst_p in all_preburst_p:
            df = all_df.loc[[preburst_p]]
        else:
            # x = all_df.index
            all_df.loc[preburst_p] = np.nan
            all_df.sort_index(axis=0, inplace=True)
            all_df.interpolate(method='index', inplace=True)
            # print(all_df)
            df = all_df.loc[[preburst_p]]
        df = df.transpose()
        # print(df)

        # Interpolate rain percentile
        # print(f'\nInterpolating rain {rain_z}')
        rain_min_z = all_rain_z[0]
        rain_max_z = all_rain_z[-1]
        if rain_z > rain_max_z:
            rain_z = rain_max_z
        elif rain_z < rain_min_z:
            rain_z = rain_min_z
        if rain_z in all_rain_z:
            preburst_proportion = df.loc[rain_z]
        else:
            df.loc[rain_z] = np.nan
            df.sort_index(axis=0, inplace=True)
            df.interpolate(method='cubic', inplace=True)
            # print(df)
            preburst_proportion = df.loc[rain_z]
            # print(f'Preburst proportion interpolation of {preburst_proportion.iloc[0]}')

        return max(0, preburst_proportion.iloc[0])

    def rearrange_preburst_data(self):
        # Rearrange the preburst rainfall with rows of AEP and columns of percentile.
        # One dataframe for each storm duration
        durations = self.prebursts_by_percentile[self.percentiles[0]]
        for duration in durations:
            duration_df = pd.DataFrame()
            for percentile in self.percentiles:
                df = self.prebursts_by_percentile[percentile]
                local = df[duration]
                duration_df = pd.concat([duration_df, local], axis=1)
            # percentile_z = ndtri(np.array(self.percentiles).astype(float) / 100)
            # print(percentile_z)
            duration_df.columns = self.percentiles  # was percentile_z
            # if not duration == 1.5:
            if duration.is_integer(): 
                duration = int(duration)  # convert the duration from float to int if possible
            duration_df['z'] = ndtri(1 - 1 / duration_df.index)
            duration_df.set_index('z', inplace=True)
            self.prebursts_by_duration[duration] = duration_df
            print (f'\nPre-burst proportion for {duration}h duration')
            print(duration_df)

    def read_preburst_proportions(self, arr_datahub_file, percentile):
        lines = self.get_prebursts_lines(arr_datahub_file, percentile)
        durations = []
        df = pd.DataFrame()
        for i, line in enumerate(lines):
            split_line = line.split()
            if i == 0:
                # Get the header
                aeps_percent = ((1 / np.array(split_line[1].split(',')[1:]).astype(int)) * 100).astype(int)
                # print(aeps_percent)
                df = pd.DataFrame(columns=aeps_percent)
            else:
                duration = float(split_line.pop(0)) / 60
                durations.append(duration)
                # print(split_line)
                proportions = []
                for element in split_line[1:]:
                    proportion = element.split(',')[0]
                    proportion = float(proportion[1:-1])
                    proportions.append(proportion)
                # print(proportions)
                df.loc[duration] = proportions
        df = df.transpose()
        return df

    def get_prebursts_lines(self, arr_datahub_file, percentile=50):
        if percentile == 50:
            percentile = ''
        preburst_lines = []
        with open(arr_datahub_file) as infile:
            copy = False
            for line in infile:
                if line.strip() == f"[PREBURST{percentile}]":
                    copy = True
                    continue
                elif line.strip() == f"[PREBURST{percentile}_META]":
                    copy = False
                    continue
                elif copy:
                    preburst_lines.append(line)
        # print(preburst_lines)
        return preburst_lines

    def read_preburst_proportions_from_csv(self, filepath, percentile):
        key_term = f'{percentile} percentile'
        data_lines = []
        with open(filepath) as f:
            switch = False
            row_counter = 0
            for line in f:
                if switch == False:
                    row_counter += 1
                    if key_term in line.strip().lower():
                        switch = True       # start reading on next loop
                        skiprows = row_counter
                        row_counter = 0
                else:
                    if 'percentile' in line.lower():
                        break               # stop reading
                    else:
                        row_counter += 1
        
        df = pd.read_csv(filepath, index_col=0, skiprows=skiprows, nrows=row_counter-1, dtype='float')
        df.dropna(axis=1, how='all', inplace=True)
        df.columns = np.array(df.columns, dtype='float')        # Convert columns from string to float
        
        return df


class ArealReduction:
    def __init__(self, region='', area=0.0, arr_datahub_file=None):
        self.region = region
        self.area = area
        if arr_datahub_file:
            self.get_arf_data_from_file(arr_datahub_file)

    def get_arf_data_from_file(self, arr_datahub_file):
        all_lines = []
        with open(arr_datahub_file) as infile:
            copy = False
            for line in infile:
                if line.strip() == "[LONGARF]":
                    copy = True
                    continue
                elif line.strip() == "[LONGARF_META]":
                    copy = False
                    continue
                elif copy:
                    all_lines.append(line)
        header = all_lines[0]
        split_header = header.split(',')
        self.region = split_header[1].strip()

    def get_areal_reduction_factor(self, duration, aep):
        # area is in km²
        if self.area <= 1:
            return 1.0
        
        duration = duration * 60  # convert from hours to minutes
        if aep == 1:
            aep = 1.582
        aep = 1 / aep  # convert from '1 in X' to %
        
        if self.area < 10:
            if duration <= 720:
                arf10 = self.short_duration_arf(duration, aep, area = 10)
            elif duration >= 1440:
                arf10 = self.long_duration_arf(duration, aep, area = 10)
            else:
                arf10 = self.mid_duration_arf(duration, aep, area = 10)
            arf = self.small_area_arf(arf10, self.area)
            return arf
        elif self.area <= 1000:
            if duration <= 720:
                arf = self.short_duration_arf(duration, aep, self.area)
            elif duration >= 1440:
                arf = self.long_duration_arf(duration, aep, self.area)
            else:
                arf = self.mid_duration_arf(duration, aep, self.area)
            return arf
        elif self.area <= 30000:
            if duration <= 720:
                print('WARNING: Generalised ARF method not applicable for catchment > 1,000 km² and duration < 12 hours')
                print('         Applying ARF for catchment area of 1000 km ²')
                print('         Be wary if the critical sorm duration is less than 12 hours!')
                # arf = np.nan
                arf = self.short_duration_arf(duration, aep, 1000)
            elif duration >= 1440:
                arf = self.long_duration_arf(duration, aep, self.area)
            else:
                arf = self.mid_duration_arf(duration, aep, self.area)
            return arf
        else:
            print('Generalised ARF method not applicable for catchment > 30,000 km²')
            return np.nan

    def mid_duration_arf(self, duration, aep, area):
        arf_12h = self.short_duration_arf(12 * 60, aep, area)
        arf_24h = self.long_duration_arf(24 * 60, aep, area)
        arf = arf_12h + (arf_24h - arf_12h) * (duration - 720) / 720
        return arf

    def long_duration_arf(self, duration, aep, area):
        coeff = self.get_coefficients()
        exp = coeff['i'] * area * duration / 1440
        arf_part_1 = 1 - coeff['a'] * (area**coeff['b'] - coeff['c'] * np.log10(duration)) * duration**-coeff['d']
        arf_part_2 = coeff['e'] * area**coeff['f'] * duration**coeff['g'] * (0.3 + np.log10(aep))
        arf_part_3 = coeff['h'] * 10**exp * (0.3 + np.log10(aep))
        arf = arf_part_1 + arf_part_2 + arf_part_3

        if arf > 1.0:
            arf = 1.0
        return arf

    def short_duration_arf(self, duration, aep, area):
        exp = -0.021 * (duration - 180)**2 / 1440
        arf_part_1 = 1-0.287*(area**0.265 - 0.439 * np.log10(duration)) * duration**-0.36
        arf_part_2 = 0.00226 * area**0.226 * duration**0.125 * (0.3 + np.log10(aep))
        arf_part_3 = 0.0141 * area**0.213 * 10**exp * (0.3 + np.log10(aep))
        arf = arf_part_1 + arf_part_2 + arf_part_3
        # print('\nARF part:')
        # print([arf_part_1, arf_part_2, arf_part_3])
        if arf > 1.0:
            arf = 1.0
        return arf
    
    def small_area_arf(arf10, area):
        # Equation 2.4.4 in ARR Book 2
        arf = 1 - 0.6614 * (1 - arf10) * (area**0.4 - 1)
        return arf
    
    def get_coefficients(self):
        if self.region == 'East Coast North':
            coefficients = {'a': 0.327,
                            'b': 0.241,
                            'c': 0.448,
                            'd': 0.36,
                            'e': 0.00096,
                            'f': 0.48,
                            'g': -0.21,
                            'h': 0.012,
                            'i': -0.0013}
        elif self.region =='Semi-arid Inland Queensland':
            coefficients = {'a': 0.159,
                            'b': 0.283,
                            'c': 0.25,
                            'd': 0.308,
                            'e': 0.00000073,
                            'f': 1.0,
                            'g': 0.039,
                            'h': 0.0,
                            'i': 0.0}
        elif self.region =='Northern Coastal':
            coefficients = {'a': 0.326,
                            'b': 0.223,
                            'c': 0.442,
                            'd': 0.323,
                            'e': 0.0013,
                            'f': 0.58,
                            'g': -0.374,
                            'h': 0.013,
                            'i': -0.0015}
        else:
            coefficients = {}
            print('Available regions are:')
            print('    East Coast North')
            print('    Semi-arid Inland Queensland')
            print('    Northern Coastal')
            raise Exception(f'The {self.region} region does not exist in the "areal reduction" class in "RainfallScheme.py"')
        return coefficients

