"""
This file contains two classes:
- ifdCurves: this class is used to manage the collection of ifd data, for which ifd objects are stored
  in the 'duration_curves' dictionary.
- DurationCurve: This is an object containing the ifd depths for a specific storm duration (depth vs aep)
"""
import json
import os
import numpy as np
import pandas as pd
from scipy.special import ndtri, ndtr
from lib.InterpolationCurves import *
import warnings


class ifdCurves:
    """
    This class is a container for the IFD curve objects. There is a separate IFD curve for each storm duration
    stored in dictionaries that use the duration as the key.
    """
    def __init__(self, json_file, areal_reduction=None):
        self.areal_reduction = None
        # open the config file and get contents
        f = open(json_file)
        self.config_data = json.load(f)
        f.close()
        self.folder = os.path.dirname(json_file)
        self.arr_ifd_curves = {}
        self.pmp = None
        self.gsdm_ifd_curves = {}
        self.gtsmr_ifd_curves = {}
        self.durations = []
        self.rare_aeps = []
        self.aep_upper = 0
        self.aep_lower = 0
        self.z_upper = 0
        self.z_lower = 0
        if 'extreme_rainfall_interpolation_method' in self.config_data.keys():
            self.interpolation_method = self.config_data['extreme_rainfall_interpolation_method']
        else:
            print('WARNING: The interpolation method to use for extreme rainfall has not ben provided.')
            print('         Defaulting to GEV')
            self.interpolation_method = "GEV"
        self.duration_range = []  # this is the range that will be simulated [shortest, longest]
        self.skip_method = 'none'
        if 'pmp_exceedance_sampling' in self.config_data.keys():
            self.pmp_exceedance_sampling = self.config_data['pmp_exceedance_sampling']
            print('Method for sampling rainfall exceeding the AEP of the PMP:', self.pmp_exceedance_sampling)
            options = ['cap', 'extrapolate']
            if self.pmp_exceedance_sampling not in options:
                print(f'WARNING: {self.pmp_exceedance_sampling} is not an option.')
                print(f'          Options are: ', options)
                raise Exception()
        else:
            self.pmp_exceedance_sampling = 'extrapolate'

    def get_aep_of_pmp(self):
        return self.config_data["AEP_of_PMP"]

    def get_rare_rainfall_depths(self, area=0):
        # Get the IFD depths for the rare events (more frequent than 1 in 2000 AEP)
        duration_data_set = self.config_data["durations"]
        for duration_data in duration_data_set:
            duration = duration_data['duration']

            # Do not import short duration rainfall (<= 12hrs) for large catchments (>1000km2).
            if duration <= 12 and area > 1000.0:
                print(f'WARNING: Catchment size is {area} km².')
                print('          No ARF method for catchment area > 1000 km² and duration <= 12 hours')
                print(f'         Importing {duration} hour duration regardless')
                # continue


            # duration_path = os.path.join(self.config_data['folder'], duration_data['filename'])
            duration_path = os.path.normpath(os.path.join(self.folder, duration_data['filename']))
            print(f'Importing {duration} hour duration: {duration_path}')
            duration_obj = DurationCurve(filepath=duration_path,
                                         duration=duration)
            duration_obj.interpret_data(self.config_data["col_suffix"])
            self.arr_ifd_curves[duration] = duration_obj

        self.durations = list(self.arr_ifd_curves.keys())
        print('\nFound the following durations:')
        print(self.durations)

        self.rare_aeps = np.sort(np.array(self.arr_ifd_curves[self.durations[0]].aeps).astype(float))
        print('\nFound the following AEPs:')
        print(self.rare_aeps)

        self.aep_upper = self.rare_aeps[-1]
        self.aep_lower = self.rare_aeps[0]

        self.z_upper = ndtri(1-1/self.aep_upper)
        self.z_lower = ndtri(1-1/self.aep_lower)

    def apply_areal_reduction(self, areal_reduction):
        # Apply the areal reduction - done before creating extreme depths!
        for duration in self.arr_ifd_curves:
            ifd_curve = self.arr_ifd_curves[duration]
            ifd_curve.apply_areal_reduction(areal_reduction)

    def set_up_extreme_rain(self):
        # Set up the PMP depths
        print('Getting the PMP depths')
        # pmp_files = [os.path.join(self.config_data['folder'], self.config_data['PMP_depths']),
        #              os.path.join(self.config_data['folder'], self.config_data['PMP_scaling'])]
        pmp_files = [
            os.path.normpath(
                os.path.join(
                    self.folder, 
                    self.config_data['PMP_depths'])),
            os.path.normpath(
                os.path.join(
                    self.folder, 
                    self.config_data['PMP_scaling']))
            ]
        pmp_obj = PMP(pmp_files, self.config_data["AEP_of_PMP"], self.skip_method)
        self.pmp = pmp_obj

        # Set up the extreme rainfall depths
        for duration in self.durations:
            if not self.skip_method == 'gsdm':
                if duration not in self.pmp.df_gsdm.index:
                    print(f'ERROR: {duration} hr duration PMP depth not found!')
                    Exception('Check the list of PMP depths and durations')
                print('GSDM rainfall for duration:', duration)
                # For GSDM:
                # gsdm_df = self.arr_ifd_curves[duration].df.loc[1000:2000].copy()
                gsdm_df = self.arr_ifd_curves[duration].df.loc[[1000, 2000]].copy()
                gsdm_df.loc[self.pmp.aep_of_pmp] = self.pmp.df_gsdm.loc[duration]
                gsdm_duration_obj = ExtremeDurationCurve(gsdm_df,
                                                         duration=duration,
                                                         interpolation_method=self.interpolation_method)
                self.gsdm_ifd_curves[duration] = gsdm_duration_obj
                print(gsdm_df)

            if not self.skip_method == 'gtsmr':
                # for GTSMR
                if duration not in self.pmp.df_gtsmr.index:
                    print(f'ERROR: {duration} hr duration PMP depth not found!')
                    Exception('Check the list of PMP depths and durations')
                print('GTSMR rainfall for duration:', duration)
                # gtsmr_df = self.arr_ifd_curves[duration].df.loc[1000:2000].copy()
                gtsmr_df = self.arr_ifd_curves[duration].df.loc[[1000, 2000]].copy()
                gtsmr_df.loc[self.pmp.aep_of_pmp] = self.pmp.df_gtsmr.loc[duration]
                gtsmr_duration_obj = ExtremeDurationCurve(gtsmr_df,
                                                          duration=duration,
                                                          interpolation_method=self.interpolation_method)
                self.gtsmr_ifd_curves[duration] = gtsmr_duration_obj
                print(gtsmr_df)

    def get_depth_z(self, z, duration, storm_method, print_msg=True):
        if self.pmp.z is not None:
            if z > self.pmp.z:
                if self.pmp_exceedance_sampling == 'cap':
                    print('Capping rainfall to the PMP')
                    z = self.pmp.z

        aep = np.around(1 / (1 - ndtr(z)), 3)
        z_2000 = ndtri(1 - 1 / 2000.0)
        z_pmp = ndtri(1 - 1 / self.get_aep_of_pmp())
        if print_msg: 
            print(f'\nInterpolating rainfall for 1 in {aep} and storm duration of {duration} hours')
        # Use ARR depths and spatial pattern if AEP rarer than 1 in 2,000
        # if aep <= 2000.0:
        if z <= z_2000:
            if duration in self.arr_ifd_curves.keys():
                duration_obj = self.arr_ifd_curves[duration]
            else:
                print(self.arr_ifd_curves.keys())
                raise Exception(f'Duration {duration} hours not found in the ARR ifd curves')
        # Otherwise use spatial pattern for whichever storm method was sampled
        elif storm_method == 'GSDM':
            if duration in self.gsdm_ifd_curves.keys():
                duration_obj = self.gsdm_ifd_curves[duration]
            else:
                raise Exception(f'Duration {duration} hours not found in the GSDM ifd curves')
        elif storm_method == 'GTSMR':
            if duration in self.gtsmr_ifd_curves.keys():
                duration_obj = self.gtsmr_ifd_curves[duration]
            else:
                raise Exception(f'Duration {duration} hours not found in the GTSMR ifd curves')
        else:
            raise Exception(f'The {storm_method} storm method is not an option!')
        depth = duration_obj.get_depth_z(z)
        return np.around(depth, 2)

    def get_depth_aep(self, aep, duration, storm_method, print_msg=True):
        if print_msg:
            print(f'\nInterpolating rainfall for 1 in {aep} and storm duration of {duration} hours')
        # Use ARR depths and spatial pattern regardless of whether AEP is rarer than 1 in 2,000
        # This is mainly for downstream hydrology generation of storm files
        if duration in self.arr_ifd_curves.keys():
            duration_obj = self.arr_ifd_curves[duration]
        else:
            print(self.arr_ifd_curves.keys())
            raise Exception(f'Duration {duration} hours not found in the ARR ifd curves')
        z = ndtri(1 - 1 / aep)
        depth = duration_obj.get_depth_z(z)
        return np.around(depth, 2)


class DurationCurve:
    """
    This class contains the IFD curves
    """
    def __init__(self, filepath, duration=0):
        self.df = pd.read_csv(filepath, index_col=0)
        self.df_average = None
        self.df_spatial_pattern = None
        self.aeps = []
        self.std_norm_vars = []
        self.log_norm_df = pd.DataFrame()
        # self.areal_reduction = areal_reduction
        self.duration = duration

    def apply_areal_reduction(self, areal_reduction):
        # apply areal reduction factors
        print(f'\nAreal reduction factors for {self.duration} hour duration:')
        arf_df = pd.DataFrame(index=self.df.index, columns=self.df.columns, dtype=float)
        for aep in self.aeps:
            arf_df.loc[aep] = areal_reduction.get_areal_reduction_factor(aep=aep, duration=self.duration)
        print(arf_df)
        self.df = self.df * arf_df
        print('\nReduced rainfall depths:')
        print(self.df)

        # recreate the log-normal transform
        log_df = np.log10(self.df)
        log_df['std_norm_var'] = self.std_norm_vars
        self.log_norm_df = log_df.set_index('std_norm_var')

    def interpret_data(self, col_suffix):
        # Get the AEPs from the column headers
        # and get the standard normal variates
        col_names = self.df.columns
        print('Rainfall column names:')
        print(col_names)
        for col_name in col_names:
            aep = float(col_name[:-len(col_suffix)])
            z = ndtri(1-1/aep)
            self.aeps.append(aep)
            self.std_norm_vars.append(z)

        # Set up the log-normal form of the rainfall for interpolation
        self.df = self.df.transpose()
        self.df['AEP'] = self.aeps
        self.df.set_index('AEP', inplace=True)

        print(f'\nRainfall depths for {self.duration} h duration:')
        print(self.df)

        log_df = np.log10(self.df)
        log_df['std_norm_var'] = self.std_norm_vars
        self.log_norm_df = log_df.set_index('std_norm_var')

        # Set up the catchment average rainfall and spatial pattern
        self.df_average = self.df.mean(axis=1)
        # number_of_catchments = self.df.shape[1]
        divisor = pd.DataFrame(index=self.df.index, columns=self.df.columns)
        for aep in self.df.index:
            divisor.loc[aep] = self.df_average[aep]  # * number_of_catchments
        self.df_spatial_pattern = self.df / divisor
        # print(f'\nSpatial pattern for {self.duration} h duration:')
        # print(self.df_spatial_pattern)

    def get_depth_z(self, z):
        all_depths = self.log_norm_df.copy()
        exactmatch = all_depths[all_depths.index == z]
        if exactmatch.empty:
            all_depths.loc[z] = np.nan
            all_depths.sort_index(inplace=True)
            all_depths.interpolate(method='cubic',
                                    inplace=True)
        depth = np.around(10 ** all_depths.loc[z], 1)
        return depth


class ExtremeDurationCurve:
    def __init__(self, df, duration=0, interpolation_method='SiriwardenaWeinmann1998'):
        print(f'Setting up extreme duration curve for {duration} hour duration using {interpolation_method} method')
        self.df = df
        self.aeps = list(df.index)
        self.std_norm_vars = []
        self.df_average = None
        # TODO: get the spatial patterns as parameters of the object
        self.df_gsdm_spatial_pattern = None
        self.df_gtsmr_spatial_pattern = None
        self.interpolation_curve = {}
        for i, catchment in enumerate(self.df.columns):
            if interpolation_method == 'SiriwardenaWeinmann1998' or interpolation_method == 'HillAndOthers2000':
                curve = CoercedQuadratic(method=interpolation_method)
            elif interpolation_method == 'GEV':
                curve = GEV()
            else:
                raise Exception(f'A {interpolation_method} interpolation method is not an option.')

            curve.setup_rainfall_boundaries(aep_2000=self.df.loc[2000, catchment],
                                            aep_1000=self.df.loc[1000, catchment],
                                            pmp_depth=self.df.iloc[-1, i],
                                            pmp_aep=self.df.index[-1])
            curve.fit_curve()
            if not curve.test_validity():
                # raise Exception(f'The {interpolation_method} is outside range of application for {duration}-hour duration!')
                warnings.warn(f'The {interpolation_method} is outside range of application for {duration}-hour duration!')
                # input('See above warning. Are you sure you want to continue? If so, hit enter.')

            self.interpolation_curve[catchment] = curve

        for aep in self.aeps:
            z = ndtri(1 - 1 / aep)
            self.std_norm_vars.append(z)

        self.duration = duration

        log_df = np.log10(self.df)
        log_df['std_norm_var'] = self.std_norm_vars
        self.log_norm_df = log_df.set_index('std_norm_var')

    def get_depth_z(self, z):
        aep = 1 / (1 - ndtr(z))
        # print(f'\nInterpolating for AEP of 1 in {aep}')
        catchments = self.df.columns
        # depths = pd.DataFrame(columns=catchments, index=['aep'])
        depths = {}         # Note: depths collated in dict, then converted to pd.Series. This returns pd.Series of type float. Previous code returned type object, which causes numpy issues
        for catchment in catchments:
            curve = self.interpolation_curve[catchment]
            depth = curve.get_quantile(aep)
            # depths.loc['aep', catchment] = depth
            depths[catchment] = depth
        # print(depths)
        # return depths.iloc[0]
        return pd.Series(depths)

    def get_depth_z_old(self, z):
        aep_of_pmp = self.aeps[-1]
        all_depths = self.log_norm_df.copy()
        all_depths.reset_index(inplace=True, drop=True)  # Integer index – easier to index

        Y = 1 / (1 - ndtr(z))
        logYpmp = np.log10(aep_of_pmp)

        depth = all_depths.apply(lambda x: interp_extreme_rain(x[0], x[1], x[2], logYpmp=logYpmp, Y=Y), axis=0)
        depth.name = z
        return depth


def interp_extreme_rain(logXy1, logXy2, logXpmp, logYpmp, Y):  # Method from Siriwardena and Weinmann (1998)

    logY1 = np.log10(1000)  # hard coded
    logY2 = np.log10(2000)  # hard coded
    logY = np.log10(Y)

    Zd = logYpmp - logY2
    Gy = (logY - logY2) / Zd

    Sgc = (logXy1 - logXy2) / (logXy2 * (logY1 - logY2))
    Sgap = (logXpmp / logXy2 - 1) / Zd

    Ry = 1 + Sgc * Zd * Gy + (Sgap - Sgc) * Zd * Gy * Gy
    exp = Ry * logXy2
    depth = round(10 ** exp, 1)
    return depth


class PMP:
    def __init__(self, ls_filepath=None, aep=None, skip_method='none'):
        self.df_gsdm = pd.DataFrame()
        self.df_gtsmr = pd.DataFrame()
        self.aep = aep
        if aep is not None:
            self.z = ndtri(1 - 1 / aep)
        else:
            self.z = None
        
        if ls_filepath:
            # pmp_by_duration = pd.read_csv(ls_filepath[0], header=9)  # TODO: header row hard-coded for sample data
            pmp_by_duration = pd.read_csv(ls_filepath[0])  # changed the format to exclude redundant headers
            spatial_pattern = pd.read_csv(ls_filepath[1], index_col=0)
            # self.aep_of_pmp = 0
            # self.std_norm_vars = []
    
            # tmp = pmp_by_duration['Duration'].apply(duration_to_int)  # Convert text duration labels to int
            # pmp_by_duration['Duration'] = tmp
            # pmp_by_duration = pmp_by_duration[['Duration', 'PMP']]
    
            # tmp = pmp_by_duration.apply(np.product, axis=1).apply(round)  # Convert intensities to depth
            # pmp_by_duration['PMP'] = tmp
            pmp_by_duration.set_index('Duration', inplace=True)
            pmp_by_duration = pmp_by_duration['PMP']  # Convert to series
    
            for duration in pmp_by_duration.index:
                if not skip_method == 'gsdm':
                    self.df_gsdm[duration] = pmp_by_duration[duration] * spatial_pattern['GSDM']
                if not skip_method == 'gtsmr':
                    self.df_gtsmr[duration] = pmp_by_duration[duration] * spatial_pattern['GTSMR']
            self.df_gsdm = self.df_gsdm.T
            self.df_gtsmr = self.df_gtsmr.T
        # Get the AEP if not provided
        if aep:
            self.aep_of_pmp = aep
        else:
            # Need the area then can use the set_aep_of_pmp function
            pass

    def set_aep_of_pmp(self, area):
        if area < 100:  # 100 km²
            aep_of_pmp = 1e-7
        elif area > 1e5:  # 100,000 km²
            aep_of_pmp = 1e-4
        else:
            aep_of_pmp = 10 ** (np.log10(area) - 9)

        aep_of_pmp = 1 / aep_of_pmp  # Express as 1 in Y
        self.aep_of_pmp = aep_of_pmp


def duration_to_int(dur):  # e.g. '24 hour'
    tmp = dur.split(' ')
    if tmp[1] != 'hour':
        raise Exception('PMP file not correctly formatted')
    value = float(tmp[0])
    if value.is_integer():
        return int(value)  # duration as int
    else:
        return float(value)


