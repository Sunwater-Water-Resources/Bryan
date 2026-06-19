"""
This class handles the antecedent lake conditions
"""
import os
import json
import math
from scipy import interpolate
import numpy as np
import pandas as pd


class LakeConditions:
    def __init__(self, parameters):
        # open the config file and get contents
        # Get the lake volume from the sims list
        self.is_correlated = False
        self.antecedent_type = None
        self.antecedent_volume = None
        self.correlation = None
        self.antecedent_volume_curves = None
        self.storage_curve = pd.DataFrame()
        self.full_supply_volume = None
        self.volume_cap = 'fsv'  # cap the lake level to the full supply volume by default

        if 'ADV' in parameters.keys():
            adv_type = parameters['ADV']
            # print(f'Using a ADV type of {adv_type}')
            if parameters.isna()['ADV']:
                print()
                print(parameters)
                raise Exception('ADV column in the sims list cannot be blank')
            else: 
                if isinstance(adv_type, (int, float, complex)):
                    self.antecedent_type = 'fixed'
                    self.antecedent_volume = adv_type
                    print(f'The antecedent dam volume is fixed at {self.antecedent_volume} ML!')
                elif str(adv_type).lower() == 'fsv':
                    print(f'The antecedent dam volume will be set to FSV')
                    self.antecedent_type = 'fsv'
                elif str(adv_type).lower() == 'varying':
                    if parameters['Method'] == 'ensemble':
                        raise Exception('Cannot use varying ADV with the ensemble method - only for Monte Carlo!')
                    else:
                        print(f'The antecedent dam volume will be varying')
                        self.antecedent_type = 'varying'
                        lake_config = parameters.get('Lake config')
                        if lake_config:
                            self.get_config_info(lake_config)
                        else:
                            raise Exception('It looks like a lake config filepath has not been given - check the sims list!')
                else:
                    raise Exception(
                        'The ADV of {} is not numeric or the keyword "fsv" or "varying"'.format(adv_type))
            # else:
            #     raise Exception('An ADV has not been specified in the Simulation list!')

    def get_config_info(self, json_file):
        # open the config file and get contents
        print(f'\nSetting up the lake conditions: {json_file}')
        f = open(json_file)
        config_data = json.load(f)
        f.close()

        # check for and set up any correlations
        if self.antecedent_type == 'varying':
            if 'correlation_layer_info' in config_data.keys():
                self.is_correlated = True
            if self.is_correlated:
                self.correlation = Correlator(config_data['correlation_layer_info'])

            # set up the volume exceedance curve
            self.antecedent_volume_curves = VolumeExceedanceCurve(config_data['exceedance_layer_info'],
                                                                  os.path.dirname(json_file))

            # set up the volume cap (either fsv, none, or ceiling)
            # fsv -- caps the ADV to  the FSV
            # none -- no capping applied
            # ceiling -- capped at the ceiling in the sigmoid curve
            if 'volume_cap' in config_data.keys():
                cap_options = ['fsv', 'none', 'ceiling']
                volume_cap = config_data['volume_cap'].lower()
                print(f'User has specified {volume_cap} method for capping of lake volume')
                if volume_cap in cap_options:
                    self.volume_cap = volume_cap
                else:
                    raise Exception(f'{volume_cap} is not an option: {cap_options}')

        # set up the storage elevation curve - this is not used anymore\\
        if "lake_storage_curve" in config_data.keys():
            csv_file = os.path.join(os.path.dirname(json_file),
                                    config_data['lake_storage_curve'])
            csv_file = os.path.normpath(csv_file)
            self.storage_curve = StorageCurve(csv_file)

    def set_full_supply_volume(self, fsv):
        self.full_supply_volume = fsv
        if self.antecedent_type == 'fsv':
            self.antecedent_volume = fsv

    def get_lake_volume(self, lake_z):
        # this is not used anymore in the URBS model
        if self.antecedent_type == 'varying':
            volume = self.antecedent_volume_curves.get_lake_volume(lake_z, self.volume_cap)
        else:
            volume = self.antecedent_volume
        # cap the volume if specified
        if self.volume_cap == 'fsv':
            if self.full_supply_volume:
                if volume > self.full_supply_volume:
                    volume = self.full_supply_volume
        return volume

    def get_correlated_z(self, rain_z, lake_z):
        print('\nApplying correlation for lake initial conditions...')
        new_z = self.correlation.apply_correlations(rain_z, lake_z)
        return new_z

    def get_level_from_volume(self, volume):
        level = self.storage_curve.get_level_from_volume(volume)
        return level

    def get_volume_below_fsl(self):
        if self.antecedent_type == 'fsv':
            volume_below_fsl = 0.0
        else:
            if self.full_supply_volume > 0:
                volume_below_fsl = self.full_supply_volume - self.antecedent_volume
                # if volume_below_fsl < 0:
                #     raise Exception(f'Volume below full supply level is negative: {volume_below_fsl} ML')
            else:
                print('Full supply volume not given or negative')
                print('Setting Volume below full supply level to zero')
                volume_below_fsl = 0
        return volume_below_fsl


class StorageCurve:
    def __init__(self, storage_file):
        print('\nReading storage file: ', storage_file)
        # open the data but also sort out some formatting issues - numbers in as text
        df = pd.read_csv(storage_file, dtype=float, thousands=',', index_col=0)
        # get rid of trailing empty lines
        df.dropna(inplace=True, axis=0, how='any')
        # get rid of any duplicate data
        df.drop_duplicates(keep='last', inplace=True)
        print(df)
        x = df[df.columns[0]]
        y = df.index
        self.curve = interpolate.interp1d(x=x, y=y, kind='cubic')

    def get_level_from_volume(self, volume):
        level = self.curve(volume)
        return np.around(level, 3)


class VolumeExceedanceCurve:
    def __init__(self, layer_info, folder):
        if isinstance(layer_info, list):
            self.number_of_layers = len(layer_info)
            self.layers = layer_info
        else:
            self.number_of_layers = 1
            self.layers = [layer_info]
        print(f'Found {self.number_of_layers} volume exceedance curve layers.')
        self.curves = []
        for layer in range(self.number_of_layers):
            self.curves.append(ExceedanceCurveLayer(self.layers[layer], folder))

    def get_lake_volume(self, lake_z, volume_cap='none'):
        print(f'\nFor lake z of {lake_z} seeking curve type')
        volume = 0.0
        for i in range(self.number_of_layers):
            layer = self.curves[i]
            if layer.test_z_bounds(lake_z):
                print(f'Success with {layer.type} volume exceedance curve')
                volume = layer.get_lake_volume(lake_z, volume_cap)
                print(f'Using ADV of {volume} ML')
                break
        return volume


class ExceedanceCurveLayer:
    def __init__(self, layer_info, folder):
        if 'lower_z' in layer_info:
            self.lower_z = layer_info['lower_z']
        else:
            self.lower_z = -99
        if 'upper_z' in layer_info:
            self.upper_z = layer_info['upper_z']
        else:
            self.upper_z = 99
        self.type = layer_info['type']
        if self.type == 'uniform':
            if 'value_ML' in layer_info.keys():
                self.value_ML = layer_info['value_ML']
                print(f'Using uniform ADV of {self.value_ML} ML')
            else:
                raise Exception('"value_ML" not given in lake config file for type: ', self.type)
        elif self.type == 'sigmoid':
            if 'coefficients' in layer_info.keys():
                self.coefficients = layer_info['coefficients']
                print('Using sigmoid curve for ADV with coefficients:')
                print(self.coefficients)
            else:
                raise Exception('"coefficients" not given in lake config file for type: ', self.type)
        elif self.type == 'empirical':
            if 'filename' in layer_info.keys():
                filepath = os.path.join(folder, layer_info['filename'])
                print('Opening empirical ADV curve:', filepath)
                df = pd.read_csv(filepath)
                print(df)
                try:
                    df['log_volume'] = np.log10(df['volume'])
                except Exception as e:
                    e.add_note('Issue with Log of ADV - perhaps there is a negative or zero volume')
                    raise
                print(df)
                x = df['z'].to_numpy()
                y = df['log_volume'].to_numpy()
                try:
                    self.curve = interpolate.interp1d(x=x, y=y, kind='linear')
                except Exception as e:
                    e.add_note('Could not set up the ADV interpolation object!')
                    e.add_note('Perhaps the x axis (std normal variate) is not increasing or there is a duplicate.')
                    raise
            else:
                raise Exception('"filename" not given in lake config file for type: ', self.type)

    def test_z_bounds(self, lake_z):
        condition = (self.lower_z < lake_z) & (lake_z < self.upper_z)
        return condition

    def get_lake_volume(self, lake_z, volume_cap):
        if self.type == 'uniform':
            return self.value_ML

        elif self.type == 'sigmoid':
            # set the coefficients
            k = self.coefficients['k']
            Vf = self.coefficients['Vf']
            Vc = self.coefficients['Vc']
            H = self.coefficients['H']
            z0 = self.coefficients['z0']
            log_vol = np.log10(Vc) / (H + np.exp(-k * (lake_z - z0))) + np.log10(Vf)
            volume = np.around(10 ** log_vol, 2)
            if volume_cap == 'ceiling':
                if volume > Vc:
                    volume = Vc
            return volume

        elif self.type == 'empirical':
            try:
                log_volume = self.curve(lake_z)
                adv_volume = 10 ** log_volume
                print(f'Empirical ADV interpolation of {log_volume} log volume or {adv_volume} ML')
            except Exception as e:
                e.add_note(f'Trying to interpolate the lake z value of {lake_z}')
                e.add_note('Something went wrong with the interpolation!')
                raise
            return adv_volume
        else:
            raise Exception(f'Antecedent dam volume of type {self.type} is not supported.')


class Correlator:
    def __init__(self, layer_info):
        if isinstance(layer_info, list):
            self.number_of_layers = len(layer_info)
            self.layers = layer_info
        else:
            self.number_of_layers = 1
            self.layers = [layer_info]
        print(f'Found {self.number_of_layers} correlation layers.')

    def apply_correlations(self, rain_z, lake_z):
        for layer in range(self.number_of_layers):
            layer_info = self.layers[layer]
            min = layer_info['lower_z']
            max = layer_info['upper_z']
            correlation = layer_info['correlation']
            # print(f'Layer {layer + 1} with min z of {min} and max z of {max}')
            condition = (min < rain_z) & (rain_z < max)
            if condition:
                if correlation > 0.01:
                    # do the correlation
                    lake_z = correlation * rain_z + lake_z * math.sqrt(1 - correlation**2)

        return lake_z


