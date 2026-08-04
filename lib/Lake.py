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
    def __init__(self, parameters=None):
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

        if parameters is None:
            # Built straight from a lake config file - see from_lake_config()
            return

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
                elif str(adv_type).lower() == 'mav':
                    print('The antecedent dam volume will be the MAV (the volume at z = 0)')
                    self.antecedent_type = 'mav'
                    self.get_config_info(self._lake_config_path(parameters, 'mav'))
                elif str(adv_type).lower() == 'varying':
                    if parameters['Method'] == 'ensemble':
                        raise Exception('Cannot use varying ADV with the ensemble method - only for Monte Carlo!')
                    else:
                        print(f'The antecedent dam volume will be varying')
                        self.antecedent_type = 'varying'
                        self.get_config_info(self._lake_config_path(parameters, 'varying'))
                else:
                    raise Exception(
                        'The ADV of {} is not numeric or one of the keywords '
                        '"fsv", "mav" or "varying"'.format(adv_type))
            # else:
            #     raise Exception('An ADV has not been specified in the Simulation list!')

    @staticmethod
    def _lake_config_path(parameters, adv_type):
        """The 'Lake config' filepath from the sims list row, for the ADV keywords
        that need the exceedance curve.

        A blank cell in the sims list comes back from pandas as NaN, which is
        truthy - so a plain 'if lake_config:' lets an empty cell through and the
        failure surfaces much later as open() being handed a float.
        """
        lake_config = parameters.get('Lake config')
        if lake_config is None or pd.isna(lake_config) or not str(lake_config).strip():
            raise Exception(f'An ADV of "{adv_type}" needs a "Lake config" filepath in the '
                            f'sims list - the volume comes from its exceedance curve.')
        return str(lake_config)

    @classmethod
    def from_lake_config(cls, json_file):
        """Build the varying-ADV machinery straight from a lake config file.

        The normal constructor takes its instructions from the sims list row (the
        ADV keyword). This alternate constructor is for callers that already know
        they want the varying distribution - e.g. the reservoir routing method
        re-sampling the ADV from the lake_z column of an existing mcdf file.
        """
        lake = cls()
        lake.antecedent_type = 'varying'
        lake.get_config_info(json_file)
        return lake

    def get_config_info(self, json_file):
        # open the config file and get contents
        print(f'\nSetting up the lake conditions: {json_file}')
        f = open(json_file)
        config_data = json.load(f)
        f.close()

        # check for and set up any correlations
        # 'mav' needs the exceedance curve and the cap, but not the correlation:
        # it reads one point off the curve rather than sampling it, and a
        # correlation against the rainfall has nothing to act on.
        if self.antecedent_type in ('varying', 'mav'):
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
        elif self.antecedent_type == 'mav':
            self.antecedent_volume = self._resolve_mav()

    def _resolve_mav(self):
        """The MAV: the antecedent volume the lake config's exceedance curve gives
        at z = 0, i.e. the median of the distribution a 'varying' ADV samples from.

        Resolved here rather than in the constructor because the fsv cap has to be
        applied against the full supply volume of the dam actually being run, which
        is not known until the model (or the rating curve, for reservoir routing)
        has been loaded. Capping and layer selection go through the same code a
        sampled z would, so the MAV is exactly the z = 0 member of that sample.
        """
        volume = self.antecedent_volume_curves.get_lake_volume(0.0, self.volume_cap)
        if self.volume_cap == 'fsv' and self.full_supply_volume:
            volume = min(volume, self.full_supply_volume)
        print(f'The MAV (volume at z = 0) is {volume:.1f} ML')
        return volume

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

    def get_lake_volumes(self, lake_z):
        """Array version of get_lake_volume - one volume (ML) per sampled z.

        Used where the whole sample is available up front (reservoir routing),
        rather than one realisation at a time as in the Monte Carlo loop.
        """
        lake_z = np.asarray(lake_z, dtype=float)
        if self.antecedent_type == 'varying':
            volumes = self.antecedent_volume_curves.get_lake_volumes(lake_z, self.volume_cap)
        else:
            volumes = np.full(lake_z.shape, float(self.antecedent_volume))
        # cap the volume if specified
        if self.volume_cap == 'fsv':
            if self.full_supply_volume:
                volumes = np.minimum(volumes, self.full_supply_volume)
        return volumes

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

    def get_lake_volumes(self, lake_z, volume_cap='none'):
        """Array version of get_lake_volume. Each z is taken by the first layer
        whose bounds contain it, matching the scalar version's break."""
        lake_z = np.asarray(lake_z, dtype=float)
        volumes = np.full(lake_z.shape, np.nan)
        unassigned = np.ones(lake_z.shape, dtype=bool)
        for i in range(self.number_of_layers):
            layer = self.curves[i]
            in_layer = np.asarray(layer.test_z_bounds(lake_z)) & unassigned
            if not in_layer.any():
                continue
            print(f'{in_layer.sum()} z values taken by the {layer.type} '
                  f'volume exceedance curve ({layer.lower_z} < z < {layer.upper_z})')
            volumes[in_layer] = layer.get_lake_volumes(lake_z[in_layer], volume_cap)
            unassigned &= ~in_layer
        if unassigned.any():
            # The scalar version quietly returns 0 ML here, which is far too easy
            # to miss when it happens across a whole sample.
            missed = lake_z[unassigned]
            raise Exception(
                f'{unassigned.sum()} lake z values fall outside every exceedance '
                f'layer in the lake config file (z from {missed.min()} to '
                f'{missed.max()}) - check the lower_z/upper_z bounds.'
            )
        return volumes


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
                self._resolve_sigmoid()
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

    def _resolve_sigmoid(self):
        """Work out the numerator A and height factor H of the sigmoid ADV curve

            log10 V(z) = A / (H + exp(-k (z - z0))) + log10(Vf)

        so the curve runs from Vf at low z up to 10**(A/H + log10(Vf)) at high z.

        TWO CONVENTIONS FOR A ARE IN CIRCULATION. A is not usually written in the
        config file, so H decides which one applies:

        H = 1, or H not given (the standard logistic - newer, fitted curves)
            A = log10(Vc) - log10(Vf), computed from the Vf and Vc already in the
            file. The upper asymptote is then exactly Vc and the curve sits
            halfway between Vf and Vc, in log space, at z0.

        H > 1 (the original ADV assessment workbooks - hand-tuned curves)
            A = log10(Vc). The "height adjustment factor" H in those workbooks was
            turned by hand until the curve passed through a chosen midpoint volume,
            which converges on H = log10(Vc)/(log10(Vc) - log10(Vf)). In these
            files Vc is really 10**A, not the ceiling: the two differ by a few
            percent. Read this way, such a model keeps exactly the curve it has
            always had.

        The split is safe because H = 1 cannot be a deliberate hand-tuned curve:
        in that convention it puts the upper asymptote at Vc * Vf, which needs a
        floor of 1 ML to mean anything. Any realistic floor gives an H well above
        1 (a 10,000 ML floor under a 140,000 ML ceiling gives 4.6). At Callide an
        H = 1 file read the old way asymptoted at 1.5 billion ML against a full
        supply volume of 129,041 ML.

        An explicit "A" in the coefficients overrides all of this, for a curve
        that is neither convention.
        """
        Vf = float(self.coefficients['Vf'])
        Vc = float(self.coefficients['Vc'])
        # H is the height adjustment factor of the older workbooks; a curve fitted
        # as a standard logistic has no use for it, so it may be left out.
        H = float(self.coefficients.get('H', 1.0))
        explicit_A = 'A' in self.coefficients
        if explicit_A:
            A = float(self.coefficients['A'])
            print(f'  A = {A:.6f} taken from the config file')
        elif abs(H - 1.0) < 1e-6:
            A = float(np.log10(Vc) - np.log10(Vf))
            print(f'  A = log10(Vc) - log10(Vf) = {A:.6f} - H is 1, so the curve is '
                  f'read as the standard logistic between Vf and Vc')
        else:
            A = float(np.log10(Vc))
            print(f'  A = log10(Vc) = {A:.6f} - H is not 1, so the curve is read in '
                  f'the original (hand-tuned H) convention')
        ceiling = 10.0 ** (A / H + np.log10(Vf))
        print(f'  ADV asymptotes: {Vf:,.0f} ML at low z, {ceiling:,.0f} ML at high z')
        if ceiling > 2.0 * Vc:
            logistic_A = np.log10(Vc) - np.log10(Vf)
            remedy = (f'the "A" of {A:.6f} in the config file is not the {logistic_A:.6f} '
                      f'(= log10(Vc) - log10(Vf))\n      a standard logistic needs - remove '
                      f'"A" to have it computed'
                      if explicit_A else
                      f'check H = {H}. For a curve fitted as a standard logistic '
                      f'(ceiling = Vc) set H to 1')
            print(f'  *** WARNING: the ADV ceiling of {ceiling:,.0f} ML is far above '
                  f'Vc = {Vc:,.0f} ML.\n'
                  f'      If this curve was meant to asymptote to Vc, {remedy}.')
        self.A, self.H = A, H

    def test_z_bounds(self, lake_z):
        condition = (self.lower_z < lake_z) & (lake_z < self.upper_z)
        return condition

    def get_lake_volume(self, lake_z, volume_cap):
        if self.type == 'uniform':
            return self.value_ML

        volume = float(self.get_lake_volumes(np.array([lake_z], dtype=float), volume_cap)[0])
        if self.type == 'empirical':
            print(f'Empirical ADV interpolation of {np.log10(volume)} log volume or {volume} ML')
        return volume

    def get_lake_volumes(self, lake_z, volume_cap):
        """Array version of get_lake_volume - the volume (ML) for each z."""
        lake_z = np.asarray(lake_z, dtype=float)
        if self.type == 'uniform':
            return np.full(lake_z.shape, float(self.value_ML))

        elif self.type == 'sigmoid':
            # set the coefficients. A and H are resolved once at load time and
            # carry the file's convention with them - see _resolve_sigmoid.
            k = self.coefficients['k']
            Vf = self.coefficients['Vf']
            Vc = self.coefficients['Vc']
            z0 = self.coefficients['z0']
            log_vol = self.A / (self.H + np.exp(-k * (lake_z - z0))) + np.log10(Vf)
            volume = np.around(10 ** log_vol, 2)
            if volume_cap == 'ceiling':
                volume = np.minimum(volume, Vc)
            return volume

        elif self.type == 'empirical':
            try:
                log_volume = self.curve(lake_z)
                adv_volume = 10 ** log_volume
            except Exception as e:
                e.add_note(f'Trying to interpolate lake z values from '
                           f'{lake_z.min()} to {lake_z.max()}')
                e.add_note('Something went wrong with the interpolation!')
                e.add_note('Perhaps a z value falls outside the range in the empirical curve.')
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

    def apply_correlations_to_array(self, rain_z, lake_z):
        """Array version of apply_correlations - correlates each lake z with its
        own rainfall z, using the layer that contains that rainfall z."""
        rain_z = np.asarray(rain_z, dtype=float)
        lake_z = np.array(lake_z, dtype=float)  # copy - not modified in place
        for layer in range(self.number_of_layers):
            layer_info = self.layers[layer]
            min = layer_info['lower_z']
            max = layer_info['upper_z']
            correlation = layer_info['correlation']
            condition = (min < rain_z) & (rain_z < max)
            if correlation > 0.01:
                lake_z[condition] = (correlation * rain_z[condition]
                                     + lake_z[condition] * math.sqrt(1 - correlation ** 2))
                print(f'Correlation of {correlation} applied to {condition.sum()} samples '
                      f'with rainfall z between {min} and {max}')
        return lake_z


