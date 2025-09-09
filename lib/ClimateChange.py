"""
Uplift rainfall to account for warmer climate projections

"""

import json
import os.path
import pandas as pd
import numpy as np


class ClimateAdjustment:
    def __init__(self, config_file, method, year=1965, ssp='SSP2-4.5', gwl=0.0):
        # Configure adjustment parameters
        self.method = method
        self.uplift_rate = 8  # % per °C
        ssp_scenarios = ['SSP1-1.9', 'SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
        if not ssp in ssp_scenarios:
            raise Exception(f'Climate SSP named {ssp} is not recognised!')
        self.config_data = {}

        # open the config file and get contents
        f = open(config_file)
        self.config_data = json.load(f)
        f.close()
        self.folder = os.path.dirname(config_file)

        # get the gwl (i.e. the temperature rise)
        self.adjust_climate = False
        self.temperature_rise = 0.0
        if method == 'ssp':
            self.get_gwl_for_ssp_method(year, ssp)
        else:
            if gwl > 0.0:
                print(f'A global warming level of {gwl}°C has been applied.')
                self.adjust_climate = True
            else:
                print(f'No global warming has been applied.')
            self.temperature_rise = gwl

        # set up change in rainfall losses
        self.nrm_cluster = self.config_data['NRM cluster']
        filepath = os.path.join(self.folder, self.config_data['loss rates file'])
        print('\nReading rainfall loss rates of change file: ', filepath)
        loss_rates_all = pd.read_csv(filepath, index_col=0)
        self.loss_rates = loss_rates_all.loc[self.nrm_cluster]

        # compute change in temporal pattern
        self.kg_classification = self.config_data['KG classification']
        filepath = os.path.join(self.folder, self.config_data['temporal pattern scaling'])
        print('\nReading temporal pattern D50 scaling file: ', filepath)
        self.tp_d50_scaling_df = pd.read_csv(filepath, index_col=0)
        self.tp_d50_scaling_df.fillna(0.0, inplace=True)
        if not self.kg_classification in self.tp_d50_scaling_df.columns:
            raise Exception(f'Koppen-Geiger classification {self.kg_classification} not found!')

    def get_gwl_for_ssp_method(self, year, ssp):
        # get the file paths for the temperature data
        temperature_files = self.config_data['files']

        # determine if climate adjustments needed
        base_period = self.config_data['base period']
        if year > base_period[1]:
            self.adjust_climate = True
        else:
            self.adjust_climate = False

        # read historic temperatures
        filepath = os.path.join(self.folder, temperature_files['historical'])
        print('\nReading historic temperatures: ', filepath)
        historic_temps = pd.read_csv(filepath, index_col=0)
        historic_temps.index = historic_temps.index.astype(int)

        # get the base temperature and adjust historic temps
        start_year = base_period[0]
        end_year = base_period[1]
        base_temp_df = historic_temps.loc[start_year: end_year]
        base_temperature = np.around(base_temp_df['Mean'].mean(), 3)
        print(f'\nThe temperature rise for the base period from {start_year} to {end_year} is {base_temperature}°C')
        historic_temps['adjusted'] = historic_temps['Mean'] - base_temperature

        # read projected temperatures
        filepath = os.path.join(self.folder, temperature_files[ssp])
        print(f'\nReading temperatures for {ssp}: ', filepath)
        projected_temps = pd.read_csv(filepath, index_col=0)
        projected_temps.index = projected_temps.index.astype(int)
        projected_temps['adjusted'] = projected_temps['Mean'] - base_temperature
        print(projected_temps)

        # get the temperature rise
        self.temperature_rise = 0
        if self.adjust_climate:
            historic_end = historic_temps.index[-1]
            projected_end = projected_temps.index[-1]
            if year <= historic_end:
                print(f'The year {year} falls within the historic record (till {historic_end})')
                print('The ten year rolling average temperature rise will be applied!')
                ave_historic = historic_temps.rolling(10).mean()
                year_temperature = ave_historic.loc[year, 'adjusted']
            else:
                print(f'The year {year} falls within the projected record (till {projected_end})')
                if year > projected_end:
                    raise Exception('The provided year is larger than in the projected record!')
                year_temperature = projected_temps.loc[year, 'adjusted']
            self.temperature_rise = np.around(year_temperature, 3)
            if self.temperature_rise < 0.1:
                print(f'\nNo climate adjustment needed for year: {year} | temperature rise is {self.temperature_rise}°C')
                self.adjust_climate = False
            else:
                print(f'\nThe temperature rise is {self.temperature_rise}°C')
        else:
            print(f'\nNo climate adjustment needed for year: {year}')


    def check_tp_d50_scaling(self, duration):
        durations = self.tp_d50_scaling_df.index
        if not duration in durations:
            # print(self.tp_d50_scaling_df)
            print('Temporal pattern scaling of D50 not provided in the standard table - interpolating')
            self.tp_d50_scaling_df.loc[duration] = np.nan
            self.tp_d50_scaling_df.sort_index(inplace=True)
            self.tp_d50_scaling_df.interpolate(method='slinear',
                                               inplace=True)
            print(self.tp_d50_scaling_df)

    def get_d50_shift_weightings(self, duration, delta_d50, pattern_area):
        all_weightings = {}
        weighting_files = self.config_data['weighting files']
        for key, file in weighting_files.items():
            if delta_d50 == 0:
                if key =='ARR point':
                    weighting = {'weightings': {'frequent': 0.1, 'intermediate': 0.1, 'rare': 0.1},
                                 'front_patterns': [1, 2, 3]}
                else:
                    weighting = {'weightings': 0.1, 'front_patterns': [1, 2, 3]}
            else:
                filepath = os.path.join(self.folder, file)
                weighting_obj = D50Weighting(filepath=filepath,
                                             pattern_type=key,
                                             duration=duration,
                                             pattern_area=pattern_area)
                weighting = weighting_obj.get_d50_weightings(delta_d50)
            all_weightings[key] = weighting
        return all_weightings

    def get_delta_d50(self, duration):
        print(f'Computing temporal pattern D50 scaling for duration of {duration} hours')
        if self.adjust_climate:
            self.check_tp_d50_scaling(duration)
            tp_d50_scaling = self.tp_d50_scaling_df.loc[duration, self.kg_classification]
            if pd.isna(tp_d50_scaling):
                print(f'WARNING: D50 scaling not available for duration of {duration} hours.')
                if duration > 72:
                    apply_duration = 72
                else:
                    apply_duration = 1
                print(f'         Applying scaling for {apply_duration} hours!')
                tp_d50_scaling = self.tp_d50_scaling_df.loc[apply_duration, self.kg_classification]
            print(f'D50 scaling is {tp_d50_scaling}/°C')
            delta_d50 = np.around((self.climate_adjustment_factor(tp_d50_scaling) - 1.0) * 100, 0)
        else:
            delta_d50 = 0
        print('d50 loading is (%): ', delta_d50)
        return delta_d50

    def get_loss_uplift_factor(self, loss_type):
        print(f'Computing rainfall loss uplift factors')
        loss_type = loss_type.upper()
        print('Loss type is: ', loss_type)
        loss_types = ['IL', 'CL']
        if not loss_type in loss_types:
            raise Exception(f'Loss type {loss_type} not recognised. Use il or cl!')
        column_header = f'{loss_type} mean'
        if self.adjust_climate:
            rate_of_change = self.loss_rates[column_header]
            loss_uplift_factor = self.climate_adjustment_factor(rate_of_change)
        else:
            loss_uplift_factor = 1.0
        print(f'{loss_type} loss scaling for climate change is {loss_uplift_factor}')
        return loss_uplift_factor

    def get_rainfall_uplift_factor(self, duration):
        print(f'Computing rainfall uplift factor for duration of {duration} hours')
        if self.adjust_climate:
            # get the uplift rate
            cap_scaling = False
            if 'cap rainfall scaling' in self.config_data.keys():
                cap_scaling = True
                scaling_cap_percent = self.config_data['cap rainfall scaling']
                print(f'Capping the rainfall uplift to no more than {scaling_cap_percent}%')
            a, b, c = [0.01177526, -0.235206, 2.70886]
            uplift_rate = np.exp(a * np.log(duration)**2 + b * np.log(duration) + c)
            if uplift_rate < 8.0 or duration == 24:
                uplift_rate = 8.0
            elif uplift_rate > 15.0 or duration == 1:
                uplift_rate = 15.0
            else:
                uplift_rate = np.around(uplift_rate, 1)
            if cap_scaling and uplift_rate > scaling_cap_percent:
                uplift_rate = float(scaling_cap_percent)
            print('Uplift rate is (%): ', uplift_rate)
            # get the uplift factor
            rainfall_uplift_factor = self.climate_adjustment_factor(uplift_rate)
        else:
            rainfall_uplift_factor = 1.0
        print('Rainfall uplift factor is:', rainfall_uplift_factor)
        return rainfall_uplift_factor

    def climate_adjustment_factor(self, rate_of_change):
        factor = np.around((1 + rate_of_change / 100) ** self.temperature_rise, 5)
        return factor


class D50Weighting:

    def __init__(self, filepath, pattern_type, duration, pattern_area=0):
        self.filepath = filepath
        self.pattern_type = pattern_type
        self.pattern_area = pattern_area
        self.duration = int(duration * 60)  # convert to minutes

    def get_d50_weightings(self, delta_d50):
        print(f'Opening d50 weighting file for {self.pattern_type}: ', self.filepath)
        df = pd.read_csv(self.filepath)

        if self.pattern_type == 'ARR areal' and self.duration < 12 * 60:
            print(f'Adjusting {self.duration} to 12 hours for D50 weighting for {self.pattern_type}')
            df = df[df['Duration'] == 12 * 60]

        elif self.pattern_type == 'GTSMR' and self.duration < 24 * 60:
            print(f'Adjusting {self.duration} to 24 hours for D50 weighting for {self.pattern_type}')
            df = df[df['Duration'] == 24 * 60]

        else:
            df = df[df['Duration'] == self.duration]

        print(df)

        if self.pattern_type == 'ARR areal' or self.pattern_type == 'GTSMR':
            df = df[df['Area'] == self.pattern_area]
            front_patters = df['Front Patterns'].iloc[0]
            df.drop(['Duration', 'Area', 'Front Patterns'], inplace=True, axis=1)
            weightings = self.format_area_weightings(df, delta_d50)
        elif self.pattern_type == 'ARR point':
            front_patters = df['Front Patterns'].iloc[0]
            df.drop(['Duration', 'Front Patterns'], inplace=True, axis=1)
            weightings = self.format_point_weightings(df, delta_d50)
        else:
            front_patters = df['Front Patterns'].iloc[0]
            df.drop(['Duration', 'Front Patterns'], inplace=True, axis=1)
            weightings = self.format_area_weightings(df, delta_d50)

        # Convert the front pattern string to a list
        front_patters = front_patters.split(',')
        front_patters = list(map(int, front_patters))  # convert the items to integers
        return {'weightings': weightings, 'front_patterns': front_patters}

    def format_point_weightings(self, df, delta_d50):
        # Set up the frequent dataframe
        frequent = df[df['AEP'] == 'frequent']
        frequent.drop(['AEP'], inplace=True, axis=1)

        # get the list of d50s
        d50s = frequent.columns
        d50s = list(map(lambda x: int(x.lstrip('DeltaD50')), d50s))
        print(d50s)

        # Continue setting up the frequent dataframe
        frequent = frequent.transpose()
        frequent['DeltaD50'] = d50s
        frequent.set_index('DeltaD50', inplace=True)
        frequent_w = frequent[frequent.columns[0]].loc[delta_d50]
        print('\nWeightings for point temporal patterns for frequent events:')
        print(frequent_w)

        # Set up the intermediate dataframe
        intermediate = df[df['AEP'] == 'intermediate']
        intermediate.drop(['AEP'], inplace=True, axis=1)
        intermediate = intermediate.transpose()
        intermediate['DeltaD50'] = d50s
        intermediate.set_index('DeltaD50', inplace=True)
        intermediate_w = intermediate[intermediate.columns[0]].loc[delta_d50]
        print('\nWeightings for point temporal patterns for intermediate events:')
        print(intermediate_w)

        # Set up the rare dataframe
        rare = df[df['AEP'] == 'rare']
        rare.drop(['AEP'], inplace=True, axis=1)
        rare = rare.transpose()
        rare['DeltaD50'] = d50s
        rare.set_index('DeltaD50', inplace=True)
        rare_w = rare[rare.columns[0]].loc[delta_d50]
        print('\nWeightings for point temporal patterns for rare events:')
        print(rare_w)

        return {'frequent': frequent_w,
                'intermediate': intermediate_w,
                'rare': rare_w}

    def format_area_weightings(self, df, delta_d50):
        d50s = df.columns
        d50s = list(map(lambda x: int(x.lstrip('DeltaD50')), d50s))
        df = df.transpose()
        df['DeltaD50'] = d50s
        df.set_index('DeltaD50', inplace=True)
        weighting = df[df.columns[0]].loc[delta_d50]
        print(f'\nWeightings for {self.pattern_type} temporal patterns:')
        print(weighting)
        return weighting

