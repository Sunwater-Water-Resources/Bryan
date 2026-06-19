import numpy as np
import pandas as pd
import json
import math

class TemporalPatterns:
    def __init__(self, filepath, storm_durations):
        self.pattern_type = ''
        self.temporal_patterns = {}
        self.timesteps = {}
        self.filepath = filepath
        self.storm_durations = storm_durations
        self.pattern_area = 0

    def get_temporal_pattern(self, duration):
        temporal_pattern = self.temporal_patterns[duration]
        return temporal_pattern

    def parse_pmp_pattern_lines(self, lines):
        header = lines[0].split(',')
        increments = int(header[0])
        timestep = int(header[1])
        times = np.linspace(timestep, increments * timestep, increments)
        times = times / 60  # convert to hours
        df = pd.DataFrame(index=times, columns=range(10))
        for i in range(10):
            row = lines[i + 1]
            divs = np.array(row.split(','), dtype=float)
            df[i] = divs
        # df.loc[0.0] = np.zeros(10)
        # df.sort_index(inplace=True)
        # print(df)
        return df

    def get_pmp_pattern_lines(self, filepath, duration):
        duration_min = int(duration * 60)
        key_term = f'{duration_min} minutes'
        pattern_lines = []
        with open(filepath) as f:
            copy = False
            for line in f:
                if line.strip() == key_term:
                    copy = True
                    continue
                elif 'minutes' in line:
                    copy = False
                    continue
                elif copy:
                    pattern_lines.append(line)
        return pattern_lines


class GtsmrPatterns(TemporalPatterns):
    def __init__(self, filepath, storm_durations, area):
        super().__init__(filepath, storm_durations)
        self.pattern_type = 'GTSMR'
        self.area = area

    def import_patterns(self):
        self.pattern_area = self.get_gtsmr_pattern_area(self.area)
        print(f'GTSMR pattern area bin is {self.pattern_area}')
        self.filepath = self.filepath.replace('~AREA~', str(self.pattern_area))
        print('Opening GTSMR temporal pattern file: ', self.filepath)
        for duration in self.storm_durations:
            duration_lines = self.get_pmp_pattern_lines(self.filepath, duration)
            temporal_pattern = self.parse_pmp_pattern_lines(duration_lines)
            print(f'\nGTSMR temporal pattern for duration {duration} hrs')
            print(temporal_pattern)
            self.temporal_patterns[duration] = temporal_pattern
            self.timesteps[duration] = temporal_pattern.index[0]

    def get_gtsmr_pattern_area(self, area):
        if area < 300.0:
            return 100
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


class GsdmPatterns(TemporalPatterns):
    def __init__(self, filepath, storm_durations):
        super().__init__(filepath, storm_durations)
        self.pattern_type = 'GSDM'

    def import_patterns(self):
        print('Opening GSDM temporal pattern file: ', self.filepath)
        for duration in self.storm_durations:
            duration_lines = self.get_pmp_pattern_lines(self.filepath, duration)
            temporal_pattern = self.parse_pmp_pattern_lines(duration_lines)
            print(f'\nGSDM temporal pattern for duration {duration} hrs')
            print(temporal_pattern)
            self.temporal_patterns[duration] = temporal_pattern
            self.timesteps[duration] = temporal_pattern.index[0]


class PointPatterns(TemporalPatterns):
    def __init__(self, filepath, storm_durations):
        super().__init__(filepath, storm_durations)
        self.pattern_type = 'point'

    def import_patterns(self):
        # TODO: code up  the importing of point temporal patterns
        print('Opening ARR point temporal pattern file: ', self.filepath)
        f = open(self.filepath)
        line = f.readline()
        headers = line.split(',')
        headers = list(map(lambda x: x.strip(), headers))  # strip whitespace from elements
        headers = [x for x in headers if x]  # remove any empty items in the header
        del headers[-1]  # remove the last header 'Increments'
        num_headers = len(headers)
        df = pd.DataFrame(columns=headers)
        increment_list = []
        for row, line in enumerate(f):
            line = line.replace('\n', '')
            line = f'{line},'
            lst = line.split(',')
            lst = [x for x in lst if x]  # remove any empty items in the split line
            data = lst[:num_headers]
            df.loc[row] = data
            increment = lst[num_headers:]
            increment_list.append(np.array(increment, dtype=float))

        df = df.astype({"EventID": int, "Duration": int, "TimeStep": int})
        df.reset_index(inplace=True)
        print(df)
        local_increments = {}
        for duration in self.storm_durations:
            converted_duration = duration * 60  # convert to minutes
            print(f'\nExtracting pattern for duration of {duration} hours')
            local_df = df[df.Duration == converted_duration]
            frequency_types = local_df['AEP'].unique()
            for frequency in frequency_types:
                print(f'Extracting {frequency} arr point temporal patterns')
                dff = local_df[local_df['AEP'] == frequency]
                print(dff)
                timestep = dff['TimeStep'].iloc[0]
                self.timesteps[duration] = timestep / 60
                print(f'Timestep is {timestep} minutes')
                idxs = dff.index.tolist()
                local_increment = pd.DataFrame(increment_list[idxs[0]: idxs[-1] + 1]).transpose()
                times = []
                current_time = 0
                first_increment = local_increment[0]
                for i in range(first_increment.shape[0]):
                    current_time += timestep / 60
                    times.append(current_time)
                local_increment['Time'] = times
                local_increment.set_index('Time', inplace=True)
                print(f'\nARR point temporal pattern for duration {duration} hrs')
                print(local_increment)
                local_increments[frequency] = local_increment.copy()
            self.temporal_patterns[duration] = local_increments.copy()


class ArealPatterns(TemporalPatterns):
    def __init__(self, filepath, storm_durations, area):
        super().__init__(filepath, storm_durations)
        self.pattern_type = 'areal'
        self.area = area

    def import_patterns(self):
        self.pattern_area = self.get_areal_pattern_area(self.area)
        if self.pattern_area == 0:
            print('Catchment area is < 75km²')
            print('ARR areal pattern not imported')
            return
        
        print(f'Importing areal temporal patterns')
        print('For storm durations:', self.storm_durations)
        print(f'ARR areal pattern area bin is {self.pattern_area}')
        print('Opening ARR Areal temporal pattern file: ', self.filepath)
        f = open(self.filepath)
        line = f.readline()
        headers = line.split(',')
        last = headers.pop()
        num_headers = len(headers)
        df = pd.DataFrame(columns=headers)
        increment_list = []
        for row, line in enumerate(f):
            lst = line.split(',')
            data = lst[:num_headers]
            if int(data[-1]) == self.pattern_area:
                df.loc[row] = data
                increment = lst[num_headers:]
                increment[-1] = increment[-1].strip()
                increment_list.append(np.array(increment, dtype=float))

        df = df.astype({"EventID": int, "Duration": int, "TimeStep": int})
        df.reset_index(inplace=True)
        print(df)
        for duration in self.storm_durations:
            duration_label = duration
            if duration < 12:
                print(f'Changing duration because it is less than 12 hours, but keeping duration label of {duration_label} hours')
                duration = 12
            converted_duration = duration * 60  # convert to minutes
            print(f'\nExtracting pattern for duration of {duration} hours')
            local_df = df[df.Duration == converted_duration]
            print(local_df)
            timestep = local_df['TimeStep'].iloc[0]
            print(f'Timestep is {timestep} minutes')
            idxs = local_df.index.tolist()
            local_increment = pd.DataFrame(increment_list[idxs[0]: idxs[-1] + 1]).transpose()
            times = []
            current_time = 0
            for i in range(local_increment.shape[0]):
                current_time += timestep / 60
                times.append(current_time)
            local_increment['Time'] = times
            local_increment.set_index('Time', inplace=True)
            print(f'\nARR areal temporal pattern for duration {duration} hrs')
            print(local_increment)
            self.temporal_patterns[duration_label] = local_increment.copy()
            self.timesteps[duration_label] = timestep / 60

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


class PreburstPatterns:
    def __init__(self, filepath):
        self.filepath = filepath
        self.preburst_patterns = {}
        self.sampling_key = {}
    
    def import_patterns(self, median_only = False):
        with open(self.filepath, 'r') as f:
            data = json.load(f)
        
        preburst_patterns = {}
        for duration, dur_df in data.items():
            duration = int(duration)

            # temppats = {}
            # for sampling_limit, bin_df in dur_df:
            #     sampling_limit = float(sampling_limit)
            #     df = pd.DataFrame.from_dict(bin_df)
            #     df.index = df.index.to_series().apply(lambda x: int(x))             # Convert index names to int
            #     df.columns = df.columns.to_series().apply(lambda x: int(x))         # Convert columns names to float
            #     temppats[sampling_limit] = df
            # preburst_patterns[duration] = temppats
            
            df = pd.DataFrame.from_dict(dur_df)
            df.index = df.index.to_series().apply(lambda x: int(x))             # Convert index names to int
            df.columns = df.columns.to_series().apply(lambda x: int(x))         # Convert columns names to float
            
            if median_only:                                         # For ensemble method
                df1 = df.cumsum() / df.sum()                        # Normalise
                df1.fillna(0, inplace = True)                       # Replace NA with zeroes, so median can be calculated
                median_pattern = df1.median(axis = 1).diff()
                
                first_idx = median_pattern.index[median_pattern.notna() & 
                                                 median_pattern.ne(0)] [0]      # Get index of first non-zero, notna 
                median_pattern = median_pattern[first_idx:]                     # Remove leading zeroes from series
                if median_pattern.sum() != 1.0:
                    raise Exception('Median pre-burst pattern error')
                
                preburst_patterns[duration] = median_pattern
                continue
            
            preburst_patterns[duration] = df
        
        self.preburst_patterns = preburst_patterns
        return
    
    def get_preburst_pattern(self, duration, preburst_proportion, timestep, sample_int = None): #, exclusion=None):
        df = self.preburst_patterns[duration]
        
        if sample_int == 'median':
            pattern = df
        else:
            if sample_int is None:
                sample_proportions = df.sum(axis = 0)
                distance = abs(sample_proportions - preburst_proportion)
                sample_int = distance[distance.rank() <= 3].index.to_list()[np.random.randint(3)]   # Randomly select pattern from 3 nearest
                # if exclusion is None:
                #     sample_int = distance[distance.rank() <= 3].index.to_list()[np.random.randint(3)]   # Randomly select pattern from 3 nearest
                # else:
                #     nearest_three_idx = distance[distance.rank() <= 3].index.to_list()
                #     lst = [idx for idx in nearest_three_idx if idx not in exclusion]
                #     if len(lst) > 0:
                #         sample_int = lst[np.random.randint(len(lst))]
                #     else:
                #         remainder = distance.loc[~distance.index.isin(exclusion)]
                #         if len(remainder) >= 1:
                #             sample_int = remainder[remainder.rank() == 1].index[0]
                #         else:
                #             return None, np.nan
                    
            pattern = df.iloc[:, sample_int]
            
            first_idx = pattern.index[pattern.notna()] [0]      # Get index of first notna 
            pattern = pattern[first_idx:]                       # Remove leading NAs from series
        
        pattern = pattern * preburst_proportion / pattern.sum()         # Scale to preburst_proportion
        pattern = pattern * 100                                         # Convert to percentage
        
        if not timestep.is_integer():
            if timestep < 1.0 and (1.0/timestep).is_integer():
                num_steps = int(len(pattern) / timestep)
                new_index = np.linspace(pattern.index[0], 0, num_steps + 1)[:-1]
                pattern = pattern.reindex(new_index).ffill() * timestep
            else:
                raise Exception('Non-integer pre-burst timestep not coded yet')
        
        if timestep > 1.0:
            resampled_df = pattern.to_frame()
            resampled_index = pattern.index.to_series() / timestep
            resampled_df['resampled timestep'] = resampled_index.apply(math.floor) * timestep
            pattern = resampled_df.groupby('resampled timestep').sum().iloc[:, 0]
        
        pattern = np.around(pattern, 3)
        
        return pattern, sample_int