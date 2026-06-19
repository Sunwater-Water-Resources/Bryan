import pandas as pd
import numpy as np
from scipy.special import ndtri
import matplotlib.pyplot as plt
import os
from scipy import interpolate


class MonteCarloSimulation:
    def __init__(self, filepath, result_type, duration, drop_aep=None):
        self.duration = duration
        self.warnings = []
        self.drop_aep = drop_aep
        self.quantiles = self.read_quantiles(filepath, result_type)
        self.filepath = filepath
        if self.quantiles is not None:
            self.std_aeps = self.quantiles.index.to_numpy()
            self.std_z = ndtri(1 - 1 / self.std_aeps)
        else:
            self.std_aeps, self.std_z = [None, None]
        self.upper_raw, self.lower_raw = self.read_raw_percentiles(filepath.replace('.csv', '_perc.csv'))
        self.upper_smooth, self.lower_smooth = self.read_smooth_percentiles(filepath.replace('.csv', '_perc_smooth.csv'))
        self.result_type = result_type
        self.raw_percentile_z = None

    def read_quantiles(self, filepath, result_type):
        print('\nOpening URBS result file:', filepath)
        try:
            df = pd.read_csv(filepath)
            df.set_index('aep (1 in x)', inplace=True)
            print(f'Result type is:', result_type)
            df = df[result_type]
            if self.drop_aep is not None:
                df.drop(self.drop_aep, inplace=True)

        except Exception:
            warning = f'WARNING: failed to open: {filepath}'
            self.warnings.append(warning)
            print(warning)
            df = None
        return df

    def read_raw_percentiles(self, filepath):
        print('\nOpening URBS result file:', filepath)
        try:
            df = pd.read_csv(filepath)
            df.set_index('AEP', inplace=True)
            headers = df.columns
            upper_perc_header = headers[-1]
            lower_perc_header = headers[-2]
            upper = df[upper_perc_header]
            lower = df[lower_perc_header]
            # self.raw_percentile_z = df['z'].to_numpy()

        except Exception:
            warning = f'WARNING: failed to open: {filepath}'
            self.warnings.append(warning)
            print(warning)
            upper, lower = [None, None]
        return upper, lower

    def read_smooth_percentiles(self, filepath):
        print('\nOpening URBS result file:', filepath)
        try:
            df = pd.read_csv(filepath)
            df.set_index('AEP', inplace=True)
            headers = df.columns
            upper_perc_header = headers[-1]
            lower_perc_header = headers[-2]
            upper = df[upper_perc_header]
            lower = df[lower_perc_header]

        except Exception:
            warning = f'WARNING: failed to open: {filepath}'
            self.warnings.append(warning)
            print(warning)
            upper, lower = [None, None]
        return upper, lower
    
    def create_volume_results(self, storage_curve):
        if self.quantiles is not None and self.result_type == 'level':
            # print(self.quantiles)
            volume_filepath = self.filepath.replace('level', 'volume')
            volume = storage_curve(self.quantiles)
            df = pd.DataFrame(index=self.quantiles.index, columns=['volume'], data=volume)
            df.to_csv(volume_filepath)
        if self.upper_raw is not None and self.result_type == 'level':
            df = pd.DataFrame(index=self.upper_raw.index)
            volume = storage_curve(self.upper_raw)
            df[self.upper_raw.name] = volume
            volume = storage_curve(self.lower_raw)
            df[self.lower_raw.name] = volume
            volume_filepath = self.filepath.replace('level', 'volume_perc')
            df.to_csv(volume_filepath)
        if self.upper_smooth is not None and self.result_type == 'level':
            df = pd.DataFrame(index=self.upper_smooth.index)
            volume = storage_curve(self.upper_smooth)
            df[self.upper_smooth.name] = volume
            volume = storage_curve(self.lower_smooth)
            df[self.lower_smooth.name] = volume
            volume_filepath = self.filepath.replace('level', 'volume_perc_smooth')
            df.to_csv(volume_filepath)


class MonteCarloSimulationGroup:
    def __init__(self, output_folder, output_name=None, drop_aeps=[]):
        self.simulations = {}
        self.std_aeps = None
        self.std_z = None
        self.upper_header = None
        self.lower_header = None
        self.critical_durations = None
        self.output_folder = output_folder
        self.output_name = output_name
        self.result_type = None
        self.drop_aeps = drop_aeps
        self.storage_curve = None

    def set_volume_curve(self, storage_filepath):
        print('Reading the storage file:', storage_filepath)
        df = pd.read_csv(storage_filepath)
        self.storage_curve = interpolate.interp1d(df['EL'], df['V'])

    def add_simulation(self, simulation, duration):
        self.simulations[duration] = simulation
        if self.result_type == 'level':
            if self.storage_curve is not None:
                simulation.create_volume_results(self.storage_curve)
        if simulation.quantiles is not None:
            self.std_aeps = simulation.std_aeps
            self.std_z = simulation.std_z
            self.result_type = simulation.result_type
        if simulation.upper_raw is not None:
            self.upper_header = simulation.upper_smooth.name
            self.lower_header = simulation.lower_smooth.name

    def compute_critical_durations(self, plot=True):
        critical_durations = self.analyse_critical_duration(plot)
        if critical_durations is not None:
            percentiles = self.analyse_percentiles(critical_durations['critical_duration'])
            critical_durations = pd.concat([critical_durations, percentiles], axis=1)
            critical_durations.index.name = 'aep (1 in x)'
            print(critical_durations.to_string())
            self.critical_durations = critical_durations
            out_file = os.path.join(self.output_folder,
                                    '{}.csv'.format(self.output_name))
            print('\nWriting critical duration file:', out_file)
            critical_durations.to_csv(out_file)
        else:
            print('WARNING: no files have been found!')

    def analyse_critical_duration(self, plot=False):
        all_quantiles = []
        for duration, simulation in self.simulations.items():
            if simulation.quantiles is not None:
                sim_df = simulation.quantiles.copy()
                sim_df.name = '{}h'.format(duration)
                all_quantiles.append(sim_df)
        if len(all_quantiles) > 0:
            df = pd.concat(all_quantiles, axis=1)
            if plot:
                self.plot_durations(df)
            crit_durations = df.idxmax(axis=1)
            df['max'] = df.max(axis=1)
            df['critical_duration'] = crit_durations
        else:
            df = None
        return df

    def plot_durations(self, df, dpi=200):
        plt_df = df.copy()
        if self.drop_aeps:
            plt_df.drop(index=self.drop_aeps, inplace=True)
        std_aeps = plt_df.index
        std_z = ndtri(1 - 1 / std_aeps)
        plt_df['z'] = std_z
        plt_df.set_index('z', inplace=True)

        if self.result_type == 'level':
            plt_df.plot()  # , title=self.output_name)
        else:
            plt_df.plot(logy=True)  # , title=self.output_name)
        plt_file = os.path.join(self.output_folder,
                                '{}_durations.png'.format(self.output_name))
        print('\nCreating plot of the durations:', plt_file)
        ylabels = {'inflow': 'Flow (m³/s)',
                   'level': 'Lake level (m AHD)',
                   'outflow': 'Flow (m³/s)',
                   'volume': 'Volume (ML)'}
        plt.ylabel(ylabels[self.result_type])
        plt.xlabel('AEP (1 in X)')
        plt.xticks(std_z, std_aeps, rotation=90)
        plt.tight_layout()
        plt.savefig(plt_file, dpi=dpi)

    def analyse_percentiles(self, critical_durations):
        all_upper = []
        all_lower = []
        for duration, simulation in self.simulations.items():
            if simulation.upper_smooth is not None:
                upper = simulation.upper_smooth.copy()
                upper.name = '{}h'.format(duration)
                all_upper.append(upper)
            if simulation.lower_smooth is not None:
                lower = simulation.lower_smooth.copy()
                lower.name = '{}h'.format(duration)
                all_lower.append(lower)
        upper_df = pd.concat(all_upper, axis=1)
        lower_df = pd.concat(all_lower, axis=1)
        for aep in self.std_aeps:
            critical_duration = critical_durations[aep]
            upper_df.loc[aep, self.upper_header] = upper_df.loc[aep, critical_duration]
            lower_df.loc[aep, self.lower_header] = lower_df.loc[aep, critical_duration]
        df = pd.concat([lower_df[self.lower_header], upper_df[self.upper_header]], axis=1)
        df = self.smoothen_percentiles(df)
        return df

    def smoothen_percentiles(self, df):
        upper_title = self.upper_header
        lower_title = self.lower_header
        # split the data and remove zero/NA flows
        df['z'] = ndtri(1 - 1 / df.index)
        upper = df[['z', upper_title]]
        lower = df[['z', lower_title]]
        upper[upper_title] = np.log10(upper[upper_title])
        lower[lower_title] = np.log10(lower[lower_title])
        upper[upper_title].replace([np.inf, -np.inf], np.nan, inplace=True)
        lower[lower_title].replace([np.inf, -np.inf], np.nan, inplace=True)
        upper.dropna(inplace=True)
        lower.dropna(inplace=True)
        lower = lower.loc[lower[lower_title] > 0]

        # Make monotonic
        lower = lower.iloc[::-1]
        x1 = upper['z']
        y1 = np.maximum.accumulate(upper[upper_title])
        x2 = lower['z']
        y2 = np.minimum.accumulate(lower[lower_title])

        # Build models of the data
        degrees = 5
        f1 = np.poly1d(np.polyfit(x1, y1, degrees))
        f2 = np.poly1d(np.polyfit(x2, y2, degrees))
        m1 = 10 ** (f1(x1))
        m2 = 10 ** (f2(x2))

        # build the smoothed dataframe
        new_lower_title = f'{lower_title}_smooth'
        lower[new_lower_title] = m2
        new_upper_title = f'{upper_title}_smooth'
        upper[new_upper_title] = m1
        df = pd.concat([df, lower[new_lower_title], upper[new_upper_title]], axis=1)
        return df
