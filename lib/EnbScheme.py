import os.path
import scipy.stats
from scipy.special import ndtri, ndtr
from scipy import stats
from scipy import interpolate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Ensemble:
    def __init__(self, aep_list, durations, number_of_temporal_patterns=10, output_folder=''):
        self.output_folder = output_folder
        self.number_of_temporal_patterns = number_of_temporal_patterns
        self.aep_list = aep_list
        self.durations = durations
        self.m = number_of_temporal_patterns * len(aep_list) * len(
            durations)  # number of simulations for a storm duration
        self.df = pd.DataFrame(index=range(self.m),
                               columns=['rain_z', 'rain_aep', 'mean_rain_mm', 'duration', 'tp', 'storm_method',
                                        'tp_frequency', 'preburst_p', 'preburst_proportion',
                                        'preburst_mm', 'initial_loss', 'continuing_loss', 'residual_depth',
                                        'embedded_bursts', 'ADV', 'inflow', 'level', 'outflow'])
        # set up the simulation list
        row = 0
        for aep in aep_list:
            for duration in durations:
                for tp in range(number_of_temporal_patterns):
                    self.df.loc[row, 'rain_aep'] = aep
                    self.df.loc[row, 'rain_z'] = self.get_rain_z(aep)
                    self.df.loc[row, 'duration'] = duration
                    self.df.loc[row, 'tp'] = tp
                    row += 1
        # print(self.df)

    def update_scheme_in_interim(self, aep_range):
        # Adjust the number of events to double up on ARR and extreme patterns in interim zone
        range_query = '{} < rain_aep <= {}'.format(*aep_range)
        extra_sims = self.df.query(range_query)
        extra_sims['storm_method'] = 'extreme'
        self.df = pd.concat([self.df, extra_sims], ignore_index=True)
        self.df.sort_values(by=['rain_aep', 'duration', 'storm_method'], inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        # print(self.df)
        # input()

    def get_rain_z(self, aep):
        f = 1 - 1 / aep
        z = ndtri(f)
        return z

    def store_simulations(self, filename='simulation_out'):
        # output_path = os.path.join(self.output_folder, filename)
        output_path = f'{filename}.csv'
        try:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            print('Storing the ensemble method analysis file:', output_path)
            self.df.to_csv(output_path)
        except IOError:
            input("Could not save the simulation file. The file may be open in Excel. Please close the file and press enter.")
            self.df.to_csv(output_path)
