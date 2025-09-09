from datetime import datetime
import time
import os
import numpy as np
import pandas as pd
import re
from scipy import optimize


class Router:
    def __init__(self, parameters):
        self.start = time.time()

        self.mcdf = pd.read_csv(parameters['Input MCDF'], index_col=0)
        std_norm_variates = self.mcdf['rain_z'].to_numpy()

        # What work needs to be done?
        do_dam_route = self.set_logic(parameters, 'Dam route')
        add_baseflow = self.set_logic(parameters, 'Baseflow')

        # Start with adding baseflow
        if add_baseflow:
            bf_obj = Baseflow(parameters)
            bf_obj.compute_baseflows(std_norm_variates)
            bf_obj.store_inflows()
            new_mcdf_path = parameters['Input MCDF']
            new_mcdf_path = new_mcdf_path.replace('.csv', '_tf.csv')
            bf_obj.write_max_to_mcdf(self.mcdf, new_mcdf_path)

        # The dam routing part
        if do_dam_route:
            pass
            dr_obj = DamRouter(parameters)
            dr_obj.set_dam_properties_from_urbs()
            dr_obj.set_dam_hydraulics()
            dr_obj.dam_route_iterative()
            dr_obj.store_outflows()
            dr_obj.write_max_to_mcdf()

        end = time.time()
        elapsed_time = np.around((end - self.start) / 60, 2)
        print(f'Elapsed time: {elapsed_time} minutes')

    def set_logic(self, parameters, key):
        print(f'Checking logic for {key}...', end='')
        if key in parameters:
            if parameters[key].lower() == 'yes':
                print('True')
                return True
            else:
                print('False')
                return False
        else:
            print('False')
            return False


class Baseflow:
    def __init__(self, parameters, b0=0.0, bm=1.0):
        print('\nSetting  up the baseflow computation object')
        self.inflows = None
        self.br24 = parameters['BR']
        self.b0 = b0
        self.timestep = None
        self.br = None
        self.bfvf10 = parameters['BFVF10']
        self.bc = None
        self.bm = bm
        self.baseflows = None
        self.filepaths = {}
        self.total_inflows = None
        self.inflow_path = parameters['Quickflow']
        self.set_inflows(self.inflow_path)

    def compute_baseflows(self, std_norm_variates):
        print('\nComputing the baseflows')
        bc = self.get_bc(std_norm_variates)
        for ind, timestep in enumerate(self.inflows.index):
            if ind == 0:
                self.baseflows.iloc[ind] = self.b0
            else:
                self.baseflows.loc[timestep] = (self.b0 +
                                                self.br * (self.baseflows.iloc[ind - 1] - self.b0) +
                                                bc * self.inflows.loc[timestep] ** self.bm)
        self.total_inflows = np.around(self.inflows + self.baseflows, 3)

    def store_inflows(self):
        baseflows_path = self.inflow_path.replace('.csv', '_bf.csv')
        print('\nWriting baseflow to file:', baseflows_path)
        self.baseflows.to_csv(baseflows_path)
        total_flows_path = self.inflow_path.replace('.csv', '_tf.csv')
        print('\nWriting total flow to file:', total_flows_path)
        self.total_inflows.to_csv(total_flows_path)

    def set_inflows(self, inflows_path):
        self.filepaths['inflows'] = inflows_path
        self.inflows = pd.read_csv(inflows_path, index_col=0)
        # some RORB inflows are run on different timesteps creating gaps. So, interpolate gaps if needed.
        self.inflows.interpolate(method='slinear', inplace=True)
        self.timestep = self.inflows.index[1] - self.inflows.index[0]
        print('Timestep:', self.timestep)
        self.br = self.br24 ** (self.timestep / 24)
        
        print(f'Baseflow recession (BR): {self.br}')
        self.baseflows = pd.DataFrame(index=self.inflows.index, columns=self.inflows.columns)

    def write_max_to_mcdf(self, mcdf, new_path):
        mcdf['inflow'] = self.total_inflows.max().to_numpy()
        print('\nWriting the new mcdf:', new_path)
        mcdf.to_csv(new_path)
    
    def get_bc(self, z):
        a = -0.02079
        b = -0.1375
        c = 0.20788
        r = a * z ** 2 + b * z + c
        bfvf = 10 ** r * self.bfvf10
        bc24 = (1 - self.br24) * bfvf
        bc = bc24 * self.timestep / 24  # check this conversion
        return bc


class DamRouter:
    def __init__(self, parameters):
        print('\nSetting up the dam routing object.')
        self.filepaths = {
            'els': parameters['ELS file'],
            'sq': parameters['SQ file'],
            'old_mcdf': parameters['Input MCDF'],
            'inflows': parameters['Inflow']
        }
        self.els = None
        self.sq_curve = None
        self.sims = None
        self.fsl = parameters['FSL']
        self.fsv = None
        self.timestep = None
        self.inflows = None
        self.lake_levels = None
        self.outflows = None
        self.volumes = None
        self.mcdf = pd.read_csv(self.filepaths['old_mcdf'], index_col=0)
        output_suffix = {'apply': False, 'value': 'tf'}  # defaults to outputting tf suffix for baseflow analysis
        if 'Output suffix' in parameters.index:
            if pd.notna(parameters['Output suffix']):
                output_suffix = {'apply': True, 'value': parameters['Output suffix']}
                print('Found output suffix:', output_suffix['value'])
        self.output_suffix = output_suffix
        self.filepaths['new_mcdf'] = self.filepaths['old_mcdf'].replace(
            '.csv', '_{}.csv'.format(output_suffix['value']))

    def set_dam_hydraulics(self):
        print('\nSetting up the dam hydraulics...')
        print(f'    Full supply level is {self.fsl} m AHD')
        self.fsv = self.storage_from_level(self.fsl)
        print(f'    Full supply volume is {self.fsv} ML')
        print('    Reading dam inflows:', self.filepaths['inflows'])
        self.inflows = pd.read_csv(self.filepaths['inflows'], index_col=0)
        # some RORB inflows are run on different timesteps creating gaps. So, interpolate gaps if needed.
        self.inflows.interpolate(method='slinear', inplace=True)  # some RORB inflows
        self.timestep = self.inflows.index[1] - self.inflows.index[0]
        print(f'    Found timestep of {self.timestep} hours.')
        self.sims = self.inflows.columns.tolist()
        time_series = self.inflows.index
        print('    Setting up containers for levels, outflows, and volumes.')
        self.lake_levels = pd.DataFrame(index=time_series, columns=self.sims)
        self.outflows = pd.DataFrame(index=time_series, columns=self.sims)
        self.volumes = pd.DataFrame(index=time_series, columns=self.sims)
        print('    Getting the antecedent dam volumes from the MCDF.')
        intial_volumes = self.mcdf['ADV'].to_numpy()
        self.volumes.iloc[0] = intial_volumes
        self.outflows.iloc[0] = self.outflow_from_storage(intial_volumes)
        self.lake_levels.iloc[0] = self.level_from_storage(intial_volumes)
        # self.volumes.iloc[0] = self.storage_from_level(self.lake_levels.iloc[0])

    def dam_route_iterative(self):
        print('    Finding the time at which each simulation starts spilling')
        spill_status = pd.Series(index=self.sims, dtype=str)
        initial_lake_levels = self.lake_levels.iloc[0]
        initial_lake_storages = self.volumes.iloc[0]
        delta_storage = self.inflows.rolling(window=2).mean() * 3600 * self.timestep / 1000
        delta_storage.iloc[0] = initial_lake_storages
        inflow_storage = delta_storage.cumsum()
        spill_times = pd.Series(index=self.sims, dtype=float)
        for sim in self.sims:
            if initial_lake_levels[sim] < self.fsl:
                local_inflow_storage = inflow_storage[sim].dropna()
                if local_inflow_storage.max() < self.fsv:  # dam does not spill
                    spill_times.loc[sim] = -1
                    spill_status.loc[sim] = 'Never spills'
                else:  # dam spills at some stage during the simulation
                    spill_times.loc[sim] = np.interp(self.fsv, local_inflow_storage, local_inflow_storage.index)
                    spill_status.loc[sim] = 'Filling'
            else:  # dam is spilling at start of simulation
                spill_times.loc[sim] = 0
                spill_status.loc[sim] = 'Spilling'
        #spill_times.to_csv('spill_times.csv')

        print('\nDoing the dam routing...')
        log_every_n = 20
        for sim_num, sim in enumerate(self.sims):
            if (sim_num % log_every_n) == 0:
                print(f'Simulation: {sim}')
            inflows = self.inflows[sim].dropna()
            for ind, timestep in enumerate(inflows.index):
                if ind > 0:
                    o1 = self.outflows.loc[timestep - self.timestep, sim]
                    i2 = inflows.loc[timestep]
                    i1 = inflows.loc[timestep - self.timestep]
                    s1 = self.volumes.loc[timestep - self.timestep, sim]

                    vol_in = 0.5 * (i1 + i2) * self.timestep * 3600 / 1000  # converts from m³ to ML
                    new_vol = s1 + vol_in
                    new_level = self.level_from_storage(new_vol)

                    if spill_status[sim] == 'Never spills':
                        self.outflows.loc[timestep, sim] = 0.0
                        self.lake_levels.loc[timestep, sim] = new_level
                        self.volumes.loc[timestep, sim] = new_vol

                    elif spill_status[sim] == 'Filling':
                        if new_vol < self.fsv:  # doesn't spill on this timestep
                            self.outflows.loc[timestep, sim] = 0.0
                            self.lake_levels.loc[timestep, sim] = new_level
                            self.volumes.loc[timestep, sim] = new_vol
                        else:  # starting to spill during this timestep
                            spill_time = timestep - spill_times.loc[sim]
                            m = (i2 - i1) / self.timestep
                            c = i2 - m * timestep
                            new_i1 = m * spill_times.loc[sim] + c
                            o2_1, o2_2 = self.initial_outflow_estimates(o1, new_vol, 0.2)
                            outflow, res = optimize.newton(
                                self.volume_func, x0=o2_1, x1=o2_2,
                                args=(o1, i2, new_i1, self.fsv, spill_time),
                                disp=False, full_output=True)
                            storage = self.storage_from_outflow(outflow)
                            level = self.level_from_storage(storage)
                            self.outflows.loc[timestep, sim] = outflow
                            self.lake_levels.loc[timestep, sim] = level
                            self.volumes.loc[timestep, sim] = storage
                            spill_status.loc[sim] = 'Spilling'

                    elif spill_status[sim] == 'Spilling':
                        o2_1, o2_2 = self.initial_outflow_estimates(o1, new_vol, 0.2)
                        if o2_1 == o2_2:
                            outflow, res = optimize.newton(self.volume_func, x0=o2_1, x1=o2_2,
                                                           args=(o1, i2, i1, s1, self.timestep),
                                                           disp=False, full_output=True)
                        else:
                            outflow, res = optimize.newton(self.volume_func, x0=o2_1, x1=o2_2,
                                                           args=(o1, i2, i1, s1, self.timestep),
                                                           disp=False, full_output=True)
                        storage = self.storage_from_outflow(outflow)
                        level = self.level_from_storage(storage)
                        self.outflows.loc[timestep, sim] = outflow
                        self.lake_levels.loc[timestep, sim] = level
                        self.volumes.loc[timestep, sim] = storage

    def initial_outflow_estimates(self, previous_outflow, new_volume, factor):
        gw_outflow = self.outflow_from_storage(new_volume)
        guess_1 = gw_outflow * (1 - factor) + previous_outflow * factor
        guess_2 = gw_outflow * factor + previous_outflow * (1 - factor)
        if guess_1 == guess_2:
            guess_2 = guess_2 * 0.9
        return guess_1, guess_2

    def store_outflows(self):
        results = {'outflows': self.outflows,
                   'levels': self.lake_levels,
                   'volumes': self.volumes}
        for result_type in results.keys():
            result_path = self.filepaths['inflows'].replace('inflows', result_type)
            if self.output_suffix['apply']:
                result_path = result_path.replace('.csv', '_{}.csv'.format(self.output_suffix['value']))
            print('\nWriting routed outflows to file:', result_path)
            result = results[result_type]
            result.to_csv(result_path)

    def volume_func(self, o2, o1, i2, i1, s1, dt):
        dt = dt * 3600  # converting hours to seconds
        lhs = 0.5 * (i1 + i2) - 0.5 * (o2 + o1)
        s2 = self.storage_from_outflow(o2)
        rhs = (s2 - s1) * 1000 / dt  # converting ML to m³
        return lhs - rhs

    def set_dam_properties_from_urbs(self):
        print('\nReading the dam properties from ELS and SQ files...')
        print('    Reading:', self.filepaths['els'])
        self.els = pd.read_csv(self.filepaths['els'])

        print('    Reading:', self.filepaths['sq'])
        sq_curve = pd.DataFrame(columns=['Storage', 'Flow'])
        counter = 0
        with open(self.filepaths['sq']) as f:
            # skip the first five lines
            for i in range(6):
                line = f.readline()
            while line:
                line = line.strip()
                if line:
                    s, q = re.split('\s+', line.strip())
                    sq_curve.loc[counter, 'Storage'] = float(s)
                    sq_curve.loc[counter, 'Flow'] = float(q)
                    line = f.readline()  # Read the next line
                    counter += 1
        self.sq_curve = sq_curve  # note this is the storage above spillway crest!

    def storage_from_level(self, level):
        return np.interp(level,
                         self.els['EL'].astype(float),
                         self.els['V'].astype(float))

    def level_from_storage(self, storage):
        return np.interp(storage,
                         self.els['V'].astype(float),
                         self.els['EL'].astype(float))

    def storage_from_outflow(self, outflow):
        storage = self.sq_curve['Storage'].astype(float) + self.fsv
        return np.interp(outflow,
                         self.sq_curve['Flow'].astype(float),
                         storage)

    def outflow_from_storage(self, outflow):
        storage = self.sq_curve['Storage'].astype(float) + self.fsv
        return np.interp(outflow,
                         storage,
                         self.sq_curve['Flow'].astype(float))

    def write_max_to_mcdf(self):
        self.mcdf['inflow'] = self.inflows.max().to_numpy()
        self.mcdf['outflow'] = self.outflows.max().to_numpy()
        self.mcdf['level'] = self.lake_levels.max().to_numpy()
        print('\nWriting the new mcdf:', self.filepaths['new_mcdf'])
        self.mcdf.to_csv(self.filepaths['new_mcdf'])
