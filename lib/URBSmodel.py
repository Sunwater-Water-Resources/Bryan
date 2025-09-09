import os
import subprocess
import json
import numpy as np
import glob
import shutil
from pathlib import Path
import pandas as pd
from scipy import interpolate
import time


class UrbsModel:
    def __init__(self, config_file, batch_file='_run.bat', dam_location=None, sub_folder=None):
        # open the config file and get contents
        f = open(config_file)
        config_data = json.load(f)
        f.close()

        # Sort out baseflow treatment
        self.bfvf = 0
        self.apply_baseflow = False
        self.baseflow_locations = None
        if 'baseflow' in config_data.keys():
            self.baseflow_info = config_data['baseflow']
            self.apply_baseflow = self.baseflow_info['apply']
            if 'locations' in self.baseflow_info:
                self.baseflow_locations = self.baseflow_info['locations']
                print('Only applying baseflow to the following locations:')
                print(self.baseflow_locations)

        # Sort out other parameterisation
        self.full_supply_volume = None
        if 'full_supply_volume' in config_data:
            self.full_supply_volume = config_data['full_supply_volume']
        self.urbs_exe = config_data['model_exe']
        if not os.path.exists(self.urbs_exe):
            print("WARNING: The URBS path could not be found:", path_to_check)
            print('Pausing the simulation... hit enter to ignore the warning above.')
            input()
        self.model_folder = os.path.normpath(os.path.join(os.path.dirname(config_file), config_data['model_folder']))
        self.vec_file = config_data['vec_file']
        self.baseflow_vec_file = '{}_baseflow.vec'.format(Path(self.vec_file).stem)

        # Get the vector file lines
        filepath = os.path.join(self.model_folder, self.vec_file)
        with open(filepath, 'r') as f:
            self.vec_lines = f.readlines()
        
        # Create sub-folder for output files. If sub-folder exists from previous run then delete first
        if sub_folder is not None:
            self.output_folder = os.path.join(self.model_folder, sub_folder)
            if os.path.exists(self.output_folder): shutil.rmtree(self.output_folder)
            os.makedirs(self.output_folder)
            self.copy_catchment_data_file()     # Copy to output folder
        else:
            self.output_folder = self.model_folder
        self.ratings_folder = os.path.normpath(os.path.join(os.path.dirname(config_file), config_data['ratings_folder']))

        if sub_folder is not None:
            self.storms_folder = os.path.normpath(os.path.join(os.path.dirname(config_file),
                                                               config_data['storms_folder'],
                                                               sub_folder))
            # Create sub-folder for storm files. If sub-folder exists from previous run then delete first
            if os.path.exists(self.storms_folder): shutil.rmtree(self.storms_folder)
            os.makedirs(self.storms_folder)
        else:
            self.storms_folder = os.path.normpath(os.path.join(os.path.dirname(config_file),
                                                               config_data['storms_folder']))

        # results_folder no longer used?
        if 'results_folder' in config_data.keys():
            self.results_folder = os.path.normpath(os.path.join(os.path.dirname(config_file),
                                                                config_data['results_folder']))
        else:
            self.results_folder = ''
        self.store_tuflow = config_data['store_tuflow']
        self.result_prefix = config_data['result_prefix']
        self.time_increment = config_data['time_increment']
        self.batch_file = batch_file
        self.paramter_string = ''
        self.results_list = []
        self.inflows = []
        self.outflows = []
        self.levels = []

        self.alpha = config_data['alpha']
        self.beta = config_data['beta']
        self.m_exponent = config_data['m_exponent']

        self.volume_below_fsl = 0

        self.set_run_parameters(alpha=self.alpha, beta=self.beta, m_exponent=self.m_exponent)
        self.simulation_periods = config_data['simulation_periods']

        # set up the header for the batch file
        self.header = []
        self.header.append(f'cd {self.output_folder}')
        self.header.append(f'del {self.output_folder}\\urbserr.log')
        self.header.append(f'del {self.output_folder}\\urbsout.log')
        self.header.append('set URBS_LOGF=TRUE')
        self.header.append(f'set URBS_LOGD={self.output_folder}')
        self.header.append(f'set URBS_RATS={self.ratings_folder}')
        if config_data['store_tuflow']:
            self.header.append('set URBS_TFLW=TRUE')

        # get the line number where the dam routing is done
        # this isn't neededanymore because now using volume below FSL
        # rather than lake level and not changing the vec file anymore - using environment variable
        self.dam_routing_line = None
        self.full_supply_level = None
        self.els_file = None
        self.storage_curve = None
        if dam_location:
            self.dam_routing_line = self.find_dam_routing_line(dam_location)
            self.dam_routing_method = self.get_dam_routing_method()
            if self.dam_routing_method == 'level':
                self.storage_curve = self.get_els_data()
                self.full_supply_volume = self.get_volume_from_level(self.full_supply_level)
                print(f'Found the FSV of {self.full_supply_volume} ML from the FSL of {self.full_supply_level} m AHD')
            else:
                if self.full_supply_volume is None:
                    raise Exception('It looks like the FSV has not been provided in the model config file?')
        if 'time_increment_override' in config_data.keys():
            self.time_increment_override = config_data['time_increment_override']
            print('Found time increment overrides:')
            print(self.time_increment_override)
        else:
            self.time_increment_override = None

    def insert_baseflow_into_vec(self, bfvf10, z):
        bfvf_f = self.baseflow_info['bfvf10_factor']
        bfvf = self.get_bfvf(bfvf10 * bfvf_f, z)
        br = self.baseflow_info['br']
        bm = self.baseflow_info['bm']
        b0 = self.baseflow_info['b0']
        bc = np.around((1 - br) * bfvf, 4)
        print(f'Baseflow-> BFVF10: {bfvf10} | Scaling: {bfvf_f} | BFVF: {bfvf} | Br: {br} | BC: {bc}')
        insertions = {}
        if self.baseflow_locations:
            for location in self.baseflow_locations:
                insertions[location] = {
                    "bf_line": f'PRINT.{location}: B0={b0} BR={br} BC={bc} BM={bm}\n',
                    "search_key": f'PRINT.{location.upper()}',
                    "used": False
                }
        else:
            insertions['default'] = {
                "bf_line": f'DEFAULT PARAMETERS: b0 = {b0} br = {br} bc = {bc} bm = {bm}\n',
                "search_key": 'DEFAULT PARAMETERS:',
                "used": False
            }

        # write baseflow_vec file to model output folder so simultaneous baseflow.vecs can be used
        filepath = os.path.join(self.output_folder, self.baseflow_vec_file)
        with open(filepath, 'w') as f:
            for i, line in enumerate(self.vec_lines):
                key_found = False
                for key, insertion in insertions.items():
                    search_key = insertion["search_key"]
                    bf_line = insertion['bf_line']
                    if line.strip().upper().startswith(search_key):
                        key_found = True
                        if key == 'default':
                            f.writelines(line)
                        f.writelines(bf_line)
                        insertion['used'] = True
                if key_found is not True:
                    f.writelines(line)

        # Check that the baseflow has been applied
        for key, insertion in insertions.items():
            if insertion['used'] is not True:
                Exception('Baseflow insertion not found:', key)

    def insert_baseflow_into_vec_old(self, bfvf10, z):
        bfvf_f = self.baseflow_info['bfvf10_factor']
        bfvf = self.get_bfvf(bfvf10 * bfvf_f, z)
        br = self.baseflow_info['br']
        bm = self.baseflow_info['bm']
        b0 = self.baseflow_info['b0']
        bc = np.around((1 - br) * bfvf, 4)
        print(f'Baseflow-> BFVF10: {bfvf10} | Scaling: {bfvf_f} | BFVF: {bfvf} | Br: {br} | BC: {bc}')
        bf_line = f'DEFAULT PARAMETERS: b0 = {b0} br = {br} bc = {bc} bm = {bm}\n'
        filepath = os.path.join(self.model_folder, self.vec_file)
        with open(filepath, 'r') as f:
            all_lines = f.readlines()
        filepath = os.path.join(self.output_folder, self.baseflow_vec_file) ## write baseflow_vec file to model output folder so simultaneous baseflow.vecs can be used
        with open(filepath, 'w') as f:
            for i, line in enumerate(all_lines):
                if line.upper().startswith('DEFAULT PARAMETERS:'):
                    f.writelines(line)
                    f.writelines(bf_line)
                else:
                    f.writelines(line)

    def get_bfvf(self, bfvf10, z):
        a = -0.02079
        b = -0.1375
        c = 0.20788
        r = a * z**2 + b * z + c
        self.bfvf = np.around(10**r * bfvf10, 5)
        # self.apply_baseflow = True
        return self.bfvf

    def set_volume_below_fsl(self, volume):
        self.volume_below_fsl = volume
        if self.dam_routing_method == 'level':
            self.set_initial_lake_level()
        else:
            last_element = self.header[-1]
            if 'set ADV_below=' in last_element:
                del self.header[-1]
            self.header.append(f'set ADV_below={self.volume_below_fsl}')

    def set_initial_lake_level(self):
        search_vol = self.full_supply_volume - self.volume_below_fsl
        search_level = self.get_level_from_volume(search_vol)
        last_element = self.header[-1]
        if 'set initial_lake_level=' in last_element:
            del self.header[-1]
        self.header.append(f'set initial_lake_level={search_level}')

    def get_level_from_volume(self, search_vol):
        els = self.storage_curve.copy()
        els.set_index('V', inplace=True)

        if search_vol in els.index:
            search_level = els.loc[search_vol, 'EL']
        else:
            els.loc[search_vol] = np.nan
            els.sort_index(inplace=True)
            els.interpolate(method='cubicspline', inplace=True)
            search_level = np.around(els.loc[search_vol, 'EL'], 3)
        return search_level

    def get_volume_from_level(self, search_level):
        els = self.storage_curve.copy()
        els.set_index('EL', inplace=True)
        if search_level in els.index:
            search_volume = els.loc[search_level, 'V']
        else:
            els.loc[search_level] = np.nan
            els.sort_index(inplace=True)
            els.interpolate(method='cubicspline', inplace=True)
            search_volume = np.around(els.loc[search_level, 'V'], 3)
        return search_volume

    def get_els_data(self):
        filepath = os.path.join(self.ratings_folder, self.els_file)
        print('Opening the URBS els file:', filepath)
        els = pd.read_csv(filepath)
        els.drop_duplicates(subset=['V'], keep='last', inplace=True)
        els.drop_duplicates(subset=['EL'], keep='last', inplace=True)
        # els.set_index('V', inplace=True)

        # elevations = els['EL']
        # volumes = els['V']
        # storage_curve = interpolate.interp1d(x=volumes, y=elevations, kind='cubic')
        return els

    def duration_string(self, duration):
        duration_test = float(duration)
        duration_str = str(duration)
        if not duration_test.is_integer():
            duration_str = duration_str.replace('.', 'p')
        return duration_str

    def get_simulation_period(self, duration):
        duration = str(duration)  #  json keys are always strings
        if duration in self.simulation_periods.keys():
            simulation_period = self.simulation_periods[duration]
        else:
            simulation_period = float(duration) * 2
            print(f'Duration period for {duration} hour storm not provided in the URBS config file.')
            print(f'Using twice the storm duration: {simulation_period} hours.')
        return simulation_period

    def find_dam_routing_line(self, dam_location):
        filepath = os.path.join(self.model_folder, self.vec_file)
        line_number = -99
        with open(filepath) as f:
            for i, line in enumerate(f):
                if line.startswith('DAM ROUTE'):
                    if dam_location in line:
                        line_number = i
                        routing_line = line
                        print(f'Found the dam routing line in the URBS vec file on line {line_number}')
                        break
        return {"line": routing_line, "line_number": line_number}

    def restore_initial_lake_level(self):
        if self.dam_routing_line:
            filepath = os.path.join(self.model_folder, self.vec_file)
            with open(filepath, 'r') as f:
                all_lines = f.readlines()
            with open(filepath, 'w') as f:
                for i, line in enumerate(all_lines):
                    if i == self.dam_routing_line['line_number']:
                        f.writelines(self.dam_routing_line['line'])
                    else:
                        f.writelines(line)

    def get_dam_routing_method(self):
        line = self.dam_routing_line['line']
        parts = line.split()
        method = None
        print('Checking the URBS model vec file to work out if using volume or level method...')
        for part in parts:
            if part.lower().startswith('fsl='):
                method = 'level'
                self.full_supply_level = float(part.split('=')[1])
                print('Found the full supply level in the URBS vec file:', self.full_supply_level)
            elif part.lower().startswith('datafile='):
                self.els_file = part.split('=')[1]
            elif part.lower().startswith('vbf='):
                method = 'volume'
                break
        if method is None:
            print('\nThe dam routing line in URBS is:')
            print(line)
            raise Exception('The keywords "FSL=" or "VBF=" were not found in the dam routing line')
        print(f'Found the dam routing method:', method)
        return method

    def adjust_initial_lake_level(self, lake_level):
        filepath = os.path.join(self.model_folder, self.vec_file)
        old_line = self.dam_routing_line['line']
        left_index = old_line.find(' il=') + 4
        right_index = old_line.find(" ", left_index)
        left_str = old_line[:left_index]
        right_str = old_line[right_index:]
        new_line = f'{left_str}{lake_level}{right_str}'

        with open(filepath, 'r') as f:
            all_lines = f.readlines()
        with open(filepath, 'w') as f:
            for i, line in enumerate(all_lines):
                if i == self.dam_routing_line['line_number']:
                    f.writelines(new_line)
                else:
                    f.writelines(line)

    def copy_catchment_data_file(self):
        filepath = os.path.join(self.model_folder, self.vec_file)
        with open(filepath, 'r') as f:
            for line in f:
                if line.lower().startswith('catchment data file'):
                    dat_file = line.split('=')[1].strip()
                    break
        src = os.path.join(self.model_folder, dat_file)
        dst = os.path.join(self.output_folder, dat_file)
        shutil.copy(src, dst)
        return
    
    def set_run_parameters(self, alpha=0, m_exponent=0, beta=0, initial_loss=0, continuing_loss=0):
        if alpha > 0:
            self.paramter_string = f'{alpha}'
        if m_exponent > 0:
            self.paramter_string = f'{self.paramter_string} {m_exponent}'
        if beta > 0:
            self.paramter_string = f'{self.paramter_string} {beta}'
        if initial_loss > 0:
            self.paramter_string = f'{self.paramter_string} {initial_loss}'
        if continuing_loss > 0:
            self.paramter_string = f'{self.paramter_string} {continuing_loss}'

    def write_ensemble_batch(self, aeps, storm_durations, ensembles):
        vec_path = os.path.join(self.model_folder, self.vec_file)
        batch_path = os.path.join(self.output_folder, self.batch_file)
        f = open(batch_path, 'w')
        f.writelines([txt + '\n' for txt in self.header])
        for ensemble in ensembles:
            for duration in storm_durations:
                for aep in aeps:
                    storm_file = f'ari{aep}_E{ensemble}.{duration}'
                    storm_path = os.path.join(self.storms_folder, storm_file)
                    result_naming = f'{self.result_prefix}{aep}{duration}_E{ensemble}'
                    self.results_list.append(result_naming)
                    f.write(f'"{self.urbs_exe}" {vec_path} {storm_path} {result_naming} {self.paramter_string}\n')
        f.close()

    def run_ensemble_storms(self, aeps, storm_durations, ensembles):
        batch_path = os.path.join(self.output_folder, self.batch_file)
        # print(f'\nRunning AEP: {aep} | Duration: {storm_duration}')
        self.write_ensemble_batch(aeps, storm_durations, ensembles)
        p = subprocess.Popen([batch_path], shell=True)
        p.wait()
        print()
        print(self.results_list)
        print('see h files for levels, p and pqh files for peak summary, q for flows')
        # collate h and q for specific locations across all ensembles then delete all outputs

    def run_storm(self, storm_name, result_name):
        if self.apply_baseflow:
            vec_path = os.path.join(self.output_folder, self.baseflow_vec_file)
        else:
            vec_path = os.path.join(self.model_folder, self.vec_file)
        batch_path = os.path.join(self.output_folder, self.batch_file)
        storm_path = os.path.join(self.storms_folder, storm_name)
        try:
            # Try opening the batch file
            f = open(batch_path, 'w')
        except IOError:
            # If it doesn't open, wait and try one more time
            time.sleep(0.5)
            f = open(batch_path, 'w')
        f.writelines([txt + '\n' for txt in self.header])
        f.write(f'"{self.urbs_exe}" {vec_path} {storm_path} {result_name} {self.paramter_string}\n')
        f.close()
        p = subprocess.Popen([batch_path], shell=True)
        p.wait()

    def create_storm_file(self, rainfall_depths, temporal_pattern, data_interval, filename, run_duration,
                          storm_duration, storm_duration_excl_pb, initial_loss, continuing_loss, ari=1, ensemble=0):
        filepath = os.path.join(self.storms_folder, filename)
        f = open(filepath, 'w')
        header = f'Storm duration of {storm_duration_excl_pb} hours'
        header = f'{header} for ARI {ari} Ensemble {ensemble}'
        if self.apply_baseflow:
            header = f'{header} BFVF={self.bfvf}'
        # f.write(f'Storm file: {filename}\n ensemble {ensemble}')
        f.write(header)
        f.write('\nDesign Run\n')
        time_increment = self.time_increment
        if self.time_increment_override is not None:
            if float(storm_duration_excl_pb).is_integer():
                duration_text = str(int(storm_duration_excl_pb))
            else:
                duration_text = str(storm_duration_excl_pb)
            if duration_text in self.time_increment_override.keys():
                time_increment = self.time_increment_override[duration_text]
        f.write(f'Time Increment: {time_increment}\n')
        f.write(f'Run Duration  : {float(run_duration)}\n')
        f.write(f'Storm Duration: {float(storm_duration)}\n')
        f.write('Pluviograph.\n')
        f.write(f'Data Interval: {float(data_interval)}\n')
        max_len = 8
        j = 0
        pattern_block = np.array2string(temporal_pattern.to_numpy(),
                                        precision=3,
                                        floatmode='fixed')
        pattern_block = pattern_block.replace('[', ' ')
        pattern_block = pattern_block.replace(']', ' ')
        f.write(pattern_block)
        f.write('\nRain on Subareas:\n')
        rain_block = np.array2string(rainfall_depths.to_numpy(),
                                     precision=2,
                                     floatmode='fixed')
        rain_block = rain_block.replace('[', ' ')
        rain_block = rain_block.replace(']', ' ')
        f.write(rain_block)
        f.write('\nLoss: Uniform Continuing\n')
        f.write(f'IL:  {initial_loss}\n')
        f.write(f'CL:  {continuing_loss}\n')
        f.close()

    def get_max_results(self, result_name, max_keys):
        result_path = os.path.join(self.output_folder, f'{result_name}.p')
        f = open(result_path)
        all_lines = f.readlines()
        all_lines = all_lines[9: -1]
        level = 9999
        for line in all_lines:
            if line.startswith(max_keys['inflow']):
                    inflow = float(line[46:55])
            if line.startswith(max_keys['level']):
                # if not 'N/A' in line:
                level_new = float(line[64:70])
                if level_new < level:
                    level = level_new
            if line.startswith(max_keys['outflow']):
                # if not 'N/A' in line:
                outflow = float(line[46:55])
        f.close()
        return {'inflow': inflow, 'level': level, 'outflow': outflow}

    def get_hydrographs(self, result_name, max_keys, simid):
        sim_label = 'sim_{}'.format(str(simid).zfill(5))
        # Get the current flows
        result_path = os.path.join(self.output_folder, f'{result_name}.q')
        f = open(result_path)
        all_lines = f.readlines()
        f.close()
        header_line = all_lines[8]
        headers = header_line.split()
        headers.insert(0, 'Time')
        
        skiprow = 12
        while True:
            line = all_lines[skiprow].split()
            if line[0] == 'Time' and line[1] == 'Cumecs':
                skiprow += 2
                break
            skiprow += 1
            
            if skiprow >= len(all_lines):
                skiprow = 16            # Fallback option if conditions not detected
                break
        
        # print('\nHeaders: ', headers)
        # df = pd.read_fwf(result_path, skiprows=16, skipfooter=5, names=headers)
        df = pd.read_csv(result_path, skiprows=skiprow, skipfooter=5, sep="\s+", names=headers, engine = 'python')
        df.set_index('Time', inplace=True)
        # print('\nDataframe:\n', df)
        inflow_label = max_keys['inflow'][:8]                   # URBS automatically truncates label to 8 char
        current_inflow = df[inflow_label]
        current_inflow.name = sim_label

        # self.inflows = pd.concat([self.inflows, current_inflow], axis=1)
        self.inflows.append(current_inflow)
        outflow_label = max_keys['outflow'][:8]                 # URBS automatically truncates label to 8 char
        current_outflow = df[outflow_label]
        current_outflow.name = sim_label
        # self.outflows = pd.concat([self.outflows, current_inflow], axis=1)
        self.outflows.append(current_outflow)

        # Get the current levels
        result_path = os.path.join(self.output_folder, f'{result_name}.h')
        f = open(result_path)
        all_lines = f.readlines()
        f.close()
        header_line = all_lines[8]
        headers = header_line.split()
        level_label = max_keys['level'][:8]                 # URBS automatically truncates label to 8 char
        # Deal with double label
        count_label = sum(header == level_label for header in headers)
        if count_label > 1:
            headers[-1] = '{}_1'.format(headers[-1])
        headers.insert(0, 'Time')
        
        skiprow = 12
        while True:
            line = all_lines[skiprow].split()
            if line[0] == 'Time' and line[1] == 'metres':
                skiprow += 2
                break
            skiprow += 1
            
            if skiprow >= len(all_lines):
                skiprow = 16            # Fallback option if conditions not detected
                break
            
        # read the data
        # df = pd.read_fwf(result_path, skiprows=16, skipfooter=1, names=headers)
        df = pd.read_csv(result_path, skiprows=skiprow, skipfooter=1, sep="\s+", names=headers, engine = 'python')
        df.set_index('Time', inplace=True)
        current_level = df[level_label]
        current_level.name = sim_label
        # self.levels = pd.concat([self.levels, current_level], axis=1)
        self.levels.append(current_level)

    def store_hydrographs(self, filename):
        folder = os.path.dirname(filename)
        results_folder = os.path.abspath(os.path.join(folder, '../urbs_results'))
        os.makedirs(results_folder, exist_ok=True)
        filename = os.path.basename(filename)
        
        filepath = os.path.join(results_folder, f'{filename}_inflows.csv')
        # filepath = os.path.join(folder, '../urbs_results', f'{filename}_inflows.csv')
        print('\nWriting hydrograph results to file:', os.path.abspath(filepath))
        # self.inflows.to_csv(filepath)
        pd.concat(self.inflows, axis=1).to_csv(filepath)
        
        filepath = os.path.join(results_folder, f'{filename}_outflows.csv')
        # filepath = os.path.join(folder, '../urbs_results', f'{filename}_outflows.csv')
        print('Writing hydrograph results to file:', os.path.abspath(filepath))
        # self.outflows.to_csv(filepath)
        pd.concat(self.outflows, axis=1).to_csv(filepath)
        
        filepath = os.path.join(results_folder, f'{filename}_levels.csv')
        # filepath = os.path.join(folder, '../urbs_results', f'{filename}_levels.csv')
        print('Writing hydrograph results to file:', os.path.abspath(filepath))
        # self.levels.to_csv(filepath)
        pd.concat(self.levels, axis=1).to_csv(filepath)
    
    # Reinstated
    def mop_storm_file(self, storm_name, skip=0, sim_id=0):
        if skip > 0:
            if sim_id > skip - 1:
                filepath = os.path.join(self.storms_folder, storm_name)
                if os.path.exists(filepath):
                    os.remove(filepath)
    
    # Reinstated
    def mop_results_files(self, results_name, which='all'):
        filepath = os.path.join(self.output_folder, results_name)

        if which == 'all':
            extensions = ['a', 'b', 'cc', 'csv', 'e', 'fs', 'hc', 'o', 'osd', 'pqh', 'prm', 'h', 'q', 'p', 'il']
        else:
            extensions = ['a', 'b', 'cc', 'csv', 'e', 'fs', 'hc', 'o', 'osd', 'pqh', 'prm']

        for ext in extensions:
            fullpath = f'{filepath}.{ext}'
            if os.path.exists(fullpath):
                os.remove(fullpath)

    def move_results(self, results_name):
        all_files = glob.glob(os.path.join(self.output_folder, f'{results_name}.*'), recursive=True)
        for filepath in all_files:
            dst_path = os.path.join(self.results_folder, os.path.basename(filepath))
            shutil.move(filepath, dst_path)
