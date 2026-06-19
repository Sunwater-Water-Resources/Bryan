import os
import subprocess
import json
import numpy as np
import glob
import shutil
from pathlib import Path
import pandas as pd


class RorbModel:
    def __init__(self, config_file, batch_file='_run.bat', dam_location=None, sub_folder=None):
        # open the config file and get contents
        f = open(config_file)
        config_data = json.load(f)
        f.close()
        self.full_supply_volume = config_data['full_supply_volume']
        self.rorb_exe = config_data['model_exe']
        self.model_folder = os.path.normpath(os.path.join(os.path.dirname(config_file), config_data['model_folder']))
        self.cat_file = config_data['cat_file']
        self.par_file = config_data['par_file']
        # self.ratings_folder = config_data['ratings_folder']
        
        # Create sub-folder for output files. If sub-folder exists from previous run then delete first
        if sub_folder is not None:
            self.storms_folder = os.path.normpath(os.path.join(os.path.dirname(config_file),
                                                               config_data['storms_folder'],
                                                               sub_folder))
            # Create sub-folder for storm files. If sub-folder exists from previous run then delete first
            if os.path.exists(self.storms_folder): shutil.rmtree(self.storms_folder)
            os.makedirs(self.storms_folder)
        else:
            self.storms_folder = os.path.normpath(os.path.join(os.path.dirname(config_file), config_data['storms_folder']))
            os.makedirs(self.storms_folder, exist_ok=True)
            
        # Copy catg file to storms_folder
        original_catg = os.path.join(self.model_folder, self.cat_file)
        shutil.copy(original_catg, self.storms_folder)
            
        if 'results_folder' in config_data.keys():
            self.results_folder = os.path.normpath(os.path.join(os.path.dirname(config_file), config_data['results_folder']))
        # self.store_tuflow = config_data['store_tuflow']
        # self.result_prefix = config_data['result_prefix']
        # self.time_increment = config_data['time_increment']
        self.batch_file = batch_file
        # self.paramter_string = ''
        self.results_list = []
        self.inflows = []       # List of series of hydrographs.
        self.outflows = []
        # self.levels = []

        # self.Kc = config_data['Kc']
        # self.beta = config_data['beta']
        # self.m_exponent = config_data['m_exponent']

        # self.set_run_parameters(alpha=self.alpha, beta=self.beta, m_exponent=self.m_exponent)
        self.simulation_periods = config_data['simulation_periods']

        # # set up the header for the batch file
        # self.header = []
        # self.header.append(f'cd {self.model_folder}')
        # self.header.append(f'del {self.model_folder}\\urbserr.log')
        # self.header.append(f'del {self.model_folder}\\urbsout.log')
        # self.header.append('set URBS_LOGF=TRUE')
        # self.header.append(f'set URBS_LOGD={self.model_folder}')
        # self.header.append(f'set URBS_RATS={self.ratings_folder}')
        # if config_data['store_tuflow']:
        #     self.header.append('set URBS_TFLW=TRUE')

        # get the line number where the dam routing is done
        if dam_location:
            self.dam_routing_line = self.find_dam_routing_line(dam_location)
        
        # Create batch file
        self.par_lines, self.num_isa = self.setup_par_file()
        
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
            print(f'Duration period for {duration} hour storm not provided in the model config file.')
            print(f'Using twice the storm duration: {simulation_period} hours.')
        return simulation_period

    def find_dam_routing_line(self, dam_location):
        filepath = os.path.join(self.storms_folder, self.cat_file)
        line_number = -99
        with open(filepath) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith('16'):
                    if lines[i+1].rstrip() == dam_location:
                        line_number = i+2
                        routing_line = lines[i+2]
                        print(f'Found the dam routing line in the RORB cat file on line {line_number}')
                        break
            # for i, line in enumerate(f):
            #     if line.startswith('DAM ROUTE'):
            #         if dam_location in line:
            #             line_number = i
            #             routing_line = line
            #             print(f'Found the dam routing line in the URBS vec file on line {line_number}')
            #             break
        return {"line": routing_line, "line_number": line_number}


    # Deprecated
    def restore_initial_lake_level(self):
        if self.dam_routing_line:
            filepath = os.path.join(self.storms_folder, self.cat_file)
            with open(filepath, 'r') as f:
                all_lines = f.readlines()
            with open(filepath, 'w') as f:
                for i, line in enumerate(all_lines):
                    if i == self.dam_routing_line['line_number']:
                        f.writelines(self.dam_routing_line['line'])
                    else:
                        f.writelines(line)
                        
        # restore PAR file
        par_path = os.path.join(self.storms_folder, self.par_file)
        with open(par_path, 'w') as f:
            for line in self.par_lines:
                f.write(line)

    def set_volume_below_fsl(self, volume):
        if volume == 0.0:
            drawdown = 0.0
        else:
            drawdown = -1* volume
        filepath = os.path.join(self.storms_folder, self.cat_file)
        old_line = self.dam_routing_line['line']
        line_parts = old_line.split(',')                    # Split line into parts by commas
        line_parts[1] = ' {:.1f}'.format(drawdown)          # Replace drawdown value
        new_line = ','.join(line_parts)                     # Join line parts together again

        with open(filepath, 'r') as f:
            all_lines = f.readlines()
        with open(filepath, 'w') as f:
            for i, line in enumerate(all_lines):
                if i == self.dam_routing_line['line_number']:
                    f.writelines(new_line)
                else:
                    f.writelines(line)
        
    def adjust_initial_lake_level(self, lake_level):                 # Note: method deprecated at REV 46
        filepath = os.path.join(self.storms_folder, self.cat_file)
        old_line = self.dam_routing_line['line']
        line_parts = old_line.split(',')                    # Split line into parts by commas
        line_parts[1] = ' {:.3f}'.format(lake_level)        # Replace drawdown value
        new_line = ','.join(line_parts)                     # Join line parts together again

        with open(filepath, 'r') as f:
            all_lines = f.readlines()
        with open(filepath, 'w') as f:
            for i, line in enumerate(all_lines):
                if i == self.dam_routing_line['line_number']:
                    f.writelines(new_line)
                else:
                    f.writelines(line)

    # def set_run_parameters(self, alpha=0, m_exponent=0, beta=0, initial_loss=0, continuing_loss=0):
    #     if alpha > 0:
    #         self.paramter_string = f'{alpha}'
    #     if m_exponent > 0:
    #         self.paramter_string = f'{self.paramter_string} {m_exponent}'
    #     if beta > 0:
    #         self.paramter_string = f'{self.paramter_string} {beta}'
    #     if initial_loss > 0:
    #         self.paramter_string = f'{self.paramter_string} {initial_loss}'
    #     if continuing_loss > 0:
    #         self.paramter_string = f'{self.paramter_string} {continuing_loss}'

    def write_ensemble_batch(self, aeps, storm_durations, ensembles):
        vec_path = os.path.join(self.storms_folder, self.vec_file)
        batch_path = os.path.join(self.storms_folder, self.batch_file)
        f = open(batch_path, 'w')
        f.writelines([txt + '\n' for txt in self.header])
        for ensemble in ensembles:
            for duration in storm_durations:
                for aep in aeps:
                    storm_file = f'ari{aep}_E{ensemble}.{duration}'
                    storm_path = os.path.join(self.storms_folder, storm_file)
                    result_naming = f'{self.result_prefix}{aep}{duration}_E{ensemble}'
                    self.results_list.append(result_naming)
                    f.write(f'{self.urbs_exe} {vec_path} {storm_path} {result_naming} {self.paramter_string}\n')
        f.close()

    def run_ensemble_storms(self, aeps, storm_durations, ensembles):
        batch_path = os.path.join(self.storms_folder, self.batch_file)
        # print(f'\nRunning AEP: {aep} | Duration: {storm_duration}')
        self.write_ensemble_batch(aeps, storm_durations, ensembles)
        p = subprocess.Popen([batch_path], shell=True)
        p.wait()
        print()
        print(self.results_list)
        # print('see h files for levels, p and pqh files for peak summary, q for flows')
        # collate h and q for specific locations across all ensembles then delete all outputs
        
    def setup_par_file(self):
        batch_path = os.path.join(self.storms_folder, self.batch_file)
        par_path = os.path.join(self.model_folder, self.par_file)
        
        with open(par_path, 'r') as f:
            all_lines = f.readlines()
        
        # get number of interstation area
        num_isa = all_lines[6].split(':')[1]
        num_isa = int(num_isa)
        
        par_path = os.path.join(self.storms_folder, self.par_file)
        with open(batch_path, 'w') as f:
            f.write(f'"{self.rorb_exe}" "{par_path}"')
        
        return all_lines, num_isa

    def run_storm(self, storm_name, initial_loss, continuing_loss):
        cat_path = os.path.join(self.storms_folder, self.cat_file)
        par_path = os.path.join(self.storms_folder, self.par_file)
        batch_path = os.path.join(self.storms_folder, self.batch_file)
        storm_path = os.path.join(self.storms_folder, storm_name)
        num_isa = self.num_isa
        
        with open(par_path, 'w') as f:
            f.write('# BEGIN\n')
            f.write(f'Cat file :{cat_path}\n')
            f.write(f'Stm file :{storm_path}\n')
            for i in range(3, 8 + num_isa):
                f.write(self.par_lines[i])
            for i in range(8 + num_isa, 2 * num_isa + 8):
                f.write(f'{":":>10}{initial_loss:.1f},{continuing_loss:.1f}\n')
            f.write('# END')
        
        # f = open(batch_path, 'w')
        # f.writelines([txt + '\n' for txt in self.header])
        # f.write(f'{self.urbs_exe} {vec_path} {storm_path} {result_name} {self.paramter_string}\n')
        # f.close()
        p = subprocess.Popen([batch_path], shell=True)
        p.wait()
        
        result_name = f'{Path(cat_path).stem}_{Path(storm_path).stem}'      # Reproduce as per RORB
        return result_name

    def create_storm_file(self, rainfall_depths, temporal_pattern, data_interval, filename, run_duration,
                          storm_duration = 0, initial_loss = 0, continuing_loss = 0):   # storm_duration and losses not used in RORB storm file, arguments retained for simplicity of Simulator code (i.e. don't need to add if block)
        filepath = os.path.join(self.storms_folder, filename)
        f = open(filepath, 'w')
        f.write(f'Storm file: {filename}\n')
        f.write('DESIGN\n')
        
        run_increment = int(round(run_duration / data_interval, 0))
        storm_increments = temporal_pattern.count()
        f.write(f'{data_interval}, {run_increment}, 1, 1, 1, -99.\n')
        f.write(f'0, {storm_increments}\n')
        f.write('Temporal pattern (% of depth)\n')
        
        # f.write(f'Time Increment: {self.time_increment}\n')
        # f.write(f'Run Duration  : {float(run_duration)}\n')
        # f.write(f'Storm Duration: {float(storm_duration)}\n')
        # f.write('Pluviograph.\n')
        # f.write(f'Data Interval: {float(data_interval)}\n')
        # max_len = 8
        # j = 0
        pattern_block = np.array2string(temporal_pattern.to_numpy(),
                                        separator=',',
                                        precision=3,
                                        floatmode='fixed',
                                        formatter={'float':lambda x: "%.3f" % x})
        pattern_block = pattern_block.replace('[', ' ')
        pattern_block = pattern_block.replace(']', ', -99.0')
        f.write(pattern_block)
        f.write('\nC Rain on Subareas:\n')
        
        rain_block = np.array2string(rainfall_depths.to_numpy(),
                                     separator=',',
                                     floatmode='fixed',
                                     formatter={'float':lambda x: "%.2f" % x})
        rain_block = rain_block.replace('[', ' ')
        rain_block = rain_block.replace(']', ', -99.0')
        f.write(rain_block)
        # f.write('\nLoss: Uniform Continuing\n')
        # f.write(f'IL:  {initial_loss}\n')
        # f.write(f'CL:  {continuing_loss}\n')
        f.close()

    def get_max_results(self, result_name, max_keys):
        result_path = os.path.join(self.storms_folder, f'{result_name}.out')
        with open(result_path, 'r') as f:
            all_lines = f.readlines()
        
        if 'storage_name' in max_keys:
            dam_location = max_keys['storage_name']
            level = np.nan                              # If dam does not fill, level = np.nan
            for i, line in enumerate(all_lines):
                if line.rstrip() == f' Results of routing through special storage {dam_location}':
                    data_line = all_lines[i+1]
                    level = float(data_line.split('=')[1].split()[0])
                    continue
                if line.rstrip() == f' *** Special storage :   {dam_location}':
                    data_line = all_lines[i+4]
                    data_line = data_line[20:].split()
                    outflow = float(data_line[0])
                    inflow = float(data_line[1])
                    break
            return {'inflow': inflow, 'level': level, 'outflow': outflow}
        else:
            print_name = max_keys['print_name']
            for i, line in enumerate(all_lines):
                if line.startswith(' *** Calc. hyd. for ungauged interstation site at:'):
                    hyd_name = line.split(sep=':', maxsplit=1)[1].strip()
                elif line.startswith(' *** Calculated hydrograph,'):
                    hyd_name = line.split(sep=',', maxsplit=1)[1].strip()
                else:
                    continue
                
                if hyd_name == print_name:
                    data_line = all_lines[i+4]
                    data_line = data_line[20:].strip()
                    inflow = float(data_line)
                    break
            return {'inflow': inflow}

    def get_hydrographs(self, result_name, max_keys, simid):
        sim_label = 'sim_{}'.format(str(simid).zfill(5))
        # Get the current flows
        result_path = os.path.join(self.storms_folder, f'{result_name}.out')
        f = open(result_path, encoding='unicode_escape')
        all_lines = f.readlines()
        f.close()
        
        for i, line in enumerate(all_lines):
            if line != ' Hydrograph summary\n':
                continue
            if all_lines[i+3] != ' Site  Description\n':
                raise Exception(f'error parsing file {result_path}')
            
            if 'storage_name' in max_keys:
                dam_location = max_keys['storage_name']
                j = 0
                while True:
                    j += 1
                    next_line = all_lines[i+3+j]
                    if(next_line[7:] != f'Special storage :   {dam_location} - Outflow\n'):
                        continue
                    next_line = all_lines[i+4+j]
                    if(next_line[7:] != f'Special storage :   {dam_location} - Inflow\n'):
                        raise Exception(f'error parsing file {result_path}')
                    hyd_no = j
                    break
            else:
                print_name = max_keys['print_name']
                j = 0
                while True:
                    j += 1
                    next_line = all_lines[i+3+j]
                    if(next_line[7:] == f'Calc. hyd. for ungauged interstation site at: {print_name}\n'):
                        hyd_no = j
                        break
                    elif(next_line[7:] == f'Calculated hydrograph,  {print_name}\n'):
                        hyd_no = j
                        break
            
            while True:
                j += 1
                next_line = all_lines[i+3+j]
                if(next_line != '\n'):
                    continue
                next_line = all_lines[i+4+j]
                if(next_line.split()[:3] != ['Inc', 'Time', 'Hyd001']):
                    raise Exception(f'error parsing file {result_path}')
                # skiprows = i+5+j
                skiprows = i+4+j
                break
            break
        
        # df = pd.read_fwf(result_path, skiprows=skiprows, encoding = 'unicode_escape').set_index('Time')
        df = pd.read_csv(result_path, skiprows=skiprows, encoding='unicode_escape', sep="\s+").set_index('Time')
        
        if 'storage_name' in max_keys:
            inflow_label = f'Hyd{hyd_no + 1:03}'
            current_inflow = df[inflow_label]
            current_inflow.name = sim_label
            self.inflows.append(current_inflow)
            
            outflow_label = f'Hyd{hyd_no:03}'
            current_outflow = df[outflow_label]
            current_outflow.name = sim_label
            self.outflows.append(current_outflow)
        
            # TODO: RORB does not output level hydrograph. 
            # Could be interpolated from S-Q curve using outflow hydrograph, but levels below spill level would not be produce.
            # Could also be derived from starting level and inflow hydrograph.
            
        else:
            inflow_label = f'Hyd{hyd_no:03}'
            current_inflow = df[inflow_label]
            current_inflow.name = sim_label
            self.inflows.append(current_inflow)

    def store_hydrographs(self, filename):
        folder = os.path.dirname(filename)
        results_folder = os.path.abspath(os.path.join(folder, '../rorb_results'))
        os.makedirs(results_folder, exist_ok=True)
        filename = os.path.basename(filename)
        
        filepath = os.path.join(results_folder, f'{filename}_inflows.csv')
        print('\nWriting hydrograph results to file:', os.path.abspath(filepath))
        pd.concat(self.inflows, axis = 1).to_csv(filepath)
        
        if len(self.outflows) > 0:
            filepath = os.path.join(results_folder, f'{filename}_outflows.csv')
            print('\nWriting hydrograph results to file:', os.path.abspath(filepath))
            pd.concat(self.outflows, axis = 1).to_csv(filepath)
        
        # filepath = os.path.join(results_folder, f'{filename}_levels.csv')
        # pd.concat(self.levels, axis = 1).to_csv(filepath)

    def mop_storm_file(self, storm_name, skip=0, sim_id=0):
        if skip > 0:
            if sim_id > skip - 1:
                filepath = os.path.join(self.storms_folder, storm_name)
                if os.path.exists(filepath):
                    os.remove(filepath)

    def mop_results_files(self, result_name, which='all'):
        result_path = os.path.join(self.storms_folder, f'{result_name}.out')
        if os.path.exists(result_path):
            os.remove(result_path)

    def move_results(self, results_name):
        os.makedirs(self.results_folder, exist_ok=True)
        all_files = glob.glob(os.path.join(self.storms_folder, f'{results_name}.*'), recursive=True)
        for filepath in all_files:
            dst_path = os.path.join(self.results_folder, os.path.basename(filepath))
            shutil.move(filepath, dst_path)
