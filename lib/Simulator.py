import matplotlib.pyplot as plt
from lib.MCScheme import SampleScheme
from lib.StormGenerator import StormBurst
from lib.Lake import LakeConditions
from lib.URBSmodel import UrbsModel
from lib.RORBmodel import RorbModel
from lib.ClimateChange import ClimateAdjustment
from lib.EnbScheme import Ensemble
from scipy.special import ndtr
from scipy import integrate
import pandas as pd
import os
import numpy as np
import sys
from datetime import datetime
import time
import json
from scipy import interpolate


class Simulator:
    def __init__(self, parameters, filepaths, test_iterations=0):
        # set up a testing parameter to exit early when testing the code:
        self.test = test_iterations  # exist after this many realizations. Use zero for not production.
        self.start = time.time()

        if parameters['Run models'].lower() == 'yes':
            self.do_runs = True
            # Create and set up the log file
            # log_file = self.config_data['log_file']
            log_file = parameters['Log file']
            print('Creating log file:', log_file)
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            # log_file = log_file.replace('.csv', 'analysis.csv')
            sys.stdout = Logger(log_file)
            life_of_brian()
        else:
            self.do_runs = False

        # Keep track of whether this is Monte Carlo or Ensemble method
        self.method = parameters['Method']

        # Get the filepaths
        self.filepaths = filepaths

        # Get the config data
        config_file = parameters['Config file']
        f = open(config_file)
        self.config_data = json.load(f)
        f.close()

        # Set up the runs - if runnning models
        self.duration = parameters['Duration']
        self.outputfile = parameters['Output file']
        self.apply_baseflow = None
        if 'Baseflow' in parameters.index:
            print('Found "Baseflow" key in the Sims List')
            if parameters['Baseflow'].lower() == 'yes':
                print('Baseflow switched on in Sims List.')
                self.apply_baseflow = True
            elif parameters['Baseflow'].lower() == 'no':
                print('Baseflow switched off in Sims List.')
                self.apply_baseflow = False
        else:
            print('"Baseflow" key not used in the Sims List - defaulting to the URBS config')

        # Set up the results analysis
        if parameters['Analyse results'].lower() == 'yes':
            self.do_analysis = True
        elif parameters['Analyse results'].lower() == 'no':
            self.do_analysis = False
        else:
            analyse_lst = str(parameters['Analyse results']).split(',')
            self.do_analysis = [x.strip().lower() for x in analyse_lst]
        
        self.do_volumes = False
        if 'Analyse volumes' in parameters.index:
            do_volumes = str(parameters['Analyse volumes']).lower()
            if do_volumes == 'yes':
                self.do_volumes = ['inflow', 'outflow']
            elif do_volumes in ['inflow', 'outflow']:
                self.do_volumes = [do_volumes]
                
        if 'Pre-burst method' in parameters.index:
            self.preburst_method = str(parameters['Pre-burst method']).lower().strip()
            if self.preburst_method == 'uniform':
                print('Applying uniform pre-burst')
        else:
            self.preburst_method = None

        # File management
        if parameters['Store hydrographs'].lower() == 'yes':
            self.store_hydrographs = True
        else:
            self.store_hydrographs = False

        if parameters['Mop up files'].lower() == 'yes':
            self.mop_files = True
        else:
            self.mop_files = False

        print("""!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Bryan: Sunwater's Design Flood Simulator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!""")

        # datetime object containing current date and time
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        print(f'\nSimulation started at: {dt_string}')
        
        # Initialise the hydrologic model
        model_info = self.initialise_model(self.filepaths['model_config'])
        self.model = model_info['model']
        self.model_type = model_info['model_type']
        self.max_keys = model_info['max_keys']
        self.bfvf10 = 0
        
        if not self.do_runs:
            return
        
        # Manage exclusions
        exclusions_line = str(parameters['Exclusions'])
        self.exclusions = self.set_exclusions(exclusions_line)

        # Manage the replication
        replicate_file = parameters['Replicate file']
        replicates_line = str(parameters['Replicates'])
        self.replicates = self.set_replicates(replicates_line)
        if any(x is True for x in self.replicates.values()):
            print('Replicates have been found')
            print('opening replication file', replicate_file)
            self.replicates_df = pd.read_csv(replicate_file, index_col=0)

        # Set up the lake
        self.lake = LakeConditions(parameters)
        self.lake.set_full_supply_volume(self.model.full_supply_volume)
        if self.lake.antecedent_type in ['fixed', 'fsv']:
            self.model.set_volume_below_fsl(self.lake.get_volume_below_fsl())

        # Set up the average rainfall losses for the complete storm
        self.storm_losses = {'IL': parameters['IL'],
                             'CL': parameters['CL']}

        # Set up the focal location
        focal_file = parameters.get('Focal subcatchments')
        if focal_file:
            print('Found the focal location file:', focal_file)
            self.focal_subcatchments = focal_file
        else:
            raise Exception('No focal subcatchments file found - check the sims list!')

        # set up the climate regime
        climate_config = filepaths['climate_config']
        if 'GWL' in parameters.keys():
            print('GWL key found in the simulation list - using the GWL method for climate change.')
            if parameters['GWL']:
                gwl = parameters['GWL']
                print(f'Using GWL of {gwl}°C.')
            else:
                print('No GWL given, using 0°C.')
                gwl = 0.0
            self.climate = ClimateAdjustment(config_file=climate_config, method='gwl', gwl=gwl)
        else:
            print('GWL key not found in the simulation list - using the SSP method for climate change.')
            self.climate = ClimateAdjustment(config_file=climate_config, method='ssp',
                                             year=parameters['Year'], ssp=parameters['SSP'])

        # Set up the percentiles when given
        self.given_percentiles = {}
        if 'IL percentile' in parameters:
            il_percentile = parameters['IL percentile']
            if pd.notna(il_percentile):
                print('\nThe percentile for initial losses has been provided:', il_percentile)
                self.given_percentiles['il'] = il_percentile
        if 'Preburst percentile' in parameters:
            pb_percentile = parameters['Preburst percentile']
            if pd.notna(pb_percentile):
                print('\nThe percentile for preburst proportion has been provided:', pb_percentile)
                self.given_percentiles['preburst'] = pb_percentile
        
    def set_replicates(self, replicates_line):
        # initialise
        replicates = {'rain_z': False,
                      'tp': False,
                      'storm_method': False,
                      'il_p': False,
                      'cl_p': False,
                      'preburst_p': False,
                      'preburst_tp': False,
                      'lake_z': False}

        # get the replicates
        replicates_line = replicates_line.strip()
        if replicates_line == '':
            print('No replicates found.')
        else:
            replicates_lst = replicates_line.split(',')
            replicates_lst = list(map(lambda x: x.strip(), replicates_lst))  # strip whitespace from elements
            for replicate in replicates_lst:
                if replicate == 'rz':
                    print('Replicating rainfall sampling')
                    replicates['rain_z'] = True
                elif replicate == 'tp':
                    print('Replicating temporal pattern number sampling')
                    replicates['tp'] = True
                elif replicate == 'stm':
                    print('Replicating storm method sampling')
                    replicates['storm_method'] = True
                elif replicate == 'ilp':
                    print('Replicating initial loss percentile sampling')
                    replicates['il_p'] = True
                elif replicate == 'clp':
                    print('Replicating continuing loss percentile sampling')
                    replicates['cl_p'] = True
                elif replicate == 'pbp':
                    print('Replicating preburst percentile sampling')
                    replicates['preburst_p'] = True
                elif replicate == 'pb_tp':
                    print('Replicating preburst temporal pattern sampling')
                    replicates['preburst_tp'] = True
                elif replicate == 'lz':
                    print('Replicating antecedent lake volume sampling')
                    replicates['lake_z'] = True
        return replicates

    def set_exclusions(self, exclusions_line):
        # initialise
        exclusions = {'preburst': False,
                      'embedded_burst_filter': False,
                      'pattern_d50_scaling': False,
                      'rainfall_uplift': False,
                      'loss_uplift': False,
                      'cl_sampling': False}

        # get the exclusions
        exclusions_line = exclusions_line.strip()
        if exclusions_line == '':
            print('No exclusions found.')
        else:
            exclusions_lst = exclusions_line.split(',')
            exclusions_lst = list(map(lambda x: x.strip(), exclusions_lst))  # strip whitespace from elements
            for exclusion in exclusions_lst:
                if exclusion == 'pb':
                    print('Excluding preburst temporal patterns')
                    exclusions['preburst'] = True
                elif exclusion == 'ebf':
                    print('Excluding the embedded burst filter')
                    exclusions['embedded_burst_filter'] = True
                elif exclusion == 'd50':
                    print('Excluding the temporal pattern D50 scaling for climate change')
                    exclusions['pattern_d50_scaling'] = True
                elif exclusion == 'ru':
                    print('Excluding the rainfall uplift for climate change')
                    exclusions['rainfall_uplift'] = True
                elif exclusion == 'lu':
                    print('Excluding the rainfall loss uplift for climate change')
                    exclusions['loss_uplift'] = True
                elif exclusion == 'clp':
                    print('Excluding stochastic sampling of continuing losses')
                    exclusions['cl_sampling'] = True
        return exclusions

    def print_elapsed_time(self):
        end = time.time()
        elapsed_time = np.around((end - self.start) / 60, 2)
        print(f'Elapsed time: {elapsed_time} minutes')

    def initialise_storm(self, filepath, storm_durations):
        # Create the storm and load some data - including losses
        storm = StormBurst(filepath)

        # List of burst durations to simulate 
        if self.method == 'monte carlo':
            storm_durations = [self.duration]   # Single duration as list
        elif self.method == 'ensemble':
            storm_durations = self.config_data['storm_durations']
        
        storm.load_subcatchment_areas(self.focal_subcatchments)
        storm.storm_initial_loss = self.storm_losses['IL']
        storm.continuing_loss = self.storm_losses['CL']
        print(f'\nUsing initial loss of {storm.storm_initial_loss} mm and continuing loss of {storm.continuing_loss} mm/h')
        do_preburst = not self.exclusions['preburst']

        # Set up the rainfall depths
        storm.import_rare_rainfall()  # importing IFD depths more frequent than 1 in 2,000 AEP
        storm.apply_areal_reduction()  # applying the areal reduction factors - BEFORE extreme rainfall!
        storm.skip_extreme_methods(storm_durations, do_preburst = do_preburst)
        storm.set_up_extreme_rainfall()  # importing the PMP depths and setting up extreme depths rarer than 1 in 2,000


        # Import the temporal patterns and associated D50 climate change weightings
        storm.import_arr_areal_patterns(storm_durations)
        storm.import_arr_point_patterns(storm_durations)
        storm.import_gsdm_temporal_patterns(storm_durations)
        storm.import_gtsmr_temporal_patterns(storm_durations)

        # Import the preburst patterns
        if do_preburst:
            if self.method == 'monte carlo':
                storm.import_preburst_patterns()
            elif self.method == 'ensemble':
                storm.import_preburst_patterns(median_only = True)

        # Get the BFVF10 if using baseflow
        if self.model_type == 'URBS':
            if self.apply_baseflow is not None:
                print('Baseflow command in urbs config overwritten by Sims List; apply:',
                      self.apply_baseflow)
                self.model.apply_baseflow = self.apply_baseflow
            else:
                print('Checking model config for application of baseflow')
            if self.model.apply_baseflow:
                print('Applying baseflow... getting BFVF10')
                self.bfvf10 = storm.get_bfvf10()
                print(f'Found a baseflow volume factor for 10% AEP of {self.bfvf10}')
            else:
                print('Baseflow not applied')

        return storm

    def initialise_model(self, config_file):
        # Get the config data
        print(f'\nSetting up hydrological model from config file: {config_file}')
        f = open(config_file)
        config_data = json.load(f)
        f.close()

        # Set up the URBS model
        model_type = config_data['model_type'].upper()
        max_keys = config_data['max_keys']
        if model_type == 'URBS':
            run_id = os.path.split(self.outputfile)[1]
            model = UrbsModel(config_file, dam_location=max_keys['level'], sub_folder = run_id)
        elif model_type == 'RORB':
            run_id = os.path.split(self.outputfile)[1]
            model = RorbModel(config_file, dam_location=config_data['ADV_name'], sub_folder = run_id)
        else:
            raise Exception(f'Model type must be "URBS" or "RORB". ({model_type} not recognised!)')
        return {'model': model, 'model_type': model_type, 'max_keys': max_keys}

    def get_initial_loss_scaling(self, percentile):
        # This data is taken from Table 5.3.13 in Book 5 Section 3.6 of ARR
        p_ax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        scale_ax = [3.19, 2.26, 1.71, 1.4, 1.2, 1, 0.85, 0.68, 0.53, 0.39, 0.14]
        # z_ax = ndtri(p_ax)
        curve = interpolate.interp1d(p_ax, scale_ax, kind='cubic')
        # z_min = z_ax[0]
        # z_max = z_ax[-1]
        # for i, element in enumerate(z):
        #     if element < z_min:
        #         z[i] = z_min
        #     elif element > z_max:
        #         z[i] = z_max
        # print(f'Interpolating initial loss for percentile {np.around(percentile, 1)}%')
        il_scaling = curve(percentile)
        return il_scaling

    def get_continuing_loss_scaling(self, percentile):
        # This data is taken from Table 5.3.13 in Book 5 Section 3.6 of ARR
        p_ax = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        scale_ax = [3.85, 2.48, 1.88, 1.50, 1.24, 1.00, 0.79, 0.61, 0.48, 0.35, 0.15]
        # z_ax = ndtri(p_ax)
        curve = interpolate.interp1d(p_ax, scale_ax, kind='cubic')
        # z_min = z_ax[0]
        # z_max = z_ax[-1]
        # for i, element in enumerate(z):
        #     if element < z_min:
        #         z[i] = z_min
        #     elif element > z_max:
        #         z[i] = z_max
        # print(f'Interpolating initial loss for percentile {np.around(percentile, 1)}%')
        cl_scaling = curve(percentile)
        return cl_scaling


class EnsembleSimulator(Simulator):
    def __init__(self, parameters, filepaths, test_iterations):
        super().__init__(parameters, filepaths, test_iterations)
        print('\nRunning in ensemble mode!')

        # output_folder = self.config_data['output_folder']
        # scheme_config = self.config_data['scheme_config']
        self.storm_durations = self.config_data['storm_durations']
        self.sample_both_methods = False
        self.enb = Ensemble(aep_list=self.config_data['aep_list'],
                            durations=self.storm_durations)
                            # output_folder=output_folder)

        if self.do_runs:
            print('\nRunning the models...')
            self.run_models()

        if self.do_analysis:
            print('\nAnalysing the model results...')
            self.analyse_results()
            
        if self.do_volumes:
            # for result_type in self.do_volumes:
            self.analyse_volumes(result_types=self.do_volumes)

    def run_models(self):
        # Initialise the simulation
        model = self.model
        storm_durations = self.config_data['storm_durations']
        # filepaths = self.config_data['file_paths']
        # max_keys = self.config_data['max_keys']
        # sim_filename_template = self.config_data['sim_filename_template']
        # sim_filename_template = sim_filename_template.replace('.csv', f'{self.suffix}.csv')  # add the suffix
        storm = self.initialise_storm(self.filepaths['storm_config'], storm_durations)
        storm.storm_method_config.update(self.config_data['storm_method_config'])
        if storm.storm_method_config['interim_for_ensemble'] == 'both':
            self.enb.update_scheme_in_interim(storm.storm_method_config['aep_changeover_to_extreme'])
            self.sample_both_methods = True

        # Set up the antecedent dam volume (ADV)
        # self.lake.set_full_supply_volume(model.full_supply_volume)
        self.enb.df['ADV'] = self.lake.antecedent_volume
        # model.set_volume_below_fsl(self.lake.get_volume_below_fsl())
        
        # Climate change adjustment on losses
        if not self.exclusions['loss_uplift']:
            il_climate_adjustment = self.climate.get_loss_uplift_factor(loss_type='il')
            cl_climate_adjustment = self.climate.get_loss_uplift_factor(loss_type='cl')
        else:
            il_climate_adjustment = 1.0
            cl_climate_adjustment = 1.0
        
        # do_filtering = storm.do_embedded_burst_filtering
        # if do_filtering:
        #     if self.exclusions['embedded_burst_filter']:
        #         do_filtering = False
        
        if self.exclusions['embedded_burst_filter']:
            do_filtering = False
            # storm.do_embedded_burst_filtering = False
        else:
            do_filtering = True
            # storm.do_embedded_burst_filtering = True
        
        # -------------------------------
        # Loop through simulations in the ensemble method
        for sim_id in self.enb.df.index:
            # Get the storm characteristics
            rain_sample_z = self.enb.df.loc[sim_id, 'rain_z']
            rain_sample_aep = self.enb.df.loc[sim_id, 'rain_aep']
            point_tp_frequency_bin = storm.point_tp_frequency_bins(rain_sample_aep)
            self.enb.df.loc[sim_id, 'tp_frequency'] = point_tp_frequency_bin
            duration = self.enb.df.loc[sim_id, 'duration']
            if not self.exclusions['rainfall_uplift']:
                rainfall_climate_adjustment = self.climate.get_rainfall_uplift_factor(duration=duration)
                rain_climate_adj = {}
                for dur in storm.rainfall.durations:
                    rain_climate_adj[dur] = self.climate.get_rainfall_uplift_factor(duration=dur)
            else:
                rainfall_climate_adjustment = 1.0
                rain_climate_adj = None

            # Sample the storm method to use: ARR/GSDM/GTSMR
            if self.enb.df.loc[sim_id, 'storm_method'] == 'extreme':
                df_method = 'extreme'
            else:
                df_method = None
            # storm_method = storm.sample_storm_method(rain_sample_z=rain_sample_z,
            storm_method = storm.sample_storm_method(rain_sample_aep=rain_sample_aep,
                                                     duration=duration,
                                                     sim_method='ensemble',
                                                     df_method=df_method)
            self.enb.df.loc[sim_id, 'storm_method'] = storm_method

            # Get the rainfall - don't apply climate change uplift at this stage
            # because we first need to do the temporal pattern filtering
            rain_depths = storm.rainfall.get_depth_z(z=rain_sample_z,
                                                     duration=duration,
                                                     storm_method=storm_method)
            # ave_rain = rain_depths.mean()
            ave_rain = storm.get_average_rain(rain_depths)

            # sample the temporal pattern for events
            tp_sample = self.enb.df.loc[sim_id, 'tp']
            temporal_pattern = storm.get_temporal_pattern(storm_method=storm_method,
                                                          duration=duration,
                                                          tp_sample=tp_sample,
                                                          rain_sample_z=rain_sample_z)
            timesteps = storm.timesteps
            # print('\nTemporal pattern is:')
            # print(temporal_pattern)

            # now apply climate change adjustment to the storm
            rain_depths = rain_depths * rainfall_climate_adjustment
            ave_rain = ave_rain * rainfall_climate_adjustment
            self.enb.df.loc[sim_id, 'mean_rain_mm'] = ave_rain

            # Get the rainfall losses
            # initial_loss = initial_loss * il_climate_adjustment
            storm_initial_loss = np.around(storm.storm_initial_loss * il_climate_adjustment, 1)
            if storm.apply_cl_cap:
                continuing_loss = storm.do_the_cl_capping(rain_sample_z)
            else:
                continuing_loss = storm.continuing_loss
            continuing_loss = np.around(continuing_loss * cl_climate_adjustment, 1)
            self.enb.df.loc[sim_id, 'continuing_loss'] = continuing_loss

            # filter embedded bursts
            if do_filtering:
                temporal_pattern, embedded_burst_comment = storm.filter_temppat(temporal_pattern, rain_sample_z,
                                                                                duration, ave_rain, storm_method,
                                                                                buffer = 1.1,  # hard coded buffer to reduce below burst depth
                                                                                climate_adjustment = rain_climate_adj)
                self.enb.df.loc[sim_id, 'embedded_bursts'] = embedded_burst_comment
            else:
                embedded_burst_comment = storm.check_embedded_burst(temporal_pattern, rain_sample_z,
                                                                    duration, ave_rain, storm_method,
                                                                    climate_adjustment = rain_climate_adj)
                self.enb.df.loc[sim_id, 'embedded_bursts'] = embedded_burst_comment

            # get the preburst rainfall
            preburst_p = 0.5  # using median preburst for ensembe event approach
            if 'preburst' in self.given_percentiles.keys():
                preburst_p = self.given_percentiles['preburst']
            preburst_proportion = storm.get_preburst_proportion(rain_sample_z, preburst_p, duration)
            self.enb.df.loc[sim_id, 'preburst_proportion'] = preburst_proportion
            preburst_depth = preburst_proportion * ave_rain  # preburst depth is uplifted for climate change
            self.enb.df.loc[sim_id, 'preburst_mm'] = preburst_depth

            # get the initial loss
            il_scaling = 1.0  # just using best estimate for ensemble approach
            if 'il' in self.given_percentiles.keys():
                il_scaling = self.get_initial_loss_scaling(self.given_percentiles['il'])
            initial_loss = np.around(il_scaling * storm_initial_loss, 1)
            self.enb.df.loc[sim_id, 'initial_loss'] = initial_loss    
            if self.exclusions['preburst']:
                initial_loss = np.around(il_scaling * storm_initial_loss - preburst_depth, 1)
                if initial_loss < 0:
                    residual_depth = -1 * initial_loss
                    initial_loss = 0.0
                else:
                    residual_depth = 0.0
                self.enb.df.loc[sim_id, 'residual_depth'] = residual_depth
            else:
                if self.preburst_method == 'uniform':
                    preburst_pattern = storm.uniform_preburst(duration, preburst_proportion, timesteps)
                    preburst_duration = preburst_pattern.index.min() * -1.0
                    temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)
                    
                else:
                    preburst_pattern, _ = storm.preburst_patterns.get_preburst_pattern(duration,
                                                                                       preburst_proportion, 
                                                                                       timesteps,  
                                                                                       sample_int='median')
    
                    preburst_pattern, filter_comment = storm.filter_embedded_bursts_in_preburst(preburst_pattern, temporal_pattern,
                                                                                                duration, rain_sample_z, ave_rain,
                                                                                                storm_method, buffer=0.9,
                                                                                                climate_adjustment = rain_climate_adj)
                    if filter_comment.split(':')[0].upper() == 'ERROR':
                        print('Pre-burst filter error with the following inputs:')
                        print(self.enb.df.loc[sim_id])
                        raise Exception()
    
                    preburst_duration = preburst_pattern.index.min() * -1.0
                    temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)
                
                    storm_duration = preburst_duration + duration
                    ifd_max_duration = max(storm.rainfall.durations)
                    if storm_duration > ifd_max_duration:
                        print(f'WARNING: storm duration {storm_duration}h exceeds max duration in IFD curves ({ifd_max_duration}h). IFD curve extrapolated for pre-burst filtering')
                

            # get lake level - not used anymore
            # lake_level = self.lake_level
            # self.enb.df.loc[sim_id, 'lake_level'] = lake_level
            # model.adjust_initial_lake_level(lake_level)

            # create storm file and run the URBS model
            simulation_period = model.get_simulation_period(duration)
            if self.model_type == 'URBS':
                # self.enb.df.loc[sim_id, 'ADV'] = self.lake_volume
                storm_duration = duration
                if self.model.apply_baseflow:
                    self.model.insert_baseflow_into_vec(self.bfvf10, rain_sample_z)
                if not self.exclusions['preburst']:
                    rain_depths = rain_depths * (1+preburst_proportion)
                    storm_duration = duration + preburst_duration
                    simulation_period += preburst_duration
                duration_str = model.duration_string(duration)
                storm_filename = f'sim_{str(sim_id).zfill(5)}.{duration_str}h'
                model.create_storm_file(rainfall_depths=rain_depths,
                                        temporal_pattern=temporal_pattern,
                                        data_interval=timesteps,
                                        filename=storm_filename,
                                        run_duration=simulation_period,
                                        storm_duration=storm_duration,
                                        storm_duration_excl_pb=duration,
                                        initial_loss=initial_loss,
                                        continuing_loss=continuing_loss,
                                        ari=np.round(rain_sample_aep, 0),
                                        ensemble=tp_sample)

                result_filename = f'sim_{str(sim_id).zfill(4)}_{duration_str}h'
                model.run_storm(storm_name=storm_filename,
                                result_name=result_filename)

            elif self.model_type == 'RORB':
                if not self.exclusions['preburst']:
                    simulation_period += preburst_duration
                duration_str = model.duration_string(duration)
                storm_filename = f'sim{duration_str}h_{str(sim_id).zfill(5)}.stm'
                model.create_storm_file(rainfall_depths=rain_depths,
                                        temporal_pattern=temporal_pattern,
                                        data_interval=timesteps,
                                        filename=storm_filename,
                                        run_duration=simulation_period)  # ,
                                        # storm_duration=duration,
                                        # initial_loss=initial_loss,
                                        # continuing_loss=continuing_loss)

                result_filename = model.run_storm(storm_name=storm_filename,
                                                  initial_loss=initial_loss,
                                                  continuing_loss=continuing_loss)

            # Store max flows
            max_data = model.get_max_results(result_filename, self.max_keys)
            self.enb.df.loc[sim_id, 'inflow'] = max_data['inflow']
            self.enb.df.loc[sim_id, 'level'] = max_data['level']
            self.enb.df.loc[sim_id, 'outflow'] = max_data['outflow']

            # Store the hydrographs
            if self.store_hydrographs:
                model.get_hydrographs(result_filename, self.max_keys, sim_id)

            # Mop up files and move outputs to results folder
            if self.mop_files:
                model.mop_storm_file(storm_filename, skip=10, sim_id=sim_id)
                model.mop_results_files(result_filename)
                # model.move_results(result_filename)

            # Increment the simulation id
            sim_id += 1

        # # Set the initial lake level back to what it was in the original vec file
        # model.restore_initial_lake_level()

        # Store simulation output
        #sim_filename = sim_filename_template.replace('xx', str(duration))
        self.enb.store_simulations(self.outputfile)
        print(self.enb.df.head())

        # Store the hydrographs to file
        if self.store_hydrographs:
            # result_filename = sim_filename.replace('.csv', '')
            model.store_hydrographs(self.outputfile)

        self.print_elapsed_time()

    def analyse_results(self):
        # output_folder = self.config_data['output_folder']
        result_types = ['inflow', 'level', 'outflow']
        result_labels = {'inflow': 'Flow (m³/s)',
                         'level': 'Elevation (m AHD)',
                         'outflow': 'Flow (m³/s)'}
        result_titles = {'inflow': 'Dam inflow',
                         'level': 'Lake level',
                         'outflow': 'Dam outflow'}
        sim_filename = f'{self.outputfile}.csv'
        # sim_path = os.path.join(output_folder, sim_filename)
        print('\nReading the Monte Carlo analysis file:', sim_filename)
        df = pd.read_csv(sim_filename, index_col=0)
        aep_lst = df['rain_aep'].unique()
        duration_lst = df['duration'].unique()
        positions = range(len(duration_lst))
        pattern_ids = df['tp'].unique()
        print(df)

        # Set up the critical results dataframe
        headers = []
        for result_type in result_types:
            headers.append(result_type)
            headers.append(f'{result_type}_duration')
            headers.append(f'{result_type}_tp')
        criticals = pd.DataFrame(columns=headers, index=aep_lst)
        
        # Create sub-folders for outputs
        output_path = os.path.join(os.path.dirname(self.outputfile), 'plots')
        os.makedirs(output_path, exist_ok=True)
        output_path = os.path.join(os.path.dirname(self.outputfile), 'csv')
        os.makedirs(output_path, exist_ok=True)

        # Analyse each AEP
        for aep in aep_lst:
            print(f'Working on 1 in {aep} AEP events')
            aep_df = df[df['rain_aep'] == aep]
            # storm_types = aep_df['storm_method'].unique()
            # number_of_storm_types = storm_types.shape[0]
            # pattern_lst_20 = None
            # if number_of_storm_types == 2:
            #     print('Found 20 temporal patterns')
            #     part_1 = pd.DataFrame(columns=['tp'], index=num_patterns, data=num_patterns)
            #     part_2 = pd.DataFrame(columns=['tp'], index=num_patterns, data=num_patterns)
            #     part_1['composite'] = storm_types[0] + ": " + part_1['tp'].astype(str)
            #     part_2['composite'] = storm_types[1] + ": " + part_2['tp'].astype(str)
            #     pattern_lst_20 = pd.concat([part_1, part_2])
            #     pattern_lst_20 = pattern_lst_20['composite'].tolist()

            # Set up the box plot
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(9, 4))
            fig.suptitle(f'1 in {aep} AEP')

            # Set up the median dataframe
            medians = pd.DataFrame(columns=headers, index=duration_lst)

            # Collate the data for each type: inflows, outflows, and levels
            for ind, result_type in enumerate(result_types):
                # if number_of_storm_types == 2:
                #     result = pd.DataFrame(columns=duration_lst, index=pattern_lst_20)
                # else:
                #     result = pd.DataFrame(columns=duration_lst, index=num_patterns)
                # extract the data for each duration into the result dataframe
                storm_types = aep_df['storm_method'].unique()
                number_of_storm_types = storm_types.shape[0]
                print(f'Found {number_of_storm_types} storm types:', storm_types)
                result_index = []
                for storm_type in storm_types:
                    for pattern_id in pattern_ids:
                        index_label = f'{storm_type}: {pattern_id}'
                        result_index.append(index_label)
                result = pd.DataFrame(index=result_index)
                for duration_ind, duration in enumerate(duration_lst):
                    print(f'Working on {duration} hour storm duration')
                    dur_df = aep_df[aep_df['duration'] == duration]
                    # median_id = 5
                    # if dur_df.shape[0] == 10:
                    #     print('Found ten temporal patterns')
                    #     dur_df.set_index('tp', inplace=True)
                    # elif dur_df.shape[0] == 20:
                        # pattern_lst = dur_df['composite_index'].unique()
                    dur_df['composite_index'] = dur_df['storm_method'] + ': ' + dur_df['tp'].astype(str)
                    dur_df.set_index('composite_index', inplace=True)
                    #     median_id = 10
                    print(dur_df)
                    result[duration] = dur_df[result_type]
                    print(result)

                    local = dur_df[[result_type]].sort_values(by=result_type, ascending=True)
                    median_id = int(np.around(local.shape[0] / 2, 0))
                    print('\nGetting median from position:', median_id)
                    median = local[result_type].iloc[median_id]
                    median_tp = local.index[median_id]
                    medians.loc[duration, result_type] = median
                    medians.loc[duration, f'{result_type}_tp'] = median_tp
                    medians.loc[duration, f'{result_type}_duration'] = duration
                    print(f'\n{result_type} results for 1 in {aep} AEP with {duration} hour duration:')
                    print('Median value is {} of {} for pattern {}'.format(result_labels[result_type],
                                                                           median,
                                                                           median_tp))
                    # Plot the local data points
                    print(f'\nDuration index: {duration_ind}')
                    x = np.full(result.shape[0], duration_ind)
                    print(x)
                    ax[ind].plot(x, result[duration], '.', mfc='b', markeredgewidth=0.0)

                # Plot the box data
                print(result)
                plot_data = []
                for res in result.columns:
                    local_result = result[res].dropna()
                    plot_data.append(local_result)
                # ax[ind].boxplot(result.to_numpy(), positions=positions, usermedians=medians[result_type])
                ax[ind].boxplot(plot_data, positions=positions, usermedians=medians[result_type])
                # plt.show()
                # Format the plot
                ax[ind].set_xticks(positions)
                ax[ind].set_xticklabels(duration_lst)
                # plt.ylim((-10, 10))
                # plt.xticks(rotation=90)
                ax[ind].set_xlabel("Storm duration (hours)")
                ax[ind].set_ylabel(result_labels[result_type])
                ax[ind].set_title(result_titles[result_type])
                ax[ind].grid()

            # Output the plot
            # output_file = self.outputfile.replace('.csv', f'_{aep}_bp.png')
            output_file = f'{self.outputfile}_{aep}_bp.png'
            output_path = os.path.join(os.path.dirname(output_file),
                                       'plots',
                                       os.path.basename(output_file))
            # output_name = sim_path.replace('.csv', '')
            # output_name = f'{output_name}_{result_type}'
            print('\nWriting figure to:', output_path)
            plt.tight_layout()
            plt.savefig(output_path)
            # plt.show()

            # Output the medians
            # outut_file = self.outputfile.replace('.csv', f'_{aep}_med.csv')
            output_file = f'{self.outputfile}_{aep}_med.csv'
            # output_path = os.path.join(output_folder, 'csv', output_file)
            output_path = os.path.join(os.path.dirname(output_file),
                                       'csv',
                                       os.path.basename(output_file))
            medians.index.name = 'aep (1 in x)'
            medians.to_csv(output_path)

            # Get the critical events
            for result_type in result_types:
                local_df = medians[result_type].astype(float)
                crit_index = local_df.squeeze().argmax()
                criticals.loc[aep, result_type] = medians[result_type].iloc[crit_index]
                criticals.loc[aep, f'{result_type}_tp'] = medians[f'{result_type}_tp'].iloc[crit_index]
                criticals.loc[aep, f'{result_type}_duration'] = medians[f'{result_type}_duration'].iloc[crit_index]
                # print(crit_durations)
                # Need to output the storm duration and then the final results!!!

        print('\nFinal results of critical events from median patterns:')
        print(criticals)
        # output_file = self.outputfile.replace('.csv', '_critical.csv')
        output_file = f'{self.outputfile}_critical.csv'
        # output_path = os.path.join(output_folder, 'csv', output_file)
        output_path = os.path.join(os.path.dirname(output_file),
                                   'csv',
                                   os.path.basename(output_file))
        criticals.to_csv(output_path)
        
    def analyse_volumes(self, result_types=['inflow', 'outflow']):
        enb = self.enb
        sim_filename = f'{self.outputfile}.csv'
        print('\nOpening Ensemble analysis file:', sim_filename)
        enb.df = pd.read_csv(sim_filename, index_col=0)
        aep_lst = enb.df['rain_aep'].unique()
        duration_lst = enb.df['duration'].unique()
        positions = range(len(duration_lst))
        pattern_lst = enb.df['tp'].unique()
        print(enb.df)
        
        if self.model_type == 'URBS':
            results_folder = os.path.abspath(os.path.join(self.outputfile, '../../urbs_results'))
        elif self.model_type == 'RORB':
            results_folder = os.path.abspath(os.path.join(self.outputfile, '../../rorb_results'))
            
        for result_type in result_types:
            print('\nAnalysing the ', result_type, ' volumes')
            hyd_filename = os.path.basename(self.outputfile) + f'_{result_type}s.csv'
            filepath = os.path.join(results_folder, hyd_filename)
            print('\nOpening hydrograph file:', filepath)
            try:
                hyd_df = pd.read_csv(filepath, index_col=0)
            except:
                print('Hydrograph file not found.')
                return
            
            time_inc = hyd_df.index[2] - hyd_df.index[1]
            print('Integrating volume using: ', time_inc, 'h time steps.')
            hyd_df.fillna(0, inplace=True)   # Replace NaNs with zero. For hydrographs shorter than required duration. Recommend extending simulation periods if practical.
            
            durations = duration_lst
            data = []
            for duration in durations:
                window = duration / time_inc
                if window.is_integer(): 
                    window = int(window)
                else: 
                    print(f'Unable to integrate {duration}h duration with {time_inc}h time steps')
                    continue
                
                print(f'Calculating max {duration}h volumes')
                vol = hyd_df.rolling(window).apply(integrate.trapz).max() 
                vol = vol / 1000 * 60 * 60 * time_inc           # Convert from m3/s to ML
                data.append(vol)
            
            vol_df = pd.concat(data, axis=1)
            print(vol_df)
            vol_types = [f'Vol{dur}h' for dur in durations]
            vol_df.columns = vol_types
            vol_df.index = vol_df.index.to_series().map(lambda x: int(x.split('_')[1]))
        
            vol_df = pd.concat([enb.df[['rain_z', 'rain_aep', 'mean_rain_mm', 'duration', 'tp', 'storm_method']], vol_df], axis=1)
            volume_file = self.outputfile + '_volume.csv'
            print('Writing volume results:', volume_file)
            vol_df.to_csv(volume_file)
            enb.df = vol_df
                              
            # Set up the critical results dataframe
            headers = []
            for vol_type in vol_types:
                headers.append(vol_type)
                headers.append(f'{vol_type}_duration')
                headers.append(f'{vol_type}_tp')
            criticals = pd.DataFrame(columns=headers, index=aep_lst)
            
            # Analyse each AEP
            for aep in aep_lst:
                print(f'Working on 1 in {aep} AEP events')
                aep_df = vol_df[vol_df['rain_aep'] == aep]
                storm_types = aep_df['storm_method'].unique()
                number_of_storm_types = storm_types.shape[0]
                pattern_lst_20 = None
                if number_of_storm_types == 2:
                    print('Found 20 temporal patterns')
                    part_1 = pd.DataFrame(columns=['tp'], index=pattern_lst, data=pattern_lst)
                    part_2 = pd.DataFrame(columns=['tp'], index=pattern_lst, data=pattern_lst)
                    part_1['composite'] = storm_types[0] + ": " + part_1['tp'].astype(str)
                    part_2['composite'] = storm_types[1] + ": " + part_2['tp'].astype(str)
                    pattern_lst_20 = pd.concat([part_1, part_2])
                    pattern_lst_20 = pattern_lst_20['composite'].tolist()
                    
                # Set up the median dataframe
                medians = pd.DataFrame(columns=headers, index=duration_lst)
                    
                # Collate the data for each type: inflows, outflows, and levels
                result_labels = 'ML'
                for ind, vol_type in enumerate(vol_types):
                    if number_of_storm_types == 2:
                        result = pd.DataFrame(columns=duration_lst, index=pattern_lst_20)
                    else:
                        result = pd.DataFrame(columns=duration_lst, index=pattern_lst)
                    # extract the data for each duration into the result dataframe
                    for duration_ind, duration in enumerate(duration_lst):
                        dur_df = aep_df[aep_df['duration'] == duration]
                        median_id = 5
                        if dur_df.shape[0] == 10:
                            print('Found ten temporal patterns')
                            dur_df.set_index('tp', inplace=True)
                        elif dur_df.shape[0] == 20:
                            # pattern_lst = dur_df['composite_index'].unique()
                            dur_df['composite_index'] = dur_df['storm_method'] + ': ' + dur_df['tp'].astype(str)
                            dur_df.set_index('composite_index', inplace=True)
                            median_id = 10
                        print(dur_df)
                        result[duration] = dur_df[vol_type]
                        print(result)

                        local = dur_df[[vol_type]].sort_values(by=vol_type, ascending=True)
                        print(local)
                        median = local[vol_type].iloc[median_id]
                        print(median)
                        
                        median_tp = local.index[median_id]
                        medians.loc[duration, vol_type] = median
                        medians.loc[duration, f'{vol_type}_tp'] = median_tp
                        medians.loc[duration, f'{vol_type}_duration'] = duration
                        print(f'\n{vol_type} results for 1 in {aep} AEP with {duration} hour duration:')
                        print('Median value for volume is {} of {} for pattern {}'.format(result_labels,
                                                                               median,
                                                                               median_tp))
                                                                               
                # Get the critical events
                for vol_type in vol_types:
                    local_df = medians[vol_type].astype(float)
                    crit_index = local_df.squeeze().argmax()
                    criticals.loc[aep, vol_type] = medians[vol_type].iloc[crit_index]
                    criticals.loc[aep, f'{vol_type}_tp'] = medians[f'{vol_type}_tp'].iloc[crit_index]
                    criticals.loc[aep, f'{vol_type}_duration'] = medians[f'{vol_type}_duration'].iloc[crit_index]
                    # print(crit_durations)
                    # Need to output the storm duration and then the final results!!!

            print('\nFinal results of critical events from volume median patterns:')
            print(criticals)
            # output_file = self.outputfile.replace('.csv', '_critical.csv')
            output_file = f'{self.outputfile}_{result_type}Vol_critical.csv'
            # output_path = os.path.join(output_folder, 'csv', output_file)
            output_path = os.path.join(os.path.dirname(output_file),
                                       'csv',
                                       os.path.basename(output_file))
            criticals.to_csv(output_path)

class MonteCarloSimulator(Simulator):
    def __init__(self, parameters, filepaths, test_iterations):
        super().__init__(parameters, filepaths, test_iterations)
        print('\nRunning in Monte Carlo mode!')
        # output_folder = self.config_data['output_folder']
        scheme_config = self.config_data['scheme_config']
        self.mc = SampleScheme(lower_aep=scheme_config['lower_aep'],
                               upper_aep=scheme_config['upper_aep'],
                               number_of_main_divisions=scheme_config['number_of_main_divisions'],
                               number_of_sub_divisions=scheme_config['number_of_sub_divisions'],
                               number_of_temporal_patterns=scheme_config['number_of_temporal_patterns'],
                               sample_method=scheme_config['sample_method'])
                               # output_folder=output_folder)
        
        if self.do_runs:
            self.mc.setup_rain_sample_space()       # Set up the sample space for running and/or analysing the results
            print('\nRunning the models...')
            self.run_models()

        if self.do_analysis == True:
            print('\nAnalysing the model results...')
            self.analyse_results()
        elif type(self.do_analysis) is list:
            print('\nAnalysing the model results...')
            self.analyse_results(result_types=self.do_analysis)
            
        if self.do_volumes:
            # for result_type in self.do_volumes:
            self.analyse_volumes(result_types=self.do_volumes)

    def run_models(self):
        mc = self.mc
        model = self.model
        # duration = int(self.duration)
        duration = float(self.duration)
        if duration.is_integer():
            duration = int(duration)

        # Initialise the simulation
        # storm_durations = self.config_data['storm_durations']
        # filepaths = self.config_data['file_paths']
        # max_keys = self.config_data['max_keys']
        # sim_filename_template = self.config_data['sim_filename_template']
        # sim_filename_template = sim_filename_template.replace('.csv', f'{self.suffix}.csv')  # add the suffix
        storm = self.initialise_storm(self.filepaths['storm_config'], [duration])

        # sample the initial losses
        # il_sample = mc.random_std_norm_variates()
        if self.replicates['il_p']:
            il_sample = self.replicates_df['il_p']
        else:
            il_sample = mc.sample_percentiles()
        mc.df['il_p'] = il_sample
        il_scaling = self.get_initial_loss_scaling(il_sample)
        mc.df['il_scaling'] = il_scaling

        # sample the continuing losses
        # Check if replicating prior simulation
        if self.replicates['cl_p']:
            # replicated sample
            cl_sample = self.replicates_df['cl_p']
        else:
            # new sample
            cl_sample = mc.sample_percentiles()
        # Check if excluding the stochastic sampling of continuing losses
        if self.exclusions['cl_sampling']:
            mc.df['cl_p'] = 0.5
            mc.df['cl_scaling'] = 1.0
        else:
            mc.df['cl_p'] = cl_sample
            cl_scaling = self.get_continuing_loss_scaling(cl_sample)
            mc.df['cl_scaling'] = cl_scaling

        # sample the pre-burst likelihood
        # preburst_sample = mc.random_std_norm_variates()
        if self.replicates['preburst_p']:
            preburst_sample = self.replicates_df['preburst_p']
        else:
            preburst_sample = mc.sample_percentiles(lower=0.1, upper=0.9)
        mc.df['preburst_p'] = preburst_sample
        
        if self.replicates['preburst_tp']:
            mc.df['preburst_tp'] = self.replicates_df['preburst_tp']
        else:
            mc.df['preburst_tp'] = np.nan   # Sampling done during run
            mc.df['preburst_filter'] = None
            

        # sample the antecedent lake volume - correlation done later if needed
        # lake = LakeConditions(filepaths['lake_config'])
        if self.replicates['lake_z']:
            mc.df['lake_z'] = self.replicates_df['lake_z']
        else:
            mc.df['lake_z'] = mc.random_std_norm_variates()

        # -------------------------------
        # Loop through Storm durations - no longer looping through durations: simulating one duration at a time
        # for duration in storm_durations:
        # Apply climate change adjustments
        if not self.exclusions['rainfall_uplift']:
            rainfall_climate_adjustment = self.climate.get_rainfall_uplift_factor(duration=duration)
            rain_climate_adj = {}
            for dur in storm.rainfall.durations:
                rain_climate_adj[dur] = self.climate.get_rainfall_uplift_factor(duration=dur)
            # print(rain_climate_adj)
        else:
            rainfall_climate_adjustment = 1.0
            rain_climate_adj = None

        if not self.exclusions['loss_uplift']:
            il_climate_adjustment = self.climate.get_loss_uplift_factor(loss_type='il')
            cl_climate_adjustment = self.climate.get_loss_uplift_factor(loss_type='cl')
        else:
            il_climate_adjustment = 1.0
            cl_climate_adjustment = 1.0
        if not self.exclusions['pattern_d50_scaling']:
            delta_d50 = self.climate.get_delta_d50(duration=duration)
        else:
            delta_d50 = 0
        delta_d50_info = self.climate.get_d50_shift_weightings(duration=duration,
                                                               delta_d50=delta_d50,
                                                               pattern_area=storm.gtsmr_pattern_area)
        print(f'\nDelta D50 is {delta_d50} and the associated information is:')
        print(delta_d50_info)
        
        # do_filtering = storm.do_embedded_burst_filtering
        # if do_filtering:
        #     if self.exclusions['embedded_burst_filter']:
        #         do_filtering = False
        
        if self.exclusions['embedded_burst_filter']:
            do_filtering = False
            # storm.do_embedded_burst_filtering = False
        else:
            do_filtering = True
            # storm.do_embedded_burst_filtering = True

        # Set up  the simulations
        sim_id = 0
        for m in range(mc.m):
            # Exit if testing
            if 0 < self.test < sim_id:
                print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                print('Running in test mode - exiting early')
                print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                break

            for n in range(mc.n):
                mc.df.loc[sim_id, 'm'] = m
                mc.df.loc[sim_id, 'n'] = n

                # Sample the rainfall
                if self.replicates['rain_z']:
                    rain_sample_z = self.replicates_df.loc[sim_id, 'rain_z']
                else:
                    rain_sample_z = mc.get_rain_sample(m, n)
                mc.df.loc[sim_id, 'rain_z'] = rain_sample_z
                rain_sample_aep = 1 / (1 - ndtr(rain_sample_z))  # 1 in X AEP
                mc.df.loc[sim_id, 'rain_aep'] = rain_sample_aep
                point_tp_frequency_bin = storm.point_tp_frequency_bins(rain_sample_aep)
                mc.df.loc[sim_id, 'tp_frequency'] = point_tp_frequency_bin

                # Sample the storm method to use: ARR/GSDM/GTSMR
                if self.replicates['storm_method']:
                    storm_method = self.replicates_df.loc[sim_id, 'storm_method']
                else:
                    storm_method = storm.sample_storm_method(rain_sample_z=rain_sample_z,
                                                             duration=duration)
                mc.df.loc[sim_id, 'storm_method'] = storm_method
                print(f'\nsim_{str(sim_id).zfill(5)} | 1 in {np.around(rain_sample_aep, 1)} AEP | {storm_method}')
                # Get the weighting of the temporal pattern D50 shift for climate change
                pattern_d50_info = delta_d50_info[storm_method]
                delta_d50_weighting = pattern_d50_info['weightings']
                delta_d50_patterns = pattern_d50_info['front_patterns']
                if storm_method == 'ARR point':
                    delta_d50_weighting = delta_d50_weighting[point_tp_frequency_bin]

                # Get the rainfall
                rain_depths = storm.rainfall.get_depth_z(z=rain_sample_z,
                                                         duration=duration,
                                                         storm_method=storm_method)

                #ave_rain = rain_depths.mean()
                ave_rain = storm.get_average_rain(rain_depths)

                # sample the temporal pattern for events
                if self.replicates['tp']:
                    tp_sample = self.replicates_df.loc[sim_id, 'tp']
                else:
                    tp_sample = mc.get_temporal_pattern_sample(delta_d50_weighting, delta_d50_patterns)
                mc.df.loc[sim_id, 'tp'] = tp_sample
                temporal_pattern = storm.get_temporal_pattern(storm_method=storm_method,
                                                              duration=duration,
                                                              tp_sample=tp_sample,
                                                              rain_sample_z=rain_sample_z)
                timesteps = storm.timesteps
                # print('\nTemporal pattern is:')
                # print(temporal_pattern)

                # apply climate change adjustment to the storm - BEFORE preburst proportion applied!
                # Therefore, uplifting preburst as well
                rain_depths = rain_depths * rainfall_climate_adjustment
                ave_rain = ave_rain * rainfall_climate_adjustment
                mc.df.loc[sim_id, 'mean_rain_mm'] = ave_rain

                # Compute rainfall losses
                storm_initial_loss = il_climate_adjustment * storm.storm_initial_loss
                if storm.apply_cl_cap:
                    continuing_loss = storm.do_the_cl_capping(rain_sample_z)
                else:
                    continuing_loss = storm.continuing_loss
                continuing_loss = np.around(continuing_loss * cl_climate_adjustment, 1)
                # continuing_loss = storm.continuing_loss * cl_climate_adjustment

                # filter embedded bursts
                if do_filtering:
                    filter_buffer = 1.1  # hard coded buffer to reduce below burst depth
                    original_temporal_pattern = temporal_pattern.copy()
                    temporal_pattern, embedded_burst_comment = storm.filter_temppat(temporal_pattern, rain_sample_z,
                                                                                    duration, ave_rain, storm_method,
                                                                                    buffer=filter_buffer,
                                                                                    climate_adjustment=rain_climate_adj)

                    if embedded_burst_comment.strip() != 'No embedded bursts':
                        print(embedded_burst_comment)
                        filtered_temporal_pattern = temporal_pattern.copy()
                        storm.store_hyetographs(sim_id, original_temporal_pattern, filtered_temporal_pattern)
                    mc.df.loc[sim_id, 'embedded_bursts'] = embedded_burst_comment
                else:
                    embedded_burst_comment = storm.check_embedded_burst(temporal_pattern, rain_sample_z,
                                                                        duration, ave_rain, storm_method,
                                                                        climate_adjustment=rain_climate_adj)  # ,
                    # buffer = 1.0)
                    mc.df.loc[sim_id, 'embedded_bursts'] = embedded_burst_comment
                    
                #if embedded_burst_comment.upper().startswith('ERROR:'):
                #    self.terminate_model_runs(print_msg=embedded_burst_comment, run_data=mc.df.loc[sim_id])

                # get the preburst rainfall
                preburst_p = mc.df.loc[sim_id, 'preburst_p']
                preburst_proportion = storm.get_preburst_proportion(rain_sample_z, preburst_p, duration)
                mc.df.loc[sim_id, 'preburst_proportion'] = preburst_proportion
                preburst_depth = preburst_proportion * ave_rain
                mc.df.loc[sim_id, 'preburst_mm'] = preburst_depth

                # get the initial loss
                il_scaling = mc.df.loc[sim_id, 'il_scaling']
                initial_loss = np.around(il_scaling * storm_initial_loss, 1)
                mc.df.loc[sim_id, 'initial_loss'] = initial_loss

                # get the continuing loss
                cl_scaling = mc.df.loc[sim_id, 'cl_scaling']
                continuing_loss = np.around(cl_scaling * continuing_loss, 1)
                mc.df.loc[sim_id, 'continuing_loss'] = continuing_loss
                
                # Do pre-bursts
                if self.exclusions['preburst']:             # Convert to burst initial loss
                    initial_loss = np.around(il_scaling * storm_initial_loss - preburst_depth, 1)
                    if initial_loss < 0:
                        residual_depth = -1 * initial_loss
                        initial_loss = 0.0
                    else:
                        residual_depth = 0.0
                    mc.df.loc[sim_id, 'residual_depth'] = residual_depth
                else:
                    if self.preburst_method == 'uniform':
                        preburst_pattern = storm.uniform_preburst(duration, preburst_proportion, timesteps)
                        preburst_duration = preburst_pattern.index.min() * -1.0
                        temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)      # prepend preburst pattern
                        
                    else:
                        mc.df.loc[sim_id, 'residual_depth'] = np.around(preburst_depth - il_scaling * storm_initial_loss, 1)
                        if preburst_proportion == 0.0:
                            sample_int = -1
                            preburst_duration = 0.0
                            # temporal pattern unchanged
                        else:
                            if self.replicates['preburst_tp']: 
                                sample_int = int(np.around(self.replicates_df.loc[sim_id, 'preburst_tp'], 0))
                                preburst_pattern, _ = storm.preburst_patterns.get_preburst_pattern(duration, preburst_proportion, timesteps, sample_int)
                                # preburst_duration = preburst_pattern.index.min() * -1.0
                                # temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)      # prepend preburst pattern
                                
                            else:
                                preburst_pattern, sample_int = storm.preburst_patterns.get_preburst_pattern(duration, preburst_proportion, timesteps)
                                
                            preburst_pattern, filter_comment = storm.filter_embedded_bursts_in_preburst(preburst_pattern, temporal_pattern,
                                                                                                        duration, rain_sample_z, ave_rain,
                                                                                                        storm_method, buffer=0.9,
                                                                                                        climate_adjustment = rain_climate_adj)
                            if filter_comment.split(':')[0].upper() == 'ERROR':
                                error_msg = 'Pre-burst filter error with the following inputs:'
                                self.terminate_model_runs(print_msg=error_msg, run_data=self.mc.df.loc[sim_id])
                                
                            preburst_duration = preburst_pattern.index.min() * -1.0
                            temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)      # prepend preburst pattern
                            
                            storm_duration = preburst_duration + duration
                            ifd_max_duration = max(storm.rainfall.durations)
                            if storm_duration > ifd_max_duration:
                                print(f'WARNING: storm duration {storm_duration}h exceeds max duration in IFD curves ({ifd_max_duration}h). IFD curve extrapolated for pre-burst filtering')
                                                    
                            mc.df.loc[sim_id, 'preburst_filter'] = filter_comment
                                
                        mc.df.loc[sim_id, 'preburst_tp']  = sample_int

                # correlate the lake level z if needed
                lake_z = mc.df.loc[sim_id, 'lake_z']
                if self.lake.is_correlated:
                    lake_z = self.lake.get_correlated_z(rain_sample_z, lake_z)
                    mc.df.loc[sim_id, 'lake_z'] = lake_z

                # get the antecedent lake volume
                if self.lake.antecedent_type == 'varying':
                    # self.lake.set_full_supply_volume(model.full_supply_volume)      
                    lake_vol = self.lake.get_lake_volume(lake_z)
                    mc.df.loc[sim_id, 'ADV'] = lake_vol
                    self.lake.antecedent_volume = lake_vol
                    model.set_volume_below_fsl(self.lake.get_volume_below_fsl())
                
                # Catch error
                if sum(rain_depths.isna()):
                    print('\nSub-area rain depths')
                    print(rain_depths)
                    error_msg = '\n##################################' \
                    f'\nERROR generating rain depths. \nRain sample z: {rain_sample_z} \nRain sample aep: {rain_sample_aep}'
                    self.terminate_model_runs(print_msg=error_msg, run_data=mc.df.loc[sim_id])
                    return
                

                # create storm file and run the URBS model
                simulation_period = model.get_simulation_period(duration)
                if self.model_type == 'URBS':
                    storm_duration = duration
                    if self.model.apply_baseflow:
                        self.model.insert_baseflow_into_vec(self.bfvf10, rain_sample_z)
                    if not self.exclusions['preburst']:
                        rain_depths = rain_depths * (1+preburst_proportion)             # sub-area rain depths need to be scaled up for URBS but not RORB
                        storm_duration = duration + preburst_duration                   # storm duration needs to be adjusted for preburst duration for URBS storm files (but not RORB)
                        simulation_period += preburst_duration
                    # model.set_volume_below_fsl(volume_below_fsl)
                    duration_str = model.duration_string(duration)
                    storm_filename = f'sim_{str(sim_id).zfill(5)}.{duration_str}h'
                    # print(f'DARN: {storm_filename}')
                    model.create_storm_file(rainfall_depths=rain_depths,
                                            temporal_pattern=temporal_pattern,
                                            data_interval=timesteps,
                                            filename=storm_filename,
                                            run_duration=simulation_period,
                                            storm_duration=storm_duration,
                                            storm_duration_excl_pb=self.duration,
                                            initial_loss=initial_loss,
                                            continuing_loss=continuing_loss,
                                            ari=np.round(rain_sample_aep, 0),
                                            ensemble=tp_sample)

                    result_filename = f'sim_{str(sim_id).zfill(5)}_{duration_str}h'
                    model.run_storm(storm_name=storm_filename,
                                    result_name=result_filename)

                elif self.model_type == 'RORB':
                    if not self.exclusions['preburst']:
                        simulation_period += preburst_duration
                    storm_filename = f'sim{duration}h_{str(sim_id).zfill(5)}.stm'
                    # model.set_volume_below_fsl(volume_below_fsl)
                    model.create_storm_file(rainfall_depths=rain_depths,
                                            temporal_pattern=temporal_pattern,
                                            data_interval=timesteps,
                                            filename=storm_filename,
                                            run_duration=simulation_period) #,
                                            # storm_duration=duration,
                                            # initial_loss=initial_loss,
                                            # continuing_loss=continuing_loss)

                    result_filename = model.run_storm(storm_name=storm_filename,
                                                      initial_loss=initial_loss,
                                                      continuing_loss=continuing_loss)

                # Store max flows
                try:
                    max_data = model.get_max_results(result_filename, self.max_keys)
                except:
                    error_msg = '\n##################################' \
                    f'\nERROR reading URBS/RORB output file. \nUsually indicates error running URBS/RORB. \nCheck inputs and output files for sim {sim_id}'
                    self.terminate_model_runs(print_msg=error_msg, run_data=mc.df.loc[sim_id])
                    print('The URBS result file or peak flow could not be found... perhaps the URBS model did not run to completion.')
                    print('I suggest checking the URBS batch file and that all file references exist.')
                    
                mc.df.loc[sim_id, 'inflow'] = max_data['inflow']
                if 'level' in max_data:
                    mc.df.loc[sim_id, 'level'] = max_data['level']
                if 'outflow' in max_data:
                    mc.df.loc[sim_id, 'outflow'] = max_data['outflow']

                # Store the hydrographs
                if self.store_hydrographs:
                    model.get_hydrographs(result_filename, self.max_keys, sim_id)

                # Mop up files and move outputs to results folder
                if self.mop_files:
                    model.mop_storm_file(storm_filename, skip=10, sim_id=sim_id)
                    model.mop_results_files(result_filename)
                    # model.move_results(result_filename)

                # Increment the simulation id
                sim_id += 1

                # Exit if testing
                if 0 < self.test < sim_id:
                    print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    print('Running in test mode - exiting early')
                    print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    break

        # Set the initial lake level back to what it was in the original vec file
        # model.restore_initial_lake_level()

        # Store simulation output
        # sim_filename = sim_filename_template.replace('xx', str(duration))
        mc.store_simulations(self.outputfile)
        print(mc.df.head())

        # Store the hydrographs to file
        if self.store_hydrographs:
            # result_filename = sim_filename.replace('.csv', '')
            model.store_hydrographs(self.outputfile)

        storm.output_hyetographs(self.outputfile)

        self.print_elapsed_time()
        
    def terminate_model_runs(self, print_msg=None, run_data=None):
        if print_msg:
            print()
            print(print_msg)
        if run_data is not None:
            print()
            print(run_data)
        
        # if self.model_type == 'RORB':
        #     self.model.restore_initial_lake_level
        
        # store mcdf and hydrographs
        self.mc.store_simulations(self.outputfile)
        if self.store_hydrographs:
            self.model.store_hydrographs(self.outputfile)

        self.print_elapsed_time()
        return

    def analyse_results(self, result_types=['inflow', 'level', 'outflow']):
        mc = self.mc
        # storm_durations = self.config_data['storm_durations']
        # duration = int(self.duration)
        # sim_filename_template = self.config_data['sim_filename_template']
        # sim_filename_template = sim_filename_template.replace('.csv', f'{self.suffix}.csv')  # add the suffix
        # output_folder = self.config_data['output_folder']

        # result_types = ['inflow', 'level', 'outflow']
        # quantile_analysis_set = self.config_data['tpt_quantile_analysis']
        # for duration in storm_durations:
        # sim_filename = sim_filename_template.replace('xx', str(duration))
        # sim_filename = Path(self.outputfile).stem + '__mcdf.csv'
        sim_filename = self.outputfile + '__mcdf.csv'
        print('\nOpening Monte Carlo analysis file:', os.path.abspath(sim_filename))
        # sim_path = os.path.join(output_folder, sim_filename)
        mc.df = pd.read_csv(sim_filename, index_col=0)
        print(mc.df)

        # Get the AEP of the PMP
        storm = StormBurst(self.filepaths['storm_config'], generate_storms=False)
        mc.aep_of_pmp = storm.get_aep_of_pmp()

        for result_type in result_types:
            if result_type not in ['inflow', 'level', 'outflow']:
                print(f'{result_type} is not a valid variable. Skipping analysis'); continue
            # quantile_analysis = quantile_analysis_set[result_type]
            # output_name = self.outputfile.replace('.csv', '')
            output_name = f'{self.outputfile}_{result_type}'
            # mc.compute_quantiles(start_q=quantile_analysis['lower'],
            #                      end_q=quantile_analysis['upper'],
            #                      step_q=quantile_analysis['step'],
            #                      result_type=result_type,
            #                      output_filename=f'{output_name}.csv')
            mc.compute_std_quantiles(result_type=result_type,
                                     output_filename=f'{output_name}.csv')

            percentiles = mc.main_division_percentiles(result_type=result_type,
                                                       output_filename=f'{output_name}_perc.csv')
            percentiles = mc.smooth_percentiles(percentiles, f'{output_name}_perc_smooth.csv')
            mc.plot_tpt_results_2(result_type, f'{output_name}.png', percentiles=percentiles)
        # Rewrite the mcdf with the aep estimates for the representative event analysis
        mc.df.to_csv(sim_filename)

    def analyse_volumes(self, result_types=['inflow', 'outflow']):
        mc = self.mc
        sim_filename = self.outputfile + '__mcdf.csv'
        print('\nOpening Monte Carlo analysis file:', sim_filename)
        mc.df = pd.read_csv(sim_filename, index_col=0)
        print(mc.df)
        
        if self.model_type == 'URBS':
            results_folder = os.path.abspath(os.path.join(self.outputfile, '../../urbs_results'))
        elif self.model_type == 'RORB':
            results_folder = os.path.abspath(os.path.join(self.outputfile, '../../rorb_results'))
            
        for result_type in result_types:
            print('\nAnalysing the ', result_type, ' volumes')
            hyd_filename = os.path.basename(self.outputfile) + f'_{result_type}s.csv'
            filepath = os.path.join(results_folder, hyd_filename)
            print('\nOpening hydrograph file:', filepath)
            try:
                hyd_df = pd.read_csv(filepath, index_col=0)
            except:
                print('Hydrograph file not found.')
                return
            
            if self.model_type == 'RORB':
                # RORB simulation increment is automatically set to temporal pattern timestep
                # Therefore RORB hydrograph files may contain mismatched indices
                # Linear interpolation is used to fill data gaps. Potential interpolation error at end of time series ignored.
                hyd_df.interpolate(method = 'linear', axis = 0, inplace = True)
            
            time_inc = hyd_df.index[2] - hyd_df.index[1]
            print('Integrating volume using: ', time_inc, 'h time steps.')
            hyd_df.fillna(0, inplace=True)   # Replace NaNs with zero. For hydrographs shorter than required duration. Recommend extending simulation periods if practical.
            
            durations = [24, 36, 48, 72, 96, 120]
            data = []
            for duration in durations:
                window = duration / time_inc
                if window.is_integer(): 
                    window = int(window)
                else: 
                    print(f'Unable to integrate {duration}h duration with {time_inc}h time steps')
                    continue
                
                print(f'Calculating max {duration}h volumes')
                vol = hyd_df.rolling(window).apply(integrate.trapz).max() 
                vol = vol / 1000 * 60 * 60 * time_inc           # Convert from m3/s to ML
                data.append(vol)
            
            vol_df = pd.concat(data, axis=1)
            vol_df.columns = [f'Vol{dur}h' for dur in durations]
            vol_df.index = vol_df.index.to_series().map(lambda x: int(x.split('_')[1]))
        
            vol_df = pd.concat([mc.df[['m','n', 'rain_z', 'rain_aep']], vol_df], axis=1)
            volume_file = self.outputfile + '_volume.csv'
            print('Writing volume results:', volume_file)
            vol_df.to_csv(volume_file)
            mc.df = vol_df
            
            for duration in durations:
                output_name = f'{self.outputfile}_{result_type}Vol{duration}h'
                mc.compute_std_quantiles(result_type=f'Vol{duration}h', 
                                         output_filename=f'{output_name}.csv')
                mc.plot_tpt_results_2(result_type=f'Vol{duration}h',
                                      output_filename=f'{output_name}.png')


class Logger(object):
    def __init__(self, name):
        sys.stdout = sys.__stdout__
        # now direct the stdout to terminal and log file
        self.terminal = sys.stdout
        self.log = open(name, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


def life_of_brian():
    z = np.random.randint(low=1, high=11)
    print('\n----------------------------------------------------------------------------------------')
    quotes = {1: 'Blessed are the cheesemakers',
              2: 'All right, but apart from the sanitation, medicine, education, wine, public order, irrigation, roads, the fresh water system and public health, what have the Romans ever done for us?',
              3: "He's not the Messiah. He's a very naughty boy!",
              4: "Welease Bwyan!",
              5: "Welease Wodewick!",
              6: "Welease Woger!",
              7: 'Always look on the bright side of life!',
              8: 'I may be of thome athithtanthe if there ith a thudden crithith!',
              9: "Of course they've brought forth juniper berries! They're juniper bushes! What do you expect?!",
              10: 'Apart from the sanitation, the medicine, education, wine, public order, irrigation, roads, the fresh water system, and public health ... what have the Romans ever done for us? ... Brought peace?'}
    print(quotes[z])
    print('----------------------------------------------------------------------------------------')
