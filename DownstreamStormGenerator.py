import pandas as pd
from lib.StormGenerator import StormBurst
from lib.URBSmodel import UrbsModel
from lib.InterpolationCurves import *
import os
import json
from scipy.special import ndtri
from scipy import interpolate
import numpy as np
from lib.Simulator import Logger
import sys
from lib.ClimateChange import ClimateAdjustment
import datetime


def main():
    config_folder = r"C:\PythonProjects\TFD_2024\03_Design\runs\E012\US_E012\sims_mc\representative_events\Downstream_storm_files"

    log_file = os.path.join(config_folder, '__log.txt')
    print('Creating log file:', log_file)
    sys.stdout = Logger(log_file)
    current_time = datetime.datetime.now()
    print("Time now is:", current_time)

    # pmf_events(config_folder)
    main_events(config_folder)


def pmf_events(config_folder):
    catchment_approaches = ['sub-area', 'full-area']
    correlations = ['c24h', 'c48h']
    config_file_template = '_config_-CATCHMENT-_L+0d_-COR-.json'
    for approach in catchment_approaches:
        for correlation_duration in correlations:
            print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('Correlation duration:', correlation_duration)
            print('Catchment approach:', approach)
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            filename = config_file_template.replace('-COR-', correlation_duration)
            filename = filename.replace('-CATCHMENT-', approach)
            config_path = os.path.join(config_folder, filename)
            all_depths = create_storm_files(config_path, approach)
            output_path = config_path.replace('.json', '_depths.csv')
            print('\nWriting the applied rainfall depths to file:', output_path)
            all_depths.to_csv(output_path)


def main_events(config_folder):
    lags = ['L+0d', 'L+1d']
    # lags = ['L+0d', 'L+1d', 'L-1d']
    # lags = ['L-1d']
    # catchment_approaches = ['sub-area', 'full-area']
    catchment_approaches = ['full-area']
    config_file_template = '_config_-CATCHMENT-_-LAG-.json'

    for approach in catchment_approaches:
        for lag in lags:
            print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('Rain lag:', lag)
            print('Catchment approach:', approach)
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            filename = config_file_template.replace('-LAG-', lag)
            filename = filename.replace('-CATCHMENT-', approach)
            config_path = os.path.join(config_folder, filename)
            all_depths = create_storm_files(config_path, approach)
            output_path = config_path.replace('.json', '_depths.csv')
            print('\nWriting the applied rainfall depths to file:', output_path)
            all_depths.to_csv(output_path)


def create_storm_files(config_path, sub_folder=''):

    # Set up the config info
    print('Opening the config file:', config_path)
    f = open(config_path)
    config_data = json.load(f)
    f.close()
    filepaths = config_data['filepaths']
    folder = filepaths['folder']

    # Get the event info
    event_df = pd.read_excel(os.path.join(folder, filepaths['event_info']), sheet_name='event_list')
    event_df = event_df.loc[event_df['Include'] == 'yes']
    print(event_df)

    # Get the rainfall data and write the ifd files
    rain_files = filepaths['rain']
    aep_files = filepaths['aep_files']
    output_aeps = []
    for duration_str, rain_file in rain_files.items():
        print(f'\nGetting IFD data for {duration_str} hour storm duration:')
        ifd_filepath = os.path.join(folder, rain_file)
        print(ifd_filepath)
        ifd_data = get_ifd_depths(ifd_filepath)
        print(ifd_data)

        print(f'\nGetting AEP matrix for {duration_str} hour storm duration')
        # aep_filepath = ifd_filepath.replace('rain', 'aep')
        aep_filepath = os.path.join(folder, aep_files[duration_str])
        print(aep_filepath)
        aep_matrix = pd.read_csv(aep_filepath, index_col=0)
        print(aep_matrix)

        print(f'Creating the rainfall depths for {duration_str} hour storm duration')
        rain = get_rain_depths(ifd_data, aep_matrix)
        print(rain)
        print(f'\nOutputting the rainfall for record keeping...')
        output_filepath = ifd_filepath.replace('rain', 'processed_rain')
        rain.to_csv(output_filepath)

        print(f'\nCreating the IFD files across all subcatchments for {duration_str} storm')
        output_aeps = write_rainfall_ifd(rain, duration_str, folder, event_df, filepaths)

    # Create the model object
    model = UrbsModel(os.path.join(folder, filepaths['urbs_config']))
    model.storms_folder = os.path.join(model.storms_folder, sub_folder)
    print('\nStorm folder', model.storms_folder)

    # event_df.set_index('Rain AEP', inplace=True)
    # for ind, aep in enumerate(output_aeps):
    all_depths = pd.read_csv(os.path.join(folder, filepaths['catchments']), index_col=0).drop_duplicates(subset='Subarea')
    for ind, event in event_df.iterrows():
        # Get storm data
        print(f'\nWorking on row {ind}:')
        print(event)
        #  need to loop through the events for event in events
        storm_duration = event_df.loc[ind, 'Duration']
        result_filename = event_df.loc[ind, 'Result filename']
        result_folder = event_df.loc[ind, 'Result folder']
        event_id = event_df.loc[ind, 'Event']
        print('Folder:', result_folder)
        print('Filename:', result_filename)
        result_filepath = os.path.join(result_folder, result_filename)
        parameter_df = pd.read_csv(result_filepath, index_col=0)
        print('\nChecking for embedded busts in the main burst...')
        embedded_bursts = parameter_df.loc[event_id, 'embedded_bursts']
        print(embedded_bursts)
        if embedded_bursts != 'No embedded bursts':
            print('WARNING: this temporal pattern has an embedded burst')

        # Get the climate scaling
        climate_config_file = os.path.join(folder, filepaths['climate_config'])
        print('Opening the climate config file:', climate_config_file)
        climate_obj = ClimateAdjustment(config_file=climate_config_file,
                                        method='gwl',
                                        gwl=event_df.loc[ind, 'GWL'])
        climate_scaling = climate_obj.get_rainfall_uplift_factor(event_df.loc[ind, 'Duration'])
        # log_filepath = os.path.join(result_folder, event_df.loc[ind, 'Log file'], )
        # print('\nOpening the log file:', log_filepath)
        # climate_scaling = 1
        # with open(log_filepath) as f:
            # lines = [line.rstrip() for line in file]
        #     for line in f:
        #         rain_key = 'Rainfall uplift factor is:'
        #         if line.startswith(rain_key):
        #             climate_scaling = float(line.strip(rain_key))
        #             print(line.strip())

        tp_sample = parameter_df.loc[event_id, 'tp']
        il = parameter_df.loc[event_id, 'initial_loss'] - parameter_df.loc[event_id, 'preburst_mm']
        if il < 0.0:
            il = 0.0
        # NOTE: no need to scale il and cl for climate change - done by Bryan already in the mcdf
        il = np.around(il * event_df.loc[ind, 'IL scaling'], 1)
        cl = parameter_df.loc[event_id, 'continuing_loss']
        print(f'\nContinuing loss is {cl} mm/hr')
        cl_scaling = event_df.loc[ind, 'CL scaling']
        print('scaling continuing loss by ', cl_scaling)
        cl = np.around(cl * cl_scaling, 1)
        if cl < event_df.loc[ind, 'CL limit']:
            cl = event_df.loc[ind, 'CL limit']
        print(f'Scaled continuing loss is {cl} mm/hr')
        storm_method = parameter_df.loc[event_id, 'storm_method']

        # Set up the storm
        storm = StormBurst(os.path.join(folder, filepaths['storm_config']))
        storm.area = config_data['catchment_area']
        storm.storm_initial_loss = il
        storm.continuing_loss = cl
        storm.import_rare_rainfall()
        # storm.apply_areal_reduction() --> already applied in the provided IFDs
        storm.import_arr_areal_patterns([storm_duration])
        storm.import_arr_point_patterns([storm_duration])
        storm.import_gsdm_temporal_patterns([storm_duration])
        storm.import_gtsmr_temporal_patterns([storm_duration])

        # get the rainfall for current AEP
        rain_sample_z = parameter_df.loc[event_id, 'rain_z']
        rain_sample_aep = int(np.around(parameter_df.loc[event_id, 'rain_aep'], 0))
        print(f'\nGetting the rainfall for 1 in {rain_sample_aep} AEP')
        rain_depths = storm.rainfall.get_depth_aep(aep=rain_sample_aep,
                                                   duration=storm_duration,
                                                   storm_method=storm_method)
        rain_depths = rain_depths * climate_scaling

        # Track the applied rainfall depths across all events -> output as csv at the end
        current_depths = rain_depths.loc[all_depths.index]
        all_depths[event['Output filename']] = current_depths

        # get the temporal pattern
        temporal_pattern = storm.get_temporal_pattern(storm_method=storm_method,
                                                      duration=storm_duration,
                                                      tp_sample=tp_sample,
                                                      rain_sample_z=rain_sample_z)

        # Apply the temporal pattern shift
        original_simulation_period = model.simulation_periods[str(storm_duration)]
        temporal_shift = event_df.loc[ind, 'Rain lag']
        rain_end_time = original_simulation_period
        if temporal_shift != 0:
            temporal_shift = temporal_shift * 24  # convert from days to hours
            print(f'\nApplying a rain shift of {temporal_shift} hours.')
            temporal_pattern = shift_temporal_pattern(temporal_shift, temporal_pattern)
            print(f'Unshifted rain period is {rain_end_time} hours')
            rain_end_time += np.absolute(temporal_shift)
            print(f'Shifted rain period is {rain_end_time} hours')
            model.simulation_periods[str(storm_duration)] = rain_end_time

        # write the storm file
        storm_filename = event_df.loc[ind, 'Output filename']
        timesteps = storm.timesteps
        simulation_period = model.get_simulation_period(storm_duration)
        model.create_storm_file(rainfall_depths=rain_depths,
                                temporal_pattern=temporal_pattern,
                                data_interval=timesteps,
                                filename=storm_filename,
                                run_duration=simulation_period,
                                storm_duration=storm_duration,
                                storm_duration_excl_pb=storm_duration,
                                initial_loss=il,
                                continuing_loss=cl,
                                ari=int(np.around(rain_sample_aep)),
                                ensemble=tp_sample)

        # Reset the simulation period
        model.simulation_periods[str(storm_duration)] = original_simulation_period

        # Write the URBS inflow (*.o) file
        if '_mcdf' in result_filename:
            urbs_filename = result_filename.replace('_mcdf', 'outflows')
        else:
            urbs_filename = result_filename.replace('.csv', '_outflows.csv')
        urbs_filepath = os.path.join(result_folder, event_df.loc[ind, 'URBS folder'], urbs_filename)
        print('\nOpening the URBS results:', urbs_filepath)
        outflows = pd.read_csv(urbs_filepath, index_col=0)
        event_col = 'sim_{}'.format(str(event_id).zfill(5))
        print('Getting simulation id:', event_col)
        outflow = outflows[[event_col]].dropna()
        outflow.index = outflow.index - (outflow.index[-1] - event_df.loc[ind, 'Upstream period'])
        outflow = outflow.loc[outflow.index >= 0]
        o_file_lines = [config_data['o_file_metadata']]
        o_file_lines.append('Monte Carlo simulation from Bryan: {}\n'.format(event_df.loc[ind, 'Result filename']))
        o_file_lines.append('From upstream event: {}, {} hr, and GWL of {}°C\n'.format(event_id,
                                                                                     storm_duration,
                                                                                     event_df.loc[ind, 'GWL']))
        start_time = 0
        flow_timestep = int(event_df.loc[ind, 'Flow timestep'] * 3600)
        flow_steps = outflow.shape[0]
        extra_steps = int(np.absolute(temporal_shift / flow_timestep))
        # if temporal_shift < 0:
        #     start_time = int(-1 * temporal_shift * 3600)
        o_file_lines.append(f'{start_time} {flow_timestep} {flow_steps + extra_steps}\n')
        if temporal_shift < 0:
            for i in range(extra_steps):
                o_file_lines.append('     0.000\n')
        for i in range(flow_steps):
            # if outflow.index[i] <= rain_end_time:
            flow = '{:0.3f}'.format(np.around(outflow.iloc[i, 0], 3))
            o_file_lines.append(f'{flow.rjust(10)}\n')
        # if temporal_shift > 0:
        #     for i in range(extra_steps):
        #         # if outflow.index[i] <= rain_end_time:
        #         o_file_lines.append('     0.000\n')
        print(outflow)

        if event_df.loc[ind, 'Name'] == 'PMF':
            inflow_filename = f'PMF_{str(event_id)}'
        else:
            inflow_filename = str(event_id).zfill(6)

        inflow_filename = '{}_{}h_GWL{}'.format(inflow_filename,
                                                storm_duration,
                                                str(event_df.loc[ind, 'GWL']).replace('.', 'p'))

        if 'Lag string' in event_df.columns:
            if pd.notna(event_df.loc[ind, 'Lag string']):
                inflow_filename = '{}_{}'.format(inflow_filename, event_df.loc[ind, 'Lag string'])

        if 'Correlation duration' in event_df.columns:
            if pd.notna(event_df.loc[ind, 'Correlation duration']):
                correlation_string = 'c{}h'.format(event_df.loc[ind, 'Correlation duration'])
                inflow_filename = '{}_{}'.format(inflow_filename, correlation_string)

        inflow_folder = os.path.join(folder, 'inflows', sub_folder)
        directory_name = os.path.join(inflow_folder, inflow_filename)

        try:
            os.mkdir(directory_name)
        except FileExistsError:
            print(f"Directory '{directory_name}' already exists.")
        output_filepath = os.path.join(directory_name, event_df.loc[ind, 'Inflow filename'])
        with open(output_filepath, 'w') as f:
            for line in o_file_lines:
                f.write(line)

    all_depths.set_index('Subarea', inplace=True)
    all_depths = all_depths.transpose()
    return all_depths


def shift_temporal_pattern(shift, temporal_pattern):
    print(f'\nShifting the temporal pattern by {shift} hours. Unshifted pattern:')
    print(temporal_pattern)
    timestep = temporal_pattern.index[0]
    steps = int(np.around(shift / timestep, 0))
    print(f'\nPattern timestep is {timestep} hours... adding {steps} time steps')
    if steps > 0:
        print('Padding out the front end...')
        temporal_pattern.index = temporal_pattern.index + shift
        for i in range(1, steps + 1):
            temporal_pattern.loc[i * timestep] = 0.0
    else:
        print('Padding out the back end...')
        steps = steps * -1
        end_time = temporal_pattern.index[-1]
        for i in range(1, steps + 1):
            temporal_pattern.loc[end_time + i * timestep] = 0.0
    temporal_pattern.sort_index(inplace=True)
    print('\nShifted pattern:')
    print(temporal_pattern)
    return temporal_pattern


def get_ifd_depths(rain_filepath):
    # rain_path = os.path.join(folder, rain_file)
    rain = pd.read_csv(rain_filepath, index_col=0)
    aeps = rain.index
    log_rain = np.log10(rain)
    log_rain['z'] = ndtri(1 - 1 / rain.index)
    log_rain.set_index('z', inplace=True)
    log_rain.interpolate(method='slinear', inplace=True)
    rain = 10 ** log_rain
    rain['AEP'] = aeps
    rain.set_index('AEP', inplace=True)
    return rain


def write_rainfall_ifd(rain, duration_str, folder, event_df, filepaths):
    # print(f'\nGenerating IFD file for {duration_str} storm')
    # rain_path = os.path.join(folder, rain_file)
    # rain = pd.read_csv(rain_path, index_col=0)
    largest_aep = rain.index[-1]
    # aeps = rain.index.to_numpy()
    print('\nRainfall prior to interpolation')
    print(rain)
    do_interpolation = False
    do_extreme = False
    output_aeps = []
    extreme_aeps = []
    for event_ind, event in event_df.iterrows():
        aep = event['Rain AEP']
        if aep not in output_aeps:
            output_aeps.append(aep)
            if aep not in rain.index:
                do_interpolation = True
                if aep > largest_aep:
                    do_extreme = True
                    extreme_aeps.append(aep)
                rain.loc[aep] = np.nan
            rain.sort_index(inplace=True)
    if do_interpolation:
        log_rain = np.log10(rain)
        log_rain['z'] = ndtri(1 - 1 / log_rain.index)
        log_rain.set_index('z', inplace=True)
        log_rain.interpolate(method='slinear', inplace=True)
        log_rain['Rain AEP'] = rain.index
        log_rain.set_index('Rain AEP', inplace=True)
        rain = np.around(10 ** log_rain, 1)
        if do_extreme:
            print('\nThere are some AEPs that need to be extracted - using GEV:')
            print(extreme_aeps)
            # extreme_depths = rain.dropna().iloc[-3:]
            # print(extreme_depths)
            curves = {}
            for catchment in rain.columns:
                curve = GEV()
                curve.setup_rainfall_boundaries(aep_2000=rain.loc[2000, catchment],
                                                aep_1000=rain.loc[1000, catchment],
                                                pmp_depth=rain.loc[largest_aep, catchment],
                                                pmp_aep=largest_aep)
                curve.fit_curve()
                print('Curve fitted for catchment:', catchment)
                print(f'Location: {curve.location} | Scale: {curve.scale} | Shape: {curve.shape}')
                for extreme_aep in extreme_aeps:

                    rain.loc[extreme_aep, catchment] = curve.get_quantile(extreme_aep)
                    print(rain.loc[extreme_aep, catchment])

        print('\nRainfall after to interpolation')
        print(rain)

    # Get the spatial data
    pattern = pd.read_csv(os.path.join(folder, filepaths['catchments']), index_col=0)
    # print(pattern)

    # Create the ifd files
    output_folder = filepaths['output_folder']
    aep_col = [f'{aep}_AEP' for aep in output_aeps]
    ifd = pd.DataFrame(index=pattern.index, columns=aep_col)
    for aep in output_aeps:
        for catchment in pattern.index:
            zone = pattern.loc[catchment, 'Subarea']
            current_rain = rain.loc[aep, zone]
            ifd.loc[catchment, f'{aep}_AEP'] = current_rain
    # print(ifd)
    output_path = os.path.join(folder, output_folder, f'ds_ifd_{duration_str}h.csv')
    print('Writing the IFD file to be loaded into the storm object:', output_path)
    ifd.to_csv(output_path)
    return output_aeps


def get_rain_depths(ifd_data, aep_matrix):
    rain = aep_matrix.copy()
    base_aeps = ifd_data.index.to_numpy()
    base_z = ndtri(1 - 1 / base_aeps)
    for catchment in ifd_data.columns:
        base_ifd = ifd_data[catchment].to_numpy()
        base_log_ifd = np.log10(base_ifd)
        # for aep in aep_matrix.index:
        catchment_aeps = aep_matrix[catchment].to_numpy()
        catchment_z = ndtri(1 - 1 / catchment_aeps)
        catchment_log_depths = np.interp(catchment_z, base_z, base_log_ifd)
        catchment_depths = 10 ** catchment_log_depths
        rain[catchment] = catchment_depths
    return rain


if __name__ == "__main__":
    main()
