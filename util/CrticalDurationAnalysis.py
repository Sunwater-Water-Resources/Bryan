"""
This script is used to do the critical duration analysis
GWL refers to  the global warming level used for including climate change
"""

from UtilModule import MonteCarloSimulation, MonteCarloSimulationGroup
import os


# File information
folder = r"Directory containing the Monte Carlo results (*_mcdf.csv)"
filename_template = 'example-model-name_-DUR-h_GWL-TEMP-_result'
suffix = ''  # If needed, a suffix can be added to the critical duration outputs
drop_aeps = [2, 2000000]  # AEPs to exclude from plotting, ignored if empty list.
storage_filepath = r"directory conataining URBS els file\your.els"  # used to compute storage volume results from the level results

# Simulation information
gwls = [0, 1.3, 1.7, 2.7]  # in degrees Celcius
durations = [2, 3, 4.5, 6, 9, 12, 18, 24, 36, 48, 72, 96, 120]  # in hours
result_types = ['inflow', 'level', 'outflow']

# Set up the simulation object
sims = MonteCarloSimulationGroup(folder, drop_aeps=[2, 2000000])
if storage_filepath is not None:
    sims.set_volume_curve(storage_filepath)
    result_types.append('volume')

# Assess the critical durations
for gwl in gwls:
    gwl_str = str(gwl).replace('.', 'p')
    for result_type in result_types:
        output_name = filename_template
        output_name = output_name.replace('result', result_type)
        output_name = output_name.replace('-TEMP-', gwl_str)
        output_name = output_name.replace('-DUR-', '00')
        output_name = f'{output_name}{suffix}'
        sims.output_name = output_name

        for duration in durations:
            duration_str = str(duration).zfill(2)
            this_filename = filename_template.replace('-DUR-', duration_str)
            this_filename = this_filename.replace('result', result_type)
            this_filename = this_filename.replace('-TEMP-', gwl_str)
            filepath = os.path.join(folder, f'{this_filename}.csv')
            simulation = MonteCarloSimulation(filepath, result_type, duration)
            sims.add_simulation(simulation, duration)

        sims.compute_critical_durations()
