# Utilities for post-processing
Python scripts for post-processing are located in the ```util``` folder. It might be easiest to create a copy of the scripts in the location where the analyses are being done to avoid changing the original scripts - unless fixing bugs or enhancing. To play with the size of plots search for the ```dpi``` key in the scripts.

## Creating frequency plots
Use the ```PlotFrequencyCurves.py``` script to create frequency plots. Users can add as many model results (Monte Carlo or Ensemble) and FFA plots as needed. The inputs are specified in a spreadsheet listing what should go into each plot. For an example spreadsheet, see ```_Plots_List_example.xlsx``` in the ```util``` folder. The only part of the script that should need editing is shown below:

```python
# Insert the folder and plot list filenames below
folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E001\plots'
plot_list_files = ['_Plots_List_01.xlsx']
min_aep = 2  # can be specified in the spreadsheet or here
max_aep = 5000000 # can be specified in the spreadsheet or here
```
 

## Getting representative events
Use the ```GetRepresentativeEvents.py``` script to extract summary information for selecting the representative events. It is recommended this analysis is done in a ```representative_events``` folder in the ```sims_mc``` folder. Inputs to the analysis are obtained from a spreadsheet; see the ```_analyseRepresentativeEvents.xlsx``` spreadsheet for an example. The only part of the script that will need changing is shown below (located at the top of the ```main()``` function:

```python
folder = r'C:\PythonProjects\TFD_2024\03_Design\runs\E001\sims_mc\representative_events'  # change this to your folder
filename = '_analyseRepresentativeEvents.xlsx'  # change this to your file
```
## Calibrating temporal pattern weights
Use the ```CalibrateTpWeights.py``` script to calibrate temporal pattern probability weights so that the Monte Carlo ensemble satisfies the sub-burst AEP-neutrality condition - see [the sub-burst check](sub_burst_check.md) for the background. The script works entirely from the mcdf file (ideally from an unfiltered simulation run with ```Run models``` set to *storms only*): trial weights are evaluated by re-weighting the recorded realisations in the TPT, so no model reruns are needed. Patterns whose sub-bursts exceed the same-z IFD are progressively down-weighted until the weighted sub-burst frequency curves sit at or below the IFD. Outputs are the calibrated weights, the neutrality margins before and after calibration, and (if flow results are present in the mcdf) the weighted flood quantiles as a preview of the effect on the flood frequency curve. If the weights land on the floor without achieving neutrality, weighting alone is not enough for that simulation - consider the embedded burst filter, or review the offending patterns. The only part of the script that should need editing is shown below:

```python
mc_config_file = r'C:\path\to\mc_config.json'   # the monte carlo config (for the scheme_config)
mcdf_files = {                                  # label: mcdf path - one entry per simulation
    '24h': r'C:\path\to\outputs\mc_24h__mcdf.csv',
}
output_folder = r'C:\path\to\outputs\tp_weights'

margin_target = 1.0      # calibrate until all tested margins <= this
aep_range = None         # e.g. [200, 500000] to set the AEP range for the convergence test
reduction_factor = 0.5   # per-iteration down-weighting: w <- w * reduction_factor**breach_fraction
weight_floor = 0.02      # minimum pattern weight
max_iterations = 25
```
