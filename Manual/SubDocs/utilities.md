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
### Extracting the events chosen in the launcher

```RepresentativeEvents.py``` is the other half of that job, for events chosen on the [run launcher's](ui.md) Events page rather than in a spreadsheet. It is a command line script, not a block of paths to edit:

```bat
python util\RepresentativeEvents.py --config sims_config.json ^
    --selection sims_mc\results\GWL1p3_representative_events.json
```

It reads the selection file the launcher saved, and for each chosen event writes:

| Output | Holds |
| ----------- | ----------- |
| ```<selection>_<loading>_sim<id>.png``` | a three-panel plot: the hyetograph at the top, drawn downwards from zero; the inflow and outflow; and the lake level |
| ```<selection>_events.xlsx``` | a sheet per series (```inflows```, ```levels```, ```outflows```, ```hyetographs```), one column per event, plus an ```events``` summary sheet |

**The hyetograph is rebuilt, not read.** An mcdf records what was *sampled* - the rainfall variate, the pattern number, the pre-burst proportion - not the rainfall series, so the storm generation is replayed for that one realisation. The replay is then checked against three things the run itself wrote down: the catchment average burst depth (```mean_rain_mm```), the pre-burst depth (```preburst_mm```) and the embedded burst comment. A rebuild that disagrees with any of them is reported on the plot and in the workbook rather than presented as the storm that was modelled. Use ```--no-hyetograph``` to skip the rebuild, which is much faster and works where the rainfall data is not to hand.

**Events chosen from a routed run take their storm from the run that was re-routed.** A ```reservoir routing``` row has no ```Duration``` and no ```Focal subcatchments``` - it re-routes hydrographs a previous run stored - so the script follows that row's ```Input MCDF``` back to the sims-list row that produced it and rebuilds from there, saying so in the notes. The realisation itself is never the problem: a routed database is the inherited one with the routed peaks written over it, so every draw the storm was made from is still in the row. Only those two keys are missing, and there are three ways to supply them:

- have the source run in the same simulation list, which is the usual case and needs nothing;
- point at the list it *is* in with ```--source-sims-list``` (repeatable) - routing rows commonly live in a sims list of their own;
- put a ```Duration``` and a ```Focal subcatchments``` on the routing row itself. The method ignores both, so they cost nothing and the rebuild will use them.

Failing all three the event still gets its hydrographs, and the note says which two keys would have fixed it.

Time on the plot runs **from the start of the main burst**, so the pre-burst is at negative times and events with different pre-burst durations can be compared. The stored hydrographs begin at the start of the storm file - that is, at the start of the pre-burst - so they are shifted by the pre-burst duration; where the hyetograph is not rebuilt, the shift comes from how much longer the run is than the simulation period in the model config, which Bryan lengthens by exactly that amount.

Options: ```--out``` for a different folder, ```--name``` to change the output basename, ```--no-plots``` for the workbook alone, and ```--dpi```.

## Calibrating temporal pattern weights
Use the ```CalibrateTpWeights.py``` script to calibrate temporal pattern probability weights so that the Monte Carlo ensemble satisfies the sub-burst AEP-neutrality condition - see [the sub-burst check](sub_burst_check.md) for the background. The script works entirely from the mcdf file (ideally from an unfiltered simulation run with ```Run models``` set to *storms only*): trial weights are evaluated by re-weighting the recorded realisations in the TPT, so no model reruns are needed. Patterns whose sub-bursts exceed the same-z IFD are progressively down-weighted until the weighted sub-burst frequency curves sit at or below the IFD. Outputs are the calibrated weights, the neutrality margins before and after calibration, and (if flow results are present in the mcdf) the weighted flood quantiles as a preview of the effect on the flood frequency curve. If the weights land on the floor without achieving neutrality, weighting alone is not enough for that simulation - consider the embedded burst filter, or review the offending patterns. For a production run, the calibrated weights file can be applied to the pattern sampling itself using the ```TP weights``` key in the [simulation list](sim_list.md). Note that the calibration is against **main-burst** sub-burst depths only - the pre-burst is not scanned, so weights that achieve neutrality say nothing about the pre-burst-inclusive storm. The only part of the script that should need editing is shown below:

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
