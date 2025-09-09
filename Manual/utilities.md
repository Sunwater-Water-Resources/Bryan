# Utilities for post-processing
Python scripts for post-processing are located in the ```util``` folder. It might be easiest to create a copy of the scripts in the location where the analyses are being done to avoid changing the original scripts - unless fixing bugs or enhancing. To play with the size of plots search for the ```dpi```  key in the scripts.

## Creating frequency plots
Use the ```PlotFrequencyCurves.py``` script to create frequency plots. Users can add as many model results (Monte Carlo or Ensemble) and FFA plots as needed. The inputs are specified in a spreasheet listing what should go into each plot. For an example spreadsheet, see ```_Plots_List_example.xlsx``` in the ```util``` folder. The only part of the script that should need editing is shown below:

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