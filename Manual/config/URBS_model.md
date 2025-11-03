# URBS wrapper

### Dam routing - level method
This is the preferred method. Dam routing in the URBS model can be set up using the full supply level (FSL) method with an initial level (IL) and the storage profile in a elevation-level-storage (els) file. This method will provide results when the peak lake level does not reach the full supply level. Bryan applies the IL using a label in the environment variables, called ```initial_lake_level```. An example of the dam routing command in URBS is below.
```
DAM ROUTE FSL=215.50 datafile=callide.els il=initial_lake_level location=CALLIDE FILE=CALLIDE.SQ
```

### Dam routing - volume method
Dam routing in the URBS model can be set up using the volume below full supply level (VBF) method. However, this method does not provide results when the lake level does not reach the full supply level and the full supply volume needs to be provided in the model config file using the ```full_supply_volume``` key. Bryan applies the VBF using a label in the environment variables, called ```ADV_below```. An example of the dam routing command in URBS is below.
```
DAM ROUTE VBF=ADV_below location=DAM FILE=DAM.SQ 
```

### Baseflow
Including **baseflow** is a little intricate. Do not include the baseflow parameters in URBS's catchment definition (_*.vec_) file. Instead, the parameters are inserted in the [URBS config](URBSModelConfig.md) file. An example of the baseflow parameter line in the config file is shown below.
```python
"baseflow": {"apply": true, "b0": 0.0, "br": 0.99, "bm": 1.0, "bfvf10_factor": 1.0},
```
The ```apply``` key can be used to turn the inclusion of baseflow on (```true```) or off (```false```). If the baseflow parameter line is not in the config file, the models will be set up without baseflow or whatever default parameters are in the *vec* file. The ```bfvf10_factor``` key is a calibration factor used (if needed) to scale the baseflow volume factor for an AEP of 10% (BFVF10) obtained from the ARR datahub to more closely match the expected baseflow constant (BC) obtained from calibration. If using baseflow, the text file downloaded from the ARR datahub must include the baseflow parameters for the catchment. Bryan will collect the *Volume Factor* from the ARR data file, which is the baseflow volume factor for a 10% AEP event (BFVF10), and multiply this by the ```bfvf10_factor```. The BFVF for the particular event magnitude will then be computed by Bryan as BFVF = BFVF10 x 10^r, where r = -0.02079z^2-0.1375z+0.2079 and z is the standard normal variate of the AEP. Finally, the BC parameter used in URBS is computed as BC=(1-BR)BFVF. Bryan will then create a copy of the *vec* file with the suffix *_baseflow* where the baseflow parameters are inserted into the default parameters, and use this *vec* file for running the simulations. 


### Result file management
URBS generates lots of result files, like a lot and a lot. Therefore, there are some **file management options** (deleting and storing results) controlled through the [simulation list](../sim_list.md). Keywords need to be configured in [Monte Carlo config](../MonteCarloConfig.md) and [Ensemble config](../EnsembleConfig.md) files for tracking the inflow, lake level, and outflow results printed by URBS; an example is given below. 
```python
"max_keys": {"inflow": "INFLOW", "level": "DAM", "outflow": "DAM"}
```
The peak results are then extracted for these three labels in the URBS ```*.p``` file. In addition, there is an option to store the full hydrographs of all simulations in a single csv file (```Store hydrographs``` key in [simulation list](../sim_list.md)), enabling the thousands of results files created by URBS to be deleted (```Mop up files``` key in [simulation list](../sim_list.md)). 

Regardless of whether these file management options are used, the peak iflow, lake level, and outflow are tracked and stored for each simulation in the simulation results file, which is then used to analyse the results.
