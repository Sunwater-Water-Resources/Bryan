# Simulation list
The simulation list provides a list of simulations to be performed. Note that the parent folder for the project is specified using the ```project folder``` key in the *sims_config.json* file. The fields contained in the list are as follows:

**Table 1: Fields used in the simulation list (Excel spreadsheet).** 
| Field | Description |
| ----------- | ----------- |
| ```Include```|If ```no```, this simulation is skipped. Else, ```yes``` to include this simulation. |
| ```Method```| Either ```ensemble``` or ```monte carlo```. |
| ```Duration```| The storm duration in hours if running a monte carlo simulation. Not used for ensemble simulations - specified in the [ensemble config file](.\config\EnsembleConfig.md).  |
| ```Run models```| If ```yes```, the RORB/URBS models are run for this simulation. Else, ```no``` if not wanting to run the models - only analysing existing results.|
| ```Analyse results```| If ```yes```, the results from the simulation are analysed: applying total probbaility theorem for monte carlo and creating box plots and tables for ensemble events. Else, ```no``` if not wanting to analyse the results. |
| ```Store hydrographs```| Type ```yes``` if wanting to store all the modelled hydrographs for this simulation in a single csv file. Useful if wanting to use some results as boundary conditions for a hydraulic model or to check the results. Though, storing the results for thousands of runs takes a bit more time. |
| ```Mop up files```| The monte carlo simulations produces thousands of storm files and results files from the URBS model. Type ```yes``` here to delete these files - keeping file management tidier. A few storm files are skipped in the deletion to aid checking. Note that the peak inflow, lake level and outflow for each run are stored and output seperately; so, are not deleted. |
| ```ADV```| The antecedend dam volume (ADV) for simulations. Three options are available: **1)** A fixed ADV can be set using a number (in ML). **2)** The full supply volume can be used by using the keyword ```fsv```. **3)** For the Monte Carlo method, a varying ADV can be used by using the keyword ```varying``` and supplying the lake config file (see below). |
|```Baseflow```|++Optional++: This key can be used to switch baseflow on and off - it overrides the ```apply``` switch in the *urbs_config* file. |
| ```Lake config```| Relative filepath to the [lake config file](config/LakeConfig.md), where the distribution and correlation for the ADV are set. Only used if the *varying* keyword is used in the ```ADV``` key (see above). |
|```Focal subcatchments``` | Filepath to a csv file containing a list of subcatchment areas upstream of the focal point for the hydrology. The area-weighted catchment average rainfall is computed from this (used for preburst depth computation and embedded burst filtering). The total catchment area (used for the ARF and selecting areal temporal patterns) is computed from the sum of these areas. |
| ```Config file```| Relative path to the main ensemble or monte carlo config file.|
| ```Replicates```| Used to replicate the monte carlo sampling in a previous simulation. Provide a list of components to replicate from a previous simulation; e.g. ```rz,tp,ilp,clp```. See Table 2 below for the key definitions. |
|```Replication file```| The realtive path to the file used to replictae the sample scheme for the monte carlo simulations. |
| ```Exclusions```| Used for sensitivity testing of excluding some aspects of the simulations; e.g. ```ebf.d50```. See Table 3 below for the key definitions. |
|```Log file``` | Relative path to the log file created for this simulation. |
|```Output file``` | Relative path for output files for this simulation. |
|```GWL```| ++Optional (used in Sunwater spec):++ The global warming level (in °C) relative to the hydrologic basleine period to apply to the simulation. Typically using 1.3°C for the near term, 1.8°C for the medium term, and 2.7°C for the long term. Note that if the ```GWL``` key is not used, the ```SSP``` and ```Year``` keys must be provided.|
| ```SSP```| ++Optional (not used in Sunwater spec):++ Shared Socioeconomic Pathway (SSP) to use for the global warming projection (climate chnage). Can be one of : SSP1-1.9, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5. |
| ```Year```| ++Optional (not used in Sunwater spec):++ Climate horizon for the simulation. If the horizon falls within the base period (usually 1961 to 1990), no climate change adjustments are made. If the horizon falls before the warming projection period begins (2015), historical temperature anomalies are used and the SSP doesn't influence the temperature rise. |


**Table 2: Keys used for the replication of monte carlo sampling.**
| Key | Description |
| ----------- | ----------- |
|```rz``` | The rainfall sampling - standard normal variate |
|```tp``` | The temporal pattern samping - integer from 0 to 9 |
|```stm```| The storm method sampling: ARR point, ARR areal, GSDM, and GTSMR. Use this one with care, as these are influenced by catchment size and storm duration. |
|```ilp```| The initial loss percentile (0-1) used for the scaling of the storm's initial loss.|
|```clp```| The continuing loss percentile (0-1) used for the scaling of the storm's continuing loss.|
|```pbp```| The sampling of the preburst percentile( 10% to 90%), which is then used to get the preburst proportion of the catchment average rainfall. |
|```lz```| The sampling of the antecedent lake volume - standard normal variate.|

**Table 3: Keys used for the excluding some compoentns of the modelling.**
| Key | Description |
| ----------- | ----------- |
|```pb``` | Exclude the residual preburst rainfall depth (after subtracting losses) from the storms. |
|```ebf``` | Exclude the filtering of embedded bursts |
|```d50```| Exclude the front loading shift in the temporal patterns due to climate change (probability weighted sampling). This shift is not applied in the ensemble simulations. |
|```ru```| Exclude the rainfall uplift due to climate change.|
|```lu```| Exclude the rainfall loss uplift due to climate change. |
|```clp``` | Exclude the sampling of continuing loss - a fixed CL is used. |
