### Getting started - the batch file,  main config file, and simulation list
As discussed above, the Python scripts sit in a separate location to the model files. The main script is accessed using a **batch file** that passes the main config file path for a project as a *key* to the Python script. This way, a copy of the Python code is not needed for each model (dam catchment). An example batch file is below.
```bat
:: Insert the name of the config file
set config_file="path\to\sims_config.json"

:: Insert the file of the python project
set pyfile="path\to\Bryan\Main.py"

:: Activate the Anaconda base environment
call path\to\Anaconda3\condabin\activate.bat

:: Run the model
python %pyfile% %config_file%
```
All config files are [JSON](config/json_files.md) files that set up the inputs to the Python scripts. The **main config file**, referenced in the batch file, tells Bryan which simulation list (Excel file) to use and provides filepaths to the model, storm, and climate config files. An example is shown below:

```python
{
    "simulation_list": "SimsList.xlsx",
	"filepaths": {
		"model_config": "model\\model_config.json",
		"storm_config": "storm_data\\storm_config.json",
		"climate_config": "climate_change\\climate_config.json"
	}
}
```
There is an optional key ```test run``` that is used for testing the code during development. If this key is included, the number of simulations in the monte carlo simulations will be killed early, only running the number of simulations provided to this key (e.g. ```test run: 4``` will exit the code after running four simulations). If ```test run: 0```, the simulations will not be killed early - the simulations will run through to completion as normal.  

From here, all simulations are controlled through the simulation list (Excel file, discussed below). A summary of the workflow is shown below. 

***batch file --> config file --> simulation list***

The **simulation list** is an Excel file containing a list of simulations to run. There are a number of fields used to set up each simulation. Two types of simulations can be performed:

- ***Monte Carlo*** simulations are set up for a single storm duration, which is specified in the ```Duration``` field in the simulation list. The simulation parameters are controlled using a [Monte Carlo config file](config/MonteCarloConfig.md) specified in the ```Config file``` field in the simulation list. 
- ***Ensemble*** simulations are run across a range of simulations specified in the [Ensemble config file](EnsembleConfig.md), which is specified in the ```Config file``` field in the simulation list. 

Other fields in the simulation list are explained [here](config/sim_list.md). 
