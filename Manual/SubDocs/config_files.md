# Config files
Bryan uses several config files for managing the catchment-specific inputs to the modelling. The config files use the [JSON](config/json_files.md) file format. The following config files are expected:

- [Simulation config file](config/SimsConfig.md) is the main config file, the one the batch file passes to Bryan. It names the simulation list and points at the model, storm and climate config files.
- [Monte Carlo config file](config/MonteCarloConfig.md) is used to set up the monte carlo framework, modelling method and result analysis parameters.
- [Ensemble config file](config/EnsembleConfig.md) is used to set up the ensemble method framework, modelling method and result analysis parameters.
- [IFD config file](config/IFDFilesConfig.md) is used to set up the rainfall depth information.
- [Storm config file](config/StormConfig.md) is used to set up the design storm parameters and inputs. 
- [Lake config file](config/LakeConfig.md) is used to set up the lake information including the antecedent conditions.
- [Model config file](config/hydrologic_modelling.md) is used to set up the hydrologic model details (URBS or RORB).
- [Climate change config file](config/ClimateConfig.md) is used to set up the climate change inputs. 