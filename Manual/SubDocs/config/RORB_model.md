# RORB wrapper

The RORB wrapper is set up through the [model config file](ModelConfig.md), which names the catchment file with the ```cat_file``` key and the parameter file with the ```par_file``` key. **Both are filenames, not paths** - Bryan looks for each of them inside the folder given by ```model_folder```.

Note also that ```full_supply_volume``` is **required** for RORB, whereas for URBS it is only needed when the volume method is used for dam routing.

### The model files are never modified
Bryan copies the catchment (_*.catg_) file into the ```storms_folder``` when it starts, and every edit it makes - the drawdown written into the dam routing line, for instance - is made to that copy. The file in the ```model_folder``` is only ever read. The parameter (_*.par_) file is treated the same way: the one in the model folder is read once as a template, and per-run copies are written into the storms folder.

If a sub-folder is used for the storm files, it is **deleted and recreated** at start-up. Do not point it at a folder holding anything worth keeping.

### One par and batch file per simulation
RORB reads its par file as it starts up, and ```cmd.exe``` reads a _*.bat_ file incrementally while executing it. Neither can safely be rewritten while a previous run might still be reading it, so Bryan writes a separate pair per simulation, named ```_run_<result name>.par``` and ```_run_<result name>.bat``` in the storms folder.

Both are deleted once the run succeeds. **A failed run leaves its pair behind**, which is deliberate - the batch file can be run by hand to reproduce the failure with exactly the inputs Bryan used.

The result name reproduces RORB's own convention: the stem of the catchment file, an underscore, then the stem of the storm file.

### Rainfall losses
The initial and continuing losses are not written into the storm file. Bryan writes them into the per-run par file, one line per interstation area, so every interstation area of a given simulation carries the same pair of losses. The number of interstation areas is read from the template par file rather than being configured.

### Dam routing
Bryan finds the dam routing line in the catchment file by looking for a line beginning with ```16``` whose following line matches the dam location name it was given. The drawdown - the volume below full supply level - is written into the second comma-separated field of the line after that.

This means the **dam location name must match the name in the catg file exactly**, including case.

### Waiting for results
RORB, or the operating system on a slow or network drive, can still be flushing the _*.out_ file when the batch process returns. Bryan waits until that file exists, is not empty, and is no longer held open by a writer before reading it. This is why a RORB simulation sometimes appears to pause after the model window closes.

### Simulation periods
As for URBS, the length of run for each storm duration comes from the ```simulation_periods``` key of the model config file. Where a duration is not listed, Bryan uses twice the storm duration and prints a note saying so.

Return to [Hydrologic modelling](hydrologic_modelling.md)
