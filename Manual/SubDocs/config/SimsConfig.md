# Simulation config file
This is the **main config file** - the one the [batch file](../getting_started.md) passes to
```Main.py```, usually named *sims_config.json*. It is short: it names the
[simulation list](../sim_list.md) to run and gives the filepaths of the three config files that
apply to the whole project. Everything else is set per simulation in the simulation list, or in
the [method config file](../config_files.md) that each row's ```Config file``` column points at.

```json
{
    "simulation_list": "SimsList.xlsx",
    "filepaths": {
        "model_config": "model\\urbs_config.json",
        "storm_config": "storm_data\\storm_config.json",
        "climate_config": "climate_change\\climate_config.json"
    },
    "test_runs": 0
}
```

**Table 1: Keys in the simulation config file**

| Key | Description |
| ----------- | ----------- |
| ```simulation_list```| The [simulation list](../sim_list.md) to run - an Excel file, read from its first sheet. Resolved against the **project folder** (see below). |
| ```filepaths```| A dictionary of the three project-wide config files: ```model_config``` (the [URBS](URBS_model.md) or [RORB](RORB_model.md) [model config](hydrologic_modelling.md)), ```storm_config``` (the [storm config](StormConfig.md)) and ```climate_config``` (the [climate change config](ClimateConfig.md)). Resolved against the folder holding **this config file** - see the note below, as that is not necessarily the project folder. |
| ```project_folder```| ++Optional:++ Where the ```simulation_list``` is looked for. Defaults to the folder holding this config file, which is what you want when the config sits in the project folder beside the batch file. The keyword ```default``` means the same as leaving it out. |
| ```test_runs```| ++Optional:++ Stops a monte carlo simulation after this many realisations, for testing the code without waiting for a full run. ```0```, or leaving the key out, runs every simulation through to completion. **The key is ```test_runs```**; earlier versions of this manual called it ```test run```, which Bryan does not read - a config written that way runs to completion with nothing said. |

> **Two base folders.** ```simulation_list``` resolves against the ```project_folder```, but the
> ```filepaths``` always resolve against the folder holding this config file, whatever
> ```project_folder``` says. They are the same folder unless ```project_folder``` is set, so this
> only bites when it is. Paths *inside* the simulation list are different again: they resolve
> against the working directory of the process, which is why the batch file starts with
> ```cd /D "%~dp0"```.

Everything in here is a path, and a wrong one is the usual reason a batch does not run. Bryan
prints each resolved path as it opens it, and the [run log](../getting_started.md)
records which config files each simulation used.

Return to [Config files](../config_files.md) or the [Main Manual](../../Manual.md)
