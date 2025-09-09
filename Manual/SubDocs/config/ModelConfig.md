# Model config file
This is the config file for managing the hydrologic model inputs and uses the keys in the tables below.

Firstly, the keys common to **both** the URBS and RORB models
| Config file | Description |
| ----------- | ----------- |
|```model_exe``` | Absolute path to the URBS or RORB executable. |
|```model_folder``` | Absolute path to the model folder (where the vec file is if using URBS). |
|```storms_folder``` | Absolute path to the folder containing the storm files. |
|```simulation_periods``` | Dictionary containing the length of time to run each storm duration event for; e.g. ```{"12": 48, "24": 48, "48": 72}```, where the key is the storm duration and the item is the simulation period. If a simulation period for a specific storm duration is missing, Bryan will use twice the storm duration. |
|```full_supply_volume``` | **Optional for the URBS model - depending on the dam routing method**. The full supply volume (ML) for the dam. For URBS models care should be taken to ensure that this matches the volume at the full supply level in the *els* file. If using the level--based method for dam routing in URBS, this key is not needed - the FSL is extracted from the vec file and the FSV is inferred from the els file.   

Secondly, the keys used in the **URBS** model
| Config file | Description |
| ----------- | ----------- |
|```vec_file``` | Name of the model file. |
|```ratings_folder``` | Absolute path to the URBS ratings tables.  |
|```store_tuflow``` | Use ```true``` if wanting to output files from URBS for the TUFLOW model. Otherwise, use ```false```.  |
|```result_prefix``` | Prefix to tack onto the URBS results files. |
|```vec_file``` | Name of the model file. |
|```time_increment```| Timestep (in hours) the URBS model will use to compute flows. |
|```time_increment_override```|(Optional) A dictionary that can be used to set time increments for specific storm durations. For example, if running shorter durations which have temporal patterns that use a shorter timesteps than the longer durations, the override can be used as follows: ```"time_increment_override": {"2": 0.05, "3": 0.25, "4.5": 0.25, "6": 0.25}```|
|```alpha``` | The URBS alpha parameter to apply to the model. |
|```beta``` | The URBS beta parameter to apply to the model. |
|```m_exponent``` | The catchment routing non-linearity exponent (usually 0.8). |

And lucky last, the keys used in the **RORB** model
>Graigan needs to check/correct/expand  the table below as needed.

 Config file | Description |
| ----------- | ----------- |
|```cat_file``` | Filename (or path?) to the catchment file. |
|```par_file``` | Filename (or path?) to the par file. |
|```results_folder``` | Absolute path to where the final results should be stored.  |

