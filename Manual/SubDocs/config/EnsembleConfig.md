# Ensemble config file
This is the main config file for the ensemble simulations and uses the keys in the table below.  
**Table 1: Keys for the storm config file.**
| Config file | Description |
| ----------- | ----------- |
|```storm_durations```| a list of the storm durations to simulate; e.g. ```[12, 24, 48]```
| ```aep_list```|A list of the AEPs to simulate; e.g. ```[50, 100, 1000]```|
|```storm_method_config```| A dictionary used to specify storm configurations. Keys are listed in Table 2 below.  |

**Table 2: Keys for the ```storm_method_config``` dictionary.**
| Filepath key | Description |
| ----------- | ----------- |
| ```interim_for_ensemble``` | The temporal pattern method to use in the AEP changeover range when running in ensemble mode. Use ```arr``` to specify the temporal patterns downloaded from the ARR datahub. Use ```extreme``` to specify GDSM/GTSMR (Bryan will pick based on the storm duration). Use ```both``` to use both ARR and extreme (all 20) patterns. |
|```extreme_pattern_for_ensemble``` | Should be either ```GTSMR``` or ```GSDM``` and specifies which temporal patterns to use in the ```gsdm_gtsmr_changover_duration``` range for the ensemble runs. |
