# Storm config file
This is the config file for setting up the design storms and uses the keys in the table below. 
**Table 1: Keys for the storm config file.**
| Key | Description |
| ----------- | ----------- |
|```file_paths``` | A dictionary of filepaths. Keys are listed in Table 2 below. |
|```storm_method_config```| A dictionary used to specify storm configurations. Keys are listed in Table 3 below.  |
|```CL_limit```|(Optional) A dictionary used to reduce the continuing loss in extreme storms (rarer than 1 in 100 AEP). For example, the folllowing command will reduce the CL to 3 mm per hour at the PMPF: ```"CL_limit":{"apply":true, "limit": 3}```. The ```apply``` key can be set to ```false``` to switch off this limit. If the applied CL is lower than the limit (e.g. applied CL of 1 mm and limit is 3 mm), the applied CL will be used. |


**Table 2: Keys for the ```file_paths``` dictionary.**
| Key | Description |
| ----------- | ----------- |
|```ARR_datahub_file``` | Relative path to the text (*.txt) file obtained from the ARR datahub. This is used to get the areal reduction factor region and preburst proportions.  |
|```rare_ifds```|Relative path to the IFD config file|
|```point_patterns``` | Relative path to the ARR point temporal patterns (csv file)|
|```areal_patterns``` | Relative path to the ARR areal temporal patterns (csv file)|
|```gsdm_patterns``` | Relative path to the GSDM temporal patterns (pat file)|
|```gtsmr_patterns``` | Relative path to the GTSMR temporal patterns (pat files). Note that separate files are provided for each catchment area bin. A keyword ```~AREA~``` is used in the filename here and replaced in the Python code with the relevant area based on  the catchment area provided. |
|```subcatchment_areas``` | Relative path to a table of subcatchment areas (csv file). This information is used to compute the catchment-average rainfall. |

**Table 3: Keys for the ```storm_method_config``` dictionary.**
| Key | Description |
| ----------- | ----------- |
| ```aep_changeover_to_extreme``` | A list setting the range in which the temporal patterns switch from ARR to GSDM/GTSMR (e.g. ```[100, 2000]```). The lower bound is excluded and the upper bound is included in the range. So, in the example, only ARR patterns will be used for the 1 in 100 AEP and for the 1 in 2000 AEP either ARR, GSDM, or GTSMR will be used. In this range, the method (ARR/GSDM/GTSMR) is randomly sampled for the Monte Carlo simulations. For the ensemble simulations, the method used in this range must be specified using the ```interim_for_ensemble``` key in the [ensemble config](EnsembleConfig.md.html) file. |
|```gsdm_gtsmr_changover_duration``` | Used to specify the changeover zone from GSDM to GTSMR due to storm duration in hours (e.g. ```[6, 24]```). The bounds are excluded. So, in the example, only GSDM will be used for 6 hours and only GTSMR for 24 hours. With the Monte Carlo method the GSDM and GTSMR patterns are randomly sampled in the changeover zone. For the ensemble method, the patterns to use must be specified using the ```extreme_pattern_for_ensemble``` key in the [ensemble config](EnsembleConfig.md.html) file.|
