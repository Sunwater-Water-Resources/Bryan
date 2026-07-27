# 16 August 2024 -- Richard
- Removed the ```catchment area``` key from the storm config file. Bryan now computes the catchment area from the provided list of catchment areas.
- List of catchment areas is no longer provided in the storm config. It is now provided (as a filepath) in the sims list with the ```Focal subcatchments``` key, and should only include the subcatchments upstream of the focal point. I'm not sure of where the best location for this file is. I created a filepath called ```sim_options/focal_locations/dam.csv```. But there may be a better way.
- Folder locations are no longer needed in the sims config and mc config files. All filepaths are relative to the batch file used to run Bryan - including the filepaths provided in the sims list.
- The path to the lake config file is now provided in the sims list using the ```Lake config``` key and is only needed for monte carlo simulations where the ```ADV``` key is set to ```varying```. I'm not sure of where the best location for the lake config file is. I created a filepath called ```sim_options/lake_conditions/lake_config.json```. But there may be a better way.
- The ```ADV``` key in the sims list is now the place where the ADV method is defined. Use:
    - A number to set a fixed ADV (monte carlo or ensemble method)
    - ```fsv``` to set the ADV to match the full supply volume (monte carlo or ensemble method)
    - ```varying``` to use a stochastically sampled ADV (monte carlo only). With this method, the lake config file must be provided (see above).
- The full supply volume is now provided in the model config file using a ```full_supply_volume``` key. 
- These changes mean that the lake config file has been simplified to, for example:
```python
{
  "exceedance_layer_info": {"type": "sigmoid", "coefficients":
                                {"k": 2.3895, "Vf":  20000.0, "H": 6.1502, "z0": 1.2956, "Vc": 136400}}
}
```
- The storm IL and CL have been moved from the storm config file to the sims list using the ```IL``` and ```CL``` keys.
- These changes mean that the storm config file only includes rainfall depth and pattern information, and will generally not change across different existing/design/historic scenarios. 
- Similarly, the model config file will not change unless there is a structural change to the model (e.g. a design scenario). The sims config file is now simplified to:
```python
{
    "simulation_list": "SimsList_2.xlsx",
	"filepaths": {
		"model_config": "urbs\\urbs_config.json",
		"storm_config": "storm_data\\storm_config.json",
		"climate_config": "climate_change\\climate_config.json"
	}
}
```
- Fixed up some bugs with the ```Store hydrographs``` option in the sims list.
- Removed ```results_folder``` key from the urbs config file. Not sure if this is needed in the RORB config? It was being used to provide the location where the hydrographs and/or URBS results were moved. But not moving URBS results anymore (they are not needed) and hard coded hydrographs to be stored in a ```urbs_results``` folder one level down from the *mcdf* output location (i.e. one folder down from the ```analysis``` folder).
- Added stochastic sampling of the continuing losses, which is on as default. Can be fixed using the ```Exclusions``` key in the sims list.
- Added a check on the storm durations to determine if the GSDM/GTSMR patterns (spatial and temporal) will be needed. If they are not needed, the GSDM/GTSMR inputs are skipped - so spatial scaling does not need to be provided by the user.
 
# 21 August 2024 -- Richard
- Added the baseflow component. This is a little intricate, as URBS did not seem to use the BFVF in the header of the storm file. So, as a workaround, Bryan creates a copy of the vec file and inserts the user-provided parameters in the URBS config file and gets the BFVF10 from the ARR datahub text file. This is all documented in the revised manual. 
- Changed the ```rorb_exe``` and ```urbs_exe``` keys to ```model_exe``` in the *model config* files to make the nomenclature more model ambiguous.  
- Changed the Ensemble event config file. Now only includes two keys, which are used to set the AEPS and durations for the simulation. Example shown below:

```python
{
  "storm_durations": [12, 24],
  "aep_list": [50, 100, 200]
}
```

# 22 August 2024 -- Richard
- Added the level--based dam routing method for URBS as an option. No config required. Bryan will look in the vec file to see what method is being used. If the 'FSL=' term is found on the dam routing line, the level-based method is assumed and:
    - The FSL is obtained from the vec file
    - The ELS file path is obtained from the vec file, then opened and the storage curve read into Bryan
    - The FSV is inferred from the FSL
    - The ```full_supply_volume``` key is no longer needed in the model config file. 
- The advantage of the level--based method is that when the peak lake level does not reach FSL, URBS will provide the actual peak lake level; the volume-based method just outputs the FSL level. Also, this method handles ADV levels higher than the FSV, if needed. This should probably become the default.

# 23 August 2024 -- Richard
- Moved the ensemble method specific keys in the ```storm_method_config``` key from the storm config file to the ensemble config file. Example of an ensemble config file's content is below. So, the ```interim_for_ensemble``` and ```extreme_pattern_for_ensemble``` keys can be removed from the storm config file.

```python
{
  "storm_durations": [12, 24],
  "aep_list": [50, 100, 1000, 2000, 5000],
  "storm_method_config": {
	"interim_for_ensemble": "arr",
	"extreme_pattern_for_ensemble": "GTSMR"
  }
}
```

- For the ```interim_for_ensemble``` key (see above), any keyword other than *arr* was interpreted as use GSDM or GTSMR depending on duration. I have changed this so that Bryan is now expecting either *arr* or *extreme*, and any other term will raise an error. You might be wondering... wouldn't it be better to have these parameters in the simulations list. That way it is easy to test the different options. But the ensemble config file is called from the sims list. So, if a few ensemble files are created, it is still easy to test a few options from the sims list by pointing at different ensemble config files. 
- On the same theme... added a new keyword for ```"interim_for_ensemble": "both"```. The *both* keyword invokes the simulation of 20 temporal patterns in the AEP changover zone: both the ARR and extreme temporal patterns. 
- This has not yet been written up in the manual.

# 30 August 2024
- Changed the column label in the mcdf file from *lake_volume* to *ADV* for clarity.
- Changed the bounds for the sampling of TP in the AEP changover zone. For example, if in the storm config ```"aep_changeover_to_extreme": [100, 2000]```: Will use ARR for 1 in 100 AEP and sample either ARR, GSDM, GTSMR for the 1 in 2000 AEP. Previously, would have only used GSDM/GTSMR for 1 in 2000. 
- Fixed a bug where there are duplicate elevations for a given volume in the URBS ELS file (i.e. multiple elevations given for zero volume). The ELS database now removes duplicates.  

# 27 July 2026 -- Richard
- Added sub-burst depth tracking to the Monte Carlo method, to support an ensemble AEP-neutrality check on embedded bursts (whether sub-bursts within the sampled storms occur more often than the IFD statistics imply). For each realisation, the maximum rolling-window depth of the main burst pattern (after filtering, if the filter is on) is recorded for each standard sub-duration, together with the same-z catchment-average IFD reference used by the embedded burst filter. These are stored in the mcdf file as ```subburst_<duration>h``` and ```ifd_<duration>h``` columns (mm, climate-adjusted, before preburst prepending).
- Added an optional ```Analyse sub-bursts``` column to the simulation list. If ```yes```, a frequency curve of sub-burst depths is derived by TPT for each sub-duration and compared with the IFD: outputs are a TPT quantile csv and a plot (TPT curve, realisation scatter, and the same-z IFD reference dashed, one branch per storm method) per sub-duration, plus a ```*_subburst_neutrality.csv``` summary of the TPT/IFD depth ratio at the standard AEPs. Ratios above 1.0 indicate sub-bursts occur more often than the IFD implies (evidence for filtering); ratios below 1.0 mean the unfiltered choice is not contradicted (the check is one-sided, as sub-bursts within one parent duration are only one contributor to the sub-duration IFD).
- Refactored the construction of the embedded burst filter curve into a shared ```build_temppat_filter_curve``` method used by the filter, the unfiltered check, and the new sub-burst measurement. Note this also fixes an inconsistency: ```filter_temppat``` computed the per-duration storm method changeover (GSDM/GTSMR for z rarer than 1 in 2000 in *interpolate_depths* mode) but then queried the depths with the unmodified storm method; the corrected behaviour (as already used by ```check_embedded_burst```) is now applied in both places. This can slightly change filter depths for extreme AEP storms in *interpolate_depths* mode.
- Sub-burst tracking is always on during Monte Carlo runs (its cost is negligible); only the analysis step is switched by the sims list.
- The ```Run models``` key in the sims list now accepts ```storms only``` (Monte Carlo only): the full sampling and storm generation runs and the mcdf is written, but the hydrologic model is skipped. Useful for cheap runs of the sub-burst check and for calibrating temporal pattern weights. The flow results analysis is skipped in this mode; the sub-burst analysis still runs if requested.
- The TPT analysis now recognises an optional ```tp_w``` column in the mcdf: when present, the exceedance fractions within each main division are weighted by it (temporal pattern probability weights). Without the column, behaviour is unchanged.
- Added ```util/CalibrateTpWeights.py```: calibrates temporal pattern probability weights so the ensemble satisfies the sub-burst AEP-neutrality condition. Works entirely offline from the mcdf (no model reruns): patterns whose sub-bursts exceed the same-z IFD are progressively down-weighted until the weighted sub-burst frequency curves sit at or below the IFD, with a floor so no pattern is silently excluded. Also outputs the weighted flood quantiles (from the existing runs) as a preview. Applying the weights to the sampling itself is not yet implemented.
- Refactoring note: the neutrality margin computation moved from ```Simulator``` to module functions in ```MCScheme``` (```interpolate_reference_curve```, ```neutrality_margin```) so the calibration script can share them.
- Added an optional ```TP weights``` key to the sims list (monte carlo only): a path to a temporal pattern weights file (as written by ```CalibrateTpWeights.py```) to apply probability-weighted sampling of the temporal patterns. The weights are grouped by storm method (plus frequency bin for ARR point) and compose with the D50 climate shift weightings. The applied weight is recorded in a ```tp_weight``` mcdf column for audit. Cannot be combined with the ```tp``` replicate. Note that the standard TPT analysis is already correct for a weighted-sampling run (the sampled frequencies embody the weights) - do not also attach a ```tp_w``` column to such a run.
