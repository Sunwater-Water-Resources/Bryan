# Example project: Juniper Creek Dam

A complete, self-contained Bryan project for a catchment that does not exist. It is here so
that someone setting up their first model has a working set of config files to copy, and can
run something end to end before pointing Bryan at a real study.

## Real data and invented data

The inputs split in two, and the split matters if you are going to copy anything.

**Real, and shipped as they come.** The generic ARR and BoM reference data that any study in
the region would use unchanged:

| File | What it is |
| --- | --- |
| `storm_data/patterns/ECnorth_Increments.csv` | ARR point temporal patterns, East Coast North |
| `storm_data/patterns/Areal_ECnorth_Increments.csv` | ARR areal temporal patterns, East Coast North |
| `storm_data/patterns/gsdm.pat` | BoM GSDM temporal patterns |
| `storm_data/patterns/gtsmr_coastal_500.pat` | BoM GTSMR coastal patterns, 500 km² bin |
| `climate_change/rainfall_loss_rates_of_change.csv` | ARR climate change loss rate table |
| `climate_change/temporal_pattern_scaling_factors.csv` | ARR D50 shift scaling factors |
| `climate_change/Temporal_Pattern_Weighting/*.csv` | ARR D50 pattern weighting factors |
| the `[LONGARF]` zone in the data hub file | ARR ARF zone, East Coast North |

Two further things are real whatever this folder contains, because Bryan carries them in its
own source rather than reading them from a file: the **ARF coefficients**, looked up by zone
name in `ArealReduction.get_coefficients` (which implements only `East Coast North`,
`Northern Coastal` and `Semi-arid Inland Queensland`), and the **rainfall uplift rates** in
`ClimateChange.get_rainfall_uplift_factor`, 15 %/°C at 1 hour easing to 8 %/°C at 24 hours
and beyond. The a–i coefficients written into the `[LONGARF]` block are never read.

**Invented, and prefixed `SYNTHETIC_`.** Everything specific to a catchment that does not
exist: the IFD depth tables, the PMP depths and spatial scaling, the pre-burst depths and
patterns, the storm losses, the dam curves, the model stubs, and the stored results the
routing demo re-routes. These are plausible in shape and magnitude so Bryan behaves the way
it would on a real study, but none of them came from a rainfall frequency analysis, a PMP
assessment or a dam. **Do not copy them into a real assessment.**

The real files keep their own terms - Bryan's GPL v3 licence covers the code and the invented
data, not these. The ARR patterns and climate change tables are published under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), © Commonwealth of Australia
(Geoscience Australia), and are included here unmodified; the GSDM and GTSMR patterns are
Bureau of Meteorology data, subject to the Bureau's terms of use.

So the example applies genuine ARR and BoM regional inputs to a fictional catchment. That is
the useful combination: the parts you would reuse are correct, and the parts you must replace
are obviously not.

## The catchment

| | |
| --- | --- |
| Dam | Juniper Creek Dam (fictional) |
| Catchment area | 620 km² over 10 subcatchments |
| Full supply level | 210.0 m AHD, 190,000 ML |
| Spillway | ungated ogee, 300 m wide, crest at FSL |
| Dam crest | 218.0 m AHD |
| AEP of the PMP | 1 in 1,600,000 |
| Storm durations | 2, 3, 6, 9, 12, 18, 24, 36, 48, 72 hours |

At 620 km² the catchment falls in the 500 km² bin for both the ARR areal and the GTSMR
temporal patterns, and uses ARR point patterns below 12 hours, so all four pattern sets get
exercised. The size was chosen for that bin: it is the one the real pattern files shipped here
cover.

## What runs, and what does not

Bryan drives URBS or RORB, which are licensed programs it does not ship with. Two of the
three methods do not need them, and those are the ones set up here to run:

| | Runs without URBS/RORB | What it demonstrates |
| --- | --- | --- |
| `run_storms_only` | **yes** | rainfall setup, ARF, storm method sampling, temporal patterns, the embedded burst filter, the sub-burst neutrality check, the PMP embedded burst screen, the mcdf |
| `run_routing` | **yes** | Modified Puls routing, antecedent volume handling, the total probability theorem, quantile tables and frequency curve plots |
| a full Monte Carlo or ensemble run | no | needs a licence and a calibrated model |

Between them the two cover the whole chain except the rainfall-runoff step itself.

```
./run_storms_only.sh        # or run_storms_only.bat on Windows
./run_routing.sh
```

Both write to `results/`. The launchers assume Bryan is the parent folder; override with
`BRYAN_ROOT` and `VENV_PY` if your clone sits elsewhere.

### The storms-only demo

`SimsList_storms.xlsx` holds four rows, three of them included:

1. 24 hour storms at present climate
2. the same with the embedded burst filter switched off (`ebf` exclusion), so the two can be
   compared
3. 24 hour storms at GWL 2.7 °C
4. a 6 hour duration, excluded by default — set `Include` to `yes` to see GSDM patterns
   sampled instead of GTSMR

`Run models` is `storms only`, so the sampling and storm generation run and the mcdf is
written, but no hydrologic model is called. Each log opens with the PMP embedded burst screen
and, for rows 1 and 2, the sub-burst neutrality analysis follows.

### The routing demo

`inflows/` holds a synthetic Monte Carlo results database and matching inflow hydrographs, in
the shape a run with the model already run would have left behind. `SimsList_routing.xlsx`
re-routes them three ways: through the dam as built, with the full supply level lowered by a
metre, and with the antecedent volume resampled from `lake_z` through the lake config. That
is the intended use of the method — the same hydrology assessed against several dam options.

## Layout

```
sims_config_storms.json      main config: names the sims list and the three config filepaths
sims_config_routing.json
SimsList_storms.xlsx         one row per simulation
SimsList_routing.xlsx
storm_data/
  storm_config.json          filepaths to the rainfall and pattern data, changeover zones
  ifd_config.json            durations, PMP, AEP of the PMP, interpolation methods
  SYNTHETIC_arr_datahub.txt  ARF zone (real), losses, pre-burst, baseflow (invented)
  ifd/                       one depth table per duration, plus PMP depths and scaling
  patterns/                  ARR point and areal, GSDM, GTSMR (all real), pre-burst
sim_options/
  mc_config.json             the sampling scheme and the TPT quantile ranges
  enb_config.json            durations and AEPs for the ensemble method
  focal_locations/           subcatchment areas upstream of the dam
  lake_conditions/           the antecedent volume distribution
climate_change/
  climate_config.json        NRM cluster, KG classification, weighting files
  *.csv                      the real ARR climate change tables
model/
  rorb_config.json           RORB template
  urbs_config.json           URBS template
  SYNTHETIC_juniper.els      elevation / area / storage
  SYNTHETIC_juniper.sq       storage / discharge above full supply
  rorb/, urbs/               model file stubs - see below
inflows/                     stored results and hydrographs for the routing demo
make_example_data.py         regenerates the invented data only; leaves the real files alone
make_example_configs.py      then writes the configs, sims lists and launchers
```

The `.catg`, `.par` and `.vec` files under `model/` are **stubs, not working models**. They
carry just enough structure for Bryan to start — RORB reads the `.catg` at construction to
find the dam routing line, URBS reads the `.vec` — so the storm generation side can be
exercised without a licence. Neither model will run against them.

## Things worth knowing before you copy the configs

- **Paths in every config are relative to the file that names them.** The main config's
  filepaths resolve against itself, the storm config's `file_paths` against the storm config,
  and the IFD config's duration filenames against the IFD config. Getting this wrong is the
  most common setup error.
- **The ARF zone in the data hub file must be one Bryan implements.** `ArealReduction` has
  coefficients for `East Coast North`, `Northern Coastal` and `Semi-arid Inland Queensland`
  only, and it raises if the `[LONGARF]` `Zone` is anything else. The a–i coefficients in the
  file are not used — Bryan carries its own for the named zone.
- **Subcatchment order matters.** The rows of the IFD tables and the PMP scaling file must be
  in the same order as the subcatchments in the URBS or RORB model.
- **`volume_cap` in the lake config takes a keyword**, one of `fsv`, `none` or `ceiling`, not
  a volume.
- **The routing method checks the results database against the sampling scheme** in the config
  the `Config file` column points at, and refuses to analyse if they disagree. Point it at the
  Monte Carlo config the input mcdf was actually generated with.
- **The temporal pattern files are order-sensitive.** All patterns for a given duration (and,
  for ARR point patterns, a given frequency bin) must be contiguous rows, because the readers
  slice by row position rather than filtering.

## Regenerating

```
python make_example_data.py
python make_example_configs.py
```

The generators are deterministic — the seed is fixed — so re-running them reproduces the same
files. Edit the constants at the top of `make_example_data.py` to build a differently shaped
example: change `SUBCATCHMENT_AREAS` and the area bins follow, change `DURATIONS` and every
depth table follows. If you resize past a pattern area bin boundary you will need the real
ARR areal and GTSMR files for the new bin, which is why the catchment is 620 km².

Neither generator touches the real ARR and BoM files, so re-running them cannot quietly
replace real reference data with invented numbers.
