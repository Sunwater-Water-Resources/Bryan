# Downstream storm generation

[Manual contents](../Manual.md) | [Run launcher](ui.md) | [Utilities](utilities.md)

A representative event chosen on the launcher's [Events page](ui.md) is an event over **one**
catchment — the dam catchment that produced it. The regional model downstream runs many JPA
rainfall regions at once, and each of them needs its own rainfall for that same event. This is
what turns the one into the other.

The essential point, and the reason the inputs look the way they do: **a quantity that is
scalar upstream is a vector here.** One AEP becomes one AEP per region, conditional on how rare
the storm was over the dam catchment. For each region, per event:

```
conditional AEP -> region-average areal depth -> within-region spatial pattern -> per-subarea depth
```

The temporal pattern is the one the upstream event drew, rebuilt from the mcdf.

## What the launcher asks for

The **Downstream** page has two tabs, **Storms** and **Results**, under the two things both are
about: the **Event selection** and the **Regional model URBS config**. A line beneath them says
how each dam will start for that kind of event - from the event's ```lake_z``` for a Monte Carlo
selection, at full supply (or the model config's ```ensemble_start```) for an ensemble one.

The **Storms** tab writes the storm files. Its inputs are in the **Storm inputs** card at the top,
above the plan:

| Field | What it is |
| --- | --- |
| **Event selection** | the `<group>_representative_events.json` the Events page saved. It is found by searching the project folder, so it appears in the list on its own, an ensemble one marked ```(ensemble)```. |
| **Duration (h)** | left blank, it is read from the source database's name. Fill it in only when the name does not carry it. Not shown for an ensemble selection, which records each event's duration. |
| **GWL (°C)** | the same, read from `GWL<n>` in the database name. |
| **Downstream storm config** | the JSON described below. This is the file to get right. |
| **Regional model URBS config** | a normal Bryan [model config](config/ModelConfig.md) for the *regional* model — the one the storm files are being written for. |
| **File name suffix** | added to each storm file's name, before the extension. See *Two groups that share a realisation* below. Blank unless typed; the page suggests one when it is needed. |
| **Run URBS on each storm** | off, the storm files and a batch file each are written and URBS is left to be run by hand. On, each storm is run as it is written. |

While the generator runs, a line under the buttons says how long it has been going and the
last thing it wrote - ```wrote 005044_72h_GWL1p3_RFSL.72``` - with **Cancel** beside it. When
it stops the line says how it ended: the number of storm files and how many need a look, or
the error and the log to read. The full output is in ```<selection>.downstream.log``` beside the
selection, written as it goes. That line is the only thing on the page that updates by itself.

It also passes the generator the sims config the events came from, so an event whose embedded
bursts were filtered keeps its filtered pattern (see below). The selection file records which
one that was, so it does not matter which sims config is open. A selection saved before
22 October 2026 records none, and the open one is used: **open the sims config it was saved
under**, since a selection from another sims list has no row there to rebuild from, and a
filtered event in it is refused.

The two config paths are remembered **per project**, keyed on the sims config's resolved
path, so they survive leaving the page and reloading the project. They are kept in
```~/.bryan_ui.json``` beside the other launcher settings. The duration and warming level are
deliberately not kept: blank means *read it from the database name*, and a stale override is
worse than retyping one.

Everything else it needs — which realisation, which database, the sampled temporal pattern, the
losses, the pre-burst — comes out of the selection file and the mcdf behind it. The page shows
the storm filenames it would write before it writes any of them.

The same run is available as a CLI, which is what the page shells out to:

```
python DownstreamStorms.py --config <downstream config> --selection <selection json>
                           --model <regional URBS config> [--storm-config ...]
                           [--climate-config ...] [--duration H] [--gwl G]
                           [--rain-lag H] [--suffix S] [--dry-run] [--project FOLDER]
                           [--sims-config <upstream sims config>]
```

`--dry-run` reports without writing, and is the quickest way to see whether the config resolves.

The selection records each event's database relative to the project folder. The command finds
it from the selection file itself, which is saved inside the project, so it can be run from any
folder. `--project` names the project folder for a selection that has been copied out of it.

### Two groups that share a realisation

A storm file is named for its realisation, duration and warming level. Two groups that share
all three write the same name - the same Monte Carlo run re-routed at two full supply levels
is the usual case - and the second replaces the first: the storm, and the URBS results named
after it. The suffix keeps them apart. For a selection that has been run, the box shows the
suffix that run used, read off the names it wrote, so a rerun keeps its names; for one never run
it starts blank. Otherwise it is only ever what is typed. The page compares the selection with the others saved in the project and, where they
would clash, says so under the plan and suggests one - the words of the group name no clashing
group has: ```_RFSL``` for ```RFSL_Design_GWL1p3``` against ```FSL_Design_GWL1p3```, and for
```CLD_enb_PMF_..._GWL1p3_RFSL``` against ```..._FSL```, or the whole name where no word is its
own. **Use** beside the warning puts it in the box.

The suffix is part of the name of the run, which is how a dam failure project finds the
regional results for an event - so a group keeps the suffix it was first run with.

### Events whose embedded bursts were filtered

The database records what a realisation *sampled* - which temporal pattern, at what rainfall -
and not the series that came of it. Where the upstream run filtered an embedded burst the two
differ: the same burst depth, spread so that no sub-duration exceeds its own IFD depth. A
regional storm written from the sampled pattern carries the event's depth in a sharper burst
than the event had, and the dam catchment nested in the regional model then floods harder than
the event that was chosen. On Callide's E013 events that was a dam inflow peak 9 to 30 % high
and a lake level up to 1.1 m high, with the volume right.

```--sims-config``` names the upstream run's sims config - the one the Events page had open.
With it each event's storm is rebuilt the way [```util/RepresentativeEvents.py```](utilities.md)
rebuilds it for plotting, checked against the depths the run recorded, and the pattern written
is the one the run applied. The launcher passes it for you. Without it, an event recorded as
```Embedded burst filtered``` is refused rather than written with the wrong pattern.

The record of what went into each storm (```..._downstream_storms.csv```) says which was used,
in its ```temporal_pattern``` column.

### Running URBS

Each storm file is written with a **batch file** beside it, carrying the URBS command line
and a ```set``` statement for every dam's starting lake level. That file is always written,
because it is the record of what the run was: a run done by hand a week later is then the
run that was described here.

```--run-urbs``` runs URBS on each storm as it is written. Without it, nothing is executed
and the batch files are left to be run by hand or collected into a batch of your own - which
is usually what you want, because URBS on a regional model is minutes per event.
```--dry-run``` and ```--run-urbs``` together are refused: the first writes no storm for the
second to run.

**```store_tuflow``` must be ```true``` in the model config.** It is what puts
```set URBS_TFLW=TRUE``` in the batch file, and so what makes URBS write the subcatchment
runoff and the flows at every print location to csv - the upstream boundaries this run exists
to produce. The run warns if it is off rather than leaving you to find the missing files.

### The boundaries for a hydraulic model

With ```store_tuflow``` on, URBS writes ```<run>_i.csv```: the flow at every print location
whose ```PRINT``` statement in the vec carries an asterisk after the name. These are a
hydraulic model's upstream boundaries, and ```<run>_r.csv``` beside it, the routed runoff of
each subarea, is its inflows along the reaches. One of the print locations is the dam, and
its column is the regional model's own routing of the dam.

Where the dam is routed in a dam failure model instead, that column is replaced by the dam
failure project, not by Bryan: it runs after Bryan, finds each event's regional run by its
name, and writes the boundaries itself. At Callide that is ```tuflow_handoff.py``` in
```callide-judyth```. **Keep the suffix the same from one run of the page to the next**, since
it is part of the run name that project looks for.

Bryan used to do this itself, from the page's lower section and ```DownstreamBoundaries.py```.
It was taken out on 22 October 2026 so that one tool writes the boundaries, not two.

### More than one dam

A regional model can route through several dams, and each needs its own antecedent storage.
Everything about a dam is read off its ```DAM ROUTE``` line in the vec - the full supply
level, the storage curve, and the ```il=``` token naming the variable URBS takes the starting
level from - so the [model config](config/ModelConfig.md) only names the extra dams and their
lake configs:

```json
"additional_dams": [
  {"location": "KROOMBIT", "lake_config": "model/kroombit_lake_config.json"}
]
```

**Each ```DAM ROUTE``` line needs its own ```il=``` name.** Two dams sharing one means a
single ```set``` statement drives both, and the second starts every event at the first one's
level - the model runs and nothing in the output says so. Bryan refuses that arrangement by
name. The regional Callide model shipped with ```il=initial_lake_level``` on both lines.

**Every dam takes its antecedent storage from the same ```lake_z```.** A representative event
carries one - the standardised variate the upstream run sampled its storage at - and each dam
maps it through *its own* volume exceedance curve, so one regional wetness gives a different
volume, and a different level, at each dam. That is an assumption: the defensible treatment is
a joint probability analysis of the two storages. The main dam is the exception - where the
realisation recorded an ```ADV```, that is the storage it actually had, and it is used in
preference to re-reading the curve.

A lake config whose ceiling does not match the dam's full supply volume is refused. A curve
fitted to another storage returns a volume on the wrong scale, and the level that comes back
is just somewhere else on the ```.els``` - plausible, and wrong.

### Ensemble events: the PMF

A selection saved from an **ensemble** group on the Events page (```"method": "ensemble"```)
goes through the same page and the same generator, with three differences.

- **The duration is the event's own**, recorded in the selection: an ensemble database holds
  every duration, so its name says none.
- **The starting storage.** An ensemble event has no ```lake_z``` - the run held one starting
  volume for every event. The main dam starts from that ```ADV```; every other dam starts at
  its ```ensemble_start```, which is full supply unless the model config gives a volume:

  ```json
  "additional_dams": [
    {"location": "KROOMBIT", "lake_config": "model/kroombit_lake_config.json",
     "ensemble_start": "fsl"}
  ]
  ```

- **The run name** starts with ```e``` and five digits - ```e00042_9h_GWL1p3_RFSL``` - because
  an ensemble event is a row of its own database, not a realisation, and the two must never
  share a name.

**The concurrent rainfall is set by the event's rainfall AEP**, as for any event - for the PMF
that is the AEP of the PMP, 1 in 1.9 million at Callide, the top of JPA's conditional table. The
PMF's notional AEP (1 in 8 million at Callide) is a label, written to the log and the record,
and is never used for depths: it is the AEP of the flood, which includes the full starting
storage and the choice of the highest pattern, not of the rain.

**An ensemble event whose embedded bursts were filtered is refused for now.** The pattern
rebuild replays a Monte Carlo row's storm, and an ensemble row's sampled pattern would carry
the same depth in a sharper burst than the event had. The PMF ensembles at Callide have no
filtered events.

### Durations the downstream inputs do not cover

JPA's conditional AEPs and areal IFDs are tabulated for a set of durations - **6 to 120 h at
Callide**. An event of any other length cannot be taken downstream: the run reports it as
**SKIPPED** with its duration, carries on with the others, and says how many it skipped at the
end. At Callide that is any 4.5 h event - the PMF at GWL 2.7 among
them - until the tables are extended. **Check without writing** reports it too.

## The Results tab

The **Results** tab reads the runs back - when it is opened, on **Reload**, and once when a run
ends; it never polls. It finds each run by the name the last run wrote, from the selection's
```_downstream_storms.csv```, whatever is in the suffix box now, so the results are there when
the page is opened again later. The page's address names the selection and the tab, so a
bookmark or a refresh comes back to the same view. One row per storm the plan names, with its status and what went into it:
the losses, the mean depth, the pattern, each dam's starting storage and any flags.

| Status | Means |
| --- | --- |
| **finished** | URBS's results are in the results folder. |
| **running** | the generator this page started is running URBS on this event now. |
| **waiting** | the generator has not reached this event yet. |
| **failed** | URBS was asked to run it and stopped without results. Its batch file is left to run by hand; ```urbserr.log``` in the results folder says why. |
| **not run** | the storm and its batch file are written, and URBS was not asked to run it. |
| **not finished** | the launcher was restarted during the run, so it cannot tell running from failed. Reload once URBS has finished. |
| **not written** | no storm file yet. |

Nothing here waits on a clock: a run is *running* because the page's own generator is still
going, and *failed* because it has stopped and left the batch file behind.

### The dams against their levels

Where the regional model has more than one dam, or its config gives ```dam_levels```, the tab
holds each dam's peak lake level in every finished run against reference levels:

```json
"dam_levels": {
  "_note": "Where each level comes from.",
  "CALLIDE":  {"DCL": 219.13, "PMF": {"1.3": 221.50, "1.7": 221.60, "2.7": 221.86}},
  "KROOMBIT": {"DCL": 270.70}
}
```

Each key is the ```location=``` on a dam's ```DAM ROUTE``` line. A level is a number, or a number
per warming level - the PMF moves with the GWL and one regional config serves several - read at
the run's own GWL. Any name can be given; **```DCL```**, the dam crest level, is the one the
chart and the scatter are drawn against. Bryan refuses a location with no ```DAM ROUTE``` line,
or a level that is not a number, when it opens the model.

- **The table** gives, per event and dam, the starting level, the peak, its time and the margin
  to each level - the peak less the level, positive above it. A level the config does not give,
  or not for that GWL, says **not given**. Above it: how many events take every dam to its crest.
  **Copy for Word** copies it.
- **The chart** draws one event's lake levels as **height above each dam's crest**, so dams tens
  of metres apart in AHD share one axis and the crest is zero; a dam's PMF level is a dotted line
  on the same basis.
- **The scatter** takes every finished run of every selection in the project, coloured by group:
  the first dam's margin to its crest across, the second's up. The quadrant above and right of
  zero holds the runs that took both dams over their crests.

Read the dams section with care, as the tab says: each loading is one event, so this is what the chosen events
did and not a probability that the dams reach their crests together; the other dams' starting
storage is an assumption (```lake_z``` through each dam's curve for a Monte Carlo event, full
supply for an ensemble one); and the rain on the other catchments is conditional on the event
over the main dam.

### The print locations

Below the dams, the flow at the print locations marked ```*``` in the vec - the columns of
```<run>_i.csv```, the hydraulic model's boundaries - with every finished event overlaid. The
dam's outflow (```max_keys``` ```outflow```) is shown first, or the last location where that is
not one of them; any others can be added.

Above them, **the end-of-run check**: a run that stops while a location still carries more than
a tenth of its peak, or peaks in the last twentieth of its period, has stopped before its flood
passed, and is named with the locations. The regional model's ```simulation_periods``` are its
upstream periods plus 48 h, which is provisional; this is the check that confirms them. All 42
Callide runs of 21 October 2026 pass it.

Last, **the peaks**: events down the side and every ```*``` location across, as the peak flow in
m³/s or, switched, the time it came in hours - with **Copy for Word** and **Copy as text** for the
report.

## The downstream storm config

All paths are relative to **this config file**, and are resolved through the same helper the
rest of Bryan uses, so Windows backslashes are fine.

```json
{
  "driver": "CLD",
  "driver_area_km2": 1180.0,
  "pmp_aep_by_region": { "CLD": 1600000, "KRO": 2100000 },
  "storm_method_config": { "aep_changeover_to_extreme": [2000, 1600000] },
  "file_paths": {
    "crosswalk":         "regions/subarea_crosswalk.csv",
    "conditional_aep":   "regions/conditional_aeps.csv",
    "ifd_folder":        "ifd/point",
    "areal_ifd_folder":  "ifd/areal",
    "pmp_scaling":       "pmp/pmp_spatial_scaling.csv",
    "pmp_by_region":     "pmp/pmp_depths_by_region.csv",
    "pmp_anchors":       "pmp/pmp_gev_anchors.csv",
    "arr_datahub_file":  "arr/datahub.txt",
    "storm_config":      "../storm_data/storm_config.json",
    "climate_config":    "../climate_change/climate_config.json"
  },
  "mass_balance": {
    "KRO_lower": { "whole": "KRO", "subtract": "KRO_dam",
                   "whole_aep_from": "KRO", "area_whole": 449.0, "area_subtract": 329.0 }
  }
}
```

| Key | Meaning |
| --- | --- |
| ```driver``` | the region the event was generated over — the dam catchment. |
| ```driver_area_km2``` | **the driver catchment's area, not the regional model's.** The ARR areal patterns are banded by area and the GTSMR pattern file is named for one, so this selects which pattern set is read. Give it the regional model's area and the event's own temporal pattern is quietly replaced by a plausible stranger. |
| ```pmp_aep_by_region``` | the AEP of the PMP for each region, as *1 in X*. |
| ```storm_method_config.aep_changeover_to_extreme``` | the two AEPs the spatial pattern is blended across, exactly as in the [storm config](config/StormConfig.md). |
| ```mass_balance``` | optional; see below. |

### The input files

| Path key | Format |
| --- | --- |
| ```crosswalk``` | one row per regional model subarea. Columns: ```new_id``` (the index, the subarea id in the regional model), ```area_km2```, ```cond_target``` (which region the subarea belongs to), ```arf_area_km2``` (the area the region's published ARF was computed over — see the note below). |
| ```conditional_aep``` | JPA's conditional analysis. Indexed by ```Duration_h``` and ```Driver_AEP_1in```, one column per region, each holding that region's AEP as *1 in X* given the driver's. Interpolated on the standard normal variate, not on 1 in X. |
| ```ifd_folder``` | ```average_ifd_<duration>.csv``` per duration: a row per subarea (```name```) and a column per AEP, named ```<aep>_AEP```. Point depths — they set the *shape* within a region. |
| ```areal_ifd_folder``` | ```ArealIFD_<duration:03.0f>h.csv``` per duration: index of ```1 in X``` labels, one column per region. Region-average areal depths — they set the *level*. |
| ```pmp_scaling``` | spatial scaling factors, indexed by ```ID_Number``` (matching ```new_id```), with a ```GTSMR``` column and either ```GSDM_<duration>``` columns or a single ```GSDM```. |
| ```pmp_by_region``` | PMP depths, indexed by region. |
| ```pmp_anchors``` | the reverse-fitted GEV used above 1 in 2,000. Columns ```Duration_h```, ```Region```, ```shape```, ```scale```, ```location```. |
| ```arr_datahub_file``` | the ARR data hub export, read for the ARF zone. |
| ```storm_config``` | the **upstream** storm config, for the temporal patterns. It has to be the one the event was generated with. |
| ```climate_config``` | optional; only needed if the events carry a warming level. |

### Two things that are easy to get wrong

**The level and the shape come from different files.** The region's depth is read off the
*areal* IFD at the conditional AEP, because that curve already has the region's own ARF in it,
computed on the right area, and is tabulated to 1 in 2,000,000. The *point* tables supply only
the within-region pattern. Building the level from the point tables instead means un-doing an
ARF computed over a different area, and extrapolating our own last two points at the PMP
undershot the published depth by 12 % at 120 h. So ```arf_area_km2``` on the crosswalk records
provenance and drives the *reported* ARF — it is not applied to the level.

**```mass_balance``` is for a region derived as a residual.** Some regions have no marginal
curve and no PMP of their own: the piece of a catchment below a dam inherits its AEP from the
whole purely as bookkeeping, and its depth is what is left once the dam catchment's share is
removed:

```
D_lower = (A_whole * D_whole - A_subtract * D_subtract) / (A_whole - A_subtract)
```

Reading the whole region's curve at the whole region's AEP instead gives a different quantity,
and was 6 % out. Name such a region under ```mass_balance``` with the region it comes out of
(```whole```), the region taken off it (```subtract```), the two areas, and which column its
inherited AEP is read from (```whole_aep_from```). A region absent from this block is read off
its curve in the normal way.

**Which file each of the three names has to appear in is not the same**, and this is the part
that catches people:

| Name | Needs a ```cond_target``` in the crosswalk | Needs a column in ```conditional_aep``` | Needs a column in the areal IFD |
| --- | --- | --- | --- |
| the residual (```KRO_lower```) | yes - it owns subareas | yes - that column *is* its AEP | no - it has no curve, that is why it is a residual |
| the parent (```whole```, ```KRO```) | no | **no** | yes - its curve is what gets read |
| the piece taken off (```subtract```, ```KRO_dam```) | yes | yes | yes |

The parent is a bookkeeping entity: no subarea belongs to it and it never gets an AEP of its
own, so a conditional column for it would have nothing to say. What it needs is the curve, read
at the AEP named by ```whole_aep_from``` - which is normally the residual's own column, since
the residual's AEP is inherited from the parent in the first place. Setting
```whole_aep_from``` to the parent instead is the 6 % error above, written down in the config.

The parent still needs an entry in ```pmp_aep_by_region``` and a row in ```pmp_anchors```,
because reading its curve above 1 in 2,000 goes through the fitted GEV like any other region.

## Outputs

One URBS storm file per chosen event, named
```<realisation>_<duration>h_GWL<gwl><suffix>.<duration>``` with the decimal point written as
`p` (```e<event>_...``` for an ensemble event), plus a CSV beside the selection recording what
went into each: the depths, losses and pattern, and the method, event, loading, run name,
notional AEP and each dam's starting storage. The launcher shows the names first;
if one is reported that never appears, the naming in the page and the naming in the generator
have drifted apart — they are written twice deliberately (the page cannot import scipy) and
```tests/test_downstream_naming.py``` exists to stop exactly that.
