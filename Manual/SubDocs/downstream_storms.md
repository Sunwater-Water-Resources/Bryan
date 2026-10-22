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

The **Downstream** page takes these to write the storm files:

| Field | What it is |
| --- | --- |
| **Event selection** | the `<group>_representative_events.json` the Events page saved. It is found by searching the project folder, so it appears in the list on its own. |
| **Duration (h)** | left blank, it is read from the source database's name. Fill it in only when the name does not carry it. |
| **GWL (°C)** | the same, read from `GWL<n>` in the database name. |
| **Downstream storm config** | the JSON described below. This is the file to get right. |
| **Regional model URBS config** | a normal Bryan [model config](config/ModelConfig.md) for the *regional* model — the one the storm files are being written for. |
| **File name suffix** | added to each storm file's name, before the extension. See *Two groups that share a realisation* below. Blank unless typed; the page suggests one when it is needed. |
| **Run URBS on each storm** | off, the storm files and a batch file each are written and URBS is left to be run by hand. On, each storm is run as it is written. |

It also passes the open project's sims config to the generator, so an event whose embedded
bursts were filtered keeps its filtered pattern (see below). **Open the sims config the
selection was saved under**: a selection from another sims list has no row here to rebuild
from, and a filtered event in it is refused.

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
after it. The suffix keeps them apart. It starts blank and is only ever what is typed in the
box. The page compares the selection with the others saved in the project and, where they
would clash, says so under the plan and suggests one: the leading word of the group name
(```_RFSL``` for ```RFSL_Design_GWL1p3```), or the whole group name where that word is shared.

The same suffix has to be given when the dam's outflow is put into the inflow files, since
it is part of the name of the run.

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

### The dam's outflow, from a dam failure model

With ```store_tuflow``` on, URBS writes ```<run>_i.csv```: the flow at every print location
whose ```PRINT``` statement in the vec carries an asterisk after the name. These are a
hydraulic model's upstream boundaries, and ```<run>_r.csv``` beside it, the routed runoff of
each subarea, is its inflows along the reaches. One of the print locations is the dam, and
its column is the regional model's own routing of the dam.

Where the dam is routed in a dam failure model instead - Judyth - that column has to come
from there, so that the no-failure boundary and the failure cases are routed by one model.
The page's lower section does that, and is the command

```
python DownstreamBoundaries.py --selection <selection json> --model <regional URBS config>
                               --dam-outflow <dam failure simulation list .xlsx> --out <folder>
                               [--suffix S] [--gate-scenario NAME] [--name-contains TEXT]
                               [--duration H] [--gwl G]
```

| Field | What it is |
| --- | --- |
| **Dam failure simulation list** | the dam failure model's simulation list. It is how a run is matched to an event. |
| **Write the inflow files to** | a folder of its own. URBS's files are left as they are. |
| **Gate scenario** | the ```Gate_scenario``` of the dam runs to use. Blank is the gates working. |
| **Only simulations whose name contains** | narrows the list where it holds the same event more than once. |

For each chosen event a copy of ```<run>_i.csv``` is written with the dam's column - the
location named as ```max_keys``` ```outflow``` in the model config - replaced by the dam
model's outflow at the same times, and a ```<group>_boundaries.csv``` records which run went
into which file. With a gate scenario the files are named ```<run>_<scenario>_i.csv```.

The run for an event is the row of the simulation list whose ```Inflow_id``` ends in the
realisation (```sim08029```), whose ```Bryan_group``` is the selection's group, and whose
```Gate_scenario``` is the one asked for. Exactly one row has to match; none or several is
reported with what was looked for. ```Bryan_group``` is a column the dam failure project adds
to its list for this purpose - without it, two groups sharing a realisation cannot be told
apart, and the name filter has to do it.

**The dam run has to reach the end of the regional run.** The regional model runs on after
the upstream one stops, for the flood to pass down the valley, and the dam is still
releasing over that time. A dam run that stops early is refused rather than padded.

Only the dam's column changes. A print location on the main stem below the dam still holds
URBS's total, which includes URBS's routing of the dam; it is not a boundary for a hydraulic
model that starts at the dam.

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
`p`, plus a CSV of the applied depths beside the selection. The launcher shows the names first;
if one is reported that never appears, the naming in the page and the naming in the generator
have drifted apart — they are written twice deliberately (the page cannot import scipy) and
```tests/test_downstream_naming.py``` exists to stop exactly that.
