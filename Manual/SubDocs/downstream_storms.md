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

The **Downstream** page takes five things:

| Field | What it is |
| --- | --- |
| **Event selection** | the `<group>_representative_events.json` the Events page saved. It is found by searching the project folder, so it appears in the list on its own. |
| **Duration (h)** | left blank, it is read from the source database's name. Fill it in only when the name does not carry it. |
| **GWL (°C)** | the same, read from `GWL<n>` in the database name. |
| **Downstream storm config** | the JSON described below. This is the file to get right. |
| **Regional model URBS config** | a normal Bryan [model config](config/ModelConfig.md) for the *regional* model — the one the storm files are being written for. |

Everything else it needs — which realisation, which database, the sampled temporal pattern, the
losses, the pre-burst — comes out of the selection file and the mcdf behind it. The page shows
the storm filenames it would write before it writes any of them.

The same run is available as a CLI, which is what the page shells out to:

```
python DownstreamStorms.py --config <downstream config> --selection <selection json>
                           --model <regional URBS config> [--storm-config ...]
                           [--climate-config ...] [--duration H] [--gwl G]
                           [--rain-lag H] [--suffix S] [--dry-run]
```

`--dry-run` reports without writing, and is the quickest way to see whether the config resolves.

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
