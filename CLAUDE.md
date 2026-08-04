# CLAUDE.md

Guidance for working in this repository.

## What this is

**Bryan** is Sunwater's Python platform for design flood hydrology simulation. It drives
external hydrologic models (**URBS** and **RORB**) through either a **Monte Carlo** or
**Ensemble** scheme to estimate design flood quantiles (peak inflow, lake level, outflow),
then post-processes results using the Total Probability Theorem (TPT). It implements
Sunwater's design flood hydrology specification.

Developers: Richard Sharpe and Graigan Panosot.

Units convention throughout: AEP is expressed as **"1 in X"**, storm durations in **hours**,
catchment area in **km²**, volumes in **ML**.

## Running

Bryan is run from a **batch file** (Windows) that passes a config file path to a top-level
script. The Python code lives separately from the model/project data files, so one copy of
Bryan serves many dam catchments.

- `python Main.py <sims_config.json>` — main entry point. Reads the config, opens the
  simulation list (Excel), and dispatches each row to a simulator by its `Method`:
  `monte carlo`, `ensemble`, or `reservoir routing`.
- `python RouteFlows.py <routing_sheet.xlsx>` — standalone: adds baseflow to quickflow
  hydrographs (for RORB) and/or routes flows through the dam. Driven by an Excel sheet.
- `IFD_export.py`, `DownstreamStormGenerator.py`, `StormInstance.py` — top-level helper
  scripts with hard-coded paths in their `main()`/module header (legacy/ad-hoc tooling for
  generating IFD tables, downstream storm files, and single storm instances). Not part of
  the main batch workflow.

### Config / control flow

```
batch file --> sims_config.json --> simulation list (Excel) --> per-row simulators
```

- **sims_config.json** (main config): names the `simulation_list` Excel file and holds
  `filepaths` to the `model_config`, `storm_config`, and `climate_config` JSON files. All
  filepaths are **relative to the batch/config file location**. Optional keys: `project_folder`,
  `test_runs` (limits run count for testing; `0` = run to completion).
- **simulation list** (Excel, sheet 0): one row per simulation. Only rows with `Include == 'yes'`
  run. Key columns include `Method`, `Duration`, `Run models`, `Analyse results`,
  `Store hydrographs`, `Mop up files`, `ADV`, `Baseflow`, `Focal subcatchments`, `Config file`
  (the MC/ensemble method config), `GWL`/`SSP`/`Year` (climate), `Replicates`, `Exclusions`.
  See `Manual/SubDocs/sim_list.md` for the full field reference and the replication/exclusion key tables.
- Per-method config files (Monte Carlo, Ensemble, IFD, Storm, Lake, Model, Climate) — documented
  under `Manual/SubDocs/`.

## Architecture

All core logic lives in `lib/`. The top-level scripts are thin dispatchers.

### Simulators (`lib/Simulator.py`)
- `Simulator` — base class: reads config, sets up logging, model, analysis flags.
- `MonteCarloSimulator(Simulator)` — stratified TPT sampling across many realizations.
- `EnsembleSimulator(Simulator)` — runs a grid of AEPs × durations from the ensemble config.
- `lib/ReservoirRouting.py` — `ReservoirRoutingSimulator`, `FastTPT` (the third method).
- `lib/EnbAnalysis.py` — `analyse_ensemble(df, outputfile)`: median pattern per duration,
  box plots, critical duration per AEP. Shared by `EnsembleSimulator.analyse_results` and
  the reservoir routing method, so the two cannot drift apart.
- `Logger` — tees stdout to a log file.

### Storm generation (`lib/StormGenerator.py`)
- `StormBurst` — central object assembling a design storm: imports rare + extreme rainfall,
  applies areal reduction, loads temporal patterns, computes catchment-average rain, filters
  embedded bursts, and prepends preburst patterns.
- `PreBursts`, `ArealReduction` — supporting components.

### Rainfall (`lib/RainfallScheme.py`)
- `ifdCurves`, `DurationCurve`, `ExtremeDurationCurve`, `PMP`, `WeightInterpolatedSpatialPattern`.

### Temporal patterns (`lib/TemporalPatterns.py`)
- `TemporalPatterns` base with `GtsmrPatterns`, `GsdmPatterns`, `PointPatterns` (ARR point),
  `ArealPatterns` (ARR areal) subclasses, plus `PreburstPatterns`.

### Sampling (`lib/MCScheme.py`)
- `SampleScheme` — Monte Carlo sampling of rainfall (truncated normal in standard-normal-variate
  space), temporal patterns (int 0–9), storm method, losses, preburst percentile, lake level.
- `TotalProbTheorem` — TPT result analysis.
- `lib/EnbScheme.py` — `Ensemble` scheme.

### Hydrologic models (external executables wrapped by Bryan)
- `lib/URBSmodel.py` — `UrbsModel`: writes storm files, runs URBS, parses results. Supports
  volume-based and level-based dam routing (auto-detected from the `.vec` file).
- `reservoir routing` accepts the output of **either** scheme, detected from the input database's
  columns (`_detect_scheme`): `m`/`n` means Monte Carlo, `duration`/`tp` means ensemble. That
  choice drives both the ADV and the re-analysis (TPT vs `analyse_ensemble`), so there is no
  sims-list column to get out of step with the data.
- Monte Carlo input takes the antecedent dam volume from the `ADV` column of the input mcdf by
  default; the optional `ADV source` sims-list column (`lake_z` / `lake_z correlated`) instead
  resamples it from the mcdf `lake_z` column via the lake config distribution, so one set of
  inflows can be re-routed under different antecedent storage distributions.
- Ensemble input holds one starting volume for the whole run, so the ADV comes from the sims-list
  `ADV` column via `LakeConditions` — a number, `fsv`, `mav`, or `database`. **Use `fsv` or `mav`
  when re-routing under a different dam**: both resolve against the curve being routed, whereas
  `database` (and the MC-style default) would start every event at the *old* dam's full supply
  volume, which fails silently. `varying` is rejected.
- Both the inflow hydrographs and the input database for `reservoir routing` may be `.csv` or
  `.parquet`; both go through `_read_indexed`. Parquet written with `index=False` carries the
  index as an ordinary column, so it is promoted on read — do not assume a positional
  `index_col` works for both formats. Parquet also returns `storm_method` as a `Categorical`,
  which will not concatenate with a string — `analyse_ensemble` casts it back.
- An ensemble inflow file is **ragged**: one file spans every duration, each padded with NaN
  after its own simulation period. This works without special handling — `interpolate('slinear')`
  does not fill trailing NaN, the routing propagates them, and the `np.nanmax` in `_write_mcdf`
  takes the peak off the finite prefix. Do not "fix" the NaNs by forward-filling: holding the last
  inflow value out to the end of the longest duration fills the dam and corrupts every short storm.
- `lib/RORBmodel.py` — `RorbModel`: equivalent wrapper for RORB.

### Routing & baseflow (`lib/Routing.py`)
- `Router`, `Baseflow`, `DamRouter`.

### Lake / antecedent conditions (`lib/Lake.py`)
- `LakeConditions`, `StorageCurve`, `VolumeExceedanceCurve`, `ExceedanceCurveLayer`, `Correlator`.
- Sims-list `ADV` keywords: a number (ML), `fsv`, `mav`, `varying`. `mav` reads the exceedance
  curve at `z = 0` — the median of what `varying` samples — so it needs a `Lake config` too.
  `fsv` and `mav` are both resolved in `set_full_supply_volume`, not the constructor, because
  the fsv cap has to be applied against the dam actually being run.
- The `sigmoid` ADV curve is `log10 V(z) = A / (H + exp(-k(z - z0))) + log10(Vf)`. **`A` has two
  conventions and `H` is what tells them apart** — see `_resolve_sigmoid` and
  `Manual/SubDocs/config/LakeConfig.md`. `H > 1` is an old hand-tuned workbook curve:
  `A = log10(Vc)`, and `Vc` is `10**A` rather than the ceiling. `H == 1` (or absent) is a fitted
  standard logistic: `A = log10(Vc) - log10(Vf)`, and `Vc` is the exact ceiling. Bryan derives
  `A` rather than requiring it; an explicit `A` coefficient overrides the rule but is rarely
  needed. Old models are untouched because their `H` is above 1, and `H == 1` can never be a
  deliberate old-style curve (it would asymptote to `Vc * Vf`, needing a 1 ML floor to mean
  anything). Getting it wrong fails *silently*: at Callide the ceiling came out at 1.5e9 ML,
  every realisation started above the top of the `.els`, and only the level frequency curve was
  blank — the inflow curve still looked fine. Bryan prints the resolved `A`, its convention and
  both asymptotes on load, and warns when the ceiling lands far above `Vc`.

### Climate change (`lib/ClimateChange.py`)
- `ClimateAdjustment` — rainfall/loss uplift and temporal-pattern shift per the 2023 draft ARR
  Climate Change Considerations update (GWL- or SSP/Year-based). `D50Weighting` — front-loading shift.

### Curve fitting (`lib/InterpolationCurves.py`)
- `Curve`, `CoercedQuadratic`, `GEV` — used to extrapolate rainfall to rare/extreme AEPs.

### Post-processing utilities (`util/`)
Standalone scripts with editable paths at the top of `main()`, e.g. `PlotFrequencyCurves.py`
(frequency plots), `GetRepresentativeEvents.py` (representative event selection),
`DesignFloodInterpolation.py`, `MaxQuantiles.py`, `ReportCollation.py`. See
`Manual/SubDocs/utilities.md`.

## Environment

- Python 3.12. Dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`, `pyarrow`.
- `pyarrow` is needed only to read `.parquet` inflow hydrographs in `ReservoirRouting.py`
  (`_read_inflows`). Everything else works without it.
- `requirements.txt` — pip pins (note: newer/looser than the conda env).
- `_env_bryan.yml` — conda environment `bryan29` (the as-used Windows environment;
  numpy 1.26, pandas 2.2, scipy 1.13). Prefer this for reproducing study results.
  **It has no `pyarrow`** — it is a `conda env export` with exact build strings, so add the
  package to the env and re-export rather than hand-editing the file.
- Designed for **Windows** (batch files, `COMPUTERNAME` env var, backslash paths in configs,
  external URBS/RORB `.exe`). Running the full pipeline on Linux requires the model executables
  and will hit path/env assumptions.

## Conventions & gotchas

- **Paths in configs are relative to the batch/config file**, converted to absolute at load time.
- `.gitignore` excludes all data files (`*.csv`, `*.xlsx`, `*.nc`, `*.tif`, etc.) and
  `outputs/`, `results/`, `figures/`, `plots/`. Models and study data live outside the repo.
- Several helper scripts (`StormInstance.py`, `DownstreamStormGenerator.py`, `IFD_export.py`,
  most of `util/`) carry **hard-coded absolute Windows paths** in their `main()`. These are
  per-study scratch tools — expect to edit paths before use; don't assume they run as-is.
- Logging is done by redirecting `sys.stdout` to a `Logger` that writes to the per-sim log file;
  expect heavy `print()` usage rather than the `logging` module. All three methods write one,
  reservoir routing included — it defaults to `<Results folder>/<Output file><suffix>_log.txt`
  when the sims list leaves `Log file` blank, and opens with a `SIMULATION INPUTS` block (the
  whole sims row, plus each input file's absolute path, size and mtime) for QA.
- **The `Logger` owns its file and Main.py closes it** in the per-simulation teardown. Do not go
  back to letting garbage collection do it: an unclosed log still holds buffered output, and the
  sims lists give several rows the same `Log file` (one per GWL group, e.g. 92 Callide rows over
  7 paths), so whichever buffer flushed last decided what the file ended up containing. That is
  what made a later run's log come out holding an earlier run's output. `lib/LogFiles.py`
  separately renames duplicate log paths within a batch so the runs cannot collide at all.
  `Logger.flush()` must stay a real flush — it was a no-op, which silently lost the tail of any
  run that crashed.
- `lib/ConsoleTitle.py` puts the running simulation in the console window title (`Bryan 3/21:
  <output name> - <sims list>`). Windows uses `SetConsoleTitleW`, which works even when a batch
  file redirects the output; elsewhere it writes the xterm OSC sequence, and only when
  `sys.__stdout__` is a tty. Write to `sys.__stdout__`, never `sys.stdout` — the latter is the
  `Logger`, and the escape codes would end up in every simulation log. Best-effort throughout:
  a title is never worth failing a run over.
- The authoritative technical reference is `Manual/Bryan_Technical_Reference_v1.pdf` (and `.docx`).
  `Manual/Manual.md` + `Manual/SubDocs/` are the user guide. `Manual/change_log.md` records
  design decisions and config-format changes — read it when config keys seem inconsistent.
