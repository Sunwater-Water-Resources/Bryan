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
- `python MainMulti.py <sims_config.json>` — same as `Main.py` but splits the simulation
  list across processes (`multiprocessing` key in config). Does **not** support the
  `reservoir routing` method.
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
  `multiprocessing`, `test_runs` (limits run count for testing; `0` = run to completion).
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
- `lib/RORBmodel.py` — `RorbModel`: equivalent wrapper for RORB.

### Routing & baseflow (`lib/Routing.py`)
- `Router`, `Baseflow`, `DamRouter`.

### Lake / antecedent conditions (`lib/Lake.py`)
- `LakeConditions`, `StorageCurve`, `VolumeExceedanceCurve`, `ExceedanceCurveLayer`, `Correlator`.

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

- Python 3.12. Dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`.
- `requirements.txt` — pip pins (note: newer/looser than the conda env).
- `_env_bryan.yml` — conda environment `bryan29` (the as-used Windows environment;
  numpy 1.26, pandas 2.2, scipy 1.13). Prefer this for reproducing study results.
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
  expect heavy `print()` usage rather than the `logging` module.
- The authoritative technical reference is `Manual/Bryan_Technical_Reference_v1.pdf` (and `.docx`).
  `Manual/Manual.md` + `Manual/SubDocs/` are the user guide. `Manual/change_log.md` records
  design decisions and config-format changes — read it when config keys seem inconsistent.
