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
- `python ui/main.py [<sims_config.json>]` — the run launcher (browser UI). Picks rows out
  of a sims list, checks them, and drives `Main.py` for you. Its own environment
  (`ui/requirements-ui.txt`); see the `ui/` section below.
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
- **Temporal pattern weights are applied in the sampling *or* in the analysis, never both**, and
  the column name is what says which. The sims-list `TP weights` column composes the weights into
  `get_temporal_pattern_sample`, so the realisations are *drawn* in proportion to them and the TPT
  then counts them plainly — a weighted sample already is the weighted probability. That run
  records what it used as **`tp_weight`**: provenance, never an input. `util/CalibrateTpWeights.py`
  works the other way, attaching **`tp_w`** to a *uniformly* sampled mcdf so a weight set can be
  tried without re-running. Renaming `tp_weight` to `tp_w` (it looks like an obvious typo, and it
  is not) would weight the patterns twice and every curve would still look plausible.
  `compute_std_quantiles` and `attach_weights` now refuse an mcdf carrying both.
- `lib/EnbScheme.py` — `Ensemble` scheme.

### Hydrologic models (external executables wrapped by Bryan)
- `lib/URBSmodel.py` — `UrbsModel`: writes storm files, runs URBS, parses results. Supports
  volume-based and level-based dam routing (auto-detected from the `.vec` file).
- `reservoir routing` accepts the output of **either** scheme, detected from the input database's
  columns (`_detect_scheme`): `m`/`n` means Monte Carlo, `duration`/`tp` means ensemble. That
  choice drives both the ADV and the re-analysis (TPT vs `analyse_ensemble`), so there is no
  sims-list column to get out of step with the data.
- **Every reservoir routing output ends with `Output suffix`**, the results database included:
  `<Output file>__mcdf<suffix>.csv` (Monte Carlo input) or `<Output file><suffix>.csv`
  (ensemble). One sims list re-routes one set of inflows under several rating curves into one
  results folder, so a name that drops the suffix is a silent overwrite — which the mcdf was
  until 26 August 2026. `_ensure_mcdf_loaded` falls back to the old unsuffixed name, warning
  that every suffix over that `Output file` wrote it; do not remove that fallback without
  re-routing the studies that depend on it.
- Monte Carlo input takes the antecedent dam volume from the `ADV` column of the input mcdf by
  default; the optional `ADV source` sims-list column (`lake_z` / `lake_z correlated`) instead
  resamples it from the mcdf `lake_z` column via the lake config distribution, so one set of
  inflows can be re-routed under different antecedent storage distributions. `ADV source =
  sims list` holds every realisation at one volume from the sims-list `ADV` column instead
  (`_sims_list_adv`, shared with the ensemble path), for testing dam operation against a
  nominated antecedent storage — the quantiles are then conditional on it, not design flood
  quantiles, so the run says so in the log. **The sims-list `ADV` is not read at all under the
  default `mcdf` source**, which is what made an `ADV` of `fsv` there look like it was being
  ignored; it now prints a note saying the key was never used and naming the source that would
  use it.
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
- `Analyse volumes` (`yes`/`inflow`) adds the inflow volume analysis: the peak volume in a moving
  window of each duration from `volume_durations` in the row's `Config file`, through the same
  analysis as the peaks (TPT or `analyse_ensemble`). Inflow only, deliberately — outflow and
  storage volumes follow from the peak level and the rating curve. The windows and the reading of
  the sims-list column both come from `lib/Volumes.py`, shared with `Simulator.analyse_volumes`.
  Output naming follows the Monte Carlo method's, because `util/MaxQuantiles.py` and
  `util/PlotFrequencyCurves.py` read both: the **files** are tagged `inflowVol24h`, the **column**
  inside them is `Vol24h`. Neither is free to change.
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

### Flood volumes (`lib/Volumes.py`)
- `rolling_max_volumes(flows, dt, durations)` — the peak volume (ML) in a moving window of each
  analysis duration, one pass per duration over the whole array. Shared by
  `Simulator.analyse_volumes` (Monte Carlo and ensemble) and `ReservoirRoutingSimulator`, which
  is the point: it was written twice, and the two windows differed. The window is
  `duration/dt + 1` samples so it spans the **full** duration — the Monte Carlo copy used
  `duration/dt`, a timestep short, and biased every volume low by about `dt/duration` (fixed
  21 August 2026, so volumes from earlier runs are slightly lower than they should be). Durations
  that are not a whole number of timesteps, or longer than the hydrographs, are reported and
  skipped; both used to be silent, and the first crashed the caller on a column-count mismatch.

### Representative events (`lib/RepresentativeEvents.py`)
- The analysis behind the launcher's Events page: rank the realisations of an mcdf against a
  design loading (an AEP, or a lake level read off the level frequency curve), and report what
  would make each candidate indefensible — an embedded burst, a pre-burst or an antecedent
  storage far off the median.
- **Ranking is a distance in standard normal variate space**, on both the result AEP and the
  rainfall AEP at once, so an event is judged on AEP neutrality as well as on reaching the
  loading. `1 in X` is not a linear scale — at a 1 in 2,000 target being 200 out is nothing and
  at 1 in 100 it is everything — which is why `util/GetRepresentativeEvents.py`'s
  `sqrt(d_rain² + d_result²)` on the raw AEPs is kept only as `delta_aep`, for comparison.
- **The units in an mcdf are mixed and nothing warns you**: `rain_aep` is **1 in X**
  (`Simulator.py:1012`) while `level_aep`/`inflow_aep`/`outflow_aep` are **probabilities**
  (`TotalProbTheorem.assign_aep`). `prepare` converts once, into `result_aep`; a test pins it.
- **pandas and the standard library only — no scipy, no matplotlib.** The UI imports this
  module directly, and `statistics.NormalDist` covers `ndtri`/`ndtr`. Anything needing scipy
  belongs in the util script that plots the chosen events.
- `prepare` (per database and result type) is split from `score` (per loading) because the
  first costs a z per row and the second is arithmetic; the UI caches the first.
- Embedded bursts are asked of the **data**, never of the flag text: `has_embedded_burst`
  reads the `embedded_bursts` comment and the `subburst_<d>h`/`ifd_<d>h` ratio. A substring
  match on the flags dropped every high *pre-burst* event, since "pre-burst" contains "burst".
- Everything except the PMP cap and the distance limit **flags** rather than excludes. A tool
  that silently drops the event someone was looking for does not get trusted twice.

### Rebuilding one storm (`lib/EventStorm.py`)
- Replays the storm generation for a single realisation so a chosen representative event can be
  given a **hyetograph** — an mcdf records what was *sampled*, never the rainfall series.
  Deterministic, because every draw is in the mcdf row.
- Written from `Simulator.run_models`, not from `StormInstance.py`, which does the same job for
  a different purpose and has drifted: it omits the `buffer=0.9` on the pre-burst filter. The
  buffers are **not the same** — 1.1 for the main burst, 0.9 for the pre-burst — and getting one
  wrong changes the pattern without changing its total, so no check on depths would catch it.
  Both are pinned by tests.
- **Every rebuild is checked against what the run wrote down**: `mean_rain_mm`, `preburst_mm` and
  the `embedded_bursts` comment. The three are independent, so agreement with all of them means
  it is the storm that was modelled. A rebuild that disagrees is reported, never quietly plotted.
- Temporal patterns are **percentages of the burst depth** throughout (`uniform_preburst`:
  `preburst_depth = preburst_proportion * 100`), so depths are `pattern / 100 * ave_rain`. The
  index is hours with the pre-burst at negative times: t = 0 is the start of the main burst,
  which is also what the stored hydrographs are shifted onto.
- Needs scipy and matplotlib, so it is **not** UI-importable and does not try to be — the split
  from `lib/RepresentativeEvents.py` is the whole reason that one stays cheap.

### Curve fitting (`lib/InterpolationCurves.py`)
- `Curve`, `CoercedQuadratic`, `GEV` — used to extrapolate rainfall to rare/extreme AEPs.

### Run launcher (`ui/`)
A NiceGUI browser app for choosing which sims-list rows to run, checking them, launching
Bryan and following it. See `ui/README.md` and `Manual/SubDocs/ui.md`.

- **It never imports Bryan's simulators.** The simulators do all their work in `__init__`
  and reassign `sys.stdout` globally, so the UI runs `Main.py` as a **subprocess** and reads
  files for progress. It writes a filtered copy of the sims list plus a config pointing at
  it — Bryan's entry point is untouched.
- **Three binding rules, pinned by `ui/tests/test_dependency_direction.py`:** nothing in
  `lib/` may reference `ui`; nothing in `ui/core/` may import `nicegui`; and `ui/` imports
  Bryan only through `ui/core/bryan.py`, whose allow-list is `lib.RunLog`, `lib.LogFiles` and
  `lib.RepresentativeEvents`. All three import only pandas, which is what keeps the UI
  environment free of scipy and matplotlib — and that, not the list itself, is the test for
  membership: `test_the_allow_list_stays_cheap` imports every entry and fails if scipy or
  matplotlib arrives with it. Do not widen the allow-list to pull in a simulator.
- **The UI never writes a master sims list.** openpyxl cannot store a formula's cached
  value, so saving a formula-driven workbook from Python turns every formula into an
  uncached one and Bryan then reads blanks. `TFD_SimsList_LongList_02.xlsx` in the Tinaroo
  project is already in that state — 96 formula cells, no values, so its first 24 rows read
  with no `Output file` and no inputs. `ui/core/simslist.audit_formulas` detects it (the
  test must be `value in (None, "")` — that workbook stores an *empty* `<v>`, so an
  `is None` check reports it clean). `runwriter._set` forces every written cell to a string
  type, because openpyxl turns a literal starting with `=` back into an uncached formula.
- **Chunking is constrained by `Output file`.** `Simulator.initialise_model` makes it the URBS working
  sub-folder and `UrbsModel.__init__` rmtree's that folder on entry, so rows sharing one are
  indivisible — and selecting several is blocked outright, since they overwrite each other
  even sequentially (24 rows share one name in `CLD_RFSL_mc_sims_01.xlsx`). Parallelism
  defaults to 1 on purpose.
- **Per-chunk run logs are free**, because `RunLog.log_filepath` builds the path from the
  raw `simulation_list` string against the CWD. The UI writes a relative
  `_ui_runs/<id>/chunk_NN.xlsx` and launches with `cwd=project_folder`, so each chunk's log
  lands in its own run folder. No locking is needed — do not add any.
- **`stdin=DEVNULL` is deliberate.** The `input()` calls in `MCScheme`/`EnbScheme.store_simulations`
  would hang a windowless run forever; with DEVNULL they raise `EOFError`, Main.py catches it
  per row, and `progress.explain_error` says why. There was a third in `URBSmodel`, on a missing
  URBS executable; it now warns and carries on, so that case fails later and only for a
  simulation that actually runs the model.
- Progress matches run-log rows to simulations **positionally**, never by name: `Simulation`
  is `Output file` (`RunLog.py:36`) and carries no duration.
- `ui/tests/test_progress.py` greps `Main.py` and `lib/RunLog.py` for the literal strings it
  parses. If you change that wording, that test fails — update both together rather than
  letting the progress display silently go blank.
- **A reservoir routing row is judged by its suffixed outputs.** Everything the method writes
  takes `Output suffix` last, the Monte Carlo database included —
  `<Output file>__mcdf<suffix>.csv` since 26 August 2026. Before that it was `__mcdf` with no
  suffix, so every rating-curve variant of one `Output file` overwrote one file: Tinaroo's
  `TFD_SimsList_LongList_01.xlsx` is 6 output names x 6 suffixes over one results folder, 36
  rows leaving 6 databases. Those files are still out there and `_ensure_mcdf_loaded` still
  falls back to one (with a warning), so `core/outputs.py` keeps the roles apart: `primary` is
  what proves *this* row ran, `shared` is the legacy unsuffixed mcdf, `databases` is what
  `_ensure_mcdf_loaded` will read back, and only the last of those lifts `needs prior results`
  for a `Run models = no` row. Truncation counts the database, never `primary`: a quantile
  table is one row per standard AEP by design.
- **`preflight.triage` decides what can run without the rest**, for the Select page's
  "Deselect problem rows" and the run dialog's "Skip N and run the rest". It re-checks after
  each drop because clearing one issue can clear another, folds in the planner's blocking
  hazards through a `plan_for` callback (a collision that only appears once the rows are
  chunked would otherwise stop a run the user was just told was fine), and refuses when a
  blocker names no row — a missing column is not fixable by deselecting. It drops every row an
  issue names, not the fewest that would clear it. Keep it off the render path: `plan` reads
  the run logs for its timings, so `refresh()` uses the cheap `blocked_rows` instead.
- **The Results page reads only the quantile tables**, never the mcdf: `<Output file>_<type>.csv`
  (monte carlo) and `<Output file>__<type>_quantiles<suffix>.csv` (reservoir routing) hold the
  same three columns, so `core/results.py` has one reader for both. It draws with `ui.echart`
  and takes its standard normal variate from `statistics.NormalDist`, because the allow-list
  rule means the UI environment has neither matplotlib nor scipy — do not reach for
  `ui.pyplot` or `scipy.special.ndtri` here.
- **The Events page is the one that reads the mcdf**, which is unavoidable: a representative
  event is a realisation, not a quantile. `core/events.py` finds the database with
  `outputs.find_database` (so a routed row gets its suffixed mcdf), and caches the prepared
  frame per file mtime and result type — an mcdf is m x n rows and the page re-ranks on every
  control change. Leave a blank source on a loading and it takes the event from the duration
  that is **critical at that loading's AEP**, reusing `results.compare`/`analyse`; for level
  that genuinely differs across the frequency range, so one list of loadings draws from several
  runs. A level loading is read off the **envelope**, not one duration's curve, because the
  envelope is the design quantile the level was quoted from. **A level loading therefore has two places on the
  result axis and the plot marks both** (`variate_at_value`, `data_z`): the design AEP off
  that envelope, and the AEP this database's own realisations reach the level at. They
  differ by the envelope-versus-duration gap, the straight-line read between the standard
  AEPs, and the smoothing — nothing to do with rounding, and marking only the first made the
  line look offset from the level that was asked for. The chosen list is saved per group
  as `<group>_representative_events.json` beside the databases — per group because a GWL series
  usually shares one results folder.
  **`core/hydrographs.py` is the only thing that reads a stored hydrograph file**, and it does
  so on the button, never on a redraw: one column per simulation is tens of megabytes, so the
  frame is cached per file (mtime and size) and the read is pushed off the event loop. It
  answers the question the mcdf cannot — whether an event is one rise or two — with
  `shape_of`, whose prominence rule (a fifth of the peak, a tenth of it above the preceding
  trough) is what stops every step on a routed recession counting as a second flood.
  `run.io_bound` returns None on cancellation *and* when there is no live pool, which is
  indistinguishable from a result, so `_off_thread` falls back to running inline unless the
  app is stopping.
  **Which loading cards are open is view state the page has to keep** (`open_cards`): picking an
  event redraws every expansion, so a `value=position == 0` on the rebuild collapsed the card
  being worked on and sprang the first one open under it.
  **Ranking has two orders and neither is a filter** (`rank(order=...)`): `delta_z` is the
  distance on both axes at once, `result` is the distance from the loading in the result's own
  units (`delta_value`, ties broken on `delta_z`) — for the common case where hitting the lake
  level matters more than AEP neutrality. A level loading names that value; a design AEP has it
  read back off the envelope by `value_for_aep`, the inverse of `aep_for_level` on the same
  (log value, z) interpolation. Both distances are always computed, so whichever is not ranked
  on is still in the table. The `result` order is **banded and rounded**
  (`banded`), because a level is not meaningful to the millimetre and an exact sort on the
  difference is decided by noise: the *band* (default 20 mm of level) is how far off still
  counts as reaching the loading — one group, ordered by `delta_z` — and the *rounding*
  (default 10 mm) is the grid everything further out is measured on, so equally distant events
  are separated by neutrality too. They are two settings answering two questions, both per
  result type, and the UI asks for them in millimetres for a level, m³/s for a flow.
- **A critical-duration crossover is judged over the range the new duration holds, not at the
  crossing point.** The margin at a crossover is near zero by definition — the curves are equal
  there — so measuring strength there would dismiss every real crossover as noise. `Band.peak_margin`
  is the quantity that matters; a switch whose band never clears the noise floor is an
  `idxmax` hop and is reported rather than pinned. **The margin's unit is per result type**
  (`margin_scale`): percent for flows and volumes, **metres for level**, because level is an
  interval scale on an arbitrary datum — a percentage of 217 m AHD is meaningless, and on the
  real Callide E010 results a percentage floor dismissed all six level crossovers as noise
  when in metres they are 0.02–0.10 m and four of them are real. Nothing assumes which way the critical duration
  moves with AEP: it lengthens with rarity for inflow but *shortens* for lake level, which goes long
  at frequent AEPs while the storage fills and short on the rare tail as the dam becomes a conveyance.
- `ui/tests/test_results_page.py` renders the page through `nicegui.testing.user_simulation`,
  which is why `pytest-asyncio` is in `requirements-ui.txt`. It builds its own fixture instead of
  enabling the nicegui pytest plugin, so the rest of the suite still runs without nicegui
  installed — and the plugin proper needs selenium, which nothing here wants.

### Post-processing utilities (`util/`)
Standalone scripts with editable paths at the top of `main()`, e.g. `PlotFrequencyCurves.py`
(frequency plots), `GetRepresentativeEvents.py` (representative event selection),
`DesignFloodInterpolation.py`, `MaxQuantiles.py`, `ReportCollation.py`. See
`Manual/SubDocs/utilities.md`.
- `RepresentativeEvents.py` is the launcher's other half: an argparse CLI that takes the Events
  page's saved JSON and writes the hydrographs, the rebuilt hyetographs, a three-panel plot per
  event and a workbook. It shares `lib/RepresentativeEvents.py` with the page, so the event that
  gets plotted is the event that was chosen. Time runs from the **start of the main burst**; the
  stored hydrographs start at the beginning of the storm file, so they are shifted by the
  pre-burst duration — taken from the rebuild, or failing that from how much longer the run is
  than the model config's simulation period, which `Simulator.run_models` lengthens by exactly
  that amount. It reads the simulation period out of the JSON rather than through `UrbsModel`,
  whose constructor rmtree's the run's working folder.
  **The storm inputs come from the row that generated the storms, not the row the event
  came from**: a reservoir routing row has no `Duration` and no `Focal subcatchments`, so
  `storm_row` follows its `Input MCDF` back to the source run. Without that the rebuild
  reached `load_subcatchment_areas(None)` and pandas reported a NoneType buffer, which
  named neither the key nor the row. **The realisation is never what is missing**: a routed
  mcdf is the inherited one with the routed peaks written over it (`_write_mcdf`), so every
  draw survives — only those two sims-list keys do not, and they can come from the source row,
  from another sims list (`--source-sims-list`), or from the routing row itself.
- `GetRepresentativeEvents.py` picks representative events and then extracts and plots their
  hydrographs, driven by an `_analyseRepresentativeEvents.xlsx` control sheet with hard-coded
  paths. The **selection** half of it now also exists as `lib/RepresentativeEvents.py`, which
  the launcher's Events page uses; the hydrograph extraction has no equivalent there and is
  still the reason to run this script. If it is ever reworked into a CLI, take the selection
  from the shared module rather than keeping this copy of it.
- `CriticalDurationAnalysis.py` is the exception: a real CLI (argparse), because the UI's
  Results page shells out to it to export what it is showing. `CrticalDurationAnalysis.py`
  (sic — the typo is the older file) is the same analysis with hard-coded paths; both go
  through `UtilModule.MonteCarloSimulationGroup`, so keep the analysis there, not in either
  entry point.
- `UtilModule` fails whole-analysis rather than per-row, so its edge cases matter: an AEP no
  duration reached (all-NA row → `idxmax` raises; a dam that does not spill has no frequent
  outflow quantile), no `_perc_smooth` files at all (reservoir routing writes none →
  `pd.concat([])`), a volume column named `Vol24h` (not in the axis-label dict), and dropping
  an AEP the set does not have. All four are fixed and pinned by
  `ui/tests/test_critical_export.py`, which runs the real script.

## Environment

- Python 3.12. Dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`, `pyarrow`.
- `pyarrow` is needed only to read `.parquet` inflow hydrographs in `ReservoirRouting.py`
  (`_read_inflows`). Everything else works without it.
- `requirements.txt` — pip pins (note: newer/looser than the conda env).
- **Tests live in two places, in two environments.** `ui/tests` runs in the UI environment
  (`ui/requirements-ui.txt`, no scipy or matplotlib) — `cd ui && python -m pytest tests -q`.
  `tests/` runs in **Bryan's** environment, for the parts of `lib/` the UI is forbidden to
  import — `python -m pytest tests -q`. Put a test where the code it tests can be imported;
  a lib test in `ui/tests` fails `test_dependency_direction` unless the module is allow-listed.
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
