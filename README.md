# Bryan

**Bryan** is Sunwater's Python platform for design flood hydrology simulation. It drives
external hydrologic models (**URBS** and **RORB**) through either a **Monte Carlo** or
**Ensemble** scheme to estimate design flood quantiles (peak inflow, lake level, outflow),
then post-processes results using the Total Probability Theorem (TPT). It implements
Sunwater's design flood hydrology specification.

Developers: Richard Sharpe and Graigan Panosot.

## What it does

- Runs URBS and RORB models through a **Monte Carlo** scheme (stratified sampling of
  rainfall, temporal patterns, storm method, losses, preburst, and antecedent lake level)
  or an **Ensemble** scheme (a grid of AEPs x durations).
- Applies **climate change** adjustments to rainfall intensity, losses, and temporal
  patterns per the draft 2023 update to the ARR Climate Change Considerations chapter.
- Analyses results using the Total Probability Theorem (Monte Carlo) or critical-duration
  selection (Ensemble), producing csv tables and plots.
- Includes a standalone reservoir routing method (`ReservoirRouting.py`) for routing
  existing inflow hydrographs through a dam.

See [`Manual/Manual.md`](Manual/Manual.md) for the full user guide, and
`Manual/Bryan_Technical_Reference_v1.pdf` for the technical reference.

## Requirements

- Python 3.12
- `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`, `pyarrow` (see `requirements.txt`
  for pinned versions, or `_env_bryan.yml` for the conda environment used to produce study
  results)
- `pyarrow` is only needed for reading `.parquet` inflow hydrographs in reservoir routing;
  everything else works without it.
- Designed for **Windows**: run via a batch file, with external URBS/RORB executables and
  backslash paths in configs. Running the full pipeline on Linux requires the model
  executables and will hit path/environment assumptions.

## Running

Bryan is driven from a batch file that passes a config file path to a top-level script:

```
python Main.py <sims_config.json>          # main entry point
python RouteFlows.py <routing_sheet.xlsx>  # standalone baseflow addition / dam routing
```

The config file names a simulation list (Excel) and points to per-method config files
(Monte Carlo, Ensemble, Storm, Lake, Model, Climate). See
[`Manual/SubDocs/getting_started.md`](Manual/SubDocs/getting_started.md) for the full
control-flow walkthrough.

There is also a browser UI for choosing which simulations in a list to run, checking their
inputs first, and following the run:

```
python ui/main.py <sims_config.json>
```

It writes its own copy of the simulation list holding the selected rows and runs `Main.py`
on it, so Bryan itself is unchanged. It installs into its own environment -- see
[`ui/README.md`](ui/README.md) and [`Manual/SubDocs/ui.md`](Manual/SubDocs/ui.md).

The Python code and the model/project data are kept separate, so one copy of Bryan's code
serves many dam catchments; models, catchment data, and study outputs are not part of this
repository.

## Try it without a model

[`example_project/`](example_project/README.md) is a complete, self-contained project for a
catchment that does not exist -- Juniper Creek Dam, 620 km² over ten subcatchments. It is
there so a new user has working config files to copy, and can run something end to end before
pointing Bryan at a real study.

URBS and RORB are licensed programs this repository does not ship, but two of the three
methods do not need them, and those are the ones set up to run:

```
cd example_project
./run_storms_only.sh     # or run_storms_only.bat on Windows
./run_routing.sh
```

- **`run_storms_only`** is a Monte Carlo row with `Run models` set to `storms only`: rainfall
  setup, areal reduction, storm method sampling, all four temporal pattern sets, the embedded
  burst filter, the sub-burst neutrality check and the PMP embedded burst screen. It writes
  the results database without calling a hydrologic model.
- **`run_routing`** re-routes a stored set of inflow hydrographs through the dam: Modified
  Puls routing, antecedent volume handling, the Total Probability Theorem, quantile tables
  and frequency curve plots.

Between them they cover the whole chain except the rainfall-runoff step. Complete URBS and
RORB config templates are included for that step, and both construct on a machine without the
executables, so a setup can be checked before a licence is involved.

The generic ARR and BoM reference data is real -- the ARR point and areal temporal patterns,
the BoM GSDM and GTSMR patterns, the ARR climate change tables. Everything specific to the
fictional catchment is invented and prefixed `SYNTHETIC_`: the IFD depths, PMP depths, dam
curves and stored results. Copy the config structure, not the numbers. See
[`example_project/README.md`](example_project/README.md) for the full split and a list of
setup gotchas worth knowing before building a real project.

## Repository layout

- `Main.py`, `RouteFlows.py` -- top-level entry points (thin dispatchers)
- `lib/` -- core logic: simulators, storm generation, rainfall, temporal patterns,
  sampling, hydrologic model wrappers, routing, lake conditions, climate change,
  interpolation curves
- `ui/` -- the run launcher (NiceGUI). `ui/core/` is Bryan-facing logic with no UI
  imports; `ui/pages/` is the interface; `ui/tests/` runs without URBS or RORB
- `util/` -- standalone post-processing scripts (frequency plots, representative events,
  design flood interpolation, etc.)
- `example_project/` -- a runnable, self-contained example project for a fictional catchment
- `Manual/` -- user guide and technical reference; run `Manual/render_manual.py` to render
  the Markdown to browsable HTML (Markdeep)
- `IFD_export.py`, `DownstreamStormGenerator.py`, `StormInstance.py` -- legacy/ad-hoc
  tooling with hard-coded paths, not part of the main batch workflow

## Documentation

- [`Manual/Manual.md`](Manual/Manual.md) -- user guide
- `Manual/Bryan_Technical_Reference_v1.pdf` -- authoritative technical reference
- [`Manual/change_log.md`](Manual/change_log.md) -- design decisions and config-format
  changes
- [`example_project/README.md`](example_project/README.md) -- the worked example, and what
  to watch for when setting up a real project
