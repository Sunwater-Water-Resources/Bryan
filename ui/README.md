# Bryan run launcher

A browser UI for choosing which simulations in a sims list to run, checking
their inputs, launching Bryan, and watching it go.

Bryan itself is unchanged. The UI writes its own copy of the sims list holding
only the selected rows, writes a config pointing at it, and runs
`python Main.py <that config>` — the same command a batch file runs.

## Install

The UI lives in its own environment. **Do not add it to `bryan29`** — that env
is frozen so study results reproduce, and `_env_bryan.yml` is a
`conda env export` with exact build strings.

```
py -m venv C:\PythonProjects\Bryan_ui\.venv
C:\PythonProjects\Bryan_ui\.venv\Scripts\pip install -r ui\requirements-ui.txt
```

## Run

```
python ui\main.py [sims_config.json] [--port 8081]
```

Or drop a `run_ui.bat` beside the model's existing `*_sims_*.bat` files — see
`Manual/SubDocs/ui.md` for the template. On the Project page, set the
interpreter and `Main.py` that Bryan should be run with: the same pair the
batch files set as `VENV_PY` and `PYFILE`.

## What it does

**Select** — the sims list as a table, filtered by group. Where there is no
`Group` column, one is derived by stripping the duration from the output name,
so the durations of a case collapse into one row of the filter.

Each row carries a status worked out from the files on disk:

| | |
|---|---|
| not run | no results database |
| up to date | results exist and are newer than every input |
| stale | results exist but an input changed since |
| incomplete | the results database is short of the expected realisations |
| needs prior results | `Run models` is no, and there is nothing to analyse |

**stale** is the one worth having. Re-routing under an edited `.sq`, or
re-running after the storm config changed, leaves results that look fine.

Pre-flight then checks what Bryan would trip over: missing columns, missing
input files, blank required values, a `Method` whose capitalisation Bryan will
not accept, and rows that would write over each other.

**Run** — one card per process: which simulation is going, how far through,
elapsed time, the per-simulation log growing, and the console output. Two stop
buttons, because they mean different things:

- *Stop after the current chunk* drops the queue and signals nothing.
- *Stop now* kills Bryan and the URBS/RORB processes underneath. The simulation
  in flight never reaches the run log, so it gets **no entry there at all**,
  and its working folder and log are left part-written.

**Results** — two tabs over the same files, and neither runs anything.

*Durations* plots the analysed frequency curves of one group on top of each
other so the critical duration can be read off. Pick the group and the result
type (inflow, level, outflow, or a volume window), tick the durations, and
optionally add the maximum envelope over them. With the mark-up on, the AEP
axis is shaded by which duration owns the envelope and each crossover is
pinned.

What it is really for is the question of whether the durations that were run
bracket the critical one, so that is said in words rather than left to be read
off the picture: whenever the shortest or longest duration you ran is critical
at some AEP, the page says so and at which AEPs.

Two things to know when reading it:

- **The direction is not the same for every result type.** For lake level the
  critical duration tends to be long at frequent AEPs — it takes rainfall
  volume to fill and charge the storage — and shorter on the rare tail as the
  dam behaves as a conveyance rather than a volume system. Peak inflow follows
  the catchment and does no such thing. Nothing in the page assumes a
  direction; every AEP is checked on its own.
- **The margin column is what separates a crossover from noise**, and it is
  measured in the units of the result type: percent for flows and volumes,
  **metres for level**. Level is an interval scale on an arbitrary datum, so a
  percentage of 217 m AHD says nothing — at Callide the durations separate by
  0.01–0.10 m, which as a percentage is 0.005–0.05%, and a percentage floor
  dismisses every real level crossover as noise. The floor itself (1%, or
  0.05 m) is on the page and can be changed per plot.
- **A crossover is judged over its range, not at the crossing point.** The two
  curves are equal where they cross, so the margin there is near zero whatever
  the crossover means. A switch is judged on how convincingly the new duration
  wins over the range it then holds; one that never gets clear of the next
  duration is reported as noise and left unpinned.

*Groups* plots one line per **group** instead — its envelope over the durations
it ran, which is the design quantile — so climate scenarios and dam options can
be compared with each other. There is deliberately no envelope across the
groups and no mark-up: groups are scenarios, not alternatives to be enveloped.
What goes underneath is the change from a chosen **baseline**, in metres for
level and percent for flows and volumes, and the critical duration is overlaid
per group so it is visible whether that moves between scenarios too. Because an
overlay hides how each envelope was built, the tab reports uneven duration
coverage, an envelope pinned to the end of its own range, a group mixing two
methods, and groups that stop at different AEPs. Groups come from the open sims
list; to compare against another one, use `util/PlotFrequencyCurves.py`.

**Export critical durations** writes the analysis out, one run per result type:

```
<name>_<type>_critical.csv           the table: every duration, the maximum,
                                     the critical duration, and the confidence
                                     limits at it
<name>_<type>_critical_durations.png the duration curves
```

The dialog shows every file before it writes any, and marks the ones already
there. It runs `util/CriticalDurationAnalysis.py` with **Bryan's** interpreter
rather than recomputing anything here: the exported table also carries the
smoothed confidence percentiles, which is a polyfit in log space, and a second
copy of that would be a second thing to drift. The files are then identical in
format to what the study post-processing already produces.

Reservoir-routing results export without the confidence columns, because that
method writes quantiles but no `_perc_smooth` files.

The ensemble method does not appear here: `lib/EnbAnalysis.py` computes the
critical duration per AEP inside the run and writes its own plots, and a second
implementation would be somewhere for the two to disagree.

**Events** — picking a representative event for each design flood loading.

Give it a list of loadings, as an AEP or as a lake level, and it ranks the
realisations of the Monte Carlo database against each one. The rank is distance
from the loading on two axes at once: the flood as rare as asked, and the
rainfall about as rare as the flood. That is AEP neutrality — an event reaching
the 1 in 2,000 level off 1 in 200 rainfall got there by coincidence and will
not behave like a 1 in 2,000 event when the dam is changed.

Distance is in **standard normal variate space**, not `1 in X`: at 1 in 2,000
being 200 out is nothing, at 1 in 100 it is everything, so a `1 in X` distance
ranks the rare end nearly at random. (`delta_aep`, the measure
`util/GetRepresentativeEvents.py` sorts on, is still computed for comparison.)

Closeness alone is not enough, so each candidate carries what would make it
indefensible anyway: the worst `subburst_<d>h / ifd_<d>h` ratio in the storm
(above 1.0 is an embedded burst, measured rather than described), the
`embedded_bursts` comment, the pre-burst percentile, and the antecedent storage
as a z. **Nothing is dropped quietly** — those all flag and the event is still
offered, because the trade-off between the closest match and the cleanest storm
belongs to whoever defends the event. Dropping them is opt-in. The one
automatic filter is the rainfall cap at 1.1x the AEP of the PMP, which is the
edge of the sampling scheme rather than a real event; the number comes from the
IFD files config the storm config points at, as a Monte Carlo run gets it.

Leave the source blank and the event comes from the duration that is **critical
at that loading's AEP**, off the same envelope the Results page draws — so a
list of level loadings spanning the frequency range legitimately draws its
events from different runs. A lake level is converted to an AEP off the level
envelope, and a level above the top of the curve is reported rather than
extrapolated, the page then offering the highest events instead.

Save writes `<group>_representative_events.json` beside the databases and reads
it back next time; Export csv writes the list with the chosen event's metrics.
The `hydrograph` column names the column to pull out of the stored flows —
simulation 42 is `sim_00042`. Extracting and plotting those hydrographs stays
with `util/GetRepresentativeEvents.py`.

This is the only page that reads an mcdf rather than the quantile tables, which
is unavoidable: a representative event is a realisation. Ensemble rows do not
appear — an ensemble run ran every combination, so it has no realisations to
choose between.

**Edit** — change cells and save to a **new** file. The master workbook is
never written; see below.

**History** — past runs, when each output last ran, and deleting old run
folders.

## Things worth knowing

**Uncached formulas.** Sims lists are formula-driven. Excel stores the result
of each formula when it saves; openpyxl throws those results away when it
writes a workbook. A list rewritten by a Python tool therefore has formulas
with no values, and pandas — and so Bryan — reads them as blank.

This is real: `TFD_SimsList_LongList_02.xlsx` in the Tinaroo project has 96
formula cells and no cached values, so its first 24 rows read with no output
name and no input files. Every other sims list in Tinaroo and Callide is
intact. The UI flags this on load and refuses to run the affected rows. The
fix is to open the workbook in Excel, press Ctrl+Alt+F9 and save.

It is also why the UI never writes your master workbook. Run copies and
save-as files hold literal values, with a banner saying so.

**Why parallelism defaults to 1.** Reservoir routing takes seconds; a Monte
Carlo simulation already drives thousands of model runs one after another. So N
processes mostly multiply peak memory and storm-file count rather than saving
wall clock. The mechanism is there — judge it by the per-chunk time estimate,
which comes from the measured durations in your own run logs where they exist.

The hard limit on splitting is `Output file`: `Simulator.initialise_model` makes it
the URBS working sub-folder and `UrbsModel.__init__'s rmtree` deletes that folder on
entry, so two simulations sharing one would destroy each other. Rows sharing an
`Output file` are never split, and selecting several of them is blocked
outright — they overwrite each other even run one at a time. In
`CLD_RFSL_mc_sims_01.xlsx` twenty-four rows share one output name; running that
list one row at a time is what avoids the problem today.

**Console questions.** Bryan calls `input()` in three places — when an output
CSV is open in Excel (`MCScheme`/`EnbScheme.store_simulations`) and when
the URBS executable is missing (`lib/URBSmodel.py:41`). A background process
has no console to answer on, so the UI gives it no stdin: those calls raise
`EOFError` immediately and the simulation is recorded as failed instead of
hanging overnight. The Run page explains it when it happens.

**A trap with `Duration`.** An xlsx cell cannot carry the int/float
distinction, and pandas infers it per column. A list with a fractional duration
anywhere reads the whole column as float, so `URBSmodel.duration_string` gives
`'120.0'`; a selection without that row reads int and gives `'120'`. That
changes every storm filename, and the simulation-period lookup misses and
silently falls back to twice the storm duration (`URBSmodel.py:287`). The UI
warns when a selection would do this. `CLD_RFSL_mc_sims_01.xlsx` has a 4.5 hour
duration, so it is affected — though only for rows that run a model, since
reservoir routing takes its durations from the stored results.

## Tests

```
cd ui && python -m pytest tests -q
```

No URBS or RORB needed. `tests/fake_main.py` stands in for `Main.py`, taking
the same command line and importing the **real** `lib.RunLog` and
`lib.LogFiles` so the run-log format cannot drift from what the UI parses. It
covers the parallel cap, tree-kill of grandchildren, reattaching after a UI
restart, and the `EOFError` path.

`tests/test_progress.py` greps `Main.py` and `lib/RunLog.py` for the literal
strings the progress display is parsed from, so a change to Bryan's wording
fails loudly instead of quietly showing 0 of N forever.

`tests/test_representative_events.py` covers the representative event analysis
in `lib/RepresentativeEvents.py` against a synthetic mcdf carrying exactly the
columns `lib/MCScheme.py` writes — including the mixed units, `rain_aep` in
"1 in X" and the TPT columns as probabilities, which is the failure that module
is most exposed to and would not look wrong on a plot.

`tests/test_real_bryan.py` runs **real Bryan** on a miniature reservoir-routing
model built by `tests/make_rr_fixture.py` — an `.els`, two `.sq` curves, a small
ensemble database and its inflows. It is the only test that proves the run
copies and configs the UI writes are acceptable to Bryan itself rather than
just to a stand-in that shares the UI's assumptions. It needs an interpreter
with Bryan's dependencies, so point it at one:

```
BRYAN_TEST_PYTHON=C:\PythonProjects\Bryan_dev\env\python.exe python -m pytest tests -q
```

Without that it skips, and the rest of the suite still runs.

## Windows acceptance checklist

The tests cover the machinery; this is what shows the UI runs Bryan correctly.

1. Open a real project. Confirm the row count, groups and statuses look right.
2. Pick one reservoir-routing group. Run it from the UI with **one** chunk.
3. Run the same rows from their existing `.bat` into a separate results folder.
4. Diff the two sets of outputs. They must be identical.
5. Repeat with three chunks and diff again.
6. Stop a run mid-way with *Stop now*; confirm no `urbs32.exe` survives in Task
   Manager.
7. Open **Events** on a real Monte Carlo group. Check the chosen event for a
   loading against what `util/GetRepresentativeEvents.py` picks for the same
   one: the two rank on different measures, so they need not agree, but a large
   disagreement is worth understanding before trusting either.

Reservoir routing makes steps 2–5 cheap — seconds per simulation, and no model
executable in the loop for step 4's comparison to be muddied by.
