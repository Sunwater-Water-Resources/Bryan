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

Reservoir routing makes steps 2–5 cheap — seconds per simulation, and no model
executable in the loop for step 4's comparison to be muddied by.
