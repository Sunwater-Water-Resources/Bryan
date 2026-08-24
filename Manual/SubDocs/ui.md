# The run launcher

The **run launcher** is a browser interface for choosing which simulations in a simulation list to run, checking their inputs before committing to them, launching Bryan, and following it to the end. It is an alternative to editing the ```Include``` column in Excel and double-clicking a batch file; it is not a replacement for either, and Bryan itself is unchanged.

The code lives in ```ui/``` beside ```Main.py```, so one copy serves every dam catchment, exactly as the rest of Bryan does.

## What it does

The launcher writes its own copy of the simulation list containing only the rows you selected, writes a config file pointing at it, and runs ```python Main.py <that config>``` — the same command a batch file runs. Nothing is passed to Bryan that a batch file could not have passed.

Each run leaves a folder under ```_ui_runs/``` in the project holding the simulation list that ran, the config, the run log, the console output, and a ```launch.bat``` that repeats the run without the launcher. That folder is the record of what was actually run.

## Setting it up

The launcher runs in its **own** Python environment. Do not add it to the conda environment used to produce study results — see [installation](installation.md) and the note in ```_env_bryan.yml```.

```bat
py -m venv C:\PythonProjects\Bryan_ui\.venv
C:\PythonProjects\Bryan_ui\.venv\Scripts\pip install -r C:\PythonProjects\Bryan_dev\ui\requirements-ui.txt
```

Then put a ```run_ui.bat``` in the model folder beside the existing simulation batch files:

```bat
@echo off
setlocal
cd /D "%~dp0"
set "BRYAN=C:\PythonProjects\Bryan_dev"
set "UI_PY=C:\PythonProjects\Bryan_ui\.venv\Scripts\python.exe"
"%UI_PY%" "%BRYAN%\ui\main.py" "TFD_sims_config_LongList_02.json" ^
    --bryan-python "%BRYAN%\env\python.exe" --bryan-main "%BRYAN%\Main.py"
if errorlevel 1 pause
```

The two ```--bryan-*``` values are the same pair the ordinary batch files set as ```VENV_PY``` and ```PYFILE```. They are saved, so they only need giving once.

## Choosing simulations

The **Select** page lists the simulation list as a table, filtered by group.

A *group* is a set of rows that differ only by storm duration — the set a critical duration is picked from. Where the simulation list has a ```Group``` column it is used as-is. Where it does not, one is derived by removing the duration from the row's output name, and only when the number matches that row's own ```Duration``` — so a fixed ```_24h``` in a model name, or a ```GWL1p7```, is never mistaken for the varying part. The **Edit** page can write the derived keys into a real ```Group``` column in a new copy of the list.

Every row also carries a state worked out from the files on disk:

| State | Meaning |
| ----------- | ----------- |
| ```not run``` | No results database exists for this row. |
| ```up to date``` | Results exist and are newer than every input the row reads. |
| ```stale``` | Results exist, but an input has changed since they were written. |
| ```incomplete``` | The results database is shorter than the Monte Carlo config's main × sub divisions — a run stopped early by ```test_runs```, or one that crashed. |
| ```needs prior results``` | ```Run models``` is ```no```, so the row only re-analyses — and there is nothing to analyse. |

```stale``` is the state worth having. Re-routing under an edited ```.sq```, or re-running after the storm or climate config changed, otherwise leaves results that look perfectly fine.

Selecting a row whose results are ```up to date``` needs an explicit confirmation before launching: Bryan rewrites the results files and deletes the model working folder on entry, and there is no undo.

## What is checked before a run

- Columns Bryan will read that the list does not have, and required values left blank. A missing column is a ```KeyError``` inside the simulator, which loses the whole simulation.
- Input files that are not where the row says. This is the usual reason a simulation fails, which is why the run log flags it separately.
- Rows that would **write over each other** — two rows sharing an ```Output file```, an ```Output suffix``` and a ```Results folder``` produce the same result files, so only the last survives. This is unsafe whether they run at the same time or one after another.
- A ```Method``` whose spelling Bryan will not accept. The comparison is exact and lower case, so ```Ensemble``` fails.
- Formulas whose results were never cached — see below.

## Running several at once

The launcher can split a selection across several Bryan processes. It defaults to **one**, and that is usually right: reservoir routing takes seconds, and a Monte Carlo simulation already drives thousands of model runs one after another, so more processes mostly multiply peak memory and the number of storm files on disk. The per-chunk time estimate shown before launching — taken from the measured durations in your own run logs where they exist — is the thing to judge it by.

Rows sharing an ```Output file``` are never split across processes. Bryan derives the model working sub-folder from that name and deletes the folder when the simulation starts, so two simulations sharing one would destroy each other's working files.

Each process writes its own run log inside the run folder, because the log path is derived from the simulation list's filename. The launcher reads them together.

## Stopping a run

Two different things:

- **Stop after the current chunk** drops anything not yet started and signals nothing. Nothing is left half-written.
- **Stop now** kills Bryan and the URBS or RORB processes underneath it. The simulation in flight never reaches the run log, so it has **no entry there at all**; its model working folder and its log file are left part-written.

## Viewing results

The **Results** page plots the frequency curves a group has already produced, one line per storm duration, so the critical duration can be read off them. Pick the group, pick the result type — ```inflow```, ```level```, ```outflow```, or one of the inflow volume windows where the volume analysis has run — and tick the durations to compare.

The files it reads are the standard-AEP quantile tables the analysis writes: ```<Output file>_level.csv``` from the Monte Carlo method, and ```<Output file>__level_quantiles<suffix>.csv``` from reservoir routing. Both hold the same three columns, so a re-routed result plots beside an original one. Nothing is computed from the raw databases and nothing is re-run; a row that has not been analysed simply does not appear.

The **maximum envelope** is the largest value at each AEP across the ticked durations. With **mark up critical durations** on, the AEP axis is shaded by which duration owns the envelope, each crossover is pinned, and the table below gives the critical duration at every standard AEP.

### Whether the durations bracket the critical one

This is what the page is for, so it is said in words rather than left in the picture. Whenever the shortest or the longest duration you ran is the critical one at some AEP, the page says which and at what AEPs — the critical duration may lie beyond the range that was run, and only running further out will show it.

The direction to expect is **not the same for every result type**. For lake level on a dam the critical duration tends to be **longer at frequent AEPs**, because it takes rainfall volume to fill and charge the storage, and **shorter on the rare tail**, as the dam starts to act more like a conveyance system than a volume system. Peak inflow is set by the catchment's response instead and does not behave this way. The page therefore assumes no direction: each AEP is judged on its own, and the *critical duration against AEP* plot under the main chart draws the shape directly so a departure from it is visible.

### Crossovers and sampling noise

Where two duration curves are nearly coincident — the rare tail of a level curve especially, once everything is spilling — the winning duration changes from one AEP to the next on Monte Carlo sampling noise alone, and that is indistinguishable from a real crossover by eye.

The margin column is how far the critical duration beat the runner-up, **in the units of the result type**: percent for inflow, outflow and volumes, and **metres for level**. The distinction matters. Flows and volumes are ratio scales, where a percentage is the natural measure; a lake level is an interval scale on an arbitrary datum, and a percentage of a number like 217 m AHD says nothing about it. At Callide the durations separate by 0.01 to 0.10 m across the frequency range — 0.005 to 0.05% — so judging level on a percentage floor would report every crossover it has as noise.

The floor a switch has to clear is shown beside the plot and can be changed: it defaults to 1% for flows and volumes and 0.05 m for level, which is a rule of thumb rather than a physical constant.

Note also that the margin *at* a crossover is near zero whatever the crossover means, since the two curves are equal there; a switch is therefore judged on how convincingly the incoming duration wins over the range it then holds. One that never gets clear of the next duration is reported as noise and is not pinned on the plot.

### Exporting the analysis

**Export critical durations** writes the analysis to file, one run per result type:

| File | Holds |
| ----------- | ----------- |
| ```<name>_<type>_critical.csv``` | every duration's quantiles, the maximum at each AEP, the critical duration, and the confidence limits at that duration |
| ```<name>_<type>_critical_durations.png``` | the duration curves |

The dialog lists every file before writing any of them and marks the ones that already exist, since it overwrites without asking afterwards. The output name defaults to the group's name with the duration removed, and the folder to wherever the quantile files are.

Behind the button is ```util/CriticalDurationAnalysis.py```, run with **Bryan's** interpreter rather than the launcher's — the exported table also carries the smoothed confidence percentiles, which is a fifth-order polynomial fit in log space, and a second implementation of that in the launcher would be a second thing to drift. The files that come out are the same ones the study post-processing has always produced. That script can also be run by hand:

```bat
python util\CriticalDurationAnalysis.py --result-type level ^
    --output-folder ...\results --output-name CLD_mc_E010_level_critical ^
    --sim 24 ...\CLD_mc_24h_E010_level.csv ^
    --sim 36 ...\CLD_mc_36h_E010_level.csv
```

```CrticalDurationAnalysis.py``` beside it is the same analysis driven by a hard-coded block of paths, for a study you are already sitting in front of.

Results routed by the reservoir routing method export without the confidence columns: that method writes the quantile files but not the ```_perc_smooth``` files the limits come from.

### What is not shown

The ensemble method does not appear. ```lib/EnbAnalysis.py``` already computes the median pattern per duration and the critical duration per AEP inside the run, and writes its own plots beside the results; a second implementation in the launcher would only be somewhere for the two to disagree.

Only inflow volumes are offered, for the reason the volume analysis itself gives: outflow and storage volumes follow from the peak level and the rating curve.

## Things worth knowing

### Formulas without cached values

Simulation lists are formula-driven — ```Output file```, ```Input MCDF```, ```Inflow```, ```SQ file```, ```Log file``` and ```Basename``` are usually concatenation formulas. Excel stores each formula's *result* when it saves the file, and that stored result is what Bryan reads. A workbook rewritten by a Python tool keeps the formulas but loses the results, and Bryan then reads those cells as **blank**.

This is not hypothetical: ```TFD_SimsList_LongList_02.xlsx``` in the Tinaroo project is in this state, so its first 24 rows read with no output name and no input files. Every other simulation list in that project and in Callide is intact.

The launcher checks for this when it opens a list, says which columns and rows are affected, and refuses to run them. The fix is to open the workbook in Excel, press **Ctrl+Alt+F9** to recalculate, and save.

It is also why the launcher **never writes your master workbook**. Edits are saved to a new file you name, and both those copies and the run copies hold literal values rather than formulas.

### Console questions

Bryan calls ```input()``` in three places — twice when an output csv cannot be written because it is open in Excel, and once when the URBS executable is not where the model config says. A background process has no console to answer on, so the launcher gives Bryan no standard input: those calls fail immediately and the simulation is recorded as ```FAILED``` rather than hanging overnight. The Run page explains this when it happens.

### Whole and fractional durations

A spreadsheet cell cannot distinguish 120 from 120.0, and the type is inferred for the whole column. A list with a fractional duration anywhere — Callide has a 4.5 hour storm — is read as fractional throughout, so Bryan builds storm filenames with ```120.0```; a selection that leaves that row out is read as whole numbers and builds them with ```120```. That also changes the ```simulation_periods``` lookup in the URBS config, which falls back to twice the storm duration when it misses. The launcher warns when a selection would do this. Reservoir routing is unaffected, as it takes its durations from the stored results.

Return to [Main Manual](../Manual.md)
