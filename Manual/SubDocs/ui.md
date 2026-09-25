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
| ```not run``` | None of this row's own results exist. |
| ```up to date``` | Results exist and are newer than every input the row reads. |
| ```stale``` | Results exist, but an input has changed since they were written. |
| ```incomplete``` | The results database is shorter than the Monte Carlo config's main × sub divisions — a run stopped early by ```test_runs```, or one that crashed. |
| ```needs prior results``` | ```Run models``` is ```no```, so the row only re-analyses — and there is nothing to analyse. |

For ```reservoir routing``` the state is read from the outputs that carry the row's ```Output suffix``` — the routed database, the quantile tables, or the stored hydrographs where the row does not analyse. A results database written before 26 August 2026 has no suffix on it (see the change log), and every rating curve run over one ```Output file``` wrote that same file: judging by it marked all six Tinaroo curves of a duration ```up to date``` as soon as the first one finished, so it is not read as evidence. It still counts as *prior results* for a ```Run models = no``` row, because it is a file the re-analysis can still read.

```stale``` is the state worth having. Re-routing under an edited ```.sq```, or re-running after the storm or climate config changed, otherwise leaves results that look perfectly fine.

Selecting a row whose results are ```up to date``` needs an explicit confirmation before launching: Bryan rewrites the results files and deletes the model working folder on entry, and there is no undo.

## What is checked before a run

- Columns Bryan will read that the list does not have, and required values left blank. A missing column is a ```KeyError``` inside the simulator, which loses the whole simulation.
- Input files that are not where the row says. This is the usual reason a simulation fails, which is why the run log flags it separately.
- Rows that would **write over each other** — two rows sharing an ```Output file```, an ```Output suffix``` and a ```Results folder``` produce the same result files, so only the last survives. This is unsafe whether they run at the same time or one after another.
- A ```Method``` whose spelling Bryan will not accept. The comparison is exact and lower case, so ```Ensemble``` fails.
- Formulas whose results were never cached — see below.

Anything that would stop the run is counted in the line above the **Check and run** button. Two ways past it, for the common case where most of the selection is fine and a handful of rows point at a file that is not there:

- **Deselect problem rows** unticks every row a check would stop the run over, and says which and why.
- **Check and run** offers *Skip N and run the rest* in its dialog when the rows that are left could run on their own. Those rows are then deselected as well, so what is ticked and what is running stay the same thing, and the skipped rows are still there to fix and run afterwards.

Both drop **every** row named by a problem rather than the fewest that would clear it: where two rows write to the same output name, which of them you meant is not the launcher's to guess. Neither can help with a problem that names no particular row — a column the list does not have, or the ```Replicate file``` spelling — because deselecting rows does not fix the workbook. It says so rather than offering the button.

## Running several at once

The launcher can split a selection across several Bryan processes. It defaults to **one**, and that is usually right: reservoir routing takes seconds, and a Monte Carlo simulation already drives thousands of model runs one after another, so more processes mostly multiply peak memory and the number of storm files on disk. The per-chunk time estimate shown before launching — taken from the measured durations in your own run logs where they exist — is the thing to judge it by.

Rows sharing an ```Output file``` are never split across processes. Bryan derives the model working sub-folder from that name and deletes the folder when the simulation starts, so two simulations sharing one would destroy each other's working files.

Each process writes its own run log inside the run folder, because the log path is derived from the simulation list's filename. The launcher reads them together.

## Stopping a run

Two different things:

- **Stop after the current chunk** drops anything not yet started and signals nothing. Nothing is left half-written.
- **Stop now** kills Bryan and the URBS or RORB processes underneath it. The simulation in flight never reaches the run log, so it has **no entry there at all**; its model working folder and its log file are left part-written.

## Viewing results

The **Results** page has two tabs. **Durations** plots the frequency curves a group has already produced, one line per storm duration, so the critical duration can be read off them; **Groups** plots one line per *group* instead, so scenarios can be compared with each other. Both read the same files and neither runs anything.

On the **Durations** tab, pick the group, pick the result type — ```inflow```, ```level```, ```outflow```, or one of the inflow volume windows where the volume analysis has run — and tick the durations to compare.

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

### Comparing groups with each other

The **Groups** tab answers the other question: what a warmer climate, a raised full supply level, or a different antecedent storage assumption does to the design flood. Each group contributes **one** line — its maximum envelope over the durations it ran — because that envelope *is* the design quantile. Tick the groups to overlay, and the legend trims each group key down to the part that tells them apart, so ```sims_mc\results\TFD_mc_GWL1p3``` appears as ```GWL1p3```. The full key is on the checkbox tooltip and in *Files read*.

There is deliberately **no envelope over the groups** and no mark-up: groups are scenarios, not alternatives to be enveloped, so a maximum across them would not be a quantity. What goes underneath instead is the change from a **baseline** group, in the same units the margin column uses — metres for level, percent for inflow, outflow and volumes. Where the baseline is zero, as an outflow quantile is on a dam that does not spill at the frequent end, no change is shown and the page says so.

The *critical duration against AEP* plot under the chart is drawn per group, which is the one place it can be seen whether the critical duration itself moves between scenarios.

Overlaying envelopes hides how each one was built, so the tab says out loud what the picture cannot:

- **Envelopes over different numbers of durations.** An envelope over fewer durations is a lower bound, so the comparison is biased towards the better covered group. The duration count is shown beside every group.
- **An envelope pinned to the end of its own duration range**, for either group, for the reason above — a difference measured against a lower bound is not the difference.
- **A group drawing on more than one method**, which grouping normally separates but does not guarantee.
- **Groups that do not reach the same AEP**, which is ordinary rather than wrong: the AEP of the PMP comes from the storm config.

Groups come from the sims list that is open. Comparing against a run held in a *different* sims list is not offered here; use ```util/PlotFrequencyCurves.py```, whose plot list takes folders and filenames directly.

### What is not shown

The ensemble method does not appear here; it has its own page, below.

Only inflow volumes are offered, for the reason the volume analysis itself gives: outflow and storage volumes follow from the peak level and the rating curve.

## Viewing ensemble results

The **Ensemble** page is the Results page for ensemble runs. An ensemble run routes every temporal pattern of every storm duration at each standard AEP, so there are no quantile tables to read: the page reads the run's database itself. Pick a **group** - an ensemble row, or a reservoir routing row that re-routed an ensemble - and a **result** (level, inflow or outflow).

- **The durations.** Each duration's **median pattern** against AEP, with the envelope over them, marked up by which duration is critical, as the Results page draws the Monte Carlo curves. The median is Bryan's own pick (```lib/EnbAnalysis.py```): the pattern at position ```int(np.around(n / 2))``` of the ascending sort, the sixth of ten. The critical duration is the one whose median is largest, and a tie goes to the first duration the database lists, as Bryan's does. Each result type has its own critical duration: the inflow's is not the level's.
- **Critical durations.** Per AEP: every duration's median, the critical duration, its **margin** over the runner-up (in **metres for lake level**, percent for flows), which pattern gave the median, and the single **highest** event over every pattern and duration, with its duration. The same warnings as the Results page follow it: the shortest or longest duration being critical, and crossovers too small to mean anything.
- **The patterns at one AEP.** The PMF page's box plot: one box per duration over the patterns, every pattern a dot, the middle line Bryan's median pick and the diamond the highest event, with a table of each duration's median and highest pattern. It opens on the rarest AEP in the run.

Where Bryan's ```csv/<name>_critical.csv``` is beside the database, the page checks itself against it - the critical duration and its median at every AEP - and says it **agrees**, or lists where it differs. A difference usually means the csv predates the last run, and the page says so when the csv is the older file. Two small differences from Bryan's analysis never move a median's value: where several patterns give the same value (a lake held at full supply) the pattern named may differ, and an event with no result is left out of the count where Bryan counts it - the page reports any.

## Choosing representative events

A design flood quantile is a statistic over thousands of realisations, but the work downstream of it — a gate operation, a dambreak run, an emergency action plan trigger — needs a single **event**: one hydrograph, one storm, one starting lake level. The **Events** page picks that event.

Give it a list of loadings, either as a design AEP or as a lake level, and it ranks the realisations of the Monte Carlo database against each one. The list of chosen events is what comes out.

The saved selection is also the input to the **Downstream** page, which turns each chosen event into rainfall for a regional model - see [Downstream storm generation](downstream_storms.md) for the inputs that needs.

### What it ranks on

The rank is the distance from the loading, measured on two axes at once: the flood should be as rare as the loading asks, and the **rainfall that produced it should be about as rare as the flood**. An event that reaches the 1 in 2,000 lake level off 1 in 200 rainfall got there through a coincidence — a full lake, a large pre-burst, a spike in the pattern — and will not behave like a 1 in 2,000 event when anything about the dam is changed. That is *AEP neutrality*, and it is the diagonal on the plot.

Distance is measured in **standard normal variate space**, not in ```1 in X```. At a 1 in 2,000 loading, being 200 out is nothing; at 1 in 100 it is everything, and a distance in ```1 in X``` would rank the rare end almost at random. The ```1 in X``` distance the older ```util/GetRepresentativeEvents.py``` sorts on is still computed, as ```delta_aep```, so a selection made with that script can be checked.

**Rank by** chooses between the two, and it is the control to reach for when reaching the loading matters more than being neutral about it:

| Rank by | Sorts on | Use it when |
| ----------- | ----------- | ----------- |
| ```Δz (AEP neutral)``` | both axes at once — the default | the loading is a design AEP and the event has to be defensible as that AEP |
| ```Closest result``` | the result alone: how far the event is from the loading **in its own units** (metres for a lake level, m³/s for a flow), with Δz breaking the ties | the event exists to put the lake at a particular level, and the rarity of the rainfall is something to check afterwards rather than to rank on |

Two settings make ```Closest result``` usable, because a lake level is not meaningful to the millimetre — the rating curve, the routing timestep and the model itself are nowhere near that precision, so an exact sort on the difference is decided by noise and neutrality never gets a look in:

| Setting | What it does |
| ----------- | ----------- |
| **Same result within** | how far off still counts as **reaching** the loading. Everything inside it is one group, ordered by AEP neutrality — so of every event that gets there, the most neutral is offered first. Defaults to 20 mm of lake level |
| **Round differences to** | the grid everything further out is measured on. Two events that round to the same difference are the same distance away as far as anyone can defend, so they tie and neutrality separates them too. Defaults to 10 mm |

Both are asked for in **millimetres** for a lake level and in **m³/s** for a flow, and both are kept per result type — 20 mm of level is not 20 m³/s. Either can be set to zero on its own: with no band the closest event still leads, and with no rounding the events outside the band keep their exact order.

A lake level loading names the target outright. A **design AEP** loading is converted the other way — the design value at that AEP is read off the same envelope curve, and the page says what it read — so ```Closest result``` works for both. Where the curve cannot answer (an AEP off the end of it), the ranking falls back to the result AEP and says so. The distance is reported as **Δ target** in the candidate table and in the exported list whichever order is in force, so the measure not being ranked on can still be read.

### What it warns about

Closeness is not enough on its own, so every candidate carries the things that would make it indefensible:

| Metric | What to look for |
| ----------- | ----------- |
| **Sub-burst** | the worst ```subburst_<d>h / ifd_<d>h``` ratio in the storm. Above 1.0 the pattern contains a window rarer than the storm around it — an embedded burst, measured rather than described. Blank on a database written before those columns existed |
| **Flags** | the ```embedded_bursts``` comment the run wrote, and any pre-burst comment reporting an error |
| **Pre-burst p** | the sampled pre-burst percentile. Far from 0.5 means the flood was helped along by antecedent rainfall |
| **Lake z** | the antecedent storage as a standard normal variate. A large positive value means the event started with the dam already full |

Nothing is dropped quietly. Everything above **flags** by default and is still offered, because the choice between the closest match and the cleanest storm belongs to whoever has to defend the event. Two filters do drop candidates, and both are opt-in: *drop events with an embedded burst* and *drop anything flagged*. The rainfall cap is the exception — rainfall rarer than about the **AEP of the PMP** is the edge of the sampling scheme rather than a real event, so it is capped at 1.1 times that AEP, exactly as the old script does. The number is read from the IFD files config the storm config points at, which is where a Monte Carlo run gets it, and can be overridden on the page.

### Is it a simple event?

Closeness is not the only thing that makes an event usable. A representative event is usually wanted **simple** — one rise, one peak, one recession — because a double-peaked outflow makes a gate operation ambiguous and a dambreak run arguable. The mcdf records the peak and says nothing about the shape either side of it, so **Preview hydrographs** reads the run's stored hydrographs:

- it fills a **Shape** column for every candidate — ```single peaked```, or ```2 peaks``` — counting only the peaks worth arguing about: at least a fifth of the flood, and standing at least a tenth of it above the trough before them. The steps and shoulders every routed outflow has on its recession are not counted, because counting them would call almost everything multi-peaked and the column would say nothing;
- it draws the chosen event — inflow and outflow on the left axis, lake level dashed on the right;
- and after that, **clicking any candidate row** draws that one, so the shapes can be compared without choosing between them.

It is a button rather than something the page does by itself because the file is one column per simulation and tens of megabytes: it takes a few seconds the first time, and nothing after that. **This only works where the run stored its hydrographs** (```Store hydrographs``` = ```yes``` in the [simulation list](sim_list.md)); where it did not, the panel says so.

### The two marks on a level loading

The plot marks the loading with a line on each axis. For a **lake level** loading there are two honest answers to *where is this level on the result axis*, and both are drawn:

- the **design AEP** — where the frequency curve puts the level. This is the loading, and it is the one the ranking uses;
- **the level in this run** (the dotted line) — the AEP at which this database's own realisations reach it.

They rarely coincide, and the gap is not a rounding error. The design AEP is read off the **envelope** over the durations, as a straight line between the standard AEPs of the quantile table, from a curve that has been through the analysis; the dots are one duration's raw realisations. Where the two are more than about 0.05 in z apart the page says so in the loading's notes — *"This run reaches 220.50 at 1 in 1,340; the design curve puts it at 1 in 1,000"* — which is worth reading, because a large gap means the run being drawn is not the one the loading was quoted from. A level nothing in the run reached gets no second mark.

### Which run the event comes from

Leave **From** blank and the page takes the event from the duration that is **critical at that loading's AEP**, using the same envelope the Results page draws, and says so. This matters for lake level, where the critical duration is long at frequent AEPs and short on the rare tail: a list of loadings spanning the frequency range will legitimately draw its events from different runs. Naming a duration explicitly overrides it.

A **lake level** loading is converted to an AEP by reading it off the level frequency curve — the **envelope** over the durations, since that is the design curve the level was quoted from. A level above the top of the curve is reported rather than extrapolated, and the page then offers the events that got highest instead of the closest ones.

### What comes out

**Save** writes ```<group>_representative_events.json``` beside the databases: the loadings, the settings and the chosen simulation for each. It is reloaded when the page is next opened on that group, so the selection survives and can go into version control. There is a file per group because several groups commonly share one results folder — a GWL series usually does.

**Export csv** writes the same list as a table, one row per loading, with the metrics of the chosen event. The ```hydrograph``` column is the name of the column to pull out of the stored ```_inflows```, ```_levels``` and ```_outflows``` files: simulation 42 is ```sim_00042```.

### Extracting the hydrographs

The page chooses events; pulling out their hydrographs and plotting them is one command, shown on the page under *Extract the hydrographs and plot them* once the selection is saved. Press **Run** beside it: the page runs the command with Bryan's interpreter (set on the Project page), and shows each plot as it is written, any warning about the rebuilt hyetograph, and the log. The command itself is shown with absolute paths and Bryan's interpreter, so the **copy button** gives one that also works pasted into any console:

```bat
C:\...\env\python.exe -u C:\...\Bryan\util\RepresentativeEvents.py ^
    --config C:\...\sims_config.json ^
    --selection C:\...\sims_mc\results\GWL1p3_representative_events.json
```

It runs with **Bryan's** interpreter rather than the launcher's — it rebuilds the rainfall, which drives the storm generator. For each chosen event it writes a three-panel plot (the hyetograph drawn downwards from zero at the top, the inflow and outflow, and the lake level, all on one time axis measured from the start of the main burst) and a workbook holding every series. The hyetograph is **rebuilt** from what the run sampled and then checked against the depths the run recorded; a rebuild that does not match is said so rather than shown as fact. Events chosen from a **reservoir routing** row need one thing more: the storms belong to the run whose inflows were inherited, so either that run is in the same sims list, or its list is named with ```--source-sims-list```, or the routing row carries a ```Duration``` and a ```Focal subcatchments``` of its own. See [the utilities](utilities.md) for the detail.

### What it needs

A **monte carlo** or **reservoir routing** row that has been analysed, because the page reads the mcdf and the ```<type>_aep``` columns the analysis adds to it. This is the only page that reads a database rather than the quantile tables; there is no way around it, since a representative event *is* a realisation. Ensemble rows do not appear: an ensemble run has no realisations to choose between, having run every combination by design.

## Lake level frequency

The **Lake levels** page puts the dam's recorded annual maximum lake levels on a frequency axis beside the Monte Carlo design floods, and exports the report figure and the series.

### The record

Name the headwater level exports one per line, **earliest gauge first**. Both layouts work: Hydstra's own export, and the WMIP web export with its quality column. Where a gauge was replaced, each file owns the record from its first reading to the next file's first reading, and a jump of more than 50 mm at the handover is reported. A Hydstra export has **no quality code against any value** - the list of codes in its third column is only a legend - so nothing can be screened out of one, and the page says so.

For a dam whose full supply level or operation changed during the record, the recorded maxima are not one population. Homogenise them first and give the page the resulting **annual maximum CSV** instead; it reads the columns this page writes, and the ```WaterYear```/```Level_max```/```Level_max_at``` style (with or without ```New_```).

### Annual maxima

Water years are labelled by the year they end in and **start in the month you choose** (October by default). A maximum within the **carry-over window** of the start of its water year is the level the year opened at - left from the wet season before - and is drawn hollow, unless the lake fell away from it and came back to it after the window. Water years covered for less than the **minimum coverage** are flagged and can be left out. The **Annual maxima** table shows each year's maximum, when it fell, its day of the year and the interval between readings around it: a maximum read off once-a-day readings can miss the peak, and the page warns when any were.

### The curves

**Fit curves** runs ```util/LakeLevelFrequency.py``` with Bryan's interpreter - of the order of a minute for 400 resamples - and the curves appear when it finishes. Two forms are offered:

- **Shouldered plateau**: a polynomial shoulder up to full supply, a plateau at full supply over the extent the record gives, and a straight line through the maxima above it. The plateau starts at the most frequent maximum within the **plateau tolerance** of full supply and extends through steps smaller than the **plateau gap**. It cannot be fitted without maxima sitting on full supply and at least three above it, and the page says which is missing.
- **Logistic**, with its ceiling free - the fallback for a record without a plateau.

The shouldered form has four choices, and they should be made from how the dam works rather than from the rmse, which more terms always lower:

- **Plateau tolerance, 0 = none.** A plateau resting on one or two years is fragile: every resample that happens to leave those years out cannot be fitted. At Kroombit the 25 mm plateau rests on a single year and fewer than half the resamples fit; 100 mm, or no plateau, fits far more with the curve barely moved.
- **Degree below FSL** - the shoulder. Callide needs 4; 2 doubles the rmse.
- **Degree above FSL.** A straight line unless the dam spills in most years. A higher degree needs 4 maxima above the plateau per coefficient and is refused with a message otherwise - which a dry-belt dam like Callide or Kroombit will be, and a tropical one like Tinaroo need not.
- **Curve above FSL starts** free, leaving a step up from the plateau - right for a gated dam, where the step is the gates ceasing to hold the lake (0.57 m at Callide) - or at full supply, continuous, as an uncontrolled spillway gives.

When fewer than 70% of the resamples can be fitted the page warns that the band is likely too narrow.

Each is fitted to all the maxima and, separately, to the storm-driven ones placed on the whole record's probability scale. The 90% bands resample whole years and refit the same form. The shouldered curve and its band stop at the rarest maximum they were fitted to, because the upper limb is a straight line.

With **Compare with the Monte Carlo results** on, the level curves of the chosen durations of one group are drawn with their envelope. Use the present-day warming level and the lake configuration the record was measured under.

### Exports

**Export figure** writes ```<name>_validation.png``` (with the design floods) or ```<name>_record.png``` (without), each with a ```.json``` beside it recording the curve settings and the fit statistics they gave: 6.3 x 4.0 in at 300 dpi, for an A4 page, reading EY at the frequent end and "1 in X" past 50% AEP. **Export AMS CSV** writes ```<name>_ams.csv```, with the source files, water year and carry-over rule in ```#``` lines above the table.

### Choosing the water year

**Water year options** scores all twelve start months on the level record alone: the carried-over maxima (in all, and above a level you give), how close any storm-driven maximum comes to a year boundary, boundaries cutting through a rise or with the lake above full supply, maxima that are one only because of where the year was cut, and the empirical level at 50, 20 and 10% AEP. The month in use is highlighted. Below it, the median daily level by calendar month does not depend on the choice at all, and is the evidence for where the dry season bottoms out.

### What it keeps

Everything on the page is saved to ```lake_frequency.json``` beside the sims_config.json as it changes, with paths inside the project stored relative to it, so the analysis can be handed on with the project. Fitted results go to ```_lake_frequency/``` and are reused whenever the settings and the input files are unchanged.

## Report tables

The **Report** page fills the design flood report's result tables from the runs, so a revised rating means re-running the groups and copying the tables again rather than re-typing them.

### The study file

A report draws on several runs - the RFSL and FSL sims lists are separate sims_config.json files, and a sensitivity or a PMF ensemble is another - so the page works from a **study file** one level above them: ```bryan_study.json```, kept at the top of the study folder. Open it, or create it with **New study**. It holds:

- **Runs** - a name and a sims_config.json for each. Paths are stored relative to the study file, so a study copied to another disk still opens. Rename a run and every table that reads it follows; point the name at a re-run (E013 for E012) and every table moves with it.
- **Tables** - one entry per report table: its kind, the group or groups it is drawn from, and its layout options.

Everything is saved as it changes.

### The kinds of table

| Kind | Report table | What it computes |
|---|---|---|
| Design flood estimates | 1, 26-31 | One group. Per AEP: the peak lake level (the envelope over the durations) and its critical duration, with the peak inflow and outflow **of that same duration** - the report's "Peak inflow for lake level critical duration". Where no duration spills, the outflow is 0. The AEP of the PMP is a row of its own, labelled (PMPF). |
| AEP of given lake levels | 32 | Per group, the AEP at which each level you list is reached - the dam crest, each embankment crest - rounded to 10. By default the critical duration is taken from the design curve and the level is then read off **that duration's own realisations** (the mcdf), linear in log level against the standard normal variate; the alternative reads it off the design curve itself, between the standard AEPs. Table 1's dam crest flood row uses the same reading. |
| Peaks at one AEP | 33 | Per group, the lake level, inflow and outflow at one AEP - the AEP of the PMP for the PMPF. |
| Ensemble peak | 34 | Per ensemble group, the event with the highest lake level: its level, inflow, outflow and duration - the PMF. |
| Representative events | 35-36 | One section per group, from the events **saved on the Events page** for it: each loading's AEP (a level loading's AEP read as above), lake level, trigger, the chosen simulation and the duration it came from, in AEP order, with the PMF's event from an ensemble group last. A level loading takes its trigger from its own comment on the Events page, or from the table's list of named levels (DCF, each embankment crest). |
| Frequent levels | 37 | A grid: one column per group (the climate horizons), one row per frequency, in sections (RFSL, FSL). A standard AEP the quantile tables carry (1 in 2) is the design curve's value; a more frequent one (1 EY, 1 in 1.582) is read off each duration's realisations and the highest taken. The last column is the critical duration, as a range where the horizons differ. **Durations to consider** limits which runs count: Callide's Table 37 was made from 6-96 h, and a 120 h run changes the 1 in 2 level where it governs. |

The last three hold **sections** (the RFSL and FSL halves of a table, each under a heading row) of **rows**: a row either reads a group, or holds fixed values - for a previous study's numbers, such as the Sunwater 2020 baseline, which are not a run in this one.

**Where the inflow and outflow are read** is a choice for the last two kinds. The default is the storm that gives the peak level, consistent with the "Level critical duration" column and with Tables 1 and 26-31. The alternative is each result's own maximum, which is what the scripts that first filled Tables 33 and 34 did - the peak inflow of the PMPF then comes from a shorter storm than the peak level beside it, and differs from the PMPF row of Table 26.

### Copying

**Copy for Word** puts the table on the clipboard twice over: as a formatted table in the report's style (white bold header on the brand cyan, section rows on the cyan tint, Rubik Light 10 pt, units with a superscript), and as tab-separated text. Word pastes the first; Excel, or Word's *Keep Text Only*, takes the second. The caption is not copied - keep Word's own, so its numbering and cross-references survive. **Copy as text** gives the tab-separated form alone. A study can change the Word styling under a ```"word"``` key (```font```, ```size_pt```, ```header_fill```, ```header_text```, ```section_fill```, ```rule```).

Anything a table cannot fill - a group with no results yet, an AEP the run did not produce, a level above the top of the curve - is shown as a dash and listed above the preview, never invented.

## The PMF

The **PMF** page reads the PMF's ensemble runs and gives the PMF a notional AEP from the Monte Carlo realisations. It keeps its settings in the study file (open it on the Report page).

### Which event is the PMF

An ensemble run routes every temporal pattern of every duration, so a result has to be picked from the spread. The PMF is the **highest event** - the largest lake level of any pattern and duration, with that event's own inflow and outflow. Every other ensemble result is taken as the **median pattern**: for each duration the pattern at position round(n/2) of the ascending sort (the sixth of ten), then the duration whose median is largest - Bryan's own pick, the one in ```csv/<name>_critical.csv```. The page shows both, and the box plots show where they sit: the middle line of each box is that median pattern, the diamond is the highest event. The report's PMF table takes the highest event by default and the median on request.

### The notional AEP

The PMF has no AEP of its own, but the Monte Carlo realisations extend beyond the AEP of the PMP as far as the storm config extrapolates the rainfall, and the PMF level usually lies within them. The page reads it there:

- **Duration** - the realisations of the PMF event's own duration by default (the nearest one run), or any other.
- **Window** - the realisations from a lower AEP (1 in 500,000 by default) to the top of the sample. An upper bound can be set, but one below the PMF leaves the PMF above the window, and the page says so.
- **Degree** - of the polynomial in log10(level) fitted to the standard normal variate. A straight line (1) is the safe choice; a higher degree can turn over, and the page says when it does.

The chart shows the realisations, those in the window, the fitted curve, the PMF level, where it lands, and the AEP of the PMP. The grid below it gives the answer for three windows (0.4, 1 and 2 times the lower bound) by the three degrees. It is the sensitivity of the **fit** only: how far the realisations reach, and so where the PMF sits in them, is set by the rainfall extrapolation beyond the PMP and by the rating and storage curves, which no fit setting touches.

A straight read off the realisations without a fit is too noisy up there - a dozen or so events near the PMF level - which is why there is a fit at all. With fewer than 30 realisations in the window the page refuses rather than fits.

**Adopted** is the value carried into the report: the analyst's judgement across the horizons and the grid, typed in with a note of how it was reached. **Summary of all** estimates every PMF group with its current settings.

## Report figures

The **Figures** page replaces the plot-list workbook of ```PlotFrequencyCurves_v03.py```. Figures are kept in the study file (open it on the Report page) and name their curves by run and group, so a re-run moves every figure with it and nothing is retyped.

- **Curves**: a group's design curve (the envelope over its durations, as the Results page draws it); a curve from a file with the AEP ('1 in X') in its first column and the result as a column, for a previous study's adopted curve; or an **FFA** - an RMC Bestfit export, drawn as v03 drew it: the posterior mode (solid) and/or mean (dashed) in black, the 90% credible interval in faint grey when there is only one FFA on the figure, the annual maxima as dots and any paleoflood lower bounds as upward triangles.
- **Labels** belong to the figure. The same group can be 'URBS' on the FFA comparison and 'GWL 1.3 °C' among the horizons; edit a label in the box beside the preview and the figure redraws.
- **Reference levels** are one ```label = level``` or ```label = level = colour``` per line (full supply, the dam crest, an embankment crest). The **AEP of the PMP** draws a vertical line.
- **Duplicate** a figure and swap one curve for a sensitivity figure.

**Export PNG** runs ```util/ReportFigure.py``` with Bryan's interpreter and writes ```<name>.png``` to the figure's folder (```figures/``` beside the study file by default) at 6.3 x 4.0 in and 300 dpi, the size and style of the figures already in the reports, with ```<name>.json``` beside it holding every series and label it was drawn from. The preview is drawn from the same data, so what is checked on the page is what is in the PNG. **Export all** redraws every figure - after a re-run, the whole set in one go.

## Lake record

The **Lake record** page makes the antecedent storage the Monte Carlo runs sample - the `lake_config.json` files - from the dam's own lake level record, in three steps that feed each other, and a fourth, the **inflow record**, off the same record. Everything is kept in the study file (open it on the Report page), with paths relative to it, and each analysis leaves its job file beside its outputs so a run can be repeated from a console exactly as the page ran it.

### 1. Catchment rainfall

The daily catchment average of the AWAP / AWRA-L grids. Give the **catchment shapefile** (a field and value pick one catchment out of a regions file; blank takes every polygon) and the **folder of daily grids on this computer** as downloaded - one netCDF file per year, the rainfall a daily grid on latitude and longitude (```rain_day``` in AWRA-L; named if the file holds more than one). A shapefile not in latitude and longitude is reprojected from its ```.prj```.

**Area-weighted** counts each grid cell by the share of it the catchment covers, and by the cell's own area, which shrinks towards the pole. **Cell centres** counts every cell whose centre is in the catchment equally - a plain mask average, as the yearly ```Average_<year>_<mask>.csv``` files made by the mask-and-average tool are. A cell with no value on a day drops out of that day and the rest are reweighted.

The dates are the grids' own: a day D is the 24 hours to 9 am on D, which is what the antecedent search expects. The series is written as ```date,rain_mm``` with what produced it in ```#``` lines above.

**The series is what the study keeps.** It is written into the study folder (```lake_record/rainfall.csv``` by default) and ships with the model, with the catchment shapefile beside it. The folder of grids is **yours, not the study's** - it is saved in your own launcher settings, because the grids are tens of gigabytes not every user has and their path differs from one computer to the next. Someone who opens the study without them never needs them: the homogenisation and the antecedent storage read only the stored series.

On Callide, the area-weighted series from the AWRA-L grids reproduces the areal series the JPA project made (and the antecedent analysis was first run on) to 0.008 mm root mean square over all 42,126 days from 1911, on the same dates; the four lake configurations then come back unchanged.

This step runs in the launcher itself and needs ```netCDF4```, ```pyshp``` and ```pyproj``` in its environment (```ui/requirements-ui.txt```).

### 2. Homogenisation

The recorded lake levels re-routed through each **target rating**, so that a record made under several spillway configurations becomes one population. The net inflow is derived by closing the water balance backwards against the rating **in force at each step**, from the **rating register** (an xlsx: a ```Register``` sheet of ```Rating, from, to, FSL, ...``` and one ```level, flow``` sheet per rating), then routed through the target. Its inputs:

- **Gauge exports** (WMIP or Hydstra), in the order the gauges operated; each owns the record from its first reading to the next one's.
- An optional **overlay gauge** that replaces the chain wherever it reads below a level, its last value held until the chain climbs clear by the reconnect margin - Callide's intake gauge, reading the working storage below the sediment bar that partitions the pool at 200.90 m.
- The **storage table** (```.els```: ```EL, A, V```), the **SILO evaporation** and the monthly **pan factors**.
- **Target ratings**: a URBS ```.sq``` (storage above full supply against outflow, tied to the FSL its header declares - never re-based to another) or a ```level,flow``` csv starting at the FSL, each with a name and its FSL.

The record is routed at the gauge's own resolution, gaps longer than the longest step filled, and clipped to where the evaporation and the register both cover it; the page says where. The chart is each water year's recorded maximum against the homogenised ones.

### 3. Antecedent storage

For each water year's homogenised maximum, the storm that produced it: the rarest burst of 1-5 days in the **search window** before the peak, scored against the **IFD** (```duration_h``` against ```1 in X``` columns, with a ```1 in 2``` column) after the **restriction factors** that correct a fixed 9 am day to an unrestricted one. A year qualifies when a burst exceeds the **significance** fraction of its 1 in 2 depth. The homogenised lake volume is read at 9 am where the burst started (**burst**) and, walking back through days wetter than the **pre-burst edge**, where the pre-burst started (**storm**). Each series gets a Cunnane plotting position and a logistic S-curve in the standard normal variate, with the floor from the data (rounded) and the ceiling at the full supply volume, and is written as a sigmoid lake configuration - ```burst``` for a design storm simulated without its pre-burst, ```storm``` for one with it. Point the sims list's ```Lake config``` at them.

Run on Callide's inputs this reproduces the four lake configurations delivered for the design floods, byte for byte.

### 4. Inflow record

The inflow step 2 derives - the change in storage plus the release through the rating in force at each step - for an inflow flood frequency analysis and for calibrating the hydrologic model. It uses step 2's gauges, storage table and register, but none of its target ratings, and writes under ```lake_record/inflow```:

- ```inflow_ams.csv```: each water year's **peak inflow**, averaged over the **peak averaged over** window (1 hour by default - over a one-minute interval a millimetre of gauge is hundreds of m3/s of noise), with the release and level at the peak, and the largest inflow volume over each **burst duration**. Given the **catchment area**, each volume is also a runoff depth, and with step 1's rainfall the catchment rain over the same days (and the day before) is beside it. The page charts the two for the longest duration and warns of any year with more runoff than rain: that is the one check on the inflow that the inflow did not produce. Years covering less than 90% of the year are marked as part years.
- ```hydrographs/```: one file per event - the largest annual peaks (**hydrographs of the largest**, from **days before** to **days after** each peak) and any windows listed as ```name, start, end```, one a line.
- ```inflow_intervals.csv.gz```: every interval of the record.

Under the annual maximum chart, **Download AMS CSV** sends ```inflow_ams.csv``` to the browser, and **Copy AMS for Word** puts it on the clipboard as a table in the report style - water year, peak inflow and its time, and each burst's volume, runoff depth and rain - with the part years footnoted (**Copy as text** for Excel).

**The record between two dates** plots any stretch of the inflow record - inflow, release, level, and the uncorrected inflow dashed where the release is uncertain - with a zoom slider under the chart. A long stretch is drawn with each point the largest in its bin, so no peak is lost. Leave the **time step** blank for the record's own intervals (the inflow averaged over the peak window, as the peaks are), or give one - ```15min```, ```1h```, ```1D``` - for the mean over each step, stamped at the step's end as a model hydrograph is; each step then holds exactly the volume that entered in it, and ```Release_uncertain_share``` is the share of the step the recession correction touched. **Save CSV** writes the stretch to ```extracts/inflow_<start>_<end>[_<step>].csv``` beside the inflow record, where it stays with the study, and downloads it. The dates last used are kept in the study.

**Evaporation** is left out by default, so the inflow is the net inflow before lake losses and the record is not cut off where the evaporation file ends. Tick **keep the lake evaporation** to add it back, as step 2 does.

**Recessions.** Above full supply the falling limb often derives negative inflow: the gates released more than the rating says. The **recession correction** books that as release, which is right for the water balance, but leaves the corrected inflow at zero there where the real inflow was falling away. The actual release is not recorded, so neither version is the true recession. Each hydrograph therefore carries both, ```Inflow_m3s``` (corrected) and ```Inflow_uncorrected_m3s```, with ```Release_uncertain``` marking every interval the correction touched; the chart draws the uncorrected inflow dashed over those intervals. Calibrate to the rising limb and the peak with confidence, and to the flagged part of the recession with care. Peaks and burst volumes are on rising limbs and are unaffected.

On Callide, with the evaporation left out, the annual maxima agree with the independent reverse routing in callide-fsl-reinstate to a median of 0.00% on the peak and on every volume. The floods differ by 1-4% on volume, because that routing joins the rating rows with straight lines, which overstates a convex spillway rating between rows.

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
