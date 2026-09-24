# Bryan Technical Reference — additions for v2

Draft text for the next revision of `Bryan_Technical_Reference_v1.docx` (30 April 2026). It covers the methods added to Bryan since v1, and the v1 sections those methods supersede. Each part is headed with where it goes in the Word document. Section numbers follow v1; the proposed renumbering is in the table below.

Equations are written in LaTeX (`$$ … $$`) so that pandoc converts them to Word equations. Figures are not drafted; places where one would help are marked **[Figure: …]**.

| v1 | v2 | Change |
|---|---|---|
| 2.7 Antecedent storage | 2.7 | Add the logistic form fitted by Bryan, and a pointer to 7.4 |
| 5.4 Level pool routing | 5.4 | Unchanged |
| — | **5.5 Reservoir routing method** | New |
| 5.5 Output file management | 5.6 | Renumbered |
| 6.3 Critical durations | 6.3 | Add the definition of the envelope |
| 6.4 Representative events | 6.4 | Rewritten |
| — | **6.5 Flood volumes** | New |
| — | **6.6 AEP of a given lake level** | New |
| — | **6.7 Notional AEP of the PMF** | New |
| — | **6.8 Design flood tables** | New |
| — | **7 Historical record analyses** (7.1–7.5) | New chapter |
| 7 Output files | 8 | Renumbered |
| 8 Version control | 9 | Needs updating (see the end of this document) |
| 9 Ancillary scripts | 10 | 9.4 and 9.5 superseded in part (see 6.4, 6.6, 6.7) |

---

## 2.7 Antecedent storage — additional text

*Append after the paragraph on correlation.*

Bryan fits the sigmoid as a standard logistic in log-volume between fixed asymptotes, which is the $H = 1$ case of the five-parameter form:

$$
\log_{10} V(z) = \log_{10} V_f + \frac{\log_{10} V_c - \log_{10} V_f}{1 + e^{-k (z - z_0)}}
$$

With $H = 1$ the upper asymptote is exactly $V_c$, and the curve lies halfway between $V_f$ and $V_c$ in log space at $z = z_0$. Lake configurations from the original antecedent storage workbooks carry a hand-tuned $H > 1$. In those, $10^A$ rather than $V_c$ is the ceiling, and Bryan reads them in that convention so that they keep the curve they were published with.

The distribution can be derived from the dam's own lake level record, as described in Section 7.4.

---

## 5.5 Reservoir routing method (new)

The reservoir routing method re-routes inflow hydrographs stored by a previous Monte Carlo or ensemble run through a different dam rating, without re-running the hydrologic model. The inflows, the storms that produced them and the sampled antecedent storage are unchanged. Any difference in the routed results is therefore the rating's. It is used for rating sensitivity, for alternative operating conditions such as gates failing to open, and for testing an operating rule against a nominated antecedent storage. It is selected in the simulation list with `Method` = `reservoir routing`.

### Routing

The storage-indication form of the Modified Puls method is used. Over one timestep, continuity is

$$
\frac{(S_2 - S_1)\,1000}{\Delta t} = \frac{I_1 + I_2}{2} - \frac{O_1 + O_2}{2}
$$

with storage $S$ in ML, flows in m³/s and $\Delta t$ in seconds. Defining the storage indication

$$
\psi(S) = \frac{2000\,S}{\Delta t} + O(S)
$$

gives

$$
\psi(S_2) = \psi(S_1) - 2\,O_1 + I_1 + I_2
$$

$\psi$ is evaluated on a grid of 20,000 storages from zero to 5% above the larger of the storage table's and the rating's top storage. $S_2$ is found by interpolating the grid, for every realisation at once. Outflow $O(S)$ is zero at and below the full supply volume. Above it, $O(S)$ is interpolated linearly from the URBS storage–discharge (`.sq`) table, as URBS does. Lake level is interpolated linearly from the elevation–storage (`.els`) table. The timestep is that of the stored hydrographs.

The inversion needs $\psi$ to increase with storage. A flat run of outflow keeps it increasing, since the storage term does, and routes correctly; a gate held shut above full supply is an example. An outflow that falls as storage rises can make $\psi$ fall, and the interpolation would then return wrong storages without raising an error. A rating whose outflow falls with storage is therefore refused. Evaporation, rainfall on the lake and seepage are not modelled, as in the hydrologic model's own dam routing.

### Antecedent storage

For Monte Carlo input, the starting storage of each realisation is set by the `ADV source` column:

- `mcdf` (default): the ADV the source run sampled, from its database.
- `lake_z`: recomputed from the source run's sampled lake variate and a new lake configuration, so the same inflows can be routed under a different antecedent storage distribution.
- `lake_z correlated`: as `lake_z`, with the configuration's correlation layers applied first. This is only for databases sampled without correlation.
- `sims list`: one fixed volume for every realisation. The routed quantiles are then conditional on that volume, and are not design flood estimates.

For ensemble input, one volume is used for the whole run. It is a number, `fsv` or `mav`, resolved against the rating being routed rather than the source run's, or `database` to reproduce the source run.

### Analysis

A routed Monte Carlo run is re-analysed with the same Total Probability Theorem implementation as the Monte Carlo method (Section 6.1). The sample stratification is read from the source run's Monte Carlo config and checked against the database before the analysis. A routed ensemble run is re-analysed as an ensemble (Section 3). The routed database, quantile tables and hydrographs are written with an optional suffix, so one simulation list can route the same inflows through several ratings into one results folder.

---

## 6.3 Critical durations — additional text

*Append.*

The **design curve** of a result is its **envelope** over the storm durations. At each standard AEP it is the largest quantile of any duration, and the **critical duration** at that AEP is the duration that gives it. For lake level the critical duration usually changes across the frequency range. It is long at frequent AEPs, where the storage fills over several days, and short on the rare tail. The critical duration for lake level is generally not the critical duration for inflow at the same AEP. The design flood tables (Section 6.8) depend on keeping the two apart.

---

## 6.4 Representative events (rewritten)

*Replaces v1 Sections 6.4 and 9.5. `util/GetRepresentativeEvents.py` is superseded by `lib/RepresentativeEvents.py`, used by the Events page of the launcher and by `util/RepresentativeEvents.py`, which extracts and plots the chosen events.*

A design flood quantile is a statistic over thousands of realisations. Work downstream of it, such as gate operations, dambreak modelling and emergency action plan triggers, needs a single event: one storm, one antecedent storage, one set of hydrographs. A representative event is a realisation chosen from the Monte Carlo database to stand for a **loading**, which is either a design AEP or a lake level.

### Distance from the loading

Realisations are ranked on two axes at once:

1. how close the event's result is to the loading; and
2. how close the event's **rainfall** AEP is to the loading's AEP (**AEP neutrality**).

An event that reaches the 1 in 2,000 AEP lake level off 1 in 200 AEP rainfall got there through a coincidence, such as a full lake, a large pre-burst or a spike in the temporal pattern. It will not behave as a 1 in 2,000 AEP event when anything about the dam is changed.

Both distances are measured in standard normal variate space, where $z = \Phi^{-1}(1 - 1/X)$ for an AEP of 1 in $X$. A distance measured in 1 in $X$ is dominated by the rare end: at 1 in 2,000 AEP a difference of 200 is immaterial, while at 1 in 100 AEP it is not. The distance is

$$
\Delta z = \sqrt{\left(z_{\text{result}} - z_{\text{target}}\right)^2 + \left(z_{\text{rain}} - z_{\text{target}}\right)^2}
$$

where $z_{\text{result}}$ comes from the result AEP the Total Probability Theorem assigned to the event, and $z_{\text{rain}}$ from its sampled rainfall AEP. The v1 distance in 1 in $X$ is still reported, as `delta_aep`, so that a selection made with the earlier script can be checked.

Where an event's purpose is to put the lake at a particular level, the realisations can instead be ranked on the result alone: its difference from the loading in the result's own units (metres of level, m³/s of flow). Two settings keep that ranking from being decided by noise below the precision of the rating and the model:

- a **band** (default 20 mm of level), within which every event counts as reaching the loading and the events are ordered by $\Delta z$; and
- a **rounding** (default 10 mm), on which events outside the band are compared, with ties again broken by $\Delta z$.

A lake level loading is converted to an AEP by reading it off the design curve (Section 6.3). A design AEP loading is converted to a target result in the same way, so either ranking works for either kind of loading.

### Which run the event is taken from

By default the event is taken from the duration that is critical at the loading's AEP, found from the same envelope. For lake level, a set of loadings across the frequency range will therefore draw events from different durations.

### Checks on each candidate

Closeness alone does not make an event defensible. Each candidate is reported with the following, and flagged where it exceeds the default limit:

| Check | Flagged when |
|---|---|
| Embedded burst: the largest ratio of a sub-period's depth to the design depth for that sub-duration | ratio > 1.0 |
| Pre-burst proportion percentile | more than 0.25 from the median |
| Antecedent storage variate | $\lvert z \rvert > 1.0$ (about 1 in 6 either way) |
| Initial and continuing loss percentiles | more than 0.35 from the median |

Flagged events are still offered unless the user chooses to exclude them. The choice between the closest match and the cleanest storm belongs to whoever has to defend the event. Rainfall rarer than 1.1 times the AEP of the PMP is excluded, because it lies at the edge of the sampling scheme rather than being a plausible storm.

Where the run stored its hydrographs, the shape of each candidate can be checked as well. A peak is counted if it is at least one fifth of the flood's peak and stands at least one tenth of it above the preceding trough, which separates double-peaked events from the minor steps on a routed recession.

### What the event carries

The chosen event's inflow, outflow and level hydrographs, and its storm, are those of **one realisation**. The event is matched on the loading, usually a lake level, so its peak inflow and outflow will not in general equal the inflow and outflow quantiles reported for the same AEP in the design flood tables (Section 6.8). Those are quantiles of each result, ranked separately. Reports that pass representative events downstream should say so.

**[Figure: candidate events against the loading on the result–rainfall z plane, with the AEP-neutral diagonal.]**

---

## 6.5 Flood volumes (new)

With `Analyse volumes` set, the largest volume within a moving window of each analysis duration is computed for every realisation. It is then put through the same analysis as the peaks: the Total Probability Theorem for Monte Carlo runs, and the median pattern and critical duration for ensemble runs. The analysis durations are window lengths, not storm durations. The volume that governs storage is usually drawn from a longer window than the burst that produced it. The Monte Carlo method uses 24, 36, 48, 72, 96 and 120 hours.

For a hydrograph $Q_j$ at timestep $\Delta t$, the window volume ending at step $i$ is the trapezoidal integral over $D/\Delta t + 1$ samples, so the window spans the full duration $D$:

$$
V_D(i) = \frac{\Delta t}{1000} \left[\sum_{j=i-n}^{i} Q_j - \frac{Q_{i-n} + Q_i}{2}\right], \qquad n = D / \Delta t
$$

with $Q$ in m³/s, $\Delta t$ in seconds and $V_D$ in ML. The largest $V_D(i)$ over the hydrograph is the realisation's volume for that duration.

The sum is evaluated as a difference of cumulative sums, over all realisations at once. A duration that is not a whole number of timesteps, or that is longer than the hydrographs, is skipped and reported. Before 21 August 2026 the Monte Carlo window was one timestep short of its stated duration, which biased volumes low by about $\Delta t / D$.

The reservoir routing method analyses inflow volumes only. Outflow and storage volumes follow from the peak lake level and the rating, which are already analysed.

---

## 6.6 AEP of a given lake level (new)

*Supersedes the first part of v1 Section 9.4 (`DesignFloodInterpolation.py`).*

The AEP at which a nominated level is reached is needed for the dam crest flood, each embankment crest and emergency action plan triggers. It can be read in two ways:

- **Off the design curve.** The level is interpolated on the envelope (Section 6.3), linearly in log level against $z$ between the standard AEPs of the quantile table.
- **Off the critical duration's realisations** (default). The critical duration is taken from the design curve at the nearest standard AEP. The level is then read off that duration's realisations: the events are sorted by level, and $z$ is interpolated linearly in $\ln(\text{level})$.

The second uses thousands of events rather than a straight line between two standard AEPs a factor of 2 to 2.5 apart. At Callide Dam (E012, near-term, reduced full supply) they gave 1 in 25,350 and 1 in 27,230 for the dam crest. Where a level lies above every realisation it is reported as such, not extrapolated.

Frequent levels more frequent than the quantile tables go, such as 1 EY (1 in 1.582), are read off each duration's realisations, linearly in AEP. The highest over the durations is taken.

---

## 6.7 Notional AEP of the PMF (new)

*Supersedes the last paragraph of v1 Section 9.4 (`AEPofPMF.py`).*

### Which ensemble event is the PMF

An ensemble run routes every temporal pattern of every duration, so a single result has to be chosen from the spread. Two conventions are used:

- **Highest event:** the largest lake level of any pattern and duration, with that event's own inflow and outflow. This is taken as the PMF.
- **Median pattern:** for each duration, the pattern at position $\mathrm{round}(n/2)$ of the ascending sort (the sixth of ten), then the duration whose median is largest. This is what Bryan writes as the ensemble's critical result, and it is used for every ensemble result other than the PMF.

At Callide Dam the two differed by 0.14–0.28 m of PMF level (E012, reduced full supply).

### The notional AEP

The PMF has no AEP of its own, but a risk assessment needs one. The Monte Carlo realisations extend beyond the AEP of the PMP as far as the extreme rainfall interpolation (Section 4.1.2) extends the rainfall, and the PMF level usually lies within them. A direct read off the realisations is too noisy there, with only a dozen or so events near the PMF level, so a curve is fitted:

$$
z = \sum_{j=0}^{d} a_j \left(\log_{10} L\right)^j
$$

The fit is made to the realisations of one duration whose result AEP lies within a window at the top of the sample. The PMF level $L_{\text{PMF}}$ is then placed on it. The defaults are:

- **duration:** the PMF event's own duration, or the nearest one run;
- **window:** from 1 in 500,000 to the rarest realisation;
- **degree:** $d = 1$. A higher degree can turn over, and is flagged where it does inside the range it is read on.

At least 30 realisations are needed in the window. The result is reported with a grid of windows (0.4, 1 and 2 times the lower bound) against degrees 1 to 3. That grid is the sensitivity of the fit only. The position of the PMF within the realisations is set by the extreme rainfall interpolation and by the rating and storage curves, which no fit setting changes. The value adopted for a study is the analyst's judgement across the horizons and the grid, recorded with a note of how it was reached.

---

## 6.8 Design flood tables (new)

The design flood tables of a report are built from the analysed runs by the Report page of the launcher, which records in a study file which group of which run fills each table. The conventions below are what the tables mean.

### Flows at the level's critical duration

For each standard AEP, the design flood table reports:

- the **peak lake level**: the envelope over the durations (Section 6.3);
- the **level critical duration**: the duration giving that envelope;
- the **peak inflow** and **peak outflow**: the inflow and outflow **quantiles of the level critical duration** at that AEP.

The inflow is not the design peak inflow at that AEP, which is the envelope of the inflow quantiles and usually comes from a shorter duration. Tables carry the footnote *"Peak inflow for lake level critical duration"*. Where the critical duration produces no spill at an AEP, the outflow is reported as zero.

The three values are quantiles, each ranked separately within the critical duration. The inflow at an AEP is not the inflow of the realisations that produced the level at that AEP; the antecedent storage and the temporal pattern vary between realisations, so the two sets of storms differ. This is the intended convention. It follows that a representative event's flows differ from the table's (Section 6.4).

The table of results at the AEP of the PMP reads its flows the same way by default. The PMF row takes the level, inflow and outflow of the highest ensemble event, all from one storm (Section 6.7). A dam crest flood row takes its AEP as described in Section 6.6.

---

## 7 Historical record analyses (new chapter)

### 7.1 Overview

Bryan derives three things from a dam's recorded lake levels:

1. the **antecedent storage distribution** the Monte Carlo runs sample (Section 2.7);
2. a **homogenised annual maximum lake level series**, for comparison with the design floods; and
3. the **inflow record**: an annual maximum series of peak inflow and burst volume for flood frequency analysis, and event hydrographs for calibration.

All three rest on one reverse routing of the level record (Section 7.3). The analyses were first written for Callide Dam (the `callide-fsl-reinstate` project) and were generalised into Bryan (`lib/homogenise`, `lib/antecedent`) without changing the method. Run on Callide's inputs, Bryan reproduces that project's homogenised record and its four lake configurations exactly. Each analysis is driven by a job file that records every input and setting, and writes the job beside its outputs so that a run can be repeated from the console.

**[Figure: flow of the four steps — catchment rainfall, homogenisation, antecedent storage, inflow record — and their inputs.]**

### 7.2 Catchment rainfall

The antecedent storage analysis needs a daily catchment-average rainfall series. It is computed from the AWAP / AWRA-L daily rainfall grids (one netCDF file per year, 0.05° cells) and a catchment polygon. Two weightings are available:

- **Area-weighted** (default): each cell is weighted by the fraction of it inside the catchment, and by its own area, which scales with $\cos(\text{latitude})$. The fraction is found by testing a 10 × 10 mesh of points in the cell against the polygon.
- **Cell centres:** every cell whose centre lies in the catchment is weighted equally, as a plain mask average.

A cell with no value on a day is left out of that day, and the remaining weights are renormalised. Dates are the grids' own: day $D$ is the 24 hours to 9 am on $D$. The series is stored with the study, and the grids themselves are not needed again.

At Callide Dam the area-weighted series reproduced the areal series previously prepared for the joint probability analysis to 0.008 mm root mean square over 42,126 days.

### 7.3 Reverse routing and homogenisation

A lake level record made under several spillway configurations, such as a change of full supply level or gate operation, contains several populations. Homogenisation puts the record onto one configuration in two stages. It first derives the net inflow by closing the water balance backwards against the rating in force at each step. It then routes that inflow forward through a single **target rating**.

#### Inputs

- **Level record.** Gauge exports (WMIP or Hydstra) spliced in operating order, each owning the record from its first reading to the next gauge's. An optional overlay gauge replaces the record wherever it reads below a nominated level, for a pool that partitions at low storage. At Callide the intake gauge reads the working storage below a sediment bar at 200.90 m.
- **Storage table.** Elevation, area and volume (`.els`), interpolated with a shape-preserving monotone cubic (PCHIP). Linear interpolation would give a piecewise-constant $dV/dL$, which puts a step into the derived inflow each time the lake crosses a table row. An unconstrained cubic spline can overshoot a rating table badly.
- **Rating register.** Each historical rating as a level–discharge table, with the dates it was in force and its full supply level. Release is taken as zero at and below the full supply level in force. Ratings are interpolated with PCHIP.
- **Evaporation.** SILO Data Drill pan evaporation with monthly pan factors. The SILO value in row $D$ is the total over 09:00 $D$ to 09:00 $D+1$. It is applied as a step function on 9 am boundaries and integrated exactly over each timestep, so a day's steps sum to that day's total however the day is subdivided.

#### Timestep

Both stages run at the **native resolution of the gauge**. Every reading is kept, and gaps longer than a maximum step (default 1 hour) are filled by interpolating the level linearly, with the filled points flagged. Resampling onto a fixed grid would straddle and under-read observed peaks. A daily timestep would smear a rise of a few hours across a day. At Callide, deriving the inflow from daily maxima attenuated the peak by a median factor of 2.3 in the years that reach full supply. The record is clipped to the period covered by both the evaporation and the rating register.

#### Stage 1: net inflow

For each interval $i$ from $t_{i-1}$ to $t_i$:

$$
I_i = \left[V(L_i) - V(L_{i-1})\right] + \frac{Q_h(L_{i-1}) + Q_h(L_i)}{2}\cdot\frac{\Delta t_i}{1000} + E_i\,\frac{A(L_{i-1}) + A(L_i)}{200}
$$

where:

| Symbol | Meaning |
|---|---|
| $I_i$ | net inflow over the interval (ML) |
| $V(L)$, $A(L)$ | storage (ML) and surface area (ha) at level $L$ |
| $Q_h(L)$ | release (m³/s) from the historical rating in force at the time of each reading |
| $\Delta t_i$ | interval length (s) |
| $E_i$ | open-water evaporation over the interval (mm): pan evaporation × pan factor |

The inflow is **net**. It absorbs rainfall on the lake, seepage and extractions, none of which are modelled separately, so it can be negative between events.

**Recession correction.** Above full supply the balance asserts the rating in force. Where that rating understates the actual release, typically because gates were opened further than the rating assumes on the falling limb, the derived inflow goes negative. That deficit is a release carrying an inflow's sign. Left in the inflow, stage 2 would charge the lake for it twice: once inside the inflow and once through the target rating. So wherever the level is above the full supply level in force and $I_i < 0$, the deficit is moved from inflow to release:

$$
R_i \leftarrow R_i - I_i, \qquad I_i \leftarrow 0
$$

The sum of the terms is unchanged, so the recorded level is still reproduced exactly. Below full supply the rating asserts nothing, and a negative net inflow there, with seepage and extraction exceeding inflow, is physical and is left alone. At Callide, 1,489 of 508,179 steps were corrected, moving 59,626 ML from inflow to release.

#### Stage 2: routing through the target rating

The net inflow is routed through the target rating with the trapezoidal level-pool equation, solved **implicitly for level** at each step:

$$
V(L_i) + \frac{\Delta t_i}{2000}\,Q_t(L_i) + \frac{E_i}{200}\,A(L_i) \;=\; V(L_{i-1}) + I_i - \frac{\Delta t_i}{2000}\,Q_t(L_{i-1}) - \frac{E_i}{200}\,A(L_{i-1})
$$

where $Q_t$ is the target rating. Storage, area and release all increase with level, so the left-hand side is monotone and has a unique root. It is found by safeguarded Newton iteration: a Newton step where it stays inside the bracket, bisection where it does not. An explicit scheme is unstable against a steep gated rating at any practical timestep. At Callide the URBS gate operation releases 142 m³/s within 120 ML (10 mm) of full supply. Solving for level keeps the scheme self-consistent: when the target rating equals the rating that was in force, the routing reproduces the recorded level to solver tolerance.

A target rating is either a URBS storage–discharge table (`.sq`), tied to the full supply level its header declares and interpolated linearly in storage as URBS does, or a level–discharge table starting at full supply. The output is the homogenised level, storage and release at every step. The annual maxima are then taken by water year (default start 1 October), with years covering less than 90% of their length flagged rather than dropped.

**[Figure: recorded against homogenised annual maximum lake levels.]**

### 7.4 Antecedent storage from the record

The antecedent storage distribution is derived from the homogenised record, following Section 8.3 of the Callide Dam design flood hydrology report (after Criste-Jones et al., 2025). For each water year's homogenised maximum:

1. **Find the burst.** Within a window of 30 days before the peak, the largest forward accumulation of catchment rainfall is formed for each duration of 1 to 5 days. Each is multiplied by a restriction factor that corrects a fixed 9 am day to an unrestricted window. The defaults, 1.15, 1.11, 1.07, 1.05 and 1.04 for 1 to 5 days, are carried from the earlier assessment. Each is given an AEP by interpolating the catchment IFD curve for that duration linearly in log depth against log AEP. The search window extends past the peak by the longest duration so that a burst still in progress at the peak forms an accumulation. The burst must, however, start no later than the day after the peak, because a rain day is labelled by the 9 am at which it ends.
2. **Test its significance.** The year yields a sample only if the scaled accumulation for some duration exceeds 0.8 of the 1 in 2 AEP depth for that duration. The **critical burst** is the one with the rarest AEP.
3. **Read the antecedent volume.** The homogenised volume is interpolated at 09:00 on each day. The **burst** volume is read at the start of the critical burst's first rain day. The **storm** volume is read where the pre-burst began: walking back from the burst start while the daily rainfall exceeds 10 mm, then taking the first drier day. The burst volume suits design storms simulated without their pre-burst rainfall; the storm volume suits storms with the pre-burst prepended.
4. **Fit the distribution.** Each series is given Cunnane plotting positions, $p = (r - 0.4)/(n + 0.2)$ with the largest volume ranked 1, and $z = \Phi^{-1}(1 - p)$. The logistic curve of Section 2.7 is fitted with the floor $V_f$ set to the smallest volume rounded down to 1,000 ML and the ceiling $V_c$ to the full supply volume of the target rating. The slope $k$ and centre $z_0$ are then fitted by least squares in $\log_{10} V$.

Each fit is written as a sigmoid lake configuration, one per basis (burst or storm) and target rating, ready for the simulation list's `Lake config`. All of the settings above can be changed in the job file.

### 7.5 Inflow record

The inflow record is the stage 1 net inflow of Section 7.3, put into the forms that flood frequency analysis and model calibration need. It uses the same gauges, storage table and rating register, but no target rating. Evaporation is left out by default, so the inflow is the net inflow before lake losses and the record is not limited to the period of the evaporation file.

#### Intervals

Each stage 1 step becomes an interval with volume $I_i$ (ML). The inflow derived over an interval is its **mean over the interval**, so it belongs at the interval's midpoint and is stamped there. Summed, the interval volumes reproduce the storage change plus the release to the megalitre.

#### Peak inflow: averaged off the cumulative volume

Differencing storage over a short interval turns gauge resolution into flow. At a surface area of 1,234 ha, a 1 mm reading step is 12.3 ML. Over a one-minute interval that is about 200 m³/s; over an hour, about 3 m³/s. The native-interval inflow therefore cannot be used for peaks. Several annual maxima at native resolution are set by one- to five-minute gauge steps rather than by floods.

The peak inflow is instead the mean inflow over a window $T$ (default 1 hour) centred on each interval. It is read off the **cumulative inflow volume**:

$$
W(t) = \sum_{t_i \le t} I_i \qquad \text{(interpolated linearly between interval ends)}
$$

$$
\bar I(t) = \frac{1000\,\left[W(t + T/2) - W(t - T/2)\right]}{T}
$$

with $W$ in ML, $T$ in seconds and $\bar I$ in m³/s. The annual peak inflow is the largest $\bar I$ at an interval midpoint in each water year.

This is the exact time average of the inflow. Averaging the interval values instead, even weighted by interval length, would give a one-minute gauge step the same standing as the six-hour interval beside it, because the neighbouring intervals' midpoints fall outside the window. Reading the average off the cumulative volume counts each interval for exactly the part of it inside the window.

#### Burst volumes

The largest inflow volume over each duration $D$ (default 24, 36, 48 and 72 hours) is read off the same cumulative volume, on a 15-minute grid:

$$
V_D = \max_t \left[W(t) - W(t - D)\right]
$$

by water year. A volume over 24 hours or more is already an integral, so it does not depend on the averaging window. Given the catchment area, each volume is also expressed as a runoff depth (1 ML over 1 km² is 1 mm). With the catchment rainfall of Section 7.2, the rainfall over the days the window spans, plus the day before, is reported beside it. A burst implying more runoff than rainfall is wrong however well the water balance closes. It is the one check on the inflow that the inflow did not produce.

#### Recessions

The recession correction of Section 7.3 is right for the water balance, and it does not affect peaks or burst volumes, which sit on rising limbs. It does, however, leave the corrected inflow at zero on a recession where the real inflow was falling gradually. The actual release is not recorded, so neither version is the true recession. Event hydrographs therefore carry both the corrected and the uncorrected inflow, and flag every interval the correction touched. The rising limb and the peak can be used for calibration with confidence; the flagged part of the recession should be used with care.

#### Limitations

- Release is a single-valued function of level. Where gates were operated differently from the rating, the error passes directly into the derived inflow. On a recession it shows as the negative inflow the correction removes; on a rising limb it would bias the inflow with the opposite sign and would not be detected.
- The gauge is assumed to read the mean pool level. Wind setup or drawdown near an operating outlet appears as inflow.
- Where the record was read coarsely (for example daily), the peak is attenuated by an amount that depends on the sampling interval at the time. The median and longest gap in each year are reported so that this can be judged.

#### Verification

On Callide Dam, with evaporation left out, the annual maxima agree with the independent reverse routing in `callide-fsl-reinstate`, which shares no code with it, to a median of 0.00% on the peak and on every burst volume over 56 water years. Individual floods differ by 1–4% in volume. The independent implementation interpolates the ratings linearly between published rows, which overstates a convex spillway rating between rows.

---

## Other v1 sections needing attention

These were noticed while drafting and are not addressed above:

- **8 Version control** describes Tortoise SVN and the 2025-05-Tinaroo-1538 release. Bryan has since moved to Git, and the development features listed there (baseflow, dam routing) are in the main line.
- **4.5 Pre-burst filtering:** Figure 4 is an unresolved cross-reference ("Error! Bookmark not defined"), and an empty sub-heading follows the section.
- **Introduction and Executive summary** cite the Design Hydrology Specification as DS PRO 030 (Sunwater, 2025a), but the reference list gives 2025a as DS PRO 029. The climate change procedure (2025b) is also cited as DS PRO 029.
- **9.4 AEP interpolation/extrapolation** and **9.5 Representative events** describe scripts superseded by Sections 6.4, 6.6 and 6.7. They could be reduced to a note of what replaced them.
- **7 Output files** does not list the flood volume outputs (Section 6.5), the routed outputs of the reservoir routing method (Section 5.5), or the outputs of Chapter 7.
- **References** to add: Criste-Jones et al. (2025); Cunnane (1978); SILO (Jeffrey et al., 2001); AWAP (Jones et al., 2009) / AWRA-L (Frost et al., 2018).
- Typo in 8: "Byran".
