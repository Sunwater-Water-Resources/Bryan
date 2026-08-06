# The sub-burst check: should embedded bursts be filtered?

## The problem

A design storm is constructed by scaling a temporal pattern to the burst depth sampled for the storm's duration and AEP. Temporal patterns taken from observed storms often contain *embedded bursts*: a window within the pattern whose depth, over some shorter duration, is rarer than the AEP of the design storm itself. For example, a 1 in 2,000 AEP 24-hour storm might contain a 2-hour window whose depth exceeds the 1 in 2,000 AEP 2-hour IFD depth.

Current Australian practice leans towards filtering these embedded bursts out (see the embedded burst filter and the ```ebf``` exclusion key in the [simulation list](sim_list.md)). The argument for filtering is consistency: a storm assigned an AEP should not contain within it a rarer event, otherwise the flood response can be driven by the embedded burst and the derived flood frequency curve is biased upwards. The argument against filtering is realism: observed storms genuinely contain embedded bursts, and flattening them removes rainfall variability that real storms possess, potentially biasing the flood estimates downwards.

For a Monte Carlo scheme, there is a third view that reframes the question. Any single realisation containing an embedded burst is not, by itself, a problem — real storms have them. What matters is whether embedded bursts occur *too often across the ensemble*: the frequency of sub-bursts of a given duration and depth, aggregated over all realisations by the total probability theorem (TPT), should not exceed what the IFD statistics for that duration imply. This is an *ensemble AEP-neutrality* condition, and it is testable.

## The method

For every Monte Carlo realisation, Bryan records the maximum rolling-window depth of the main burst temporal pattern for each standard sub-duration shorter than the storm duration (the *sub-burst* depths), together with the catchment-average IFD depth for that sub-duration at the same sampled z (the same reference curve used by the embedded burst filter, including areal reduction and climate adjustment). Both are stored in millimetres in the mcdf file as ```subburst_<duration>h``` and ```ifd_<duration>h``` columns.

At the analysis stage, the TPT is applied to the sub-burst depths for each sub-duration — exactly as it is applied to peak inflows — producing a derived frequency curve of sub-burst depth. This curve is plotted against the IFD reference, and the ratio of the two (the *neutrality margin*) is tabulated at the standard AEPs.

### Interpreting the result

The check is **one-sided**:

- **Sub-burst curve above the IFD** (margin > 1.0): sub-bursts of that duration occur more often in the simulated storms than the IFD statistics imply. This is direct evidence that the unfiltered patterns over-generate embedded bursts, and justifies filtering (or re-weighting the offending patterns).
- **Sub-burst curve below the IFD** (margin < 1.0): the choice not to filter is *not contradicted*. It is not conclusively proven either, because sub-bursts within storms of one parent duration are only one of the mechanisms contributing to the sub-duration IFD — some margin below the IFD is expected. How much margin is a matter of judgement.

Because storms are simulated one duration at a time, the check applies per parent duration. If the sub-burst curve for a single parent duration comes close to the IFD, remember that the other parent durations also contribute sub-bursts of the same sub-duration — the combined frequency will be higher than any single curve shows.

The check is also informative when run with the filter **on**: if the filtered sub-burst curve sits far below the IFD while the unfiltered curve sits near it, the filter is removing variability that the rainfall statistics could accommodate. Running a pair of simulations (filtered and unfiltered, using the ```ebf``` exclusion key) brackets the question and, together with the resulting flood frequency curves, shows whether the choice matters for the catchment at all.

## Implementation

### Recording (during model runs)

Sub-burst tracking is always on for Monte Carlo runs — its cost is negligible. For each realisation:

- With filtering on, the sub-bursts of the *filtered* pattern are measured (so the check tests the storms actually simulated).
- With filtering excluded, the measurement doubles as the embedded burst check that populates the ```embedded_bursts``` comment column.
- Measurement happens **before the pre-burst is prepended** — the check currently applies to the main burst only (see limitations below).
- Sub-durations are the standard IFD durations shorter than the storm duration whose window is a whole number of pattern timesteps.

### Cheap runs without the hydrologic model

Because the sub-burst measurements depend only on the sampled storms (not the flood response), the check can be run without URBS/RORB: set ```Run models``` in the [simulation list](sim_list.md) to *storms only*. The full sampling and storm generation runs and the mcdf is written, but the hydrologic model is skipped — a Monte Carlo "run" then takes seconds to minutes instead of hours. This is the recommended way to explore the filtering question and to calibrate pattern weights (below) before committing to full model runs.

### Analysis (after model runs)

Set the optional ```Analyse sub-bursts``` key in the [simulation list](sim_list.md) to *yes*. As with ```Analyse results```, the models do not need to be rerun — but the mcdf file must have been produced by a version of Bryan that records the sub-burst columns. For each sub-duration this produces:

| Output | Description |
| ----------- | ----------- |
| ```<output>_subburst<duration>h.csv``` | The TPT quantiles of sub-burst depth at the standard AEPs. |
| ```<output>_subburst<duration>h.png``` | Plot of the sub-burst TPT curve (solid), the individual realisations (scatter), and the same-z IFD reference (dashed, one branch per storm method). |
| ```<output>_subburst_neutrality.csv``` | Summary table: the neutrality margin (TPT sub-burst depth ÷ same-z IFD depth) at the standard AEPs, one column per sub-duration. Values above 1.0 indicate sub-bursts occur more often than the IFD implies. |

The IFD reference appears as separate branches for each storm method (ARR point/areal, GSDM, GTSMR) because the reference depth at a given z depends on which method generated the storm; the branches overlap in the AEP changeover zone. The neutrality margin adopts the lowest branch (conservative) and interpolates across the gap between branches.

### Reading the outputs — practical notes

- The TPT sub-burst curve is edge-affected near the sampling bounds, just like the flood frequency curves: it flattens towards ```lower_aep``` and the margin is not computed beyond the sampled range. Read the margins a little inside the sampled AEP range.
- The tail of the sub-burst curve is *pattern-limited*: within each sampling division only ten temporal patterns are available, so the curve reflects the embedded burst content of those ten patterns scaled by the depth distribution. A single unusual pattern can dominate the curve — the realisation scatter on the plot will show this.
- Do the comparison in the same climate scenario throughout; the recorded IFD reference already includes the climate adjustment applied to the run.

## If neutrality is breached: calibrated pattern weights

If the unfiltered ensemble breaches neutrality, the blunt remedy is the embedded burst filter. A softer alternative is to keep the observed patterns intact but assign them unequal probabilities: down-weight the patterns whose embedded bursts drive the breach until the ensemble satisfies neutrality. The ```util/CalibrateTpWeights.py``` script does this calibration — see [the utilities page](utilities.md). Because the sub-burst depths depend only on the sampled pattern, depth, and storm method (all recorded in the mcdf), the calibration is entirely offline: trial weights are evaluated by re-weighting the recorded realisations in the TPT, with no model reruns. The script also re-computes the flood quantiles under the calibrated weights (a ```tp_w``` column in the mcdf is recognised by the TPT analysis), previewing the effect on the flood frequency curve from the existing model runs.

The weights inherit the scope of the measurements they are calibrated against. Those are main-burst sub-burst depths, so **weights that achieve neutrality say nothing about the pre-burst-inclusive storm** — see the limitations below.

Notes on the weights:

- Weights are grouped by storm method (and by frequency bin for ARR point patterns), and calibrated per storm duration.
- A floor is applied (default 0.02) so that no pattern is silently excluded. If the margins cannot be brought below the target with the offending patterns at the floor, weighting alone cannot achieve neutrality for that simulation — the script says so, and the choice reverts to filtering or pattern review.
- The calibration stops as soon as neutrality is achieved, so patterns are not down-weighted further than the condition requires.

### Applying the weights to the sampling

For a production run, the calibrated weights can be applied to the pattern sampling itself: set the optional ```TP weights``` key in the [simulation list](sim_list.md) to the weights file written by the calibration script. The patterns are then sampled with the calibrated probabilities (composed with the D50 climate shift weightings), and the weight applied to each realisation is recorded in a ```tp_weight``` column of the mcdf for audit. Three things to be aware of:

- **The standard TPT analysis is already correct for a weighted-sampling run** — the sampled frequencies embody the weights, which *are* the asserted pattern probability model. Do **not** also attach a ```tp_w``` column to such a run's mcdf (that would double-count the weights).
- Rerunning the sub-burst check on the weighted-sampling run verifies the neutrality of the ensemble as actually simulated — a worthwhile closing check, since the sampled realisations will differ from the calibration run's.
- The ```TP weights``` key cannot be combined with replication of the pattern sampling (the ```tp``` replicate), because the replicated patterns were sampled under different weights. Since heavily down-weighted patterns are rarely sampled, their response is also thinly resolved — another reason the weight floor matters.

## The separate question: embedded bursts above the PMP

Everything above is a *frequency* argument. A sub-burst is acceptable if it occurs no more often across the ensemble than the sub-duration IFD implies, and if it occurs too often the remedy is to re-weight the offending patterns. That argument needs a reference curve to hold it up, and the reference curve terminates at the PMP. Above the PMP for a sub-duration there is no assignable frequency — the depth is outside the domain, not merely rare within it — so there is no target for the weight calibration to hit and the only admissible weight is zero. Filtering there is not a contradiction of the unfiltered approach, it is what the approach reduces to at the physical bound.

The admissible bound on a sub-burst of duration *s* within a burst of duration *d* is

```
bound(s) = PMP_s × max(1, burst_depth / PMP_d)
```

Where the burst depth is below the PMP this is the absolute physical bound, ```PMP_s```. Where the sampling extrapolates past the AEP of the PMP (the default ```pmp_exceedance_sampling``` of ```extrapolate```) the absolute bound would be infeasible: as *s* approaches *d* the rolling window tends to the whole burst depth, which already exceeds ```PMP_d```. The multiplier scales the bound by the factor by which the realisation exceeds the envelope, so the constraint becomes one on the storm's *shape* — its internal depth-duration profile may not be peakier than the PMP depth-duration envelope — while the magnitude is left free, which is what the extrapolation already asserts at the parent duration. One expression covers both settings of ```pmp_exceedance_sampling```, and it reads the realised burst depth rather than the setting.

### The screen

Bryan does not filter against this bound. It screens for it, in ```screen_pmp_embedded_bursts```, and reports near the top of the simulation log. The screen runs before any sampling because the breach condition does not depend on it. Like the neutrality check, it measures the **main burst only**: ```r_s``` is taken from the main burst temporal pattern as it comes from the pattern file, before any pre-burst is prepended. For a pattern holding a fraction ```r_s``` of its depth in the worst window of sub-duration *s*:

- the bound engages once the burst depth exceeds ```D_crit = PMP_s / r_s```;
- if ```D_crit >= PMP_d```, equivalently ```r_s <= PMP_s / PMP_d```, the pattern never breaches at any depth;
- between ```D_crit``` and ```PMP_d``` the exceedance ratio grows linearly with depth;
- above ```PMP_d``` it is **exactly constant** at ```PMP_d / D_crit```, because the sub-burst and the bound scale together.

So the worst case is fully realised at the PMP anchor and the extrapolation adds nothing to it. Severity is bounded, and bounded by a number computed without running anything.

The binding sub-duration is the one with the largest ```r_s / PMP_s```, which is both the first to engage and the most severe. The report gives it, along with ```r_s```, ```D_crit```, the AEP at which the burst depth reaches ```D_crit```, the exceedance ratio at and above the anchor, and the share of that pattern set's sampled realisations above ```D_crit```.

### Climate change

All the depths above are the climate-adjusted ones the run actually simulates. The rainfall uplift scales a duration's whole frequency curve by the same factor at every AEP, so the PMP anchor is uplifted with it, and ```PMP_s``` and ```PMP_d``` each carry their own duration's factor. The uplift factors used are printed above the table.

Two consequences follow, and they pull in opposite directions to what a first reading suggests:

- **A uniform uplift changes nothing.** If every duration took the same factor, ```D_crit``` and the sampled depths would scale together and both the exceedance ratio and the engagement AEP would be untouched. The same is true of the parent duration on its own: the burst depth reaches the PMP at the same AEP whatever the warming level, because both sides move by the same factor.
- **The uplift is duration dependent**, currently 15%/°C at 1 hour easing to 8%/°C at 24 hours and beyond, so a sub-duration is generally uplifted *more* than its parent. The exceedance ratio becomes ```r_s × PMP_d·u(d) / (PMP_s·u(s))```, and with ```u(s) > u(d)``` warming **eases** the bound rather than tightening it. ```r_s``` itself is a fixed property of the pattern and never moves; only which sub-duration binds can change.

So a climate change run is not the adverse case for this check. Where the parent duration and the binding sub-duration are both 24 hours or longer the factors are equal and nothing moves at all.

ARR areal and point patterns are screened as well as GSDM/GTSMR, each over the AEP window it is actually sampled in. Because the MC scheme samples GSDM/GTSMR patterns above 1 in 2,000, and above that the depth cancels out of the breach condition, most of the screen is a fixed property of the PMP patterns against the PMP depth-duration curve. For a small catchment the ARR patterns cannot get near the PMP and the screen is inert for them. The case worth watching is a catchment large enough that the AEP of the PMP is only half a decade beyond the changeover — at 10⁵ km² it is 1 in 10,000 — where the top of the ARR range does approach the PMP.

### Reading the result

- **Nothing breaches.** The exception does not arise for this catchment and no filtering against the PMP is needed. Worth recording: "we checked and it cannot happen here" is a better answer than an untested safeguard.
- **Breaches only beyond the sampled range.** Reported as a ```CHECK```. No action for that run, but it will matter if the sampling is extended.
- **Breaches inside the sampled range.** Reported as a ```WARNING```. Three things to weigh before acting:
  - *Magnitude, not count.* A few per cent is within the uncertainty of a generalised PMP estimate and a hard clip at the PMP would be spuriously crisp. Twenty per cent or more is structural.
  - *Which side is wrong.* If GSDM/GTSMR patterns breach, two BoM products disagree, and it is not obvious the pattern is the wrong one. Breaches clustered at sub-durations either side of the join in the composite PMP depth-duration curve point to the curve and its areal treatment; breaches well inside one method's own range point to the patterns. The report lists the binding sub-durations for exactly this reason.
  - *Consistency with the PMF.* GSDM/GTSMR patterns also produce the deterministic PMF. Clipping them in the Monte Carlo while leaving the PMF alone would lower the derived curve and hold the PMF fixed, pushing the notional AEP of the PMF rarer for a purely procedural reason. Whatever is decided has to apply to both.
- Before any of that, check whether the breaching sub-durations are anywhere near the reservoir's critical duration. A one-hour exceedance on a catchment with a 24-hour critical duration and meaningful flood storage does not change the routed peak, and the screen is depth-independent so that judgement can be made straight off the table.

The screen also validates the composite PMP depth-duration curve, warning if catchment-average PMP depth is not increasing with duration or PMP intensity is not decreasing. Both would make the ```PMP_s / PMP_d``` ratios meaningless.

## Limitations and possible extensions

- **The PMP bound is screened, not applied.** If a breach ever needs fixing, note that pre-built filtered patterns cannot do it: the clip depth is ```PMP_s × max(1, depth / PMP_d)```, so how much has to come out depends on the sampled depth. A single filtered pattern would have to assume the fully-clipped case and apply it everywhere above ```D_crit```, putting a step in the pattern there and a kink in the derived frequency curve. Applied on the fly at storm generation the clip goes to zero as the depth approaches ```D_crit``` from above, so the routed flood stays a continuous function of depth.
- **Pre-burst is not covered by any of this.** Everything on this page — the neutrality check, the PMP embedded burst screen, and the pattern weight calibration that follows from them — measures the main burst only, before the pre-burst is prepended. Prepending the pre-burst can create embedded bursts spanning the pre-burst/main-burst boundary, including windows *equal to the parent duration itself*: a rolling window of the storm duration positioned across the boundary can exceed the sampled burst depth, which would breach the method framework. The pre-burst filter partially addresses this (see the enveloping burst filter), but none of the checks here scan the pre-burst-inclusive series. Extending them to do so (sub-durations, the parent duration, and enveloping durations against the pre-burst envelope curve) is a natural next step. **A clean neutrality margin or PMP screen result therefore says nothing about the pre-burst**, and neither do calibrated weights. This only matters for runs that actually prepend pre-burst rainfall — where the antecedent condition is carried as burst antecedent storage instead, there is nothing prepended for the checks to miss.
- **Combining parent durations.** The per-duration curves could be combined into an approximate all-durations sub-burst frequency curve. This has not been implemented; treating the durations as independent would be conservative because the same storm population contributes to neighbouring durations.
