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

Notes on the weights:

- Weights are grouped by storm method (and by frequency bin for ARR point patterns), and calibrated per storm duration.
- A floor is applied (default 0.02) so that no pattern is silently excluded. If the margins cannot be brought below the target with the offending patterns at the floor, weighting alone cannot achieve neutrality for that simulation — the script says so, and the choice reverts to filtering or pattern review.
- The calibration stops as soon as neutrality is achieved, so patterns are not down-weighted further than the condition requires.

### Applying the weights to the sampling

For a production run, the calibrated weights can be applied to the pattern sampling itself: set the optional ```TP weights``` key in the [simulation list](sim_list.md) to the weights file written by the calibration script. The patterns are then sampled with the calibrated probabilities (composed with the D50 climate shift weightings), and the weight applied to each realisation is recorded in a ```tp_weight``` column of the mcdf for audit. Three things to be aware of:

- **The standard TPT analysis is already correct for a weighted-sampling run** — the sampled frequencies embody the weights, which *are* the asserted pattern probability model. Do **not** also attach a ```tp_w``` column to such a run's mcdf (that would double-count the weights).
- Rerunning the sub-burst check on the weighted-sampling run verifies the neutrality of the ensemble as actually simulated — a worthwhile closing check, since the sampled realisations will differ from the calibration run's.
- The ```TP weights``` key cannot be combined with replication of the pattern sampling (the ```tp``` replicate), because the replicated patterns were sampled under different weights. Since heavily down-weighted patterns are rarely sampled, their response is also thinly resolved — another reason the weight floor matters.

## Limitations and possible extensions

- **Pre-burst is not yet covered.** Prepending the pre-burst can create embedded bursts spanning the pre-burst/main-burst boundary, including windows *equal to the parent duration itself*: a rolling window of the storm duration positioned across the boundary can exceed the sampled burst depth, which would breach the method framework. The pre-burst filter partially addresses this (see the enveloping burst filter), but the neutrality check does not yet scan the pre-burst-inclusive series. Extending the check to do so (sub-durations, the parent duration, and enveloping durations against the pre-burst envelope curve) is a natural next step.
- **Combining parent durations.** The per-duration curves could be combined into an approximate all-durations sub-burst frequency curve. This has not been implemented; treating the durations as independent would be conservative because the same storm population contributes to neighbouring durations.
