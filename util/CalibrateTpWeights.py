"""
Calibrate temporal pattern probability weights so that the Monte Carlo ensemble satisfies the
sub-burst AEP-neutrality condition (see Manual/SubDocs/sub_burst_check.md).

Needs an mcdf file produced with sub-burst tracking (the subburst_<d>h / ifd_<d>h columns),
normally from an unfiltered run (the 'ebf' exclusion). The models themselves are not needed:
use 'Run models: storms only' in the simulation list for a cheap calibration run.

The calibration is entirely offline. Sub-burst depths depend only on the sampled pattern, depth,
and storm method - all recorded in the mcdf - so trial weights are evaluated by re-weighting the
recorded realisations in the TPT, with no model reruns. Patterns whose sub-bursts exceed the
same-z IFD are progressively down-weighted until the weighted sub-burst frequency curves sit at
or below the IFD.

Outputs (per mcdf file, into the output folder):
- <label>_tp_weights.csv: the calibrated weights. Rows are the weight groups (storm method, plus
  the frequency bin for ARR point patterns, e.g. 'ARR point|rare'); columns are the pattern
  indices 0-9 as sampled in the mcdf 'tp' column.
- <label>_subburst_margins.csv: the neutrality margins (TPT sub-burst depth / same-z IFD depth)
  at the standard AEPs, before (uniform) and after (weighted) calibration, per sub-duration.
- <label>_subburst<d>h_weighted.csv/.png: the weighted sub-burst TPT curves and plots.
- <label>_weighted_<inflow/level/outflow>.csv: if flow results are present in the mcdf, the
  weighted flood quantiles - a preview of the effect of the weights on the flood frequency
  curve, computed from the existing model runs.

Interpretation notes:
- If a pattern's weight lands on the floor (weight_floor) and the margins still exceed the
  target, weighting alone is not enough - consider the embedded burst filter, or review that
  pattern.
- The mcdf is per storm duration, so weights are calibrated per duration. Whether to pool or
  smooth weights across durations is a judgement call.
- The margins are edge-affected near the sampling bounds (like all TPT results), so the
  convergence test excludes the first and last standard AEPs by default (see aep_range).
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from scipy.special import ndtri

# Bryan's lib is imported for the TPT machinery (run this script from anywhere)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.MCScheme import SampleScheme, TotalProbTheorem, interpolate_reference_curve, neutrality_margin

# ------------------------------------------------------------------------------- user settings
mc_config_file = r'C:\path\to\mc_config.json'   # the monte carlo config (for the scheme_config)
mcdf_files = {                                  # label: mcdf path - one entry per simulation
    '24h': r'C:\path\to\outputs\mc_24h__mcdf.csv',
}
output_folder = r'C:\path\to\outputs\tp_weights'

margin_target = 1.0      # calibrate until all tested margins <= this
aep_range = None         # e.g. [200, 500000] to set the AEP range for the convergence test;
                         # None tests all but the first and last standard AEPs (TPT edge effects)
reduction_factor = 0.5   # per-iteration down-weighting: w <- w * reduction_factor**breach_fraction
weight_floor = 0.02      # minimum pattern weight (a weight at the floor means weighting may not be enough)
max_iterations = 25
# ----------------------------------------------------------------------------------------------


def get_sub_cols(mcdf):
    return [col for col in mcdf.columns if col.startswith('subburst_')]


def get_group_key(row):
    # Weights are grouped by storm method; ARR point patterns also differ by frequency bin
    if row['storm_method'] == 'ARR point':
        return f"ARR point|{row['tp_frequency']}"
    return row['storm_method']


def prepare_mcdf(mcdf):
    # Weight group and per-realisation severity: the worst sub-burst relative to the same-z IFD
    mcdf['tp'] = mcdf['tp'].astype(int)
    mcdf['group'] = mcdf.apply(get_group_key, axis=1)
    sub_cols = get_sub_cols(mcdf)
    ifd_cols = [col.replace('subburst_', 'ifd_') for col in sub_cols]
    ratio = pd.DataFrame(mcdf[sub_cols].to_numpy() / mcdf[ifd_cols].to_numpy(), index=mcdf.index)
    mcdf['severity'] = ratio.max(axis=1)        # skips NaN; NaN if no measurements for the row
    return mcdf


def attach_weights(mcdf, weights):
    mcdf['tp_w'] = [weights[group][tp] for group, tp in zip(mcdf['group'], mcdf['tp'])]


def probability_margins(mcdf, tpt_by_col, ifd_at_std_by_col, std_aeps):
    # For each sub-duration and standard AEP: the (weighted) TPT annual exceedance probability
    # of the same-AEP IFD depth, multiplied by the AEP. Above 1.0 means sub-bursts of that depth
    # occur more often than the IFD implies. Equivalent to the depth-space neutrality margin
    # exceeding 1.0, but much cheaper to evaluate inside the calibration loop.
    margins = {}
    for col, tpt in tpt_by_col.items():
        tpt.mcdf = mcdf[['m', col, 'tp_w']]
        ifd_at_std = ifd_at_std_by_col[col]
        ratios = []
        for aep, ifd_depth in zip(std_aeps, ifd_at_std.to_numpy()):
            if np.isnan(ifd_depth):
                ratios.append(np.nan)
            else:
                ratios.append(tpt.assign_aep(ifd_depth, col) * aep)
        margins[col.replace('subburst_', '')] = pd.Series(ratios, index=std_aeps)
    return pd.DataFrame(margins)


def calibrate(mcdf, mc, label):
    sub_cols = get_sub_cols(mcdf)
    num_patterns = mc.number_of_temporal_patterns
    std_aeps = mc.get_standard_aeps()
    std_z = ndtri(1 - 1 / np.array(std_aeps, dtype=float))

    # The AEP window for the convergence test
    if aep_range is None:
        test_aeps = std_aeps[1:-1]
    else:
        test_aeps = [aep for aep in std_aeps if aep_range[0] <= aep <= aep_range[1]]
    print('\nTesting neutrality margins over AEPs:', test_aeps)

    # Build the TPT machinery and IFD reference curves once per sub-duration
    tpt_by_col = {}
    ifd_at_std_by_col = {}
    reference_by_col = {}
    for col in sub_cols:
        ifd_col = col.replace('subburst_', 'ifd_')
        reference = mcdf[['rain_z', 'storm_method', ifd_col]].dropna()
        reference_by_col[col] = reference
        ifd_at_std_by_col[col] = interpolate_reference_curve(reference, ifd_col, std_z)
        tpt_by_col[col] = TotalProbTheorem(mc.m, mc.n, mc.main_divisions, mcdf[['m', col]])

    # Per-group, per-pattern breach fraction: how often the pattern's worst sub-burst exceeds
    # the same-z IFD (static: the measurements don't change with the weights). This is the
    # heuristic that directs the down-weighting; the margin test is the arbiter of convergence.
    breach_by_group = {}
    for group, group_df in mcdf.groupby('group'):
        breach = group_df.groupby('tp')['severity'].apply(lambda s: (s.dropna() > 1.0).mean())
        breach = breach.reindex(range(num_patterns)).fillna(0.0)
        breach_by_group[group] = breach.to_numpy()
        print(f'\nEmbedded burst fraction by pattern for {group} ({len(group_df)} realisations):')
        print(np.around(breach_by_group[group], 3))

    # Iteratively down-weight the offending patterns until the margins are neutral, the weights
    # stop changing (offenders at the floor), or the iteration limit is reached
    weights = {group: np.full(num_patterns, 1.0 / num_patterns) for group in breach_by_group}
    converged = False
    for iteration in range(max_iterations + 1):
        attach_weights(mcdf, weights)
        margins = probability_margins(mcdf, tpt_by_col, ifd_at_std_by_col, std_aeps)
        worst = margins.loc[test_aeps].max().max()
        print(f'Iteration {iteration}: worst probability margin = {np.around(worst, 3)}')
        if worst <= margin_target:
            converged = True
            break
        stagnant = True
        for group in weights:
            w = weights[group] * reduction_factor ** breach_by_group[group]
            w = np.maximum(w, weight_floor)
            w = w / w.sum()
            if not np.allclose(w, weights[group], atol=1e-6):
                stagnant = False
            weights[group] = w
        if stagnant:
            print('The weights are no longer changing (offending patterns at the floor) - stopping')
            break

    if converged:
        print(f'\nConverged: all tested margins <= {margin_target}')
    else:
        print('\nWARNING: neutrality was not achieved.')
        print('If the offending patterns are at the weight floor, weighting alone cannot achieve')
        print('neutrality for this simulation. Consider the embedded burst filter, or review the')
        print('offending patterns.')

    # Report the weights
    weights_df = pd.DataFrame(weights).transpose()
    weights_df.index.name = 'group'
    print('\nCalibrated temporal pattern weights:')
    print(np.around(weights_df, 3))
    floored = (weights_df <= weight_floor + 1e-9).sum().sum()
    if floored:
        print(f'NOTE: {floored} pattern weight(s) are at the floor of {weight_floor}')
    weights_file = os.path.join(output_folder, f'{label}_tp_weights.csv')
    print('Writing weights:', weights_file)
    weights_df.to_csv(weights_file)

    return weights, reference_by_col


def report_margins(mcdf, mc, reference_by_col, label, suffix, write_outputs=False):
    # Depth-space neutrality margins (as reported by Bryan's 'Analyse sub-bursts' step)
    mc.df = mcdf
    margins = {}
    for col, reference in reference_by_col.items():
        dur = col.replace('subburst_', '')
        ifd_col = f'ifd_{dur}'
        if write_outputs:
            output_name = os.path.join(output_folder, f'{label}_subburst{dur}_{suffix}')
            mc.compute_std_quantiles(result_type=col, output_filename=f'{output_name}.csv')
            mc.plot_tpt_results_2(col, f'{output_name}.png', reference=reference)
        else:
            mc.compute_std_quantiles(result_type=col, output_filename=None)
        margins[f'{dur} {suffix}'] = neutrality_margin(mc.quantiles[col], reference, col, ifd_col)
    return pd.DataFrame(margins)


def preview_flood_quantiles(mcdf, mc, label):
    # Preview the effect of the weights on the flood frequency curve using the existing runs
    mc.df = mcdf
    for result_type in ['inflow', 'level', 'outflow']:
        if result_type in mcdf.columns and mcdf[result_type].notna().any():
            output_name = os.path.join(output_folder, f'{label}_weighted_{result_type}')
            mc.compute_std_quantiles(result_type=result_type, output_filename=f'{output_name}.csv')
            mc.plot_tpt_results_2(result_type, f'{output_name}.png')
        else:
            print(f'No {result_type} results in the mcdf - skipping the weighted flood preview')


def main():
    os.makedirs(output_folder, exist_ok=True)
    print('Opening the monte carlo config file:', mc_config_file)
    with open(mc_config_file) as f:
        scheme_config = json.load(f)['scheme_config']

    for label, mcdf_file in mcdf_files.items():
        print(f'\n{"="*100}\nCalibrating temporal pattern weights for {label}\n{"="*100}')
        print('Opening the Monte Carlo analysis file:', mcdf_file)
        mcdf = pd.read_csv(mcdf_file, index_col=0)
        if not get_sub_cols(mcdf):
            print('No subburst_ columns found in this mcdf - was it produced with sub-burst tracking?')
            continue
        mcdf = prepare_mcdf(mcdf)

        mc = SampleScheme(lower_aep=scheme_config['lower_aep'],
                          upper_aep=scheme_config['upper_aep'],
                          number_of_main_divisions=scheme_config['number_of_main_divisions'],
                          number_of_sub_divisions=scheme_config['number_of_sub_divisions'],
                          number_of_temporal_patterns=scheme_config['number_of_temporal_patterns'])

        # Margins before calibration (uniform weights)
        mcdf_uniform = mcdf.drop(columns=['tp_w'], errors='ignore')
        sub_cols = get_sub_cols(mcdf)
        reference_by_col = {col: mcdf[['rain_z', 'storm_method', col.replace('subburst_', 'ifd_')]].dropna()
                            for col in sub_cols}
        before = report_margins(mcdf_uniform, mc, reference_by_col, label, 'uniform')

        # Calibrate and re-report
        weights, reference_by_col = calibrate(mcdf, mc, label)
        attach_weights(mcdf, weights)
        after = report_margins(mcdf, mc, reference_by_col, label, 'weighted', write_outputs=True)

        margins_df = pd.concat([before, after], axis=1).sort_index(axis=1)
        margins_df.index.name = 'aep (1 in x)'
        margins_file = os.path.join(output_folder, f'{label}_subburst_margins.csv')
        print('\nNeutrality margins before (uniform) and after (weighted) calibration:')
        print(margins_df)
        print('Writing margins:', margins_file)
        margins_df.to_csv(margins_file)

        preview_flood_quantiles(mcdf, mc, label)


if __name__ == '__main__':
    main()
