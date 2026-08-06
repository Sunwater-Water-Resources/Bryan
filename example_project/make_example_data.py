"""
Build the synthetic data files for the Juniper Creek Dam example project.

Everything this writes is invented. The numbers are plausible in shape and magnitude so
that Bryan behaves the way it would on a real study, but none of them came from a rainfall
analysis, a PMP assessment, the ARR data hub, or a dam. Files are prefixed SYNTHETIC_ so
they cannot be mistaken for study data. Do not copy them into a real assessment.

Run from this folder:  python make_example_data.py

The generator is shipped with the example so the data can be regenerated, and so it is
obvious how each file was made. It is deterministic - the seed is fixed - so re-running it
reproduces the same files byte for byte.
"""
import json
import os
import shutil

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20260806

# ----------------------------------------------------------------- the catchment
DAM = 'Juniper Creek Dam'
SUBCATCHMENT_AREAS = {  # km2, summing to 900
    'JC01': 120.0, 'JC02': 95.0, 'JC03': 140.0, 'JC04': 88.0, 'JC05': 76.0,
    'JC06': 110.0, 'JC07': 65.0, 'JC08': 82.0, 'JC09': 71.0, 'JC10': 53.0,
}
CATCHMENT_AREA = sum(SUBCATCHMENT_AREAS.values())
AREA_BIN = 1000          # both the ARR areal and GTSMR pattern bins for 750-1750 km2
AEP_OF_PMP = 1100000     # 10 ** (log10(900) - 9) inverted, rounded

DURATIONS = [2, 3, 6, 9, 12, 18, 24, 36, 48, 72]
AREAL_DURATIONS = [d for d in DURATIONS if d >= 12]   # shorter ones use the 12 h pattern
TIMESTEP_MIN = {2: 5, 3: 5, 6: 15, 9: 15, 12: 30, 18: 30, 24: 60, 36: 60, 48: 60, 72: 60}
IFD_AEPS = [1.582, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]

# The Monte Carlo sampling scheme. Defined here rather than in make_example_configs.py
# because the synthetic results database for the routing demo has to be sampled on the
# same scheme the config declares - the routing method checks one against the other
# before it will analyse the results.
SCHEME = {'lower_aep': 2, 'upper_aep': 2000000,
          'number_of_main_divisions': 20, 'number_of_sub_divisions': 10}
POINT_FREQUENCIES = ['frequent', 'intermediate', 'rare']

# Dam geometry
FSL = 210.0              # m AHD
CREST = 210.0
DAM_TOP = 218.0          # crest of the non-overflow section
ELS_TOP = 220.0          # the storage curve is carried above it, so routing never runs off the end
SPILLWAY_WIDTH = 300.0   # m, ungated ogee

REGION = 'Juniper Coast'      # label for the temporal patterns and the NRM cluster
# The ARF zone must be one Bryan implements in ArealReduction.get_coefficients:
# 'East Coast North', 'Northern Coastal' or 'Semi-arid Inland Queensland'. The a-i
# coefficients in the data hub file are ignored - Bryan uses its own for the named zone.
ARF_ZONE = 'East Coast North'


def depth_100(duration):
    """Catchment-average 1 in 100 AEP depth (mm) for a duration in hours."""
    return 250.0 * (duration / 24.0) ** 0.42


def growth(duration, aep):
    """Growth factor from the 1 in 100 AEP depth, lognormal in standard normal variate."""
    cv = 0.60 - 0.03 * np.log(duration)
    return float(np.exp(cv * (ndtri(1 - 1 / aep) - ndtri(1 - 1 / 100))))


def pmp_depth(duration):
    """Catchment-average PMP depth (mm). Monotonic in depth, decreasing in intensity."""
    return round(620.0 * (duration / 6.0) ** 0.42, -1)


def spatial_factors(rng, spread=0.09):
    """Subcatchment scaling factors with an area-weighted mean of exactly 1.0."""
    names = list(SUBCATCHMENT_AREAS)
    areas = np.array([SUBCATCHMENT_AREAS[n] for n in names])
    factors = 1.0 + rng.uniform(-spread, spread, len(names))
    factors = factors / ((factors * areas).sum() / areas.sum())
    return pd.Series(np.round(factors, 4), index=names)


def temporal_pattern(rng, n_increments, peakiness, peak_position):
    """One temporal pattern as percentages summing to exactly 100."""
    x = (np.arange(n_increments) + 0.5) / n_increments
    shape = np.exp(-0.5 * ((x - peak_position) / peakiness) ** 2)
    shape = shape * rng.uniform(0.55, 1.45, n_increments)      # storm-to-storm roughness
    shape = shape + 0.04 * shape.max()                          # no truly dry increments
    pattern = np.round(100.0 * shape / shape.sum(), 2)
    pattern[-1] = round(pattern[-1] + (100.0 - pattern.sum()), 2)
    return pattern


def pattern_set(rng, n_increments, n_patterns=10, peakiness=(0.13, 0.30)):
    """Ten patterns spanning a range of peakiness and timing."""
    out = []
    for i in range(n_patterns):
        p = peakiness[0] + (peakiness[1] - peakiness[0]) * i / (n_patterns - 1)
        out.append(temporal_pattern(rng, n_increments, p, rng.uniform(0.30, 0.70)))
    return out


# ----------------------------------------------------------------- writers
def write_ifd_tables(folder, rng):
    os.makedirs(folder, exist_ok=True)
    factors = spatial_factors(rng)
    written = []
    for duration in DURATIONS:
        data = {}
        for aep in IFD_AEPS:
            catchment_average = depth_100(duration) * growth(duration, aep)
            data[f'{aep:g}_AEP'] = np.round(catchment_average * factors, 1)
        df = pd.DataFrame(data, index=factors.index)
        df.index.name = 'name'
        name = f'SYNTHETIC_ifd_{duration:g}h.csv'
        df.to_csv(os.path.join(folder, name))
        written.append((duration, name))
    return written, factors


def write_pmp_tables(folder, rng):
    depths = pd.Series({d: pmp_depth(d) for d in DURATIONS}, name='PMP')
    depths.index.name = 'Duration'
    depths.to_csv(os.path.join(folder, 'SYNTHETIC_pmp_depths.csv'))

    # GSDM and GTSMR spatial patterns, each area-weighting to 1.0
    scaling = pd.DataFrame({'GTSMR': spatial_factors(rng, spread=0.07),
                            'GSDM': spatial_factors(rng, spread=0.11)})
    scaling.index.name = 'ID_Number'
    scaling.to_csv(os.path.join(folder, 'SYNTHETIC_pmp_scaling.csv'))
    return depths


def write_point_patterns(path, rng):
    rows, event = [], 1000
    for duration in DURATIONS:                       # grouped so each block is contiguous
        timestep = TIMESTEP_MIN[duration]
        n = int(duration * 60 / timestep)
        for frequency in POINT_FREQUENCIES:
            for pattern in pattern_set(rng, n):
                event += 1
                rows.append([event, int(duration * 60), timestep, REGION, frequency]
                            + [f'{v:g}' for v in pattern])
    width = max(len(r) for r in rows)
    with open(path, 'w') as f:
        f.write('EventID, Duration, TimeStep, Region, AEP, Increments'
                + ',' * (width - 5) + '\n')
        for row in rows:
            f.write(','.join(str(v) for v in row) + ',' * (width - len(row)) + '\n')


def write_areal_patterns(path, rng):
    rows, event = [], 5000
    for duration in AREAL_DURATIONS:
        timestep = TIMESTEP_MIN[duration]
        n = int(duration * 60 / timestep)
        for pattern in pattern_set(rng, n, peakiness=(0.16, 0.34)):   # areal is smoother
            event += 1
            rows.append([event, int(duration * 60), timestep, REGION, AREA_BIN]
                        + [f'{v:g}' for v in pattern])
    with open(path, 'w') as f:
        f.write('EventID,Duration,TimeStep,Region,Area,Increments\n')
        for row in rows:                    # no trailing commas - the areal reader keeps blanks
            f.write(','.join(str(v) for v in row) + '\n')


def write_pmp_patterns(path, rng, durations, peakiness):
    """GSDM/GTSMR .pat file: a '<minutes> minutes' marker, then 'increments,timestep', then
    ten rows of percentages."""
    lines = [f'{len(durations)},10']
    for duration in durations:
        timestep = TIMESTEP_MIN[duration]
        n = int(duration * 60 / timestep)
        lines.append(f'{int(duration * 60)} minutes')
        lines.append(f'{n},{timestep}')
        for pattern in pattern_set(rng, n, peakiness=peakiness):
            lines.append(','.join(f'{v:.2f}' for v in pattern))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def write_preburst_patterns(path, rng):
    """Pre-burst temporal patterns: ten per duration, indexed by hours before the burst."""
    out = {}
    for duration in DURATIONS:
        length = max(6, int(duration))                 # hours of pre-burst record
        index = list(range(-length, 0))
        patterns = {}
        for i in range(10):
            values = rng.uniform(0.2, 1.0, length) ** 2
            values[: int(length * rng.uniform(0.0, 0.5))] = np.nan   # ragged start
            total = np.nansum(values)
            series = {str(t): (None if np.isnan(v) else round(float(v / total), 5))
                      for t, v in zip(index, values)}
            patterns[str(i)] = series
        out[str(duration)] = patterns
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)


def write_arr_datahub(path, rng):
    """A data hub export with the blocks Bryan reads: LONGARF, LOSSES, PREBURST*, BASEFLOW."""
    aeps_percent = [50, 20, 10, 5, 2, 1]
    preburst_durations = [60, 90, 120, 180, 360, 720, 1080, 1440, 2160, 2880, 4320]

    def preburst_block(scale):
        lines = ['min (h)\\AEP(%),' + ','.join(str(a) for a in aeps_percent)]
        for minutes in preburst_durations:
            hours = minutes / 60
            cells = []
            for j, aep in enumerate(aeps_percent):
                ratio = scale * (0.02 + 0.010 * j) * (1.0 + 0.35 * np.log10(max(hours, 1.0)))
                depth = ratio * depth_100(max(hours, 1.0)) * 0.8
                cells.append(f'{depth:.1f} ({ratio:.3f})')
            lines.append(f'{minutes} ({hours:.1f}),' + ','.join(cells))
        return lines

    text = [
        'Results - ARR Data Hub (SYNTHETIC - not a real data hub export)',
        '[STARTTXT]', '',
        'Input Data Information', '[INPUTDATA]',
        'Latitude,-00.000000', 'Longitude,000.000000', '[END_INPUTDATA]', '',
        'River Region', '[RIVREG]',
        'Division,Synthetic Division', 'River Number,0',
        'River Name,Juniper Creek (synthetic)', 'Shape Intersection (%),100.0',
        '[RIVREG_META]', 'Version,SYNTHETIC', '[END_RIVREG]', '',
        'ARF Parameters', '[LONGARF]',
        f'Zone,{ARF_ZONE}',
        'a,0.327', 'b,0.241', 'c,0.448', 'd,0.36', 'e,0.00096',
        'f,0.48', 'g,-0.21', 'h,0.012', 'i,-0.0013',
        'Shape Intersection (%),100.0',
        '[LONGARF_META]', 'Version,SYNTHETIC', '[END_LONGARF]', '',
        'Storm Losses', '[LOSSES]',
        'Storm Initial Losses (mm),30.0', 'Storm Continuing Losses (mm/h),2.0',
        '[LOSSES_META]', 'Version,SYNTHETIC', '[END_LOSSES]', '',
        'Temporal Patterns', '[TP]', 'code,Juniper', f'Label,{REGION}',
        '[TP_META]', 'Version,SYNTHETIC', '[END_TP]', '',
        'Baseflow', '[BASEFLOW]',
        'Volume Factor,0.24', 'Peak Factor,0.11', 'Baseflow Region,Synthetic',
        '[BASEFLOW_META]', 'Version,SYNTHETIC', '[END_BASEFLOW]', '',
    ]
    for percentile, scale in [('', 1.0), ('10', 0.15), ('25', 0.55),
                              ('75', 1.55), ('90, ', 2.30)]:
        tag = percentile.strip().rstrip(',')
        label = 'Median' if tag == '' else f'{tag}%'
        text += [f'{label} Preburst Depths and Ratios', f'[PREBURST{tag}]']
        text += preburst_block(scale)
        text += [f'[PREBURST{tag}_META]', 'Version,SYNTHETIC',
                 f'[END_PREBURST{tag}]From preburst class', '']
    with open(path, 'w') as f:
        f.write('\n'.join(text) + '\n')


def write_climate_data(folder, rng):
    os.makedirs(os.path.join(folder, 'Temporal_Pattern_Weighting'), exist_ok=True)

    losses = pd.DataFrame(
        [['Juniper Coast', 0.6, 2.0, 4.3, 1.1, 3.8, 8.0],
         ['Juniper Inland', 0.4, 1.1, 2.2, -0.5, 2.0, 7.5]],
        columns=['Natural Resource Management Cluster', 'IL lower', 'IL mean', 'IL upper',
                 'CL lower', 'CL mean', 'CL upper'])
    losses.to_csv(os.path.join(folder, 'SYNTHETIC_rainfall_loss_rates_of_change.csv'),
                  index=False)

    scaling = pd.DataFrame({'Duration': DURATIONS,
                            'Cfa': np.round(-0.5 + 0.06 * np.arange(len(DURATIONS)), 2),
                            'Aw': np.round(-0.4 + 0.05 * np.arange(len(DURATIONS)), 2)})
    scaling.to_csv(os.path.join(folder, 'SYNTHETIC_temporal_pattern_scaling_factors.csv'),
                   index=False)

    def weights(extra_column, extra_values, durations):
        rows = []
        for duration in durations:
            for value in extra_values:
                front = ','.join(str(v) for v in rng.choice(np.arange(1, 11), 3, replace=False))
                deltas = np.round(0.12 + 0.021 * np.arange(5) + rng.uniform(-0.004, 0.004), 5)
                row = {'Duration': int(duration * 60)}
                if extra_column:
                    row[extra_column] = value
                row['Front Patterns'] = front
                row.update({f'DeltaD50-{i + 1}': deltas[i] for i in range(5)})
                rows.append(row)
        return pd.DataFrame(rows)

    base = os.path.join(folder, 'Temporal_Pattern_Weighting')
    weights('AEP', POINT_FREQUENCIES, DURATIONS).to_csv(
        os.path.join(base, 'SYNTHETIC_point.csv'), index=False)
    weights('Area', [AREA_BIN], AREAL_DURATIONS).to_csv(
        os.path.join(base, 'SYNTHETIC_areal.csv'), index=False)
    weights(None, [None], DURATIONS).to_csv(
        os.path.join(base, 'SYNTHETIC_gsdm.csv'), index=False)
    weights('Area', [AREA_BIN], DURATIONS).to_csv(
        os.path.join(base, 'SYNTHETIC_gtsmr.csv'), index=False)


def write_dam_curves(folder):
    """A URBS .els elevation/area/storage table and a .sq storage/discharge rating."""
    os.makedirs(folder, exist_ok=True)
    levels = np.round(np.arange(180.0, ELS_TOP + 0.25, 0.25), 2)
    # Storage from a power law in depth above the base, giving 190,000 ML at FSL
    depth = np.clip(levels - 180.0, 0, None)
    volume = np.round(190000.0 * (depth / (FSL - 180.0)) ** 2.35, 1)
    area = np.round(np.gradient(volume, levels) / 10.0, 1)     # ML per m -> km2 (1 ML = 1000 m3)
    els = pd.DataFrame({'EL': levels, 'A': np.clip(area, 0, None), 'V': volume})
    els.to_csv(os.path.join(folder, 'SYNTHETIC_juniper.els'), index=False)

    fsv = float(np.interp(FSL, levels, volume))
    # Ungated ogee spillway, 90 m wide, discharging above the crest at FSL
    head = np.round(np.arange(0.0, ELS_TOP - CREST + 0.05, 0.05), 2)
    flow = np.round(2.1 * SPILLWAY_WIDTH * head ** 1.5, 1)
    storage_above_fsv = np.round(
        np.interp(CREST + head, levels, volume) - fsv, 1)
    with open(os.path.join(folder, 'SYNTHETIC_juniper.sq'), 'w') as f:
        f.write(f'{DAM} (SYNTHETIC)\n*\n*\n*\n{len(head)} PAIRS:\n')
        for s, q in zip(storage_above_fsv, flow):
            f.write(f'{s}\t{q}\n')
    return fsv


def write_model_stubs(folder):
    """Skeleton RORB and URBS model files.

    These are stubs, not working hydrologic models. They carry enough structure for Bryan to
    start - RORB reads the .catg at construction and looks for the dam routing line, URBS
    reads the .vec - so the storm generation side can be exercised without a licence. Neither
    model will actually run against them.
    """
    rorb = os.path.join(folder, 'rorb')
    urbs = os.path.join(folder, 'urbs')
    ratings = os.path.join(urbs, 'ratings')
    os.makedirs(rorb, exist_ok=True)
    os.makedirs(urbs, exist_ok=True)
    os.makedirs(ratings, exist_ok=True)
    # URBS reads the elevation/storage and storage/discharge curves out of the ratings
    # folder named in the config, so put copies where the vec file's DATAFILE= points
    for curve in ['SYNTHETIC_juniper.els', 'SYNTHETIC_juniper.sq']:
        shutil.copy(os.path.join(folder, curve), ratings)

    names = list(SUBCATCHMENT_AREAS)
    catg = [f'{DAM} - SYNTHETIC STUB, not a working RORB model',
            'C This file exists so Bryan can find the dam routing line and copy a catg into',
            'C the storms folder. It will not run in RORB.',
            '1']
    for name in names:                      # type 1 = route through a reach, 2 = add subarea
        catg += ['2', name, '1,0.00']
    catg += ['16', 'JuniperDam', '0,1,190000.0,0.0', '7', '0']
    catg.append('C AREAS')
    catg.append(','.join(f'{SUBCATCHMENT_AREAS[n]:g}' for n in names))
    with open(os.path.join(rorb, 'juniper.catg'), 'w') as f:
        f.write('\n'.join(catg) + '\n')

    # Bryan reads line 7 (index 6) for the interstation area count, so the line ordering
    # here matters even though the file is a stub.
    par = [f'{DAM} - SYNTHETIC STUB RORB parameter file',
           'C Bryan reads this as a template and writes one par file per simulation with the',
           'C rainfall losses filled in. Replace with a calibrated par file for a real study.',
           'Run name: juniper',
           'Catchment file: juniper.catg',
           'Storm file: placeholder.stm',
           'Number of interstation areas: 1',
           'Kc: 2.10',
           'm: 0.80',
           'Initial loss (mm): 30.0',
           'Continuing loss (mm/h): 2.0']
    with open(os.path.join(rorb, 'juniper.par'), 'w') as f:
        f.write('\n'.join(par) + '\n')

    vec = [f'{DAM} - SYNTHETIC STUB, not a working URBS model',
           'C This file exists so the URBS config in this example has something to point at.',
           'C It will not run in URBS. Replace with a calibrated vec file for a real study.',
           'CATCHMENT DATA FILE = juniper.cat',
           'MODEL: SPLIT',
           'USES: L, Sc, I',
           'DEFAULT PARAMETERS: alpha = 0.30 m = 0.80 beta = 1.50',
           '*']
    for name in names:
        vec += [f'RAIN #{names.index(name) + 1}', f'SUB AREA #{names.index(name) + 1} = '
                f'{SUBCATCHMENT_AREAS[name]:g}', 'ROUTE THRU L = 5.0 Sc = 1.0']
    # The level-based dam routing form: Bryan reads FSL= and DATAFILE= off this line
    vec += [f'DAM ROUTE FSL={FSL} DATAFILE=SYNTHETIC_juniper.els il=190000.0 '
            f'PRINT JuniperDam.',
            'PRINT JuniperDam.', 'END OF CATCHMENT DATA.']
    with open(os.path.join(urbs, 'juniper.vec'), 'w') as f:
        f.write('\n'.join(vec) + '\n')


def write_focal_subcatchments(path):
    df = pd.DataFrame({'Name': list(SUBCATCHMENT_AREAS),
                       'Area': list(SUBCATCHMENT_AREAS.values())})
    df.to_csv(path, index=False)


def write_routing_inputs(folder, rng, fsv):
    """A Monte Carlo results database and matching inflow hydrographs, as a run with the
    hydrologic model already run would have left behind. The 'reservoir routing' demo
    re-routes these, so no rainfall-runoff model is needed to see the analysis chain work."""
    os.makedirs(folder, exist_ok=True)
    m = SCHEME['number_of_main_divisions']
    n = SCHEME['number_of_sub_divisions']
    duration = 24
    z_lower = ndtri(1 - 1 / SCHEME['lower_aep'])
    z_upper = ndtri(1 - 1 / SCHEME['upper_aep'])
    rows = []
    for i in range(m):                      # one sample per sub-division of each division
        for j in range(n):
            z = z_lower + (i + (j + 0.5) / n) / m * (z_upper - z_lower)
            rows.append({'m': i, 'n': j, 'rain_z': z})
    mcdf = pd.DataFrame(rows)
    size = len(mcdf)
    # The columns a real Monte Carlo run writes, so the routing method sees the database
    # shape it expects rather than a cut-down one
    mcdf['rain_aep'] = np.round(1 / (1 - ndtr(mcdf['rain_z'])), 3)
    mcdf['mean_rain_mm'] = np.round(depth_100(duration) * np.exp(0.55 * mcdf['rain_z']), 1)
    mcdf['tp'] = rng.integers(0, 10, size)
    mcdf['storm_method'] = np.where(mcdf['rain_aep'] > 2000,
                                    rng.choice(['GSDM', 'GTSMR'], size), 'ARR areal')
    mcdf['tp_frequency'] = np.where(mcdf['rain_aep'] > 31.25, 'rare', 'intermediate')
    mcdf['il_p'] = np.round(rng.uniform(0, 1, size), 4)
    mcdf['il_scaling'] = np.round(3.19 - 3.05 * mcdf['il_p'], 3)
    mcdf['initial_loss'] = np.round(30.0 * mcdf['il_scaling'], 1)
    mcdf['cl_p'] = np.round(rng.uniform(0, 1, size), 4)
    mcdf['cl_scaling'] = np.round(0.4 + 1.2 * mcdf['cl_p'], 3)
    mcdf['continuing_loss'] = np.round(2.0 * mcdf['cl_scaling'], 2)
    mcdf['preburst_p'] = np.round(rng.uniform(0.1, 0.9, size), 4)
    mcdf['preburst_proportion'] = np.round(rng.uniform(0.0, 0.15, size), 4)
    mcdf['preburst_mm'] = np.round(mcdf['mean_rain_mm'] * mcdf['preburst_proportion'], 1)
    mcdf['residual_depth'] = np.round(np.clip(mcdf['preburst_mm'] - mcdf['initial_loss'],
                                              0, None), 1)
    mcdf['embedded_bursts'] = 'No embedded bursts'
    mcdf['lake_z'] = np.round(rng.normal(0, 1, size), 5)
    mcdf['ADV'] = np.round(np.clip(45000 + 90000 * (1 + mcdf['lake_z'] / 3), 45000, fsv), 0)

    # Peak inflow rising with the sampled rainfall, with realisation-to-realisation scatter.
    # This stands in for the rainfall-runoff model the example cannot run.
    peak = 250.0 * np.exp(0.755 * mcdf['rain_z']) * rng.lognormal(0, 0.16, size)
    mcdf['inflow'] = np.round(peak, 1)
    mcdf.index.name = 'index'
    mcdf.to_csv(os.path.join(folder, 'SYNTHETIC_mc_24h__mcdf.csv'))

    # Gamma-shaped inflow hydrographs on a 1 hour timestep
    hours = np.arange(0.0, 96.0, 1.0)
    shape = (hours / 14.0) ** 3.2 * np.exp(-(hours / 14.0) * 3.2 + 3.2)
    shape = shape / shape.max()
    inflows = pd.DataFrame({str(k): np.round(shape * p, 2) for k, p in enumerate(peak)},
                           index=np.round(hours, 3))
    inflows.index.name = 'time'
    inflows.to_csv(os.path.join(folder, 'SYNTHETIC_mc_24h_inflows.csv'))


def main():
    rng = np.random.default_rng(SEED)
    storm = os.path.join(HERE, 'storm_data')
    patterns = os.path.join(storm, 'patterns')
    ifd = os.path.join(storm, 'ifd')
    os.makedirs(patterns, exist_ok=True)

    print(f'{DAM}: {CATCHMENT_AREA:g} km2 over {len(SUBCATCHMENT_AREAS)} subcatchments')
    ifd_files, _ = write_ifd_tables(ifd, rng)
    pmp = write_pmp_tables(ifd, rng)
    print('IFD tables and PMP depths written')
    print('  PMP (mm):', ', '.join(f'{d:g}h {v:g}' for d, v in pmp.items()))

    write_point_patterns(os.path.join(patterns, 'SYNTHETIC_point_increments.csv'), rng)
    write_areal_patterns(os.path.join(patterns, 'SYNTHETIC_areal_increments.csv'), rng)
    write_pmp_patterns(os.path.join(patterns, 'SYNTHETIC_gsdm.pat'), rng,
                       DURATIONS, peakiness=(0.15, 0.32))
    write_pmp_patterns(os.path.join(patterns, f'SYNTHETIC_gtsmr_{AREA_BIN}.pat'), rng,
                       DURATIONS, peakiness=(0.18, 0.36))
    write_preburst_patterns(os.path.join(patterns, 'SYNTHETIC_preburst_patterns.json'), rng)
    print('Temporal patterns written')

    write_arr_datahub(os.path.join(storm, 'SYNTHETIC_arr_datahub.txt'), rng)
    write_climate_data(os.path.join(HERE, 'climate_change'), rng)
    print('ARR data hub file and climate change tables written')

    os.makedirs(os.path.join(HERE, 'sim_options', 'focal_locations'), exist_ok=True)
    write_focal_subcatchments(
        os.path.join(HERE, 'sim_options', 'focal_locations',
                     'SYNTHETIC_dam_subcatchments.csv'))
    fsv = write_dam_curves(os.path.join(HERE, 'model'))
    write_model_stubs(os.path.join(HERE, 'model'))
    print(f'Dam curves and model stubs written: FSL {FSL} m AHD -> FSV {fsv:,.0f} ML')

    write_routing_inputs(os.path.join(HERE, 'inflows'), rng, fsv)
    print('Stored inflows and results database for the routing demo written')

    return ifd_files, fsv


if __name__ == '__main__':
    main()
