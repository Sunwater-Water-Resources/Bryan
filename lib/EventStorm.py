"""Rebuilding the storm behind one Monte Carlo realisation.

The mcdf records what was *sampled* - a rainfall variate, a pattern number, a
pre-burst proportion - not the rainfall series that came out of it. So a
representative event has peaks, and a hyetograph only if it is rebuilt.

Rebuilding is deterministic: every random draw is already in the mcdf row, so
this replays ``MonteCarloSimulator.run_models`` for one realisation with the
sampling taken out. It has to replay it *exactly* - the same filter buffers, the
same order of climate adjustment, the same pre-burst pattern - or the hyetograph
belongs to a storm that was never run. Two things guard that:

- the sequence below is written from ``lib/Simulator.py`` and cites it, rather
  than from ``StormInstance.py``, which does the same job for a different
  purpose and has drifted (it omits the ``buffer=0.9`` on the pre-burst filter);
- every rebuild is **checked against what the run recorded**. The catchment
  average depth must come back as ``mean_rain_mm``, the pre-burst as
  ``preburst_mm``, and the embedded burst comment as ``embedded_bursts``. Those
  three are independent of each other and were written by the run itself, so a
  rebuild that agrees with all three is the storm that was modelled. A rebuild
  that does not is reported and not quietly plotted.

Temporal patterns are **percentages of the main burst depth** throughout
(``uniform_preburst``: ``preburst_depth = preburst_proportion * 100``), so the
depths here are ``pattern / 100 * ave_rain``. The index is hours, with the
pre-burst at negative times: t = 0 is the start of the main burst, which is also
what the stored hydrographs need shifting onto.

Needs Bryan's own dependencies - it drives ``StormGenerator`` - so it is not
importable by the launcher and does not try to be.
"""

import numpy as np
import pandas as pd

from lib.ClimateChange import ClimateAdjustment
from lib.StormGenerator import StormBurst

# Simulator.run_models hard-codes both of these, and they are not the same:
# the main burst is filtered to 1.1x the sub-duration IFD depth, the pre-burst
# to 0.9x. Getting either wrong changes the pattern without changing its total.
BURST_FILTER_BUFFER = 1.1
PREBURST_FILTER_BUFFER = 0.9

# How closely a rebuild has to match what the run recorded before it is trusted.
# Bryan rounds several of these to one decimal on the way into the mcdf, so this
# is a tolerance on rounding, not on method.
DEPTH_TOLERANCE = 0.005          # 0.5%


class StormContext:
    """The run that produced a realisation - everything but the realisation.

    All of it comes from the simulation list row and the configs it names, so
    the caller reads them once per source run rather than per event.
    """

    def __init__(self, storm_config, focal_subcatchments, duration,
                 climate_config=None, gwl=None, year=None, ssp=None,
                 initial_loss=None, continuing_loss=None,
                 exclusions='', preburst_method=None,
                 lower_aep=2, upper_aep=3000):
        self.storm_config = storm_config
        self.focal_subcatchments = focal_subcatchments
        self.duration = float(duration)
        self.climate_config = climate_config
        self.gwl = gwl
        self.year = year
        self.ssp = ssp
        self.initial_loss = initial_loss
        self.continuing_loss = continuing_loss
        self.exclusions = _parse_exclusions(exclusions)
        self.preburst_method = (str(preburst_method).lower().strip()
                                if preburst_method else None)
        self.lower_aep = float(lower_aep)
        self.upper_aep = float(upper_aep)

    @property
    def filters_embedded_bursts(self):
        return not self.exclusions['embedded_burst_filter']

    @property
    def applies_preburst(self):
        return not self.exclusions['preburst']


def _parse_exclusions(line):
    """The two exclusion keys that change the storm's shape.

    ``Simulator.set_exclusions`` is the authority on the key spellings; only
    ``ebf`` and ``pb`` reach the rainfall series, so only they are read here.
    The others (climate uplifts, CL sampling) change depths that the mcdf
    already records, and those recorded values are what is used below.
    """
    keys = {part.strip() for part in str(line or '').split(',')}
    return {'embedded_burst_filter': 'ebf' in keys, 'preburst': 'pb' in keys}


class Hyetograph:
    """One realisation's catchment-average rainfall, and whether it can be trusted."""

    def __init__(self, depths, timestep, preburst_hours, checks):
        self.depths = depths              # mm per timestep, hours on the index
        self.timestep = timestep          # hours
        self.preburst_hours = preburst_hours
        self.checks = checks              # list of (what, ok, detail)

    @property
    def trustworthy(self):
        return all(ok for _, ok, _ in self.checks)

    @property
    def problems(self):
        return [f'{what}: {detail}' for what, ok, detail in self.checks if not ok]

    @property
    def burst_mm(self):
        return float(self.depths[self.depths.index >= 0].sum())

    @property
    def preburst_mm(self):
        return float(self.depths[self.depths.index < 0].sum())

    @property
    def total_mm(self):
        return float(self.depths.sum())

    @property
    def intensity(self):
        """mm/h, for a plot that should not depend on the pattern's timestep."""
        return self.depths / self.timestep


def build_storm(context: StormContext):
    """The rainfall data for a run, set up once and reused for every event.

    Mirrors ``Simulator.initialise_storm``: the imports are the expensive part
    (IFD tables, four families of temporal patterns, the pre-burst set), and
    they depend on the run rather than on the realisation.
    """
    storm = StormBurst(context.storm_config)
    storm.load_subcatchment_areas(context.focal_subcatchments)
    if context.initial_loss is not None:
        storm.storm_initial_loss = context.initial_loss
    if context.continuing_loss is not None:
        storm.continuing_loss = context.continuing_loss

    durations = [context.duration]
    storm.import_rare_rainfall()
    storm.apply_areal_reduction()                    # before the extreme rainfall
    storm.skip_extreme_methods(durations, do_preburst=context.applies_preburst)

    if context.upper_aep > 2000:
        storm.set_up_extreme_rainfall(durations)
        storm.import_gsdm_temporal_patterns(durations)
        storm.import_gtsmr_temporal_patterns(durations)

    storm.import_arr_areal_patterns(durations)
    storm.import_arr_point_patterns(durations)

    if context.applies_preburst:
        storm.import_preburst_patterns()

    climate = _climate(context)
    return storm, climate


def _climate(context: StormContext):
    """The run's climate adjustment, or None where no climate config was given."""
    if not context.climate_config:
        return None
    if context.gwl is not None:
        return ClimateAdjustment(config_file=context.climate_config,
                                 method='gwl', gwl=float(context.gwl))
    if context.year is not None and context.ssp is not None:
        return ClimateAdjustment(config_file=context.climate_config,
                                 method='ssp', year=context.year, ssp=context.ssp)
    return None


def _rainfall_adjustment(storm, climate, context):
    """(factor for this duration, factor per duration) - Simulator.run_models:944."""
    if climate is None or context.exclusions.get('rainfall_uplift'):
        return 1.0, None
    factor = climate.get_rainfall_uplift_factor(duration=context.duration)
    per_duration = {dur: climate.get_rainfall_uplift_factor(duration=dur)
                    for dur in storm.rainfall.durations}
    return factor, per_duration


def hyetograph(storm, climate, context: StormContext, sim) -> Hyetograph:
    """Rebuild one realisation's rainfall from its mcdf row.

    ``sim`` is that row - a Series indexed by the mcdf's columns.
    """
    rain_z = float(sim['rain_z'])
    storm_method = str(sim['storm_method'])
    duration = context.duration

    factor, per_duration = _rainfall_adjustment(storm, climate, context)

    rain_depths = storm.rainfall.get_depth_z(z=rain_z, duration=duration,
                                             storm_method=storm_method)
    ave_rain = storm.get_average_rain(rain_depths) * factor

    temporal_pattern = storm.get_temporal_pattern(storm_method=storm_method,
                                                  duration=duration,
                                                  tp_sample=int(sim['tp']),
                                                  rain_sample_z=rain_z)
    timestep = storm.timesteps

    comment = None
    if context.filters_embedded_bursts:
        temporal_pattern, comment = storm.filter_temppat(
            temporal_pattern, rain_z, duration, ave_rain, storm_method,
            buffer=BURST_FILTER_BUFFER, climate_adjustment=per_duration)

    preburst_hours = 0.0
    proportion = float(sim.get('preburst_proportion') or 0.0)
    if context.applies_preburst and proportion > 0:
        if context.preburst_method == 'uniform':
            preburst_pattern = storm.uniform_preburst(duration, proportion, timestep)
        else:
            sample = int(sim['preburst_tp'])
            preburst_pattern, _ = storm.preburst_patterns.get_preburst_pattern(
                duration, proportion, timestep, sample)
            preburst_pattern, _ = storm.filter_embedded_bursts_in_preburst(
                preburst_pattern, temporal_pattern, duration, rain_z, ave_rain,
                storm_method, buffer=PREBURST_FILTER_BUFFER,
                climate_adjustment=per_duration)
        preburst_hours = float(preburst_pattern.index.min()) * -1.0
        temporal_pattern = pd.concat([preburst_pattern, temporal_pattern], axis=0)

    depths = temporal_pattern / 100.0 * ave_rain      # patterns are percentages
    depths.name = 'rainfall_mm'
    return Hyetograph(depths=depths, timestep=timestep,
                      preburst_hours=preburst_hours,
                      checks=_checks(depths, ave_rain, comment, sim))


def _checks(depths, ave_rain, comment, sim):
    """Hold the rebuild against the three things the run wrote down."""
    checks = []

    recorded = _number(sim.get('mean_rain_mm'))
    if recorded is not None:
        checks.append(('catchment average burst depth',
                       _close(ave_rain, recorded),
                       f'rebuilt {ave_rain:.2f} mm, run recorded {recorded:.2f} mm'))

    recorded = _number(sim.get('preburst_mm'))
    if recorded is not None:
        rebuilt = float(depths[depths.index < 0].sum())
        # An excluded or zero pre-burst is nothing in both, which agrees.
        ok = _close(rebuilt, recorded) or (abs(rebuilt) < 0.05 and abs(recorded) < 0.05)
        checks.append(('pre-burst depth', ok,
                       f'rebuilt {rebuilt:.2f} mm, run recorded {recorded:.2f} mm'))

    recorded = sim.get('embedded_bursts')
    if comment is not None and isinstance(recorded, str) and recorded.strip():
        ok = comment.strip() == recorded.strip()
        checks.append(('embedded burst filtering', ok,
                       f'rebuilt "{comment.strip()}", run recorded "{recorded.strip()}"'))
    return checks


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(value) else value


def _close(rebuilt, recorded, tolerance=DEPTH_TOLERANCE):
    if recorded == 0:
        return abs(rebuilt) < 0.05
    return abs(rebuilt / recorded - 1.0) <= tolerance
