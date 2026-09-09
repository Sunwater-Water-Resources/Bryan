"""Rainfall depths for the regional model, region by region.

The upstream models run one catchment under one AEP. The regional model runs
eleven JPA rainfall regions at once, each at *its own* AEP conditional on how
rare the storm was over the dam catchment that produced the event - so almost
everything here exists because a quantity that is scalar upstream is a vector
here.

For one representative event, per region:

    conditional AEP  ->  point depth per subarea  ->  areal reduction on that
    region's own area  ->  spatial pattern  ->  depth

**The AEP is what the conditional analysis hands over, not the depth.** JPA
publishes both, but its depths are areal for the region's own area, and the
storm file needs per-subarea depths built from point IFD. Carrying the AEP and
reading the point tables avoids un-doing an ARF computed over a different area.
The two agree - JPA's conditional depth reproduces its own areal IFD read at the
conditional AEP to 0.4% - so nothing is lost.

**The level is read off JPA's areal IFD; the shape comes from the point
tables.** A region's conditional AEP indexes a curve JPA built from rainfall
averaged over that region's own polygon, so the areal reduction is already in
that curve, computed on the right area. Our point depths times our own ARF
reproduce it exactly below 1 in 2,000 - so reading it loses nothing there and
gains the rest of the range, which is tabulated to 1 in 2,000,000. Extrapolating
our own last two points at the PMP instead undershot the published depth by 12%
at 120 h.

``arf_area_km2`` on the crosswalk therefore records provenance and drives the
reported ARF rather than being applied to the level. It still matters, because
the Kroombit split breaks the one-area-per-region rule twice: the dam subareas
take KRO_dam's 329 km2, and the two below the dam take KRO_lower, whose AEP JPA
inherits from KRO and whose area is therefore KRO's 449 km2, not its own 120.

**Two regions are not a curve read at all.** KRO_lower has no marginal and no
PMP; JPA derives its depth as what is left of KRO once the dam catchment's share
is removed. Reading KRO's curve at KRO's AEP gives KRO's depth, a different
quantity, and was 6% out. The residual is computed here the way JPA computes it.

**Normalisation is per region.** ``WeightInterpolatedSpatialPattern`` normalises
over whatever subcatchment areas it is given; handing it the whole model would
move depth between regions, because each carries its own AEP and its own PMP.
So the spatial factors are built and applied one region at a time.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

from lib.ConfigPaths import resolve
from lib.StormGenerator import ArealReduction

def z_of(aep):
    """Standard normal variate of a '1 in X' AEP. Vectorised: this module works
    on eleven regions at once, where the upstream code works on one."""
    return ndtri(1 - 1 / np.asarray(aep, dtype=float))


def interp_in_z(aeps, values, aep, log=False):
    """Interpolate against AEP in standard normal variate space.

    '1 in X' is not a linear scale - at 1 in 2,000 a difference of 200 is
    nothing and at 1 in 100 it is everything - so every interpolation against
    AEP in this study is done on z. Depths are interpolated on log depth as
    well, which is what makes a straight line between two IFD points a
    reasonable growth curve rather than a chord across a curve.
    """
    order = np.argsort(z_of(aeps))
    x = z_of(np.asarray(aeps)[order])
    y = np.asarray(values, dtype=float)[order]
    if log:
        return 10 ** np.interp(z_of(aep), x, np.log10(y))
    return np.interp(z_of(aep), x, y)


class DownstreamRainfall:
    """Per-subarea rainfall depths for the regional model."""

    def __init__(self, config_file):
        with open(config_file) as handle:
            cfg = json.load(handle)
        folder = os.path.dirname(config_file)
        paths = {k: resolve(folder, v) for k, v in cfg['file_paths'].items()}
        self.config = cfg
        self.paths = paths

        self.subareas = pd.read_csv(paths['crosswalk']).set_index('new_id').sort_index()
        self.conditional = (pd.read_csv(paths['conditional_aep'])
                            .set_index(['Duration_h', 'Driver_AEP_1in']))
        self.pattern = pd.read_csv(paths['pmp_scaling']).set_index('ID_Number').sort_index()
        self.region_pmp = pd.read_csv(paths['pmp_by_region'], index_col=0)
        self.arr_file = paths['arr_datahub_file']
        self.ifd_folder = paths['ifd_folder']
        self.bounds = cfg['storm_method_config']['aep_changeover_to_extreme']
        self.pmp_aep = cfg['pmp_aep_by_region']
        self.mass_balance = cfg.get('mass_balance', {})
        self.areal_folder = paths['areal_ifd_folder']
        self._ifd = {}
        self._areal = {}
        self._arf = {}
        self._last_conditional = None

    # -- inputs ------------------------------------------------------------
    def ifd(self, duration):
        """Point depths, a row per subarea and a column per AEP."""
        if duration not in self._ifd:
            path = os.path.join(self.ifd_folder, f'average_ifd_{duration:g}.csv')
            table = pd.read_csv(path).set_index('name').sort_index()
            table.columns = [float(c.replace('_AEP', '')) for c in table.columns]
            self._ifd[duration] = table
        return self._ifd[duration]

    def areal(self, duration):
        """JPA's region-average areal IFD, indexed by AEP as a number."""
        if duration not in self._areal:
            table = pd.read_csv(os.path.join(self.areal_folder,
                                             f'ArealIFD_{duration:03.0f}h.csv'))
            table = table.set_index(table.columns[0])
            table.index = [float(str(i).replace('1 in ', '')) for i in table.index]
            self._areal[duration] = table
        return self._areal[duration]

    def arf(self, area, duration, aep):
        if area not in self._arf:
            self._arf[area] = ArealReduction(area=area, arr_datahub_file=self.arr_file)
        return self._arf[area].get_areal_reduction_factor(duration=duration, aep=aep)

    def conditional_aep(self, duration, driver_aep):
        """Conditional AEP per JPA target region, at an arbitrary driver AEP."""
        block = self.conditional.xs(duration, level='Duration_h').sort_index()
        zx = z_of(block.index.to_numpy())
        z = z_of(driver_aep)
        return pd.Series({target: float(1 / (1 - ndtr(np.interp(z, zx, z_of(block[target].to_numpy())))))
                          for target in block.columns})

    # -- the pattern -------------------------------------------------------
    def spatial_factors(self, members, duration, aep, storm_method):
        """Within-region factors, area-weighting to 1.0, at this AEP.

        Below the changeover the IFD's own shape is the pattern. Above it the
        pattern is the PMP one. Between, the two are blended on z - the
        ``interpolate_weights`` method, applied to one region rather than to a
        whole catchment.
        """
        area = self.subareas.loc[members, 'area_km2']
        lower, upper = self.bounds
        ifd = self.ifd(duration)

        def normalise(values):
            return values / float((values * area).sum() / area.sum())

        implied = normalise(ifd.loc[members, lower])
        if aep <= lower:
            return implied
        column = 'GTSMR' if storm_method == 'GTSMR' else f'GSDM_{duration:g}'
        if column not in self.pattern.columns:
            column = 'GSDM'
        extreme = normalise(self.pattern.loc[members, column])
        if aep >= upper:
            return extreme
        w = (z_of(aep) - z_of(lower)) / (z_of(upper) - z_of(lower))
        return normalise((1 - w) * implied + w * extreme)

    def region_depth(self, target, duration, aep):
        """Region-average **areal** depth at this AEP.

        Read off JPA's own areal IFD, which is tabulated from 1 EY to 1 in
        2,000,000 - so the extreme range comes from the fitted curve rather than
        from extrapolating our own last two points at it, which undershot the
        published depth by 12% at 120 h.

        This is the region's ARF already applied, on the region's own area:
        our point depths times our ARF reproduce this table exactly, so nothing
        is lost by reading it, and the two agree at every AEP rather than only
        below the last one we tabulate.
        """
        if target in self.mass_balance:
            return self._by_mass_balance(target, duration, aep)
        table = self.areal(duration)
        return float(interp_in_z(table.index.to_numpy(), table[target].to_numpy(),
                                 min(aep, self.pmp_aep[target]), log=True))

    def _by_mass_balance(self, target, duration, aep):
        """A region JPA derives as a residual rather than fitting a curve to.

        KRO_lower is the piece of KRO below the dam. It has no marginal and no
        PMP of its own - JPA inherits its AEP from KRO purely as bookkeeping -
        and takes its depth from what is left once the dam catchment's share is
        removed: D_lower = (A_KRO*D_KRO - A_dam*D_dam) / A_lower. Reading KRO's
        curve at KRO's AEP instead gives KRO's depth, which is a different
        quantity and was 6% out.
        """
        parts = self.mass_balance[target]
        whole, taken = parts['whole'], parts['subtract']
        conditional = self._last_conditional
        # The whole region has no conditional column of its own here - the
        # residual's column *is* its AEP, which is what "inherited" means.
        d_whole = self.region_depth(whole, duration,
                                    float(conditional[parts['whole_aep_from']]))
        d_taken = self.region_depth(taken, duration, float(conditional[taken]))
        a_whole, a_taken = parts['area_whole'], parts['area_subtract']
        return (a_whole * d_whole - a_taken * d_taken) / (a_whole - a_taken)

    def _area_for(self, members):
        areas = self.subareas.loc[members, 'arf_area_km2'].unique()
        if len(areas) != 1:
            raise ValueError(f'one ARF area per group expected, got {areas}')
        return float(areas[0])

    # -- the event ---------------------------------------------------------
    def depths(self, duration, driver_aep, storm_method='ARR Areal'):
        """Areal rainfall depth for every subarea, for one representative event."""
        conditional = self.conditional_aep(duration, driver_aep)
        self._last_conditional = conditional
        out = pd.Series(index=self.subareas.index, dtype=float)
        for target, members in self.subareas.groupby('cond_target').groups.items():
            aep = float(conditional[target])
            level = self.region_depth(target, duration, aep)
            out[members] = level * self.spatial_factors(members, duration, aep, storm_method)
        return out

    def region_summary(self, duration, driver_aep, storm_method='ARR Areal'):
        """What each region got, for the record and for checking."""
        conditional = self.conditional_aep(duration, driver_aep)
        self._last_conditional = conditional
        rows = []
        for target, members in self.subareas.groupby('cond_target').groups.items():
            aep = float(conditional[target])
            area = self.subareas.loc[members, 'area_km2']
            depth = self.depths(duration, driver_aep, storm_method)[members]
            rows.append(dict(
                cond_target=target, subareas=len(members),
                conditional_aep=round(aep, 1),
                arf_area_km2=self._area_for(members),
                arf=round(self.arf(self._area_for(members), duration, aep), 4),
                mean_depth_mm=round(float((depth * area).sum() / area.sum()), 1),
                min_mm=round(float(depth.min()), 1), max_mm=round(float(depth.max()), 1)))
        return pd.DataFrame(rows).set_index('cond_target')
