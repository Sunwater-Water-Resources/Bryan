"""Route a derived net inflow through one constant rating curve.

The method is two-stage:

1.  Derive a net inflow series from the observed record by closing the water
    balance backwards, using the rating curve that was actually in force at each
    timestep::

        inflow = dStorage + spillway release + net evaporation

    Above full supply this asserts the rating that was in force, and on the
    falling limb of every event since 2011 the assertion fails: the derived
    inflow goes strongly negative, which is a release the rating did not
    account for wearing an inflow's sign.  Those steps are reattributed to
    release before stage 2 sees them -- see
    ``_reattribute_unaccounted_release``.  Without it the lake is charged for
    the same drawdown twice and the homogenised record sits up to half a metre
    low for decades at a time.

2.  Re-route that inflow through a single, constant rating curve, solving the
    level-pool routing equation forward in time.

Both stages run at the **native resolution of the gauge**, not daily.  That
matters for the peaks.  Deriving the net inflow from daily maxima attenuates the
peak by a median of 2.3x across the water years that reach full supply, because
a rise that takes eight hours is smeared across twenty-four.  The lake rises
0.4 m in fifteen minutes at Callide; a daily timestep cannot see that.

Antecedent storage was never affected by this and the daily answer for it stands.

Stage 2 is solved *implicitly*.  Callide's gated spillway is very steep relative
to the lake's surface area, and the URBS gate operation in a ``.sq`` file is
steeper still -- it releases 142 m3/s at 120 ML above full supply, which is
10 mm of level.  An explicit forward difference is unstable against that by orders
of magnitude at any timestep worth running.
"""

from __future__ import annotations

import logging
import time

import numpy as np
import pandas as pd

log = logging.getLogger("bryan.homogenise.model")

SECONDS_PER_DAY = 86400.0

# Implicit solve tolerances.  The volume tolerance is 0.1 litre; the level
# tolerance is the floor for the bracket, needed because the sq rating's
# gate-operation segment makes dQ/dLevel reach 2.5e6 m3/s/m and the residual
# cannot then be driven to the volume tolerance in double precision.
VOLUME_TOLERANCE_ML = 1e-7
LEVEL_TOLERANCE_M = 1e-10
MAX_SOLVER_ITERATIONS = 200


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------

class LakeModel:
    def __init__(self, level_record, evaporation, register, ratings,
                 volume_of_level, area_of_level, start=None, end=None,
                 recession_correction=True):
        self.record = level_record
        self.recession_correction = recession_correction
        if start is not None or end is not None:
            self.record = self.record.loc[start:end]
            if self.record.empty:
                raise ValueError(f"no level data between {start} and {end}")
            log.info("Restricted to %s -> %s (%d steps)",
                     self.record.index.min(), self.record.index.max(),
                     len(self.record))
        self.levels = self.record["Level"]
        self.evaporation = evaporation
        self.register = register
        self.ratings = ratings
        self.volume_of_level = volume_of_level
        self.area_of_level = area_of_level

        self._check_coverage()
        self.rating_ids, self.fsls = self._assign_ratings()

    # -- setup -------------------------------------------------------------

    def _check_coverage(self):
        index = self.levels.index
        if not self.evaporation.covers(index):
            raise ValueError(
                f"evaporation record ({self.evaporation.start} to "
                f"{self.evaporation.end}) does not cover the level record "
                f"({index.min()} to {index.max()})"
            )
        if (index.min() < self.register["from"].min()
                or index.max() > self.register["to"].max()):
            raise ValueError(
                f"rating register ({self.register['from'].min().date()} to "
                f"{self.register['to'].max().date()}) does not cover the level "
                f"record ({index.min().date()} to {index.max().date()})"
            )
        lo, hi = self.volume_of_level.lo, self.volume_of_level.hi
        outside = self.levels[(self.levels < lo) | (self.levels > hi)]
        if len(outside):
            raise ValueError(
                f"{len(outside)} recorded levels fall outside the storage table "
                f"({lo:.2f} to {hi:.2f} m AHD)"
            )

    def _assign_ratings(self):
        """Map each timestep to its rating and FSL using half-open intervals.

        The register's ``to`` for one rating equals the ``from`` of the next, so
        closed intervals would assign the changeover instant twice.
        """
        edges = self.register["from"].tolist() + [self.register["to"].iloc[-1]]
        position = np.searchsorted(pd.DatetimeIndex(edges), self.levels.index,
                                   side="right") - 1
        position = np.clip(position, 0, len(self.register) - 1)
        rating_ids = pd.Series(self.register.index.to_numpy()[position],
                               index=self.levels.index, name="Rating ID")
        fsls = pd.Series(self.register["FSL"].to_numpy()[position],
                         index=self.levels.index, name="FSL", dtype=float)
        log.info("Steps per rating: %s",
                 rating_ids.value_counts().sort_index().to_dict())
        return rating_ids, fsls

    # -- stage 1: derive inflow -------------------------------------------

    def derive_inflow(self):
        """Close the water balance backwards to get net inflow.

        ``Inflow_ML`` is a *net* inflow: it absorbs rainfall on the lake,
        seepage, and any extractions, because those are not modelled separately.
        Fully vectorised -- only stage 2 needs a loop.
        """
        record = self.record.copy()
        record["Rating ID"] = self.rating_ids
        record["FSL"] = self.fsls
        record["Volume_ML"] = self.volume_of_level(record["Level"])
        record["Area_ha"] = self.area_of_level(record["Level"])

        seconds = record.index.to_series().diff().dt.total_seconds()
        record["dt_s"] = seconds

        # Historical spillway discharge, using the rating in force at the time.
        flow = pd.Series(0.0, index=record.index, name="Release_flow_m3s")
        for rating_id, group in record.groupby("Rating ID"):
            rating = self.ratings[rating_id]
            spilling = group["Level"] > group["FSL"]
            if spilling.any():
                levels = group.loc[spilling, "Level"].to_numpy()
                flow.loc[group.index[spilling]] = rating(levels)
        record["Release_flow_m3s"] = flow
        record["Release_ML"] = (
            0.5 * (flow + flow.shift(1)) * record["dt_s"] / 1000.0
        )

        # Exact integral of the 9am-stepped evaporation over each timestep.
        record["Evaporation_mm"] = self.evaporation.over(record.index)
        mean_area = 0.5 * (record["Area_ha"] + record["Area_ha"].shift(1))
        # mm * ha / 100 == ML
        record["Evaporation_ML"] = record["Evaporation_mm"] * mean_area / 100.0

        record["dStorage_ML"] = record["Volume_ML"].diff()
        record["Inflow_ML"] = (
            record["dStorage_ML"] + record["Release_ML"] + record["Evaporation_ML"]
        )

        first = record.index[0]
        record.loc[first, ["Release_ML", "Evaporation_ML", "Evaporation_mm",
                           "dStorage_ML", "Inflow_ML"]] = 0.0
        record.loc[first, "dt_s"] = 0.0

        unresolved = record["Inflow_ML"].isna()
        if unresolved.any():
            raise ValueError(
                f"{int(unresolved.sum())} steps have an unresolved water balance "
                f"(first {record.index[unresolved][0]})"
            )
        if self.recession_correction:
            self._reattribute_unaccounted_release(record)
        else:
            record["Unaccounted_Release_ML"] = 0.0
            record["Unaccounted_Release_m3s"] = 0.0
        self._report_clamps()
        return record

    def _reattribute_unaccounted_release(self, record):
        """Book an unexplained above-full-supply drawdown as release, not as
        negative inflow.

        Above full supply the balance asserts the spillway rating, so a step
        where it cannot close is evidence about that rating.  Ratings 3-6 are
        one table with the first row slid down the level axis -- rating 4 runs
        ``214.000 -> 0`` straight to ``216.300 -> 196`` -- so the whole
        operational band is unconstrained interpolation, and the gates did not
        track it.  On the falling limb of every event since 2011 the derived
        inflow goes strongly negative, to -459 m3/s in January 2011 and
        -173 m3/s in February 2015, two days after a cyclone with the creek
        still delivering.  Adding the deficit back as release resolves the
        February 2015 recession to a flat ~140 m3/s held from 216.1 m down to
        215.4 m: a gate at a fixed opening, not a taper following the level.

        Left alone, that deficit is a *release* wearing an inflow's sign, and
        stage 2 charges the lake for it twice -- once inside the inflow, once
        through the new rating -- so the homogenised lake drives through full
        supply instead of settling onto it and, with no spill afterwards to
        reset it, carries the error to the end of the record.  Across regime 6
        that cost 0.26 to 0.47 m against a true gate-rule difference of 0.01 m.

        Deliberately scoped to steps **above** the full supply in force.  Below
        it the rating asserts nothing (release is zero by construction), so a
        negative net inflow there is ordinary: seepage and extraction exceeding
        inflow, 1.34 GL of it across the record, which is physical and is left
        alone.  Correcting below full supply as well swings regime 6 to
        +0.22 m, because it credits the homogenised lake with water the dam
        released through outlets this method cannot see.

        The correction reattributes; it does not invent.  ``Release_ML`` gains
        the deficit, ``Inflow_ML`` loses it, and the sum is unchanged, so the
        recorded lake is still reproduced exactly by the terms it was derived
        from.  ``Release_flow_m3s`` is left as the tabulated rating value and
        the reconstruction is reported beside it, so the two stay separable.
        """
        unaccounted = (record["Level"] > record["FSL"]) & (record["Inflow_ML"] < 0.0)
        record["Unaccounted_Release_ML"] = 0.0
        record["Unaccounted_Release_m3s"] = 0.0
        if not unaccounted.any():
            log.info("Recession correction: no step needed one")
            return record

        deficit = -record.loc[unaccounted, "Inflow_ML"]
        record.loc[unaccounted, "Unaccounted_Release_ML"] = deficit
        record.loc[unaccounted, "Unaccounted_Release_m3s"] = (
            deficit * 1000.0 / record.loc[unaccounted, "dt_s"]
        )
        record.loc[unaccounted, "Release_ML"] += deficit
        record.loc[unaccounted, "Inflow_ML"] = 0.0

        log.info("Recession correction: %d of %d steps above the full supply in "
                 "force could not close the balance; %s ML reattributed from "
                 "inflow to release (peak %.0f m3/s)",
                 int(unaccounted.sum()), len(record), f"{deficit.sum():,.0f}",
                 record["Unaccounted_Release_m3s"].max())
        return record

    # -- stage 2: re-route through a single rating ------------------------

    def simulate(self, rating, record=None, start_level=None):
        """Route the derived inflow through one constant rating curve.

        Solves, for each step, the trapezoidal level-pool equation

            V(L[i]) = V(L[i-1]) + I[i] - E[i](L[i]) - R[i](L[i])

        implicitly for the *level* ``L[i]``.  Storage, area and discharge are all
        monotone increasing in level, so the residual is monotone and the root is
        unique.  Grouping the level-dependent terms,

            G(L) = V(L) + c1.Q(L) + c2.A(L) = K

        with ``c1 = dt/2000`` and ``c2 = E/200`` known constants for the step, so
        the solve is a monotone root-find on ``G``.  It runs safeguarded Newton:
        Newton where the step stays inside the bracket, bisection where it does
        not.  Pure bisection would need forty halvings per step across half a
        million steps; unsafeguarded Newton would stall on the sq rating's
        gate-operation segment, where dQ/dLevel reaches 2.5e6 m3/s per metre.

        Working in level space rather than volume space keeps the scheme exactly
        self-consistent: when the new rating matches the rating that was in
        force, the simulation reproduces the recorded level to solver tolerance
        rather than drifting.  ``start_level`` overrides the initial condition,
        which is useful when matching another model's spin-up state.
        """
        if record is None:
            record = self.derive_inflow()

        inflow = record["Inflow_ML"].to_numpy()
        step_seconds = np.nan_to_num(record["dt_s"].to_numpy())
        step_evaporation = np.nan_to_num(record["Evaporation_mm"].to_numpy())

        volume_of_level, area_of_level = self.volume_of_level, self.area_of_level
        volume_at, volume_slope = volume_of_level.at, volume_of_level.slope_at
        area_at, area_slope = area_of_level.at, area_of_level.slope_at
        flow_at, flow_slope = rating.at, rating.slope_at
        floor, ceiling = volume_of_level.lo, volume_of_level.hi

        count = len(record)
        levels = np.empty(count)
        volumes = np.empty(count)
        flows = np.zeros(count)
        releases = np.zeros(count)
        evaporation = np.zeros(count)

        level = float(record["Level"].iloc[0] if start_level is None else start_level)
        levels[0] = level
        volumes[0] = volume_at(level)
        flows[0] = flow_at(level)

        previous_level = level
        previous_volume = volumes[0]
        previous_flow = flows[0]
        previous_area = area_at(level)

        iterations = 0
        pinned = 0
        started = time.perf_counter()

        for i in range(1, count):
            c1 = 0.5 * step_seconds[i] / 1000.0
            c2 = step_evaporation[i] / 200.0
            target = (previous_volume + inflow[i]
                      - c1 * previous_flow - c2 * previous_area)

            if (volume_at(floor) + c1 * flow_at(floor)
                    + c2 * area_at(floor) - target) >= 0.0:
                level = floor
            elif (volume_at(ceiling) + c1 * flow_at(ceiling)
                    + c2 * area_at(ceiling) - target) <= 0.0:
                level = ceiling
                pinned += 1
            else:
                low, high = floor, ceiling
                level = previous_level
                if not low < level < high:
                    level = 0.5 * (low + high)
                for _ in range(MAX_SOLVER_ITERATIONS):
                    iterations += 1
                    residual = (volume_at(level) + c1 * flow_at(level)
                                + c2 * area_at(level) - target)
                    if residual < 0.0:
                        low = level
                    else:
                        high = level
                    if -VOLUME_TOLERANCE_ML < residual < VOLUME_TOLERANCE_ML:
                        break
                    if high - low < LEVEL_TOLERANCE_M:
                        level = 0.5 * (low + high)
                        break
                    gradient = (volume_slope(level) + c1 * flow_slope(level)
                                + c2 * area_slope(level))
                    if gradient > 0.0:
                        candidate = level - residual / gradient
                        if not low < candidate < high:
                            candidate = 0.5 * (low + high)
                    else:
                        candidate = 0.5 * (low + high)
                    if abs(candidate - level) < LEVEL_TOLERANCE_M:
                        level = candidate
                        break
                    level = candidate

            volume = volume_at(level)
            area = area_at(level)
            flow = flow_at(level)
            levels[i] = level
            volumes[i] = volume
            flows[i] = flow
            releases[i] = c1 * (previous_flow + flow)
            evaporation[i] = c2 * (previous_area + area)
            previous_level, previous_volume = level, volume
            previous_flow, previous_area = flow, area

        elapsed = time.perf_counter() - started
        log.info("Routed %d steps in %.1f s (%.1f solver iterations per step)",
                 count, elapsed, iterations / max(count - 1, 1))
        if pinned:
            log.warning("%d steps pinned at the top of the storage table "
                        "(%.2f m AHD)", pinned, ceiling)

        result = record.copy()
        result["New_FSL"] = rating.fsl
        result["New_Level"] = levels
        result["New_Volume_ML"] = volumes
        result["New_Storage_above_FSL_ML"] = volumes - float(volume_of_level(rating.fsl))
        result["New_Release_flow_m3s"] = flows
        result["New_Release_ML"] = releases
        result["New_Evaporation_ML"] = evaporation
        self._report_clamps(extra=[rating])

        # The solver probes the top of the storage table to bracket the root, so
        # the rating's own clamp counter cannot distinguish a probe from a real
        # extrapolation.  Ask the question of the answer instead.
        top_level = rating.top_level
        beyond = result["New_Level"] > top_level
        if beyond.any():
            log.warning(
                "%d steps routed above the top of %s (%.3f m AHD, %.0f m3/s); "
                "the release is held flat there and the peak is a lower bound",
                int(beyond.sum()), rating.name, top_level, rating.top_flow,
            )
        else:
            log.info("Peak %.3f m AHD stays within %s, which tops out at %.3f m",
                     result["New_Level"].max(), rating.name, top_level)
        return result

    # -- diagnostics -------------------------------------------------------

    def _report_clamps(self, extra=()):
        holders = [self.volume_of_level, self.area_of_level]
        holders += list(self.ratings.values()) + list(extra)
        for holder in holders:
            for name, count in getattr(holder, "clamps", {}).items():
                log.warning("%s: %d evaluations clamped to the table endpoints",
                            name, count)
            getattr(holder, "clamps", {}).clear()


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

TRACE_COLUMNS = [
    "Gauge", "Interpolated", "Rating ID", "FSL", "Level", "Volume_ML",
    "Release_flow_m3s", "Unaccounted_Release_m3s", "Inflow_ML",
    "New_FSL", "New_Level", "New_Volume_ML",
    "New_Storage_above_FSL_ML", "New_Release_flow_m3s",
]


def daily_summary(result):
    """Collapse the sub-daily trace to one row per day, on true calendar days."""
    daily = pd.DataFrame({
        "Level_max": result["Level"].resample("D").max(),
        "Level_mean": result["Level"].resample("D").mean(),
        "New_Level_max": result["New_Level"].resample("D").max(),
        "New_Level_mean": result["New_Level"].resample("D").mean(),
        "New_Volume_ML_max": result["New_Volume_ML"].resample("D").max(),
        "Inflow_ML": result["Inflow_ML"].resample("D").sum(),
        "New_Release_ML": result["New_Release_ML"].resample("D").sum(),
        "New_Release_flow_m3s_max": result["New_Release_flow_m3s"].resample("D").max(),
        "New_Evaporation_ML": result["New_Evaporation_ML"].resample("D").sum(),
        "Observed_steps": (~result["Interpolated"]).resample("D").sum(),
    })
    return daily.dropna(subset=["Level_max"])


def summarise(result, rating):
    spilling = result["New_Release_flow_m3s"] > 0
    new_peak = result["New_Level"].idxmax()
    log.info("Routed through %s (FSL %.3f m AHD)", rating.name, rating.fsl)
    log.info("Recorded peak   %.3f m AHD at %s",
             result["Level"].max(), result["Level"].idxmax())
    log.info("Homogenised peak %.3f m AHD at %s (%.0f ML above FSL, release "
             "%.0f m3/s)", result["New_Level"].max(), new_peak,
             result["New_Storage_above_FSL_ML"].max(),
             result["New_Release_flow_m3s"].max())
    log.info("Release above FSL on %d of %d steps (%.1f%% of the record by time)",
             int(spilling.sum()), len(result),
             100.0 * result.loc[spilling, "dt_s"].sum() / result["dt_s"].sum())
    log.info("Total net inflow %.0f GL over %.1f years",
             result["Inflow_ML"].sum() / 1000.0,
             (result.index.max() - result.index.min()).days / 365.25)




# The name it had in callide-fsl-reinstate, where it was written.
CallideDam = LakeModel
