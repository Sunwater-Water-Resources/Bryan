"""
lib/Volumes.py

The flood volume analysis: the largest volume within a moving window of each
analysis duration, in ML.

Shared by every method - Simulator.analyse_volumes for the Monte Carlo and
ensemble runs, ReservoirRoutingSimulator for the inflow volumes of a re-routed
run - so the window they integrate over cannot drift apart. It was written
twice before, and the reservoir routing copy is what showed up the timestep the
Monte Carlo one was short (see the note in rolling_max_volumes).

Analysis durations are moving-window lengths, not storm durations: the volume
that sizes the storage is usually drawn from a longer window than the burst
that produced it.
"""
import numpy as np


# The fixed list the Monte Carlo method integrates over, and what the other
# methods fall back to when nothing says otherwise.
DEFAULT_VOLUME_DURATIONS = [24, 36, 48, 72, 96, 120]


def rolling_max_volumes(flows_arr, dt_hours, durations):
    """Largest volume (ML) within a moving window of each duration, per simulation.

    flows_arr : (timesteps, simulations) array of flows in m³/s
    dt_hours  : the timestep of that array, in hours
    durations : the analysis durations, in hours

    Returns {duration: (simulations,) array}, skipping - with a note - any
    duration that is not a whole number of timesteps or that is longer than the
    hydrographs. Both used to be silent: the first dropped a column while the
    caller went on naming the columns from the full list, and the second wrote a
    column of NaN.

    The window is trapezoidal over ``duration / dt + 1`` samples, so it spans the
    full duration. Up to 21 August 2026 the Monte Carlo method used
    ``duration / dt`` samples, one timestep short of the duration it was labelled
    with, which biased every volume low by roughly ``dt / duration``.

    Vectorised over every simulation at once: the windowed integral is a
    difference of cumulative sums, so each duration is one pass over the array
    rather than a rolling().apply() per column.
    """
    # NaN -> 0: a hydrograph shorter than the analysis period is padded, and an
    # ensemble inflow file is padded after each duration's own simulation
    # period. A window reaching into the padding is genuinely carrying no flow.
    flows = np.nan_to_num(np.asarray(flows_arr, dtype=float), nan=0.0)
    n_steps = flows.shape[0]
    cumulative = np.concatenate([np.zeros((1, flows.shape[1])),
                                 np.cumsum(flows, axis=0)], axis=0)

    volumes = {}
    for duration in durations:
        steps = duration / dt_hours
        if not np.isclose(steps, round(steps)):
            print(f'  Unable to integrate the {duration} h duration with '
                  f'{dt_hours} h timesteps - skipped')
            continue
        points = int(round(steps)) + 1              # samples spanning the duration
        if points > n_steps:
            print(f'  The {duration} h duration is longer than the '
                  f'{(n_steps - 1) * dt_hours:g} h of hydrograph - skipped')
            continue
        # Trapezoidal integral of every window, in m³/s.h, then converted to ML.
        window_sum = cumulative[points:] - cumulative[:-points]
        half_ends = 0.5 * (flows[:n_steps - points + 1] + flows[points - 1:])
        window_volume = (window_sum - half_ends) * dt_hours * 3600.0 / 1000.0
        volumes[duration] = window_volume.max(axis=0)
        print(f'  {duration} h window ({points} timesteps): '
              f'volumes from {volumes[duration].min():,.0f} to '
              f'{volumes[duration].max():,.0f} ML')
    return volumes
