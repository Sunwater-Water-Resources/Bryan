"""Frequency curves through the annual maximum lake levels, and their uncertainty.

The second half of the lake level frequency analysis. ``lib/LakeLevelRecord.py``
reads the record and derives the annual maxima with pandas alone, so the run
launcher can show them; this module fits curves through them, resamples the
record for the uncertainty of each curve, and reads the design flood
realisations the curves are compared against. It needs scipy, so the launcher
never imports it - ``util/LakeLevelFrequency.py`` runs it and writes the
results the page draws.

**Why no standard distribution.** An annual maximum lake level is bounded and
lumpy in a way inflows are not. A gated or ungated storage spends a run of years
sitting at full supply - held there by the gates or by the spillway crest - so
the frequency curve carries a *plateau* no single sigmoid or GEV can follow, and
above it the level is driven by the flood rather than by the operation. Two
forms are offered:

``shouldered``
    Three pieces, each fitted to the maxima that belong to it: a polynomial
    *shoulder* below full supply, pinned to reach it; the *plateau*, at full
    supply over the extent the record itself gives; and an *upper limb* through
    the maxima above the plateau. Nothing about the plateau's extent is fitted -
    it is read off the record (:func:`plateau_span`). A plateau tolerance of zero
    drops the plateau, and the shoulder and upper limb then meet at full supply
    between the last maximum at or below it and the first above. The upper limb
    is a straight line by default; a dam that spills in most years has the
    maxima above full supply to carry a higher degree, and the degree is refused
    where it does not (:data:`MAXIMA_PER_UPPER_TERM`). It either starts free,
    leaving the step a gated dam shows where the gates stop holding the lake, or
    starts at full supply, continuous, as an uncontrolled spillway would give.
    Developed on the Callide record; the tolerances are parameters for that reason.
``logistic``
    A four-parameter logistic with the ceiling free. It drives straight through
    a plateau, but it needs no full supply level and nothing sitting on it, so it
    is the fallback when the shouldered form cannot be placed.

**The band is resampled in the form the curve is drawn in.** Each resample
draws the annual maxima with replacement, re-derives the plotting positions -
they are a property of the sample, not labels attached to the values - and
refits. A logistic band around a shouldered curve would be the sampling
uncertainty of a different model from the one on the page.

**Carried-over maxima.** A water year's maximum is not always a flood: the lake
can open the year at its highest and fall from there. Those years are drawn
hollow and a second curve is fitted to the storm-driven maxima alone, placed on
the whole record's probability scale (``LakeLevelRecord.censored_positions``).
Its band resamples years whole, storm-driven or carried over, so the number of
storm-driven years varies from resample to resample as it would in another
record of the same length.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from lib import LakeLevelRecord as record

SHOULDERED = "shouldered"
LOGISTIC = "logistic"
NO_FIT = "none"
FORMS = (SHOULDERED, LOGISTIC, NO_FIT)

SHOULDER_DEGREE = 4

# Maxima this close to full supply are taken to be sitting on it: the lake held
# at its operating level rather than driven by a flood. Zero means no plateau.
PLATEAU_TOLERANCE = 0.025

# The plateau does not end where that tolerance does. Past the maxima sitting
# exactly on full supply a record carries a few more only slightly above it, and
# then a clear jump - where the lake starts being driven by the flood rather than
# held. The plateau is extended through steps smaller than this and ended at the
# first one larger.
PLATEAU_GAP = 0.100

# The upper limb, and where it starts.
UPPER_DEGREE = 1
FREE = "free"            # its own intercept: a step up from the plateau is allowed
AT_FSL = "fsl"           # starts at full supply where the plateau ends: continuous
UPPER_JOINS = (FREE, AT_FSL)

# A curved upper limb needs this many maxima above the plateau per coefficient.
# The straight line keeps the bare minimum of one more maximum than it has
# coefficients, which is what the Callide figures were drawn with. Callide and
# Kroombit carry 6 and 4-7 maxima above full supply, and a quadratic through
# them left the rmse unchanged to the millimetre while halving the resamples it
# could be fitted to - flexibility the record cannot pay for. A dam that spills
# most years can.
MAXIMA_PER_UPPER_TERM = 4

# Below this share of resamples fitted, the band is drawn from an unrepresentative
# subset - the resamples that happened to keep the maxima the form depends on -
# and is narrower than the record justifies.
LOW_FIT_SHARE = 0.7

# Fixed so a figure is reproducible byte for byte. 400 draws place a 5th and 95th
# percentile to well inside the width of the line they are drawn with.
BOOTSTRAP_SEED = 20260826
BOOTSTRAP_DRAWS = 400
BAND_PERCENTILES = (5.0, 95.0)

# A storm-driven resample with fewer maxima than this is not refitted.
MINIMUM_STORM_YEARS = 12

EY1_AEP = float(1.0 - math.exp(-1.0))


# -- the shouldered plateau -----------------------------------------------------

def plateau_span(z, level, fsl, tol=PLATEAU_TOLERANCE, gap=PLATEAU_GAP):
    """``(z1, z2)`` bounding the plateau, both read from the record.

    ``z1`` is the most frequent maximum within ``tol`` of full supply. ``z2``
    walks on to rarer maxima while each is within ``gap`` of the one before, and
    stops at the first larger step.
    """
    z, level = np.asarray(z, float), np.asarray(level, float)
    order = np.argsort(z)
    zs, ls = z[order], level[order]
    on = np.abs(ls - fsl) <= tol
    if not on.any():
        raise ValueError(f"no annual maximum within {tol * 1000:.0f} mm of full "
                         f"supply ({fsl:.3f} m), so there is no plateau to place")
    first, last = int(np.argmax(on)), int(len(on) - 1 - np.argmax(on[::-1]))
    while last + 1 < len(ls) and (ls[last + 1] - ls[last]) <= gap:
        last += 1
    return float(zs[first]), float(zs[last])


def _monotone_fit(evaluate_with, start, zz, yy, grid, bounds, what):
    """Least squares under a non-decreasing constraint, from an unconstrained start."""
    from scipy.optimize import minimize

    def rmse(c):
        try:
            return float(np.sqrt(np.mean((evaluate_with(zz, c) - yy) ** 2)))
        except Exception:
            return 1e9

    got = minimize(rmse, list(start), bounds=bounds, method="SLSQP",
                   constraints=[{"type": "ineq",
                                 "fun": lambda c: np.diff(evaluate_with(grid, c)).min()}],
                   options={"maxiter": 800, "ftol": 1e-13})
    if not got.success:
        raise RuntimeError(f"the {what} fit did not converge")
    return got.x


def _upper(z, z2, fsl, coef, join):
    """The upper limb: a polynomial in (z - z2), its constant fixed at full supply
    when it starts there."""
    d = np.asarray(z, float) - z2
    if join == AT_FSL:
        return fsl + sum(c * d ** (i + 1) for i, c in enumerate(coef))
    return sum(c * d ** i for i, c in enumerate(coef))


def fit_shouldered_plateau(z, level, fsl, degree=SHOULDER_DEGREE,
                           tol=PLATEAU_TOLERANCE, gap=PLATEAU_GAP,
                           upper_degree=UPPER_DEGREE, upper_join=FREE):
    """Fit shoulder, plateau and upper limb, each to its own maxima.

    The plateau ends at the last maximum actually sitting on it, which matches
    the record's own extent; with ``upper_join="free"`` that leaves a step up to
    the upper limb. ``tol=0`` fits no plateau at all. Both pieces are held
    non-decreasing.

    Returns ``(params, rmse)``: ``rmse`` is over every maximum given, and
    ``params`` is the dict :func:`shouldered_plateau` evaluates.
    """
    if upper_join not in UPPER_JOINS:
        raise ValueError(f"upper_join must be one of {UPPER_JOINS}, not {upper_join!r}")
    upper_degree = int(upper_degree)
    if upper_degree < 1:
        raise ValueError("the upper limb needs a degree of at least 1")
    z, level = np.asarray(z, float), np.asarray(level, float)
    if tol > 0:
        z1, z2 = plateau_span(z, level, fsl, tol, gap)
    else:
        at_or_below, over = z[level <= fsl], z[level > fsl]
        if not len(at_or_below) or not len(over):
            raise ValueError("with no plateau the maxima have to lie on both sides of "
                             "full supply")
        z1 = z2 = float((at_or_below.max() + over.min()) / 2)
    above = z > z2
    below = z < z1

    terms = upper_degree + (1 if upper_join == FREE else 0)
    needed = terms + 1 if upper_degree == 1 else MAXIMA_PER_UPPER_TERM * terms
    if above.sum() < needed:
        raise ValueError(f"{int(above.sum())} annual maxima above the plateau; a degree "
                         f"{upper_degree} upper limb "
                         f"{'starting at full supply ' if upper_join == AT_FSL else ''}"
                         f"needs at least {needed}")

    # The upper limb: one linear least squares, constrained only if it turns back.
    za, ya = z[above], level[above]
    da = za - z2
    if upper_join == AT_FSL:
        columns, target = [da ** (i + 1) for i in range(upper_degree)], ya - fsl
    else:
        columns, target = [da ** i for i in range(upper_degree + 1)], ya
    upper_coef, *_ = np.linalg.lstsq(np.vstack(columns).T, target, rcond=None)
    upper_grid = np.linspace(z2, float(za.max()), 400)
    if np.diff(_upper(upper_grid, z2, fsl, upper_coef, upper_join)).min() < -1e-9:
        upper_coef = _monotone_fit(
            lambda zz, c: _upper(zz, z2, fsl, c, upper_join), upper_coef, za, ya,
            upper_grid, [(None, None)] * len(upper_coef), "upper limb")

    # Pinned to reach full supply at z1, the shoulder is full supply plus a
    # polynomial in (z - z1) with no constant term - one linear least squares
    # rather than a search. Only when that answer turns back on itself is the
    # monotonic constraint imposed, starting from it.
    zb, yb = z[below], level[below]
    grid = np.linspace(float(zb.min()), z1, 400)

    def shoulder(zz, coef):
        d = np.asarray(zz, float) - z1
        return fsl + sum(c * d ** (i + 1) for i, c in enumerate(coef))

    d = zb - z1
    design = np.vstack([d ** (i + 1) for i in range(degree)]).T
    coef, *_ = np.linalg.lstsq(design, yb - fsl, rcond=None)
    if np.diff(shoulder(grid, coef)).min() < -1e-9:
        coef = _monotone_fit(shoulder, coef, zb, yb, grid, [(-80, 80)] * degree,
                             "shoulder")

    params = {"coef": [float(c) for c in coef], "z1": z1, "z2": z2,
              "upper_coef": [float(c) for c in upper_coef],
              "upper_join": upper_join, "fsl": float(fsl),
              "step": float(_upper(z2, z2, fsl, upper_coef, upper_join) - fsl)}
    resid = shouldered_plateau(z, params) - level
    return params, float(np.sqrt(np.mean(resid ** 2)))


def shouldered_plateau(z, params):
    """Evaluate a fit returned by :func:`fit_shouldered_plateau`."""
    z = np.asarray(z, dtype=float)
    fsl, z1, z2 = params["fsl"], params["z1"], params["z2"]
    d = z - z1
    low = fsl + sum(c * d ** (i + 1) for i, c in enumerate(params["coef"]))
    high = _upper(z, z2, fsl, params["upper_coef"], params["upper_join"])
    return np.where(z < z1, low, np.where(z <= z2, fsl, high))


# -- the logistic ---------------------------------------------------------------

def logistic(z, floor, k, z0, ceiling):
    return floor + (ceiling - floor) / (1.0 + np.exp(-k * (z - z0)))


def fit_logistic(z, level):
    """All four parameters fitted, the ceiling included.

    Least squares tends to put the ceiling just below the largest maximum,
    because the record carries no information above its own top. That is the
    honest reading of the data and the reason the logistic is the fallback here.
    """
    from scipy.optimize import curve_fit

    z, level = np.asarray(z, float), np.asarray(level, float)
    popt, _ = curve_fit(
        logistic, z, level,
        p0=[float(level.min()) - 1.0, 1.5, 0.0, float(level.max()) + 0.5],
        maxfev=400000)
    params = {"floor": float(popt[0]), "k": float(popt[1]), "z0": float(popt[2]),
              "ceiling": float(popt[3])}
    return params, float(np.sqrt(np.mean((evaluate(LOGISTIC, z, params) - level) ** 2)))


# -- either form ------------------------------------------------------------------

def fit(form, z, level, fsl=None, degree=SHOULDER_DEGREE, tol=PLATEAU_TOLERANCE,
        gap=PLATEAU_GAP, upper_degree=UPPER_DEGREE, upper_join=FREE):
    if form == SHOULDERED:
        if fsl is None:
            raise ValueError("the shouldered form needs the full supply level")
        return fit_shouldered_plateau(z, level, fsl, degree, tol, gap,
                                      upper_degree, upper_join)
    if form == LOGISTIC:
        return fit_logistic(z, level)
    raise ValueError(f"unknown curve form {form!r}")


def evaluate(form, z, params):
    if form == SHOULDERED:
        return shouldered_plateau(z, params)
    if form == LOGISTIC:
        return logistic(np.asarray(z, float), params["floor"], params["k"],
                        params["z0"], params["ceiling"])
    raise ValueError(f"unknown curve form {form!r}")


_SKIPPED = (RuntimeError, ValueError, np.linalg.LinAlgError)


def _variate(p):
    """z for exceedance probabilities, by scipy - see LakeLevelRecord.plotting_positions."""
    from scipy.special import ndtri
    return ndtri(1.0 - np.asarray(p, dtype=float))


def bootstrap_curves(level, z_eval, form, fsl=None, draws=BOOTSTRAP_DRAWS,
                     seed=BOOTSTRAP_SEED, **options):
    """Every resampled curve, one row per draw, evaluated at ``z_eval``.

    Resamples the chosen form cannot be fitted to are skipped - a resample may
    leave too few maxima on or above the plateau to place it - so the count can
    fall short of ``draws``, and the caller should say by how much.
    """
    rng = np.random.default_rng(seed)
    level = np.asarray(level, dtype=float)
    drawn = []
    for _ in range(draws):
        sample = np.sort(rng.choice(level, size=len(level), replace=True))
        _, zz = record.plotting_positions(sample, variate=_variate)
        zz = np.sort(zz)
        try:
            params, _ = fit(form, zz, sample, fsl, **options)
        except _SKIPPED:
            continue
        curve = evaluate(form, z_eval, params)
        if np.isfinite(curve).all():
            drawn.append(curve)
    if not drawn:
        raise RuntimeError("no resample of the record could be fitted")
    return np.vstack(drawn)


def bootstrap_censored_curves(level, carried, z_eval, form, fsl=None,
                              draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED,
                              minimum=MINIMUM_STORM_YEARS, **options):
    """The storm-driven reading, refitted to resamples of whole years."""
    rng = np.random.default_rng(seed)
    level = np.asarray(level, dtype=float)
    carried = np.asarray(carried, dtype=bool)
    n = len(level)
    drawn = []
    for _ in range(draws):
        pick = rng.integers(0, n, n)
        storm = level[pick][~carried[pick]]
        if len(storm) < minimum:
            continue
        zz, values = record.censored_positions(storm, n, variate=_variate)
        try:
            params, _ = fit(form, zz, values, fsl, **options)
        except _SKIPPED:
            continue
        curve = evaluate(form, z_eval, params)
        if np.isfinite(curve).all():
            drawn.append(curve)
    if not drawn:
        raise RuntimeError("no resample of the storm-driven maxima could be fitted")
    return np.vstack(drawn)


# -- design flood realisations ----------------------------------------------------

# Past this the per-duration curves are thinned for the results file. Every
# realisation is kept at the rare end, where there are few.
MAX_CURVE_POINTS = 1500


def read_design_curve(path):
    """(aep, level) from a Monte Carlo database, ascending in AEP.

    ``level_aep`` in an mcdf is a **probability**, not 1 in X - unlike
    ``rain_aep`` beside it (see lib/RepresentativeEvents.py).
    """
    frame = pd.read_parquet(path) if str(path).lower().endswith(".parquet") \
        else pd.read_csv(path, usecols=["level", "level_aep"])
    frame = frame[["level", "level_aep"]].dropna()
    frame = frame[frame["level_aep"] > 0].sort_values("level_aep", kind="stable")
    return frame["level_aep"].to_numpy(float), frame["level"].to_numpy(float)


def design_envelope(curves, aep):
    """The design flood level: the worst duration at each probability.

    Interpolated in log AEP, since the realisations span seven orders of
    magnitude. Beyond a duration's own range its curve is held at its end value.
    """
    stacked = [np.interp(np.log(aep), np.log(x), y) for x, y in curves.values()]
    return np.max(np.vstack(stacked), axis=0)


def thinned(x, y, limit=MAX_CURVE_POINTS):
    if len(x) <= limit:
        return x, y
    keep = np.unique(np.round(np.geomspace(1, len(x), limit)).astype(int) - 1)
    return x[keep], y[keep]


# -- the whole analysis -------------------------------------------------------------

GRID_POINTS = 400


def frequency_grid(frequent_aep, rare_aep, points=GRID_POINTS):
    aep = np.logspace(np.log10(frequent_aep), np.log10(rare_aep), points)
    return aep, _variate(aep)


def _band(draws):
    lo, hi = np.percentile(draws, BAND_PERCENTILES, axis=0)
    return lo, hi


def _fit_block(form, z, level, grid_z, fsl, options, draws_of, draws):
    """One curve, its rmse and its band, as plain lists, or the reason it failed."""
    block = {"form": form, "rmse": None, "params": None, "curve": None,
             "band_lo": None, "band_hi": None, "draws_used": 0, "draws": int(draws),
             "z_min": float(np.min(z)) if len(z) else None,
             "z_max": float(np.max(z)) if len(z) else None, "error": None,
             "warning": None}
    try:
        params, rmse = fit(form, z, level, fsl, **options)
    except _SKIPPED as exc:
        block["error"] = str(exc)
        return block
    block.update(params=params, rmse=rmse,
                 curve=evaluate(form, grid_z, params).tolist())
    try:
        draws_drawn = draws_of()
    except _SKIPPED as exc:
        block["error"] = f"band: {exc}"
        return block
    lo, hi = _band(draws_drawn)
    block.update(band_lo=lo.tolist(), band_hi=hi.tolist(), draws_used=int(len(draws_drawn)))
    if len(draws_drawn) < LOW_FIT_SHARE * block["draws"]:
        block["warning"] = (
            f"only {len(draws_drawn)} of {block['draws']} resamples could be fitted, so "
            f"the band comes from the resamples that kept what this form depends on and "
            f"is likely too narrow - try a wider plateau tolerance, no plateau, or "
            f"fewer terms")
    return block


def analyse(ams: pd.DataFrame, *, fsl=None, form=SHOULDERED,
            degree=SHOULDER_DEGREE, plateau_tolerance=PLATEAU_TOLERANCE,
            plateau_gap=PLATEAU_GAP, upper_degree=UPPER_DEGREE, upper_join=FREE,
            storm_driven=True, draws=BOOTSTRAP_DRAWS,
            seed=BOOTSTRAP_SEED, design_sources=(), frequent_aep=EY1_AEP,
            rare_aep=5e-4, progress=print) -> dict:
    """Fit, resample and read the design floods, into a JSON-ready dict.

    ``ams`` is ``LakeLevelRecord.with_positions`` output - the maxima that go
    into the curve, already filtered for coverage. ``design_sources`` is a
    sequence of ``(duration hours, mcdf path)``.
    """
    grid_aep, grid_z = frequency_grid(frequent_aep, rare_aep)
    z = ams["z"].to_numpy(float)
    level = ams["level"].to_numpy(float)
    carried = ams["carried_over"].to_numpy(bool)
    options = {} if form == LOGISTIC else {"degree": int(degree),
                                           "tol": float(plateau_tolerance),
                                           "gap": float(plateau_gap),
                                           "upper_degree": int(upper_degree),
                                           "upper_join": upper_join}
    out = {"grid": {"aep": grid_aep.tolist(), "z": grid_z.tolist()},
           "fits": {}, "design": None}

    if form != NO_FIT:
        progress(f"Fitting the {form} form to {len(level)} annual maxima and "
                 f"resampling {draws} times")
        out["fits"]["all"] = _fit_block(
            form, z, level, grid_z, fsl, options,
            lambda: bootstrap_curves(level, grid_z, form, fsl, draws, seed, **options),
            draws)
        if storm_driven and carried.any():
            storm = ~carried
            storm_z, storm_level = record.censored_positions(level[storm], len(level))
            progress(f"Fitting the {int(storm.sum())} storm-driven maxima")
            out["fits"]["storm"] = _fit_block(
                form, storm_z, storm_level, grid_z, fsl, options,
                lambda: bootstrap_censored_curves(level, carried, grid_z, form, fsl,
                                                  draws, seed, **options),
                draws)

    if design_sources:
        curves, drawn = {}, {}
        for hours, path in design_sources:
            progress(f"Reading the {float(hours):g} h design floods from {path}")
            x, y = read_design_curve(path)
            if not len(x):
                continue
            curves[float(hours)] = (x, y)
            tx, ty = thinned(x, y)
            drawn[f"{float(hours):g}"] = {
                "aep": tx.tolist(), "level": ty.tolist(),
                "z": _variate(tx).tolist()}
        if curves:
            out["design"] = {"durations": drawn,
                             "envelope": design_envelope(curves, grid_aep).tolist()}
    return out
