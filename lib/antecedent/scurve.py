"""Logistic (sigmoid) frequency curve for antecedent storage.

Several Sunwater dam studies fit the antecedent-storage distribution with a
sigmoid -- a logistic function in log-volume against the standard normal variate
-- and Callide fits one too (design hydrology report, Figure 27). This is the
curve a Monte Carlo design-flood assessment samples antecedent storage from.

The fit is done in two steps, following the earlier workbook:

1. **Decide the floor and ceiling.**  A logistic has two asymptotes, and left to
   a free optimiser they wander -- the ceiling in particular drifts above the
   largest sampled storage and the curve then extrapolates higher still.  So the
   floor and ceiling are fixed first, as a decision about the plausible range of
   antecedent storage, bracketing the data.  By default the floor is the
   observed minimum rounded down and the ceiling the observed maximum rounded up
   (to ``ROUND_ML``), so the curve asymptotes just outside the data and cannot
   exceed it; either can be overridden (e.g. the ceiling set to full supply
   volume).

2. **Optimise position and slope.**  With the asymptotes fixed, only the centre
   ``z0`` (position) and the rate ``k`` (slope) are fitted, by least squares in
   log10(volume) -- the "log-normal space" the distribution is described in.

The curve is the standard logistic between the two asymptotes:

    log10 V(z) = log10(floor) + (log10(ceiling) - log10(floor)) / (1 + exp(-k (z - z0)))

so ``V`` runs from ``floor`` at low z to ``ceiling`` at high z and never beyond.

The earlier workbook wrote the same family with a fifth parameter, a height
adjustment factor ``H``; the form here is its ``H = 1`` case.  Fits carry ``H``
explicitly so they are complete in either form -- see :func:`asym_parameters`.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

CUNNANE_A = 0.4
ROUND_ML = 1000.0        # floor/ceiling are rounded to this when chosen from data


def plotting_position(volumes, weights=None):
    """Cunnane plotting position and its standard normal variate.

    Largest volume ranked 1 (rarest, highest z).  Returns (p, z) aligned to the
    input order.

    ``weights`` generalises the rank to a cumulative weight, so a sample that is
    not one-per-year can be made to carry one year's worth of weight per year:
    with unit weights the cumulative weight at rank r is r and the expression is
    the plain Cunnane formula, unchanged.  This is only needed for the
    peaks-over-threshold conditioning, where a wet year contributes four storms
    and a dry one none -- see ``selection.year_weights``.
    """
    volumes = np.asarray(volumes, dtype=float)
    n = len(volumes)
    weights = np.ones(n) if weights is None else np.asarray(weights, dtype=float)
    order = volumes.argsort()[::-1]           # descending
    cumulative = np.empty(n)
    cumulative[order] = np.cumsum(weights[order])
    p = (cumulative - CUNNANE_A * weights) / (weights.sum() + 1 - 2 * CUNNANE_A)
    z = norm.ppf(1 - p)
    return p, z


def _log10_model(z, lf, lc, k, z0):
    return lf + (lc - lf) / (1.0 + np.exp(-k * (z - z0)))


def logistic_volume(z, floor, ceiling, k, z0):
    """Antecedent storage (ML) at standard normal variate ``z``.

    Standard logistic between ``floor`` and ``ceiling`` (ML), so the value is
    bounded by the two asymptotes and never extrapolates beyond the ceiling.
    """
    return np.power(10.0, _log10_model(z, np.log10(floor), np.log10(ceiling), k, z0))


def choose_bounds(volumes, round_ml=None):
    """A reasonable floor and ceiling: observed min rounded down, max rounded up.

    The decision is deliberately simple and reproducible -- the asymptotes sit
    just outside the data, so the fitted curve spans the samples without
    exceeding them.
    """
    # Resolved when called, so the job's setting takes effect (Bryan's change).
    round_ml = ROUND_ML if round_ml is None else round_ml
    volumes = np.asarray(volumes, dtype=float)
    volumes = volumes[np.isfinite(volumes)]
    floor = np.floor(volumes.min() / round_ml) * round_ml
    ceiling = np.ceil(volumes.max() / round_ml) * round_ml
    return float(floor), float(ceiling)


def fit(volumes, floor=None, ceiling=None, weights=None):
    """Fit the logistic S-curve with fixed asymptotes.

    ``floor`` and ``ceiling`` (ML) are the decision; if omitted they are taken
    from the data by :func:`choose_bounds`.  Only the position ``z0`` and slope
    ``k`` are optimised, in log10 space.  Returns the parameters, the (fixed)
    floor and ceiling, the storage at z = 0 and the fit rmse in ML.

    ``weights`` are carried through both the plotting positions and the least
    squares, for the year-weighted peaks-over-threshold variant; omitted, every
    sample counts once and the fit is unchanged.

    ``H`` is carried at 1.0 so the fit is complete in the earlier workbook's
    five-parameter form as well; see :func:`asym_parameters`.
    """
    volumes = np.asarray(volumes, dtype=float)
    keep = np.isfinite(volumes)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)[keep]
    volumes = volumes[keep]
    if floor is None or ceiling is None:
        f, c = choose_bounds(volumes)
        floor = f if floor is None else floor
        ceiling = c if ceiling is None else ceiling

    _, z = plotting_position(volumes, weights)
    log_v = np.log10(volumes)
    lf, lc = np.log10(floor), np.log10(ceiling)

    def model(z, k, z0):
        return _log10_model(z, lf, lc, k, z0)

    # curve_fit divides the residuals by sigma, so a weight of w is a sigma of
    # 1/sqrt(w); absolute_sigma is irrelevant here because only the parameters
    # are used, not their covariance.
    sigma = None if weights is None else 1.0 / np.sqrt(weights)
    (k, z0), _ = curve_fit(model, z, log_v, p0=[1.5, 0.0], sigma=sigma,
                           maxfev=200000)
    residual = volumes - logistic_volume(z, floor, ceiling, k, float(z0))
    if weights is None:
        rmse = float(np.sqrt(np.mean(residual ** 2)))
    else:
        rmse = float(np.sqrt(np.sum(weights * residual ** 2) / weights.sum()))
    return {
        "kind": "logistic", "floor_ML": floor, "ceiling_ML": ceiling,
        "H": 1.0, "k": float(k), "z0": float(z0),
        "mean_z0_ML": float(logistic_volume(0.0, floor, ceiling, k, float(z0))),
        "rmse_ML": rmse,
        "n": int(len(volumes)),
        "effective_n": float(len(volumes) if weights is None else weights.sum()),
    }


# The curve adopted in the earlier assessment (CLD_ADV_Assessment_rev3.xlsm,
# sheet "S-Curve (REV3_burst)"), read straight from the workbook's solved cells.
# Its parameterisation was log10 V = A/(H+exp(-k(z-z0))) + C, hand-tuned with the
# floor and ceiling pinned to round values; carried verbatim for comparison
# rather than re-fitted.  mean_z0_ML = 96,250 ML, the report's 75% of FSV.  Its H
# is 5.146, not 1: the workbook's "height adjustment factor", turned by hand until
# the curve passed through a chosen midpoint volume (see :func:`asym_parameters`).
#
# CAUTION: this curve does not describe antecedent storage.  The sheet it was
# fitted on holds ADV_Selection_S column B -- the ANNUAL PEAK storage, exactly,
# for all 25 years -- not column L, the volume on the day the main burst began.
# Its floor (14,000 ML) is the smallest annual peak, not the smallest antecedent
# volume; the same workbook's own burst volumes floor at 12,000 ML and fit at
# 66,742 ML, within 1% of this model.  ``extract_antecedent`` re-establishes the
# mis-paste at build time (``old_scurve_sheet_provenance``); the constants are
# kept only so the published curve can still be drawn and attributed.
OLD_PUBLISHED = {
    "kind": "old_asym",
    "A": 5.146128035678238, "H": 5.146298595784089,
    "C": 4.146128035678238, "k": 2.3809396567863494, "z0": 0.0,
    "floor_ML": 14000.0, "ceiling_ML": 140000.0, "mean_z0_ML": 96250.0,
    "n": 25, "label": "old assessment (hand-tuned)",
}

# The old assessment's pre-burst curve (sheet "S-Curve (REV3_preburst)"): the
# antecedent storage at the start of the pre-burst rainfall, excluding the rise
# the pre-burst caused.  Unlike the burst sheet above, this one does hold the
# antecedent volumes it claims to (ADV_Selection_S column J, exactly).  It is
# still hand-tuned well above its own data, which refits at 64,319 ML.
OLD_PUBLISHED_PREBURST = {
    "kind": "old_asym",
    "A": 5.107209969647869, "H": 4.5526973330837075,
    "C": 4.0, "k": 1.7194628363249493, "z0": 0.38958583060657126,
    "floor_ML": 10000.0, "ceiling_ML": 132373.0, "mean_z0_ML": 83133.0,
    "n": 25, "label": "old assessment pre-burst (hand-tuned)",
}


def asym_parameters(fit_dict):
    """The same curve in the earlier workbook's five-parameter form.

    ``CLD_ADV_Assessment_rev3.xlsm`` wrote the S-curve as

        log10 V(z) = A / (H + exp(-k (z - z0))) + C

    over five cells: ``k``, ``z0``, ``C = log10(Vf)`` the floor, ``A``, and ``H``,
    which its legend calls the "height adjustment factor" and its author turned by
    hand ("adjust to get height correct") until the curve passed through a chosen
    midpoint volume ``Vm``.  With H free the ceiling is only approximately ``Vc``
    and the value at ``z0`` is wherever the tuning put it.

    The standard logistic fitted here is the ``H = 1`` case of that same family,
    with ``A = log10(Vc) - log10(Vf)``: the ceiling is then exactly ``Vc`` and the
    value at ``z0`` exactly halfway between the asymptotes in log space.  So a fit
    from this module drops straight into a sheet built on the old formula, with
    H set to 1 -- which is why H is carried explicitly rather than left implied.

    Returns ``{"A": ..., "H": ..., "C": ...}``; ``k`` and ``z0`` are unchanged.
    """
    if fit_dict.get("kind") == "old_asym":
        return {k: float(fit_dict[k]) for k in ("A", "H", "C")}
    lf = float(np.log10(fit_dict["floor_ML"]))
    return {"A": float(np.log10(fit_dict["ceiling_ML"])) - lf,
            "H": float(fit_dict.get("H", 1.0)), "C": lf}


def evaluate(fit_dict, z):
    """Antecedent storage (ML) at ``z`` for any fit dict (new or old published)."""
    if fit_dict.get("kind") == "old_asym":
        g = fit_dict["A"] / (fit_dict["H"] + np.exp(
            -fit_dict["k"] * (z - fit_dict["z0"]))) + fit_dict["C"]
        return np.power(10.0, g)
    return logistic_volume(z, fit_dict["floor_ML"], fit_dict["ceiling_ML"],
                           fit_dict["k"], fit_dict["z0"])


def curve(fit_dict, z=None):
    """Sample the fitted curve; returns (z, volume) for plotting."""
    if z is None:
        z = np.linspace(-2.6, 2.6, 240)
    return z, evaluate(fit_dict, z)
