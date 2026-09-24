"""Correlation between burst rainfall magnitude and antecedent storage.

The design-flood Monte Carlo samples antecedent storage independently of storm
magnitude.  That is only defensible if the two are uncorrelated for the storm
sizes of interest.  The report's finding was a weak correlation overall that
vanishes once frequent storms are excluded; this reproduces that check on the
new record.

Two views, following the report:

* raw -- burst depth (mm) against antecedent storage (ML);
* standard-normal -- the storm's severity as a standard normal variate from its
  AEP, against the storage's plotting-position variate.

Storms more frequent than 1 EY (1 exceedance/year, AEP 1 in 1.582, severity
z < -0.34) are then excluded, leaving the magnitudes that matter for dam safety.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import linregress, norm

from .scurve import plotting_position

# 1 EY = one exceedance per year = AEP 1 in 1.582; below this a storm is "frequent".
ONE_EY_1INX = 1.582
ONE_EY_Z = float(norm.ppf(1 - 1 / ONE_EY_1INX))       # about -0.337


def storm_severity_z(aep_1inx):
    """Standard normal variate of a storm's severity; larger = rarer."""
    aep_1inx = np.asarray(aep_1inx, dtype=float)
    return norm.ppf(1 - 1 / aep_1inx)


def _score(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    fit = linregress(x, y)
    return {"n": int(len(x)), "r": float(fit.rvalue), "p": float(fit.pvalue),
            "slope": float(fit.slope), "intercept": float(fit.intercept),
            "strength": _strength(abs(fit.rvalue))}


def _strength(a):
    if a > 0.7:
        return "strong"
    if a > 0.49:
        return "high"
    if a > 0.3:
        return "moderate"
    if a > 0.1:
        return "low"
    return "negligible"


def analyse(table):
    """Correlation stats on the qualified antecedent-storage rows of ``table``.

    Returns raw and standard-normal scores over all retained storms and over the
    subset rarer than 1 EY, plus the arrays for plotting.
    """
    kept = table[table["qualified"]].copy()
    depth = kept["burst_depth_mm"].to_numpy(dtype=float)
    volume = kept["adv_burst_vol_ML"].to_numpy(dtype=float)
    severity = storm_severity_z(kept["burst_aep_1inx"].to_numpy(dtype=float))
    _, storage_z = plotting_position(volume)

    rare = severity >= ONE_EY_Z

    return {
        "raw_all": _score(depth, volume),
        "z_all": _score(severity, storage_z),
        "z_rare": _score(severity[rare], storage_z[rare]),
        "one_ey_z": ONE_EY_Z,
        "n_rare": int(rare.sum()),
        "n_frequent": int((~rare).sum()),
        "arrays": {
            "depth": depth, "volume": volume,
            "severity_z": severity, "storage_z": storage_z, "rare": rare,
        },
    }
