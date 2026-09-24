"""Antecedent storage: the lake volume at the start of the flood-producing storm.

For each water year, the storm behind the annual maximum lake level is found in
the catchment-average daily rainfall (the rarest burst of one to five days in the
thirty before the peak, against the IFD), and the homogenised lake volume is read
on the day it started - and, walking back through the wet days before it, on the
day the pre-burst started. A logistic S-curve in the standard normal variate is
fitted to each, bounded by a floor from the data and the full supply volume, and
written as the ``lake_config.json`` Bryan's Monte Carlo scheme samples the
antecedent storage from.

Written for Callide Dam in callide-fsl-reinstate; the method modules here are
its ``antecedent_storage`` modules, copied, with one change each so the settings
bind at call time:

    antecedent   the peak-conditioned burst search and the antecedent series
    scurve       the S-curve fit, and the model's sigmoid parameterisation
    correlation  storm severity against antecedent storage
    job          the job file, the run, the outputs and the lake configs

Needs scipy, so it runs under Bryan's interpreter (``util/AntecedentStorage.py``).
"""
