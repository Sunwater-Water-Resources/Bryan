"""Put a dam's recorded lake level record onto one spillway configuration.

The recorded level reflects whichever spillway was in place at the time - a
dam whose full supply level or gate operation changed during the record holds
several populations in one series. This homogenises it in two stages: derive the
net inflow by closing the water balance backwards against the rating that was in
force at each step, then re-route that inflow through one target rating. The
annual maxima of the result are one population, and the lake volumes are what
the antecedent storage analysis reads.

Written for Callide Dam in callide-fsl-reinstate (the general modules here are
its ``callide/`` package, copied); made general by ``job.py``, which takes every
input and setting from a job file instead of a repository layout:

    curves       storage table (.els), rating register (xlsx), target ratings
                 (.sq storage-discharge, or a level,flow .csv)
    gauges       WMIP / Hydstra exports, spliced in operating order, and an
                 optional overlay gauge below a level (Callide's intake gauge)
    evaporation  SILO Data Drill pan evaporation, monthly pan factors
    model        the two-stage routing (``LakeModel``)
    peaks        annual maxima by water year, independent peaks
    job          the job file, the run, and the outputs

Needs scipy (PCHIP storage curves), so it runs under Bryan's interpreter -
``util/HomogeniseLakeLevels.py`` - never in the launcher's process.
"""
