"""A miniature study for the Report page: one design group, one PMF run.

The design group is two durations whose curves are chosen so the answers are
unambiguous: the **level** is critical at 72 h for the frequent AEPs and at 36 h
for the rare ones, while the **inflow** is always largest at 36 h. A table that
reads the inflow at its own critical duration instead of the level's therefore
gets a different number at every frequent AEP, and the tests can tell.

Outflow is blank for both durations at 1 in 5 - the dam does not spill - which
is what drops the AEP from the outflow comparison altogether; the table must
still say 0 there.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from conftest import reservoir_row, write_workbook

GROUP = "RFSL_Design_GWL1p3"
PMF_GROUP = "PMF_GWL1p3"
AEPS = [2, 5, 10, 100, 1000, 10000, 100000, 1900000]

# level (m AHD): 72 h wins up to 1 in 100, 36 h from 1 in 1,000
LEVEL = {72: [214.0, 215.10, 216.00, 216.80, 217.10, 218.20, 219.80, 220.80],
         36: [213.5, 214.90, 215.80, 216.70, 217.30, 218.40, 220.00, 221.00]}
# inflow (m3/s): 36 h always larger, so 'own' and 'at level' disagree at 72 h
INFLOW = {72: [300, 580, 950, 2480, 4580, 6410, 9030, 15250],
          36: [400, 700, 1100, 2900, 5000, 7000, 10000, 16000]}
# outflow (m3/s): nothing at 1 in 2 or 1 in 5 for either duration
OUTFLOW = {72: [None, None, 150, 2250, 4360, 5480, 7710, 12900],
           36: [None, None, 120, 2200, 4400, 5500, 7800, 13100]}

COLUMNS = ["Include", "Group", "Output suffix", "Duration", "Method", "Output file",
           "Run models", "Analyse results", "Store hydrographs", "Input MCDF",
           "Inflow", "ELS file", "SQ file", "FSL", "Hydrographs folder",
           "Results folder", "Config file", "Log file"]


def _quantiles(path: Path, column: str, values) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"aep (1 in x)": AEPS,
                  "probability": [1.0 / aep for aep in AEPS],
                  column: values}).to_csv(path, index=False)


def _config(folder: Path, name: str, sims_list: str) -> Path:
    path = folder / name
    path.write_text(json.dumps({"simulation_list": sims_list, "filepaths": {}}),
                    encoding="utf-8")
    return path


def build_design_run(folder: Path) -> Path:
    """A reservoir-routing sims list of two durations and its quantile files."""
    rows = []
    for duration in (36, 72):
        name = f"CLD_mc_{duration}h_GWL1p3_RFSL"
        rows.append(reservoir_row(**{
            "Group": GROUP, "Duration": duration, "Output file": name,
            "Output suffix": "", "Results folder": r"sims_mc\results"}))
        results = folder / "sims_mc" / "results"
        _quantiles(results / f"{name}__level_quantiles.csv", "level", LEVEL[duration])
        _quantiles(results / f"{name}__inflow_quantiles.csv", "inflow", INFLOW[duration])
        _quantiles(results / f"{name}__outflow_quantiles.csv", "outflow", OUTFLOW[duration])
    write_workbook(folder / "CLD_RFSL_mc_sims_01.xlsx", COLUMNS, rows)
    return _config(folder, "CLD_RFSL_mc_sims_01.json", "CLD_RFSL_mc_sims_01.xlsx")


def build_pmf_run(folder: Path) -> Path:
    """An ensemble sims list of one row, with its database of events."""
    output = r"sims_enb\results\PMF_GWL1p3"
    row = reservoir_row(**{"Group": PMF_GROUP, "Method": "ensemble",
                           "Output file": output, "Duration": None})
    write_workbook(folder / "PMF_sims.xlsx", COLUMNS, [row])
    database = folder / "sims_enb" / "results" / "PMF_GWL1p3.csv"
    database.parent.mkdir(parents=True, exist_ok=True)
    # the largest level is event 2 (9 h); the largest inflow is event 0 (4.5 h)
    pd.DataFrame({"duration": [4.5, 6, 9, 12],
                  "tp": [0, 0, 0, 0], "storm_method": ["GSDM"] * 4,
                  "inflow": [18000.0, 17000.0, 17500.0, 16000.0],
                  "level": [221.20, 221.30, 221.38, 221.10],
                  "outflow": [15000.0, 15800.0, 16310.0, 14000.0]}).to_csv(database)
    return _config(folder, "PMF_sims.json", "PMF_sims.xlsx")


def build_study(folder: Path):
    """A study file with both runs; returns the open Study."""
    from core import study as studies

    runs = folder / "runs" / "E099"
    runs.mkdir(parents=True, exist_ok=True)
    design = build_design_run(runs)
    add_realisations(runs)
    pmf = build_pmf_run(runs)
    study = studies.new_study(folder / "bryan_study.json", "Miniature study")
    study.add_run("E099 RFSL", design)
    study.add_run("E099 PMF", pmf)
    study.save()
    return study


# -- Monte Carlo realisations, for the notional AEP of the PMF -----------------
# Level linear in the standard normal variate: 214 + 1.4 z. The PMF level of the
# ensemble fixture, 221.38 m, then sits at z = 5.271 - 1 in 14.7 million - inside
# a sample that reaches z = 5.8 (1 in 300 million).
MC_BASE, MC_SLOPE, MC_TOP_Z = 214.0, 1.4, 5.8


def write_mcdf(path: Path, count: int = 5000) -> Path:
    import numpy as np
    from statistics import NormalDist

    z = np.linspace(-1.0, MC_TOP_Z, count)
    level = MC_BASE + MC_SLOPE * z
    probability = [NormalDist().cdf(-value) for value in z]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"m": range(count), "n": 0, "level": level,
                  "level_aep": probability, "inflow": level * 10.0,
                  "inflow_aep": probability}).to_csv(path)
    return path


def add_realisations(folder: Path) -> None:
    """An mcdf beside each design row's quantile files, as reservoir routing writes."""
    for duration in (36, 72):
        write_mcdf(folder / "sims_mc" / "results" / f"CLD_mc_{duration}h_GWL1p3_RFSL__mcdf.csv")
