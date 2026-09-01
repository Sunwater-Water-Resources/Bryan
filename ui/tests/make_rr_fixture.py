"""Build a miniature reservoir-routing model that real Bryan can run.

Reservoir routing is the only method that needs no model executable - it takes
stored inflow hydrographs and re-routes them through a rating curve - so this
is what lets one test drive Bryan itself end to end, on Linux, in seconds.

The formats come from lib/ReservoirRouting.py:

    .els   csv with EL and V columns                      (_read_els:39)
    .sq    five header lines, then 'storage_above_FSV  flow' pairs   (_read_sq:48)
    mcdf   csv indexed by row; ensemble input is detected from its
           'duration'/'tp' columns, Monte Carlo from 'm'/'n'         (_detect_scheme)
    inflow csv indexed by time in hours, one column per simulation

Called by conftest, and runnable directly to refresh ui/tests/data_rr.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

FSL = 100.0
N_STEPS = 60
DT_HOURS = 1.0


def write_els(path: Path, fsl: float = FSL) -> None:
    """Elevation-storage, linear enough to be obvious, wide enough not to clip."""
    lines = ["EL,V"]
    for step in range(21):
        elevation = fsl - 5.0 + step          # 95 m to 115 m
        volume = 1000.0 * step + 500.0        # 500 ML to 20500 ML
        lines.append(f"{elevation:.2f},{volume:.1f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_sq(path: Path, coefficient: float = 12.0) -> None:
    """Storage above FSV against outflow. Five header lines, then pairs.

    ``coefficient`` sizes the spillway, so two curves can differ the way the
    'base' and 'raised' options of a real study do: a smaller one passes less
    flow and holds the lake higher for the same inflow. They were identical
    until 26 August 2026, which let a routed result match under both.
    """
    header = [
        "Miniature test rating curve",
        "written by ui/tests/make_rr_fixture.py",
        "storage is ML above full supply, flow is m3/s",
        "not a real dam",
        "-" * 40,
    ]
    rows = []
    for step in range(11):
        storage = 500.0 * step                # 0 to 5000 ML above FSV
        flow = coefficient * step ** 1.5      # a plausible spillway shape
        rows.append(f"{storage:.1f}  {flow:.3f}")
    path.write_text("\n".join(header + rows) + "\n", encoding="utf-8")


def write_ensemble_inputs(folder: Path, durations=(18, 24), patterns=(0, 1)) -> None:
    """An ensemble results database and the inflow hydrographs behind it.

    Ensemble input is the simpler of the two: one ADV for the whole run, and no
    Monte Carlo scheme parameters to check against.
    """
    folder.mkdir(parents=True, exist_ok=True)

    columns, records = [], []
    for duration in durations:
        for pattern in patterns:
            name = f"d{duration}_tp{pattern}"
            columns.append(name)
            records.append({
                "": name,
                "rain_z": 3.0,
                "rain_aep": 100.0,
                "mean_rain_mm": 200.0 + 10 * pattern,
                "duration": duration,
                "tp": pattern,
                "storm_method": "ARR point",
                "tp_frequency": 0.1,
                "preburst_p": 0.5,
                "preburst_proportion": 0.1,
                "preburst_mm": 20.0,
                "initial_loss": 30.0,
                "continuing_loss": 2.5,
                "residual_depth": 0.0,
                "embedded_bursts": 0,
                "ADV": 10500.0,
                "inflow": 0.0,
                "level": 0.0,
                "outflow": 0.0,
            })

    header = list(records[0])
    lines = [",".join(header)]
    for record in records:
        lines.append(",".join(str(record[key]) for key in header))
    (folder / "ensemble_results.csv").write_text("\n".join(lines) + "\n",
                                                 encoding="utf-8")

    # Ragged on purpose: each duration is padded with blanks after its own
    # simulation period, which is how a real ensemble hydrograph file looks.
    inflow_lines = ["time," + ",".join(columns)]
    for step in range(N_STEPS):
        hour = step * DT_HOURS
        cells = []
        for name in columns:
            duration = int(name.split("_")[0][1:])
            span = duration + 12
            if hour > span:
                cells.append("")               # past this duration's period
                continue
            peak = 300.0 if duration == 24 else 220.0
            shape = math.exp(-((hour - span / 2.0) ** 2) / (2 * (span / 6.0) ** 2))
            cells.append(f"{peak * shape:.4f}")
        inflow_lines.append(f"{hour:.2f}," + ",".join(cells))
    (folder / "ensemble_inflows.csv").write_text("\n".join(inflow_lines) + "\n",
                                                 encoding="utf-8")


def write_monte_carlo_inputs(folder: Path, m_count=4, n_count=5,
                             lower_aep=2.0, upper_aep=1000.0) -> dict:
    """A Monte Carlo results database and the inflows behind it.

    The sample has to survive _validate_sample, which is the point of building
    it properly: every realisation's rain_z sits inside its own main division,
    rain_aep agrees with rain_z, and every (m, n) position appears once. The
    scheme_config returned here is what the sims list's Config file must hold.
    """
    from math import erf, exp, sqrt

    def ndtri(p):
        """The standard normal variate of p, by bisection - no scipy here."""
        def ndtr(z):
            return 0.5 * (1.0 + erf(z / sqrt(2.0)))
        low, high = -10.0, 10.0
        for _ in range(200):
            middle = (low + high) / 2.0
            if ndtr(middle) < p:
                low = middle
            else:
                high = middle
        return (low + high) / 2.0

    def ndtr(z):
        return 0.5 * (1.0 + erf(z / sqrt(2.0)))

    folder.mkdir(parents=True, exist_ok=True)

    z_low, z_up = ndtri(1.0 - 1.0 / lower_aep), ndtri(1.0 - 1.0 / upper_aep)
    edges = [z_low + (z_up - z_low) * i / m_count for i in range(m_count + 1)]

    columns, records = [], []
    for m in range(m_count):
        for n in range(n_count):
            # Spread inside the division, never on an edge - a z on a boundary
            # is legal but makes a failure ambiguous to read.
            span = edges[m + 1] - edges[m]
            z = edges[m] + span * (n + 1) / (n_count + 1)
            name = f"m{m}_n{n}"
            columns.append(name)
            records.append({
                "": name,
                "m": m,
                "n": n,
                "rain_z": z,
                "rain_aep": 1.0 / (1.0 - ndtr(z)),
                "storm_method": "ARR point",
                "initial_loss": 30.0,
                "continuing_loss": 2.5,
                # Full supply exactly - see write_els. Starting above the top
                # of the .sq instead pegs every realisation at the last point
                # of the rating curve, and both curves then agree.
                "ADV": 5500.0,
                "lake_z": 0.0,
                "inflow": 0.0,
                "level": 0.0,
                "outflow": 0.0,
            })

    header = list(records[0])
    lines = [",".join(header)]
    for record in records:
        lines.append(",".join(str(record[key]) for key in header))
    (folder / "mc_results.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # One hydrograph per realisation, growing with the sampled rainfall so the
    # routed peaks come out in a sensible order.
    span = 30.0
    inflow_lines = ["time," + ",".join(columns)]
    for step in range(N_STEPS):
        hour = step * DT_HOURS
        cells = []
        for record in records:
            # Sized so the flood volume lands between about 400 and 4000 ML,
            # which is inside the 0 to 5000 ML the .sq covers.
            peak = 5.0 + 25.0 * record["rain_z"]
            shape = exp(-((hour - span / 2.0) ** 2) / (2 * (span / 6.0) ** 2))
            cells.append(f"{max(peak, 2.0) * shape:.4f}")
        inflow_lines.append(f"{hour:.2f}," + ",".join(cells))
    (folder / "mc_inflows.csv").write_text("\n".join(inflow_lines) + "\n",
                                           encoding="utf-8")

    return {
        "number_of_main_divisions": m_count,
        "number_of_sub_divisions": n_count,
        "lower_aep": lower_aep,
        "upper_aep": upper_aep,
        "aep_of_pmp": None,
    }


SIMS_COLUMNS = [
    "Include", "Method", "Output suffix", "Duration", "Run models",
    "Analyse results", "Store hydrographs", "Input database", "Inflow",
    "ELS file", "SQ file", "FSL", "ADV", "ADV source", "Hydrographs folder",
    "Results folder", "Config file", "Log file", "Output file", "Comment",
]


def sims_rows(suffixes=("base", "raised")):
    """One row per rating-curve option - the thing the method exists for."""
    return [
        {
            "Include": "yes",
            "Method": "reservoir routing",
            "Output suffix": suffix,
            "Duration": "",
            "Run models": "yes",
            "Analyse results": "yes",
            "Store hydrographs": "yes",
            "Input database": "inputs/ensemble_results.csv",
            "Inflow": "inputs/ensemble_inflows.csv",
            "ELS file": "reservoir/dam.els",
            "SQ file": f"reservoir/dam_{suffix}.sq",
            "FSL": FSL,
            # 'fsv' resolves against the curve BEING routed, which is the whole
            # point when the dam differs from the source run's.
            "ADV": "fsv",
            "ADV source": "",
            "Hydrographs folder": "results/hydrographs",
            "Results folder": "results",
            "Config file": "",
            "Log file": "",
            "Output file": "mini_enb",
            "Comment": f"miniature fixture, {suffix} rating",
        }
        for suffix in suffixes
    ]


def monte_carlo_rows(suffixes=("base", "raised"), output_file="mini_mc",
                     adv="", adv_source=""):
    """The case the suffix exists for: one Output file, several rating curves.

    Both rows write a Monte Carlo database into one results folder, so they are
    told apart only by 'Output suffix' - which is why lib/ReservoirRouting.py
    puts it on the mcdf name too.

    ``adv`` and ``adv_source`` are left blank by default, which is the ordinary
    re-route: the antecedent storage of every realisation comes from the input
    database. Setting both is how one fixed storage is tested instead.
    """
    return [
        {
            "Include": "yes",
            "Method": "reservoir routing",
            "Output suffix": suffix,
            "Duration": "",
            "Run models": "yes",
            "Analyse results": "yes",
            "Store hydrographs": "no",
            "Input database": "inputs/mc_results.csv",
            "Inflow": "inputs/mc_inflows.csv",
            "ELS file": "reservoir/dam.els",
            "SQ file": f"reservoir/dam_{suffix}.sq",
            "FSL": FSL,
            "ADV": adv,
            "ADV source": adv_source,
            "Hydrographs folder": "results/hydrographs",
            "Results folder": "results",
            "Config file": "mc_config.json",
            "Log file": "",
            "Output file": output_file,
            "Comment": f"miniature monte carlo fixture, {suffix} rating",
        }
        for suffix in suffixes
    ]


def build(folder: Path) -> Path:
    """Write the whole project. Returns the sims_config.json path."""
    folder = Path(folder)
    (folder / "reservoir").mkdir(parents=True, exist_ok=True)
    (folder / "inputs").mkdir(parents=True, exist_ok=True)

    write_els(folder / "reservoir" / "dam.els")
    write_sq(folder / "reservoir" / "dam_base.sq")
    write_sq(folder / "reservoir" / "dam_raised.sq", coefficient=6.0)
    write_ensemble_inputs(folder / "inputs")
    scheme = write_monte_carlo_inputs(folder / "inputs")
    (folder / "mc_config.json").write_text(
        json.dumps({"scheme_config": scheme}, indent=2), encoding="utf-8")

    for name in ("model_config.json", "storm_config.json", "climate_config.json"):
        (folder / name).write_text("{}", encoding="utf-8")

    config = folder / "mini_sims_config.json"
    config.write_text(json.dumps({
        "simulation_list": "MiniSimsList.xlsx",
        "filepaths": {
            "model_config": "model_config.json",
            "storm_config": "storm_config.json",
            "climate_config": "climate_config.json",
        },
        "test_runs": 0,
    }, indent=2), encoding="utf-8")
    return config


if __name__ == "__main__":
    target = Path(__file__).parent / "data_rr"
    print("wrote", build(target))
