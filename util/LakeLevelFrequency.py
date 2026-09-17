"""Lake level frequency: fit the annual maxima, and draw the report figure.

The run launcher's Lake levels page shells out to this, because the curve fits
need scipy and the figure needs matplotlib, and the launcher's environment has
neither. It can equally be run by hand against a job file:

    python util/LakeLevelFrequency.py --job lake_job.json --results lake_results.json
    python util/LakeLevelFrequency.py --job lake_job.json --results lake_results.json \\
        --png figures/CLD_lake_levels.png [--without-design]
    python util/LakeLevelFrequency.py --job lake_job.json --ams-csv CLD_lake_level_ams.csv

The job is JSON - see ``LakeLevelRecord.JOB_DEFAULTS`` for every key - naming
the level record (one or more Hydstra or WMIP exports, in gauge order) or an
annual maximum CSV, the water year, full supply, the curve form and the Monte
Carlo databases to compare against. Paths in it are absolute.

``--results`` is reused when its fingerprint matches the job and the files it
reads, so exporting a figure after the page has fitted the curves does not
resample the record a second time. The resampling takes of the order of half a
minute for the shouldered form.

The figure is formatted for an A4 report page: 6.3 x 4.0 in at 300 dpi, base
font 8 pt, black text and a full box, no title (the caption carries it), and a
frequency axis reading in EY at the frequent end and "1 in X" past 50% AEP.
``--without-design`` leaves the Monte Carlo durations and the envelope off, for
the part of a report that presents the record before the modelling.
"""

import argparse
import json
import math
import os
import sys

BRYAN_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BRYAN_ROOT not in sys.path:
    sys.path.insert(0, BRYAN_ROOT)

import numpy as np  # noqa: E402

from lib import LakeLevelFrequency as frequency  # noqa: E402
from lib import LakeLevelRecord as record  # noqa: E402

RESULTS_VERSION = 1


# -- the analysis ----------------------------------------------------------------

def _clean(value):
    """JSON cannot hold NaN; the page reads None as a gap."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def ams_rows(positions):
    rows = []
    for row in positions.to_dict("records"):
        row["level_at"] = None if row["level_at"] is None or row["level_at"] != row["level_at"] \
            else str(row["level_at"])
        rows.append({key: _clean(value.item() if hasattr(value, "item") else value)
                     for key, value in row.items()})
    return rows


def run_analysis(job, progress=print):
    job = record.complete_job(job)
    ams, level_record = record.ams_for_job(job)
    positions = record.with_positions(ams, bool(job["include_incomplete"]))
    fit = job["fit"]
    design = job["design"]
    sources = [(source["duration"], source["path"])
               for source in (design["sources"] or [])] if design["include"] else []
    rare = 1.0 / float(job["axes"]["rare_aep_1_in_x"] or 2000)
    analysis = frequency.analyse(
        positions, fsl=job["fsl"], form=fit["form"], degree=fit["degree"],
        plateau_tolerance=fit["plateau_tolerance"], plateau_gap=fit["plateau_gap"],
        storm_driven=fit["storm_driven"], draws=int(fit["draws"]),
        seed=int(fit["seed"]), design_sources=sources,
        rare_aep=min(rare, float(positions["aep"].min())), progress=progress)
    analysis.update({
        "version": RESULTS_VERSION,
        "fingerprint": record.fingerprint(job),
        "job": job,
        "ams": ams_rows(positions),
        "notes": list(level_record.notes) if level_record is not None else [],
        "site": level_record.site_label if level_record is not None else "",
        "metadata": record.ams_metadata(job, level_record),
    })
    return analysis


def load_or_run(job, results_path, progress=print):
    """The saved results when they are for this job and these inputs."""
    if results_path and os.path.isfile(results_path):
        try:
            with open(results_path, encoding="utf-8") as stream:
                saved = json.load(stream)
            if (saved.get("version") == RESULTS_VERSION
                    and saved.get("fingerprint") == record.fingerprint(job)):
                progress(f"Reusing {results_path}")
                return saved
        except (OSError, ValueError):
            pass
    results = run_analysis(job, progress)
    if results_path:
        folder = os.path.dirname(os.path.abspath(results_path))
        os.makedirs(folder, exist_ok=True)
        temporary = results_path + ".tmp"
        with open(temporary, "w", encoding="utf-8") as stream:
            json.dump(results, stream)
        os.replace(temporary, results_path)
        progress(f"Wrote {results_path}")
    return results


# -- the report figure -------------------------------------------------------------

FIGSIZE = (6.3, 4.0)
DPI = 300
BASE_FONTSIZE = 8

INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d8d7d2"
RECORD_COLOUR = "#2a78d6"
DESIGN_COLOUR = "#eb6834"

# More frequent than 50% AEP the axis reads in EY, per ARR 2019.
EY_CROSSOVER_AEP = 0.5
EY_TICKS = (4, 3, 2, 1)
AEP_TICKS = (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 2e-3, 5e-4, 1e-4, 1e-5, 1e-6)

# The AEPs the axis may end at when the design floods are left off and the
# record alone sets the reach.
RARE_LIMITS = (0.1, 0.05, 0.02, 0.01, 5e-3, 2e-3, 1e-3, 5e-4)


def _z(aep):
    from scipy.special import ndtri
    return float(ndtri(1.0 - aep))


def frequency_ticks(frequent_aep, rare_aep):
    """Tick positions and labels: EY inside 50% AEP, "1 in X" past it."""
    marks = [(1.0 - math.exp(-ey), True) for ey in EY_TICKS]
    marks += [(aep, aep > EY_CROSSOVER_AEP) for aep in AEP_TICKS]
    positions, labels = [], []
    for aep, as_ey in sorted(marks, reverse=True):
        if not rare_aep * (1 - 1e-9) <= aep <= frequent_aep * (1 + 1e-9):
            continue
        positions.append(_z(aep))
        labels.append(f"{-math.log(1.0 - aep):.3g}EY" if as_ey
                      else f"1 in {1.0 / aep:,.0f}")
    return positions, labels


def level_text(value):
    """217.23 -> '217.23', 215.5 -> '215.5': as many places as the level carries."""
    text = f"{value:.3f}".rstrip("0")
    return text + "0" if text.endswith(".") else text


def _rare_limit(results, with_design):
    job = results["job"]
    if with_design and results.get("design"):
        return 1.0 / float(job["axes"]["rare_aep_1_in_x"] or 2000)
    rarest = min(row["aep"] for row in results["ams"])
    fitting = [limit for limit in RARE_LIMITS if limit <= rarest * 0.999]
    return fitting[0] if fitting else RARE_LIMITS[-1]


def _visible(values, keep):
    values = np.asarray([np.nan if v is None else v for v in values], float)
    return values[np.asarray(keep, bool) & np.isfinite(values)]


def figure(results, path, with_design=True):
    """The report figure. Returns the path written."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    job = results["job"]
    fsl = job["fsl"]
    with_design = with_design and bool(results.get("design"))
    frequent = frequency.EY1_AEP
    rare = _rare_limit(results, with_design)
    z_lo, z_hi = _z(frequent), _z(rare)

    grid_z = np.asarray(results["grid"]["z"], float)
    in_view = (grid_z >= z_lo - 1e-9) & (grid_z <= z_hi + 1e-9)
    ams = results["ams"]
    z = np.array([row["z"] for row in ams], float)
    level = np.array([row["level"] for row in ams], float)
    carried = np.array([bool(row["carried_over"]) for row in ams])
    storm_z = np.array([np.nan if row["storm_z"] is None else row["storm_z"]
                        for row in ams], float)

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    low_values, high_values = [], []

    fits = results.get("fits") or {}
    whole = fits.get("all") or {}
    storm = fits.get("storm") or {}
    shouldered = whole.get("form") == frequency.SHOULDERED

    def extent(block):
        """Where a curve is drawn: the shouldered upper limb is a straight line,
        so neither it nor its band is carried past the rarest maximum it was
        fitted to."""
        keep = in_view.copy()
        if block.get("form") == frequency.SHOULDERED and block.get("z_max") is not None:
            keep &= grid_z <= block["z_max"]
        return keep

    if whole.get("band_lo"):
        keep = extent(whole)
        lo = np.asarray(whole["band_lo"], float)
        hi = np.asarray(whole["band_hi"], float)
        ax.fill_between(grid_z[keep], lo[keep], hi[keep], color=RECORD_COLOUR,
                        alpha=0.13, linewidth=0, zorder=0,
                        label=f"90% band, all maxima ({whole['draws_used']:,} resamples)")
        low_values += list(lo[keep]); high_values += list(hi[keep])

    if with_design:
        design = results["design"]
        for curve in design["durations"].values():
            aep = np.asarray(curve["aep"], float)
            keep = aep <= 0.5
            ax.plot(np.asarray(curve["z"], float)[keep],
                    np.asarray(curve["level"], float)[keep], color=MUTED,
                    linewidth=0.6, alpha=0.30, zorder=1)
        ax.plot([], [], color=MUTED, linewidth=0.6, alpha=0.5,
                label="Monte Carlo durations")
        envelope = np.asarray(design["envelope"], float)
        ax.plot(grid_z, envelope, color=DESIGN_COLOUR, linewidth=1.8, zorder=4,
                label="Design flood envelope")
        low_values += list(envelope[in_view]); high_values += list(envelope[in_view])

    ax.scatter(z[~carried], level[~carried], s=14, color=RECORD_COLOUR,
               linewidth=0.8, zorder=5, label="Storm-driven maxima")
    if carried.any():
        ax.scatter(z[carried], level[carried], s=14, facecolor="white",
                   edgecolor=RECORD_COLOUR, linewidth=0.8, zorder=5,
                   label="Carried over the water year")
    shown = (z >= z_lo) & (z <= z_hi)
    low_values += list(level[shown]); high_values += list(level[shown])

    if whole.get("curve"):
        keep = extent(whole)
        ax.plot(grid_z[keep], np.asarray(whole["curve"], float)[keep],
                color=RECORD_COLOUR, linewidth=1.5 if shouldered else 1.4, zorder=3,
                label=f"Fit to all maxima (RMSE {whole['rmse']:.2f} m)")

    if carried.any() and storm:
        ax.scatter(storm_z[~carried], level[~carried], s=10, marker="x",
                   color=RECORD_COLOUR, alpha=0.55, linewidth=0.8, zorder=4,
                   label="…at censored positions")
        # The storm-driven curve starts no more frequent than its own most
        # frequent point: the frequent end of the scale is where the excluded
        # years belong, and it holds no data.
        keep = extent(storm) & (grid_z >= np.nanmin(storm_z[~carried]))
        if storm.get("band_lo"):
            for edge in (storm["band_lo"], storm["band_hi"]):
                ax.plot(grid_z[keep], np.asarray(edge, float)[keep],
                        color=RECORD_COLOUR, linewidth=0.8, linestyle=(0, (1, 1.6)),
                        alpha=0.75, zorder=2)
            ax.plot([], [], color=RECORD_COLOUR, linewidth=0.8,
                    linestyle=(0, (1, 1.6)), alpha=0.75, label="90% band, storm-driven")
        if storm.get("curve"):
            ax.plot(grid_z[keep], np.asarray(storm["curve"], float)[keep],
                    color=RECORD_COLOUR, linewidth=1.3, linestyle=(0, (2, 2)),
                    zorder=3, label=f"Fit to storm-driven (RMSE {storm['rmse']:.2f} m)")

    references = [(item["label"], float(item["level"]))
                  for item in job.get("reference_levels") or []
                  if item.get("level") is not None]
    if fsl is not None:
        references.append((job.get("fsl_label") or "FSL", float(fsl)))
    high_values += [level for _, level in references]

    y_lo = job["axes"].get("level_min")
    y_hi = job["axes"].get("level_max")
    finite_low = [v for v in low_values if v is not None and np.isfinite(v)]
    finite_high = [v for v in high_values if v is not None and np.isfinite(v)]
    if y_lo is None:
        y_lo = math.floor(min(finite_low) - 0.3) if finite_low else 0.0
    if y_hi is None:
        y_hi = math.ceil((max(finite_high) + 0.8) * 5) / 5 if finite_high else 1.0
    offset = 0.0105 * (y_hi - y_lo)

    for label, value in references:
        ax.axhline(value, color=MUTED, linewidth=0.8, linestyle=":")
        ax.text(_z(0.60), value + offset, f"{label} {level_text(value)} m",
                fontsize=BASE_FONTSIZE, color=INK)

    ax.set_xlim(z_lo, z_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_ylabel("Annual maximum lake level (m AHD)", color=INK,
                  fontsize=BASE_FONTSIZE)
    positions, labels = frequency_ticks(frequent, rare)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=BASE_FONTSIZE,
                       color=INK)
    ax.set_xlabel("Exceedance frequency (EY, then 1 in X AEP)", color=INK,
                  fontsize=BASE_FONTSIZE)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_color(INK)
    ax.tick_params(axis="both", colors=INK, labelsize=BASE_FONTSIZE, length=3)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92, facecolor="white",
              edgecolor="none", fontsize=BASE_FONTSIZE - 1, labelcolor=INK)
    fig.tight_layout()
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return path


# -- entry point ----------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="LakeLevelFrequency.py",
        description="Fit the annual maximum lake levels and draw the report figure.")
    parser.add_argument("--job", required=True, help="the analysis job (JSON)")
    parser.add_argument("--results", default=None,
                        help="where the fitted results are kept, and reused from")
    parser.add_argument("--png", default=None, help="write the report figure here")
    parser.add_argument("--without-design", action="store_true",
                        help="leave the Monte Carlo design floods off the figure")
    parser.add_argument("--ams-csv", default=None,
                        help="write the annual maximum series here")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        with open(args.job, encoding="utf-8") as stream:
            job = record.complete_job(json.load(stream))
    except (OSError, ValueError) as exc:
        print(f"ERROR: could not read the job {args.job}: {exc}")
        return 2

    missing = [path for path in record.job_inputs(job)
               if not os.path.isfile(path)
               and not (path in [s["path"] for s in job["design"]["sources"]]
                        and not job["design"]["include"])]
    if missing:
        print("ERROR: input files not found:")
        for path in missing:
            print(" ", path)
        return 2

    if args.ams_csv:
        ams, level_record = record.ams_for_job(job)
        table = record.with_positions(ams, bool(job["include_incomplete"]))
        table = table.sort_values("water_year").reset_index(drop=True)
        folder = os.path.dirname(os.path.abspath(args.ams_csv))
        os.makedirs(folder, exist_ok=True)
        with open(args.ams_csv, "w", encoding="utf-8", newline="") as stream:
            stream.write(record.ams_csv_text(table, record.ams_metadata(job, level_record)))
        print(f"Wrote {args.ams_csv}")

    if args.results or args.png:
        try:
            results = load_or_run(job, args.results)
        except ValueError as exc:
            print(f"ERROR: {exc}")
            return 1
        for name, block in (results.get("fits") or {}).items():
            if block.get("error"):
                print(f"WARNING: {name} maxima: {block['error']}")
        if args.png:
            figure(results, args.png, with_design=not args.without_design)
            print(f"Wrote {args.png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
