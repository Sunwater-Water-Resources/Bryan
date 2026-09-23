"""Draw one report frequency-curve figure from a job file.

    python util/ReportFigure.py <figure>.json

The launcher's Figures page writes the job - every series already reduced to
(AEP '1 in X', value) pairs, with its label - and runs this under Bryan's
interpreter, because the UI environment has no matplotlib. This script only
draws: it reads no results of its own, so the PNG carries exactly the numbers and
labels the page previewed, and the job file left beside the PNG records them.

The look is PlotFrequencyCurves_v03.py's, the figures already in the Callide
report: 6.3 x 4.0 in at 300 dpi for an A4 page, 8 pt text, the standard AEPs as
the x ticks on a standard-normal-variate axis, a logarithmic axis with thousands
separators for flows, and the same markers - posterior mode solid and mean dashed
in black, the 90% credible interval in faint grey, the annual maxima as dots, a
paleoflood lower bound as an upward triangle, and the AEP of the PMP as a
vertical line.
"""

from __future__ import annotations

import json
import math
import sys
from itertools import cycle
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator           # noqa: E402

FIG_WIDTH_IN, FIG_HEIGHT_IN, FIG_DPI, BASE_FONT_SIZE = 6.3, 4.0, 300, 8
LINE_STYLES = ["-", "--", "-.", ":"]
_NORMAL = NormalDist()


def variate(aep: float) -> float:
    """z for '1 in X', from the upper tail so a rare AEP keeps its digits."""
    return -_NORMAL.inv_cdf(1.0 / float(aep))


def standard_aeps(lower: float, upper: float) -> list:
    out, aep = [], float(lower)
    while aep <= float(upper):
        out.append(aep)
        aep = aep * 5 // 2 if math.log10(aep / 2) % 1 == 0 else aep * 2
    return out


def _xy(points):
    return [variate(a) for a, _ in points], [v for _, v in points]


def draw(job: dict) -> Path:
    plt.rcParams.update({"font.size": BASE_FONT_SIZE})
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN), dpi=FIG_DPI)
    ax.grid(which="both", linewidth=0.4, alpha=0.6)
    styles = cycle(LINE_STYLES)
    ffa_styles = {"ffa-mode": "-", "ffa-mean": "--"}

    for series in job["series"]:
        points = series["points"]
        if not points:
            continue
        x, y = _xy(points)
        label = series["label"] if series.get("legend", True) else None
        style = series.get("style", "line")
        if style == "line":
            ax.plot(x, y, next(styles), label=label)
        elif style in ffa_styles:
            ax.plot(x, y, ffa_styles[style], color="black", label=label)
        elif style == "ffa-ci":
            ax.plot(x, y, "--", color="k", alpha=0.2, label=label)
        elif style == "ams":
            ax.plot(x, y, "o", markeredgewidth=0, markersize=4, label=label)
        elif style == "paleo":
            ax.plot(x, y, "^", color="firebrick", markersize=5, markeredgewidth=0,
                    linestyle="none", label=label)

    low, high = float(job["min_aep"]), float(job["max_aep"])
    z_range = [variate(low), variate(high)]
    for item in job.get("reference_levels") or []:
        ax.plot(z_range, [item["level"]] * 2, "-", label=item.get("label") or None,
                color=item.get("colour") or "grey", alpha=0.5)
    if job.get("aep_of_pmp"):
        ax.axvline(variate(job["aep_of_pmp"]), color="k")

    if job.get("log_y"):
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    bottom, top = job.get("y_min"), job.get("y_max")
    if job.get("log_y") and bottom is None:
        # v03: the axis starts at the power of ten below the lowest flow drawn
        lowest = [v for s in job["series"] if s.get("style") != "ffa-ci"
                  for _, v in s["points"] if v and v > 0]
        if lowest:
            bottom = 10 ** math.floor(math.log10(min(lowest)))
    if bottom is not None or top is not None:
        ax.set_ylim(bottom=bottom, top=top)

    ticks = standard_aeps(low, high)
    ax.set_xticks([variate(aep) for aep in ticks])
    ax.set_xticklabels([f"{aep:.0f}" for aep in ticks], rotation=90)
    ax.set_xlabel("AEP (1 in X)")
    ax.set_ylabel(job.get("y_label") or "")
    ax.legend(loc="upper left", ncol=1, fontsize="x-small")
    if job.get("show_title") and job.get("title"):
        ax.set_title(job["title"])

    png = Path(job["png"])
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    return png


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    job = json.loads(Path(argv[0]).read_text(encoding="utf-8"))
    png = draw(job)
    print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
