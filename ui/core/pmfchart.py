"""ECharts options for the PMF page: the patterns by duration, and the AEP fit.

Plain dicts, no nicegui, as core/resultchart.py. The fit chart uses the same
standard-normal-variate axis labelled in AEPs as the Results page.
"""

from __future__ import annotations

import math

import numpy as np

from . import ensemble
from .palette import MUTED, OCHRE, PALETTE, RULE, WATER
from .resultchart import aep_axis, numeric
from .results import Y_LABELS

# The standard AEPs the fit chart is ticked at, from the window up.
TICKS = [10_000, 20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000, 2_000_000,
         5_000_000, 10_000_000, 20_000_000, 50_000_000, 100_000_000]
MAX_POINTS = 1500


def _duration_label(hours) -> str:
    return f"{float(hours):g} h"


def _reference_lines(levels) -> dict | None:
    lines = [{"yAxis": float(item["level"]), "name": item.get("label", ""),
              "label": {"formatter": f"{item.get('label', '')} {float(item['level']):.2f}",
                        "position": "insideEndTop", "color": MUTED},
              "lineStyle": {"color": RULE, "type": "dashed"}}
             for item in levels or [] if numeric(item.get("level")) is not None]
    if not lines:
        return None
    return {"symbol": "none", "silent": True, "data": lines}


def box_chart(stats, points, result: str = "level", *, highest=None,
              reference_levels=None) -> dict:
    """One box per duration over the temporal patterns, with every pattern as a dot.

    The box's middle line is Bryan's median pick, not a quartile computed here,
    so it marks the same event ``csv/<name>_critical.csv`` names.
    """
    if stats is None or stats.empty:
        return {"series": []}
    durations = list(stats.index)
    categories = [_duration_label(hours) for hours in durations]
    position = {hours: index for index, hours in enumerate(durations)}
    boxes = [[stats.loc[h, "min"], stats.loc[h, "q1"], stats.loc[h, "median"],
              stats.loc[h, "q3"], stats.loc[h, "max"]] for h in durations]
    dots = [[position[hours], numeric(value), pattern]
            for hours, value, pattern in points if hours in position]

    series = [
        {"name": "patterns", "type": "boxplot",
         "itemStyle": {"color": "#FFFFFF", "borderColor": WATER},
         "data": [[numeric(v) for v in box] for box in boxes],
         "tooltip": {"trigger": "item"}},
        {"name": "each pattern", "type": "scatter", "symbolSize": 5,
         "itemStyle": {"color": PALETTE[2], "opacity": 0.55},
         "data": dots},
    ]
    if highest is not None and highest.found and highest.duration in position:
        value = getattr(highest, result)
        series.append({"name": "highest", "type": "scatter", "symbol": "diamond",
                       "symbolSize": 14, "itemStyle": {"color": OCHRE},
                       "data": [[position[highest.duration], numeric(value),
                                 highest.pattern]]})
    if result == "level":
        lines = _reference_lines(reference_levels)
        if lines:
            series[0]["markLine"] = lines

    return {
        "tooltip": {"trigger": "item"},
        "legend": {"top": 0},
        "grid": {"left": 70, "right": 30, "top": 40, "bottom": 50},
        "xAxis": {"type": "category", "data": categories, "name": "Storm duration",
                  "nameLocation": "middle", "nameGap": 30},
        "yAxis": {"type": "value", "scale": True, "name": Y_LABELS.get(result, result),
                  "nameLocation": "middle", "nameGap": 50,
                  "splitLine": {"lineStyle": {"opacity": 0.25}}},
        "series": series,
    }


def fit_chart(real, estimate, *, pmp_aep=None, context_factor=20.0) -> dict:
    """The top of the Monte Carlo sample, the window, the fit, and where it lands."""
    if real is None or not len(real.z):
        return {"series": []}
    z_window = ensemble.variate(1.0 / estimate.lower_aep)
    z_upper = (ensemble.variate(1.0 / estimate.upper_aep)
               if math.isfinite(estimate.upper_aep) else math.inf)
    z_context = ensemble.variate(context_factor / estimate.lower_aep)
    shown = real.z >= z_context
    z, value = real.z[shown], real.value[shown]
    if len(z) > MAX_POINTS:                  # even thinning, keeping the order
        keep = np.linspace(0, len(z) - 1, MAX_POINTS).astype(int)
        z, value = z[keep], value[keep]
    inside = (z >= z_window) & (z <= z_upper)

    series = [
        {"name": "realisations", "type": "scatter", "symbolSize": 3,
         "itemStyle": {"color": MUTED, "opacity": 0.5},
         "data": [[numeric(a), numeric(b)] for a, b in zip(z[~inside], value[~inside])]},
        {"name": "in the fit window", "type": "scatter", "symbolSize": 4,
         "itemStyle": {"color": WATER, "opacity": 0.7},
         "data": [[numeric(a), numeric(b)] for a, b in zip(z[inside], value[inside])]},
    ]
    curve = ensemble.fitted_curve(estimate, real)
    if curve:
        series.append({"name": f"fit, degree {estimate.degree}", "type": "line",
                       "symbol": "none", "lineStyle": {"width": 2, "color": PALETTE[1]},
                       "itemStyle": {"color": PALETTE[1]},
                       "data": [[numeric(a), numeric(b)] for a, b in curve]})
    marks = [{"yAxis": estimate.target, "name": "PMF",
              "label": {"formatter": f"PMF {estimate.target:.2f}", "color": OCHRE},
              "lineStyle": {"color": OCHRE}}]
    if estimate.ok:
        marks.append({"xAxis": estimate.z,
                      "label": {"formatter": f"1 in {estimate.aep:,.0f}", "color": OCHRE},
                      "lineStyle": {"color": OCHRE}})
    if pmp_aep:
        marks.append({"xAxis": ensemble.variate(1.0 / float(pmp_aep)),
                      "label": {"formatter": "AEP of PMP", "color": MUTED},
                      "lineStyle": {"color": MUTED, "type": "dashed"}})
    series[0]["markLine"] = {"symbol": "none", "silent": True, "data": marks}

    low_z = float(z.min()) if len(z) else z_window
    high_z = max(float(real.z.max()), estimate.z if estimate.ok else -math.inf)
    ticks = [aep for aep in TICKS
             if low_z - 0.05 <= ensemble.variate(1.0 / aep) <= high_z + 0.05]
    axis = aep_axis(ticks or TICKS[:3])
    axis["min"], axis["max"] = round(low_z, 3), round(high_z + 0.05, 3)
    axis["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}
    return {
        "tooltip": {"trigger": "item"},
        "legend": {"top": 0},
        "grid": {"left": 70, "right": 30, "top": 40, "bottom": 60},
        "xAxis": axis,
        "yAxis": {"type": "value", "scale": True,
                  "name": Y_LABELS.get(real.result, real.result),
                  "nameLocation": "middle", "nameGap": 50,
                  "splitLine": {"lineStyle": {"opacity": 0.25}}},
        "series": series,
    }
