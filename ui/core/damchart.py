"""Charts for the Downstream Results tab: the dams against their crests.

Two dams' levels are metres apart in AHD - Callide's crest is 219.13 m, Kroombit's 270.70 m -
so they are drawn as **height above each dam's own crest level**: one axis, zero the crest,
and whether and when each reached it is read straight off. Where a dam has a PMF level it is
drawn on the same basis.

The scatter takes every finished run of every selection in the project: one point per run,
the first dam's margin to its crest across, the second's up. The quadrant above and right of
zero is the runs where both reached their crests.

Option dicts only, coloured from the house palette; no nicegui.
"""

from __future__ import annotations

from .palette import MUTED, OCHRE, PALETTE
from .urbsresults import CREST, level_at

LIMIT = 1500           # points per line, thinned evenly beyond this


def _thin(times, values) -> list:
    step = max(1, len(times) // LIMIT)
    return [[float(t), round(float(v), 3)] for t, v in zip(times[::step], values[::step])]


def crest_chart(series_by_dam: dict, references: dict, gwl, title="") -> dict:
    """Each dam's level less its crest level through one run. ``series_by_dam`` maps a dam
    to that series (from ``urbsresults.crest_series``); ``references`` its dam_levels."""
    drawn = []
    for number, (dam, series) in enumerate(series_by_dam.items()):
        if series is None or series.empty:
            continue
        colour = PALETTE[number % len(PALETTE)]
        line = {"name": dam, "type": "line", "showSymbol": False, "smooth": False,
                "lineStyle": {"width": 1.8}, "itemStyle": {"color": colour},
                "data": _thin(list(series.index), list(series.values))}
        crest = level_at(references.get(dam, {}).get(CREST), gwl)
        pmf = level_at(references.get(dam, {}).get("PMF"), gwl)
        if crest is not None and pmf is not None:
            line["markLine"] = {"symbol": "none", "silent": True,
                                "lineStyle": {"color": colour, "type": "dotted"},
                                "label": {"formatter": f"{dam} PMF", "color": colour},
                                "data": [{"yAxis": round(pmf - crest, 3)}]}
        drawn.append(line)
    if not drawn:
        return {"series": []}
    # the crest itself, common to every dam on this axis
    crest_line = {"name": "crest", "type": "line", "data": [], "silent": True,
                  "markLine": {"symbol": "none", "silent": True,
                               "lineStyle": {"color": OCHRE, "type": "dashed", "width": 1.5},
                               "label": {"formatter": "dam crest", "color": OCHRE},
                               "data": [{"yAxis": 0}]}}
    return {
        "title": {"text": title, "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "axis", ":valueFormatter": "(v) => v.toFixed(2) + ' m'"},
        "legend": {"top": 30},
        "grid": {"left": 70, "right": 40, "top": 64, "bottom": 50},
        "xAxis": {"type": "value", "name": "Time (hours)", "nameLocation": "middle",
                  "nameGap": 30, "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}}},
        "yAxis": {"type": "value", "name": "Above the dam crest (m)", "nameGap": 45,
                  "scale": True, "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}}},
        "series": drawn + [crest_line],
    }


def margin_scatter(points_by_group: dict, x_dam: str, y_dam: str) -> dict:
    """Each finished run's peak margins to the two dams' crests, coloured by group."""
    series = []
    for number, (group, points) in enumerate(points_by_group.items()):
        series.append({"name": group, "type": "scatter", "symbolSize": 9,
                       "itemStyle": {"color": PALETTE[number % len(PALETTE)]},
                       "data": [{"value": [round(x, 3), round(y, 3)], "name": run}
                                for x, y, run in points]})
    if not series:
        return {"series": []}
    zero = {"symbol": "none", "silent": True,
            "lineStyle": {"color": OCHRE, "type": "dashed", "width": 1.2},
            "label": {"show": False}}
    series[0]["markLine"] = {**zero, "data": [{"xAxis": 0}, {"yAxis": 0}]}
    return {
        "title": {"text": f"Peak level above the crest: {x_dam} against {y_dam}",
                  "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "item",
                    ":formatter": "(p) => `${p.seriesName}<br>${p.name}<br>"
                                 f"{x_dam} ${{p.value[0].toFixed(2)}} m, "
                                 f"{y_dam} ${{p.value[1].toFixed(2)}} m`"},
        "legend": {"type": "scroll", "top": 30},
        "grid": {"left": 70, "right": 40, "top": 64, "bottom": 55},
        # no "scale": the axes keep zero in view, which is the line that matters here
        "xAxis": {"type": "value", "name": f"{x_dam} above its crest (m)",
                  "nameLocation": "middle", "nameGap": 32,
                  "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}}},
        "yAxis": {"type": "value", "name": f"{y_dam} above its crest (m)", "nameGap": 45,
                  "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}},
                  "axisLabel": {"color": MUTED}},
        "series": series,
    }
