"""ECharts preview of a report figure - the same series the PNG is drawn from.

Plain dicts, no nicegui. The PNG's line styles and markers are matplotlib's
(util/ReportFigure.py); the preview carries the same series, labels, axis and
reference lines in the house palette, so what is checked here is what is drawn.
"""

from __future__ import annotations

from .figures import FigureData, logarithmic, result_label, standard_aeps
from .palette import MUTED, PALETTE
from .resultchart import aep_axis, numeric
from .results import normal_variate

DASHES = ["solid", "dashed", "dotted", "dashed"]


def _points(series):
    return [[numeric(normal_variate(aep)), numeric(value)] for aep, value in series.points]


def preview(data: FigureData) -> dict:
    spec = data.spec
    low, high = float(spec.get("min_aep") or 2), float(spec.get("max_aep") or 2e6)
    out, line_index = [], 0
    for series in data.series:
        if not series.points:
            continue
        base = {"name": series.label, "data": _points(series), "z": 3}
        if series.style == "line":
            base.update(type="line", symbol="none",
                        lineStyle={"width": 2, "type": DASHES[line_index % len(DASHES)]},
                        itemStyle={"color": PALETTE[line_index % len(PALETTE)]})
            line_index += 1
        elif series.style in ("ffa-mode", "ffa-mean"):
            base.update(type="line", symbol="none", itemStyle={"color": "#031E2F"},
                        lineStyle={"width": 2, "type": "solid" if series.style == "ffa-mode"
                                   else "dashed"})
        elif series.style == "ffa-ci":
            base.update(type="line", symbol="none", silent=True,
                        itemStyle={"color": MUTED},
                        lineStyle={"width": 1, "type": "dashed", "opacity": 0.4})
        elif series.style == "ams":
            base.update(type="scatter", symbolSize=5, itemStyle={"color": PALETTE[0]})
        elif series.style == "paleo":
            base.update(type="scatter", symbol="triangle", symbolSize=8,
                        itemStyle={"color": "#B3261E"})
        out.append(base)

    marks = [{"yAxis": float(item["level"]), "name": item.get("label", ""),
              "label": {"formatter": item.get("label", ""), "position": "insideEndTop"},
              "lineStyle": {"color": item.get("colour") or MUTED, "opacity": 0.6,
                            "type": "solid"}}
             for item in spec.get("reference_levels") or [] if numeric(item.get("level")) is not None]
    if spec.get("aep_of_pmp"):
        marks.append({"xAxis": numeric(normal_variate(spec["aep_of_pmp"])),
                      "label": {"formatter": "AEP of PMP"},
                      "lineStyle": {"color": "#031E2F", "type": "solid"}})
    if marks and out:
        out[0]["markLine"] = {"symbol": "none", "silent": True, "data": marks}

    axis = aep_axis(standard_aeps(low, high))
    axis["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}
    key = data.key
    y_axis = {"type": "log" if logarithmic(key) else "value", "scale": True,
              "name": result_label(key), "nameLocation": "middle", "nameGap": 50,
              "splitLine": {"lineStyle": {"opacity": 0.25}}}
    if spec.get("y_min") is not None:
        y_axis["min"] = spec["y_min"]
    if spec.get("y_max") is not None:
        y_axis["max"] = spec["y_max"]
    return {
        "tooltip": {"trigger": "item"},
        "legend": {"type": "scroll", "top": 0,
                   "data": [s.label for s in data.series if s.legend and s.points]},
        "grid": {"left": 70, "right": 30, "top": 40, "bottom": 60},
        "xAxis": axis, "yAxis": y_axis, "series": out,
    }
