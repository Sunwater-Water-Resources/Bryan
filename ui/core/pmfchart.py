"""ECharts options for the PMF page: the patterns by duration, and the AEP fit.

Plain dicts, no nicegui, as core/resultchart.py. The fit chart uses the same
standard-normal-variate axis labelled in AEPs as the Results page.
"""

from __future__ import annotations

import math

import numpy as np

from . import ensemble
from .palette import INK, MUTED, OCHRE, PALETTE, RULE, WATER
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


# -- the world's maximum observed floods ----------------------------------------

#: Log axes, a decade either side of Table 1.
WORLD_AREA_RANGE = (0.1, 10_000_000)
WORLD_PEAK_RANGE = (10, 1_000_000)
# Thousands separators on a log axis: ECharts writes 100000 otherwise.
_GROUPED = "value => value.toLocaleString('en-AU', {maximumFractionDigits: 1})"

#: The page draws the chart 32rem high. The legend stands to the right of the
#: plot, one column of up to that height, and wraps into more columns past it.
#: A legend under the plot wraps into as many rows as the groups need, and ECharts
#: does not move the plot up for them.
WORLD_CHART_HEIGHT = 512
_LEGEND_TOP = 20     # px above and below the legend
_LEGEND_ROW = 24     # px per entry: the 14 px icon and the 10 px gap
_LEGEND_ICON = 40    # px per column besides the text: the icon, its gap and air
_LEGEND_CHAR = 7.5   # px per character of a 12 px label, generously
_LEGEND_PAD = 20     # px between the plot and the legend


def legend_width(names, height: float = WORLD_CHART_HEIGHT) -> float:
    """The room a vertical legend of ``names`` takes on the right of a chart
    ``height`` px high. ``copy_chart_script`` works it out again, the same way,
    for the entries left switched on."""
    if not names:
        return _LEGEND_PAD
    rows = max(1, int((height - 2 * _LEGEND_TOP) // _LEGEND_ROW))
    columns = math.ceil(len(names) / rows)
    return _LEGEND_PAD + columns * (_LEGEND_ICON + _LEGEND_CHAR * max(map(len, names)))


def world_floods_chart(floods, area_km2=None, dam_floods=(), *, herschy=False) -> dict:
    """Peak discharge against catchment area: the world's maximum observed floods
    (gauged filled, surveyed open), the WMO envelope over them, and each PMF group's
    PMF and PMPF inflow at the dam's catchment area.

    Drawn without animation, on white, so a copy taken the moment it is drawn is
    the finished chart; see ``copy_chart_script``.
    """
    from . import worldfloods

    def flood_point(flood):
        note = " (surveyed)" if flood.surveyed else " (modified value)" if flood.modified else ""
        return {"value": [flood.area_km2, flood.peak_m3s],
                "name": f"{flood.station}, {flood.country}, {flood.date}{note}"}

    gauged = [flood_point(f) for f in floods if not f.surveyed]
    surveyed = [flood_point(f) for f in floods if f.surveyed]
    series = [
        {"name": "World maximum floods", "type": "scatter", "symbol": "circle",
         "symbolSize": 6, "itemStyle": {"color": MUTED}, "data": gauged, "z": 2},
        {"name": "World maximum floods (surveyed)", "type": "scatter", "symbol": "emptyCircle",
         "symbolSize": 6, "itemStyle": {"color": MUTED, "borderWidth": 1.2},
         "data": surveyed, "z": 2},
        {"name": worldfloods.WMO_LABEL, "type": "line", "symbol": "none", "silent": True,
         "lineStyle": {"color": INK, "width": 1.6}, "itemStyle": {"color": INK},
         "data": [[a, q] for a, q in worldfloods.envelope_curve(worldfloods.wmo_envelope)],
         "z": 1},
    ]
    if herschy:
        series.append(
            {"name": worldfloods.HERSCHY_LABEL, "type": "line", "symbol": "none",
             "silent": True, "lineStyle": {"color": MUTED, "width": 1.2, "type": "dashed"},
             "itemStyle": {"color": MUTED},
             "data": [[a, q] for a, q in worldfloods.envelope_curve(
                 worldfloods.herschy_envelope)], "z": 1})

    area = numeric(area_km2)
    if area and area > 0:
        series[2]["markLine"] = {
            "symbol": "none", "silent": True,
            "lineStyle": {"color": RULE, "type": "dashed"},
            "label": {"formatter": f"{area:,.0f} km²", "color": MUTED,
                      "position": "insideEndTop"},
            "data": [{"xAxis": area}]}
        for index, flood in enumerate(dam_floods):
            colour = PALETTE[(index + 1) % len(PALETTE)]
            for name, value, symbol in (("PMF", flood.pmf, "diamond"),
                                        ("PMPF", flood.pmpf, "triangle")):
                peak = numeric(value)
                if peak is None or peak <= 0:
                    continue
                label = f"{flood.label} {name}".strip() if len(dam_floods) > 1 else name
                series.append({"name": label, "type": "scatter", "symbol": symbol,
                               "symbolSize": 13, "z": 5,
                               "itemStyle": {"color": colour, "borderColor": "#FFFFFF",
                                             "borderWidth": 1},
                               "data": [{"value": [area, peak],
                                         "name": f"{label}: {peak:,.0f} m³/s"}]})

    axis = {"type": "log", "logBase": 10, "nameLocation": "middle",
            "axisLabel": {":formatter": _GROUPED},
            "splitLine": {"show": True, "lineStyle": {"opacity": 0.35}},
            "minorTick": {"show": True}, "minorSplitLine": {"show": True,
                                                           "lineStyle": {"opacity": 0.12}}}
    # every entry named, so a copy can take the switched-off ones out and put them back
    names = [s["name"] for s in series]
    return {
        "animation": False,
        "backgroundColor": "#FFFFFF",
        "tooltip": {"trigger": "item"},
        "legend": {"data": names, "orient": "vertical", "right": 10,
                   "top": _LEGEND_TOP, "bottom": _LEGEND_TOP,
                   "itemGap": _LEGEND_ROW - 14, "itemHeight": 14},
        "grid": {"left": 80, "right": legend_width(names), "top": 20, "bottom": 60},
        "xAxis": {**axis, "name": "Catchment area (km²)", "nameGap": 30,
                  "min": WORLD_AREA_RANGE[0], "max": WORLD_AREA_RANGE[1]},
        "yAxis": {**axis, "name": "Peak discharge (m³/s)", "nameGap": 62,
                  "min": WORLD_PEAK_RANGE[0], "max": WORLD_PEAK_RANGE[1]},
        "series": series,
    }


def copy_chart_script(element_id: int, filename: str = "", *, pixel_ratio: int = 3) -> str:
    """JavaScript that puts an ``ui.echart`` on the clipboard as a PNG, or, given a
    ``filename``, downloads it. Returns 'image', 'saved' or 'failed: <why>'.

    The picture is ECharts' own render at ``pixel_ratio`` times the size on
    screen, on white, so it pastes into Word sharp at the width of a page. The
    clipboard needs ``ClipboardItem`` and a secure context, which localhost is; a
    browser without it is told to use Save instead.

    Legend entries switched off on screen are left out of the picture rather than
    drawn greyed: the legend is cut to the entries still on for the render and put
    back after it, and a vertical legend's room is worked out again for what is left
    (``legend_width``). ECharts redraws synchronously on ``setOption``, so the render
    in between is of the cut legend.
    """
    import json
    return (
        "(async () => {"
        f"const element = getElement({int(element_id)});"
        "const chart = element && element.chart;"
        "if (!chart) return 'failed: the chart has not been drawn yet';"
        "const option = chart.getOption();"
        "const legend = (option.legend || [])[0];"
        "let restore = null;"
        "if (legend && Array.isArray(legend.data) && legend.data.length) {"
        "  const nameOf = item => typeof item === 'string' ? item : item.name;"
        "  const selected = legend.selected || {};"
        "  const shown = legend.data.filter(item => selected[nameOf(item)] !== false);"
        "  if (shown.length < legend.data.length) {"
        "    restore = {legend: {data: legend.data}};"
        "    const change = {legend: {data: shown}};"
        "    if (legend.orient === 'vertical' && option.grid && option.grid.length) {"
        f"      const rows = Math.max(1, Math.floor((chart.getHeight() - {2 * _LEGEND_TOP})"
        f" / {_LEGEND_ROW}));"
        "      const columns = Math.ceil(shown.length / rows);"
        "      const longest = Math.max(0, ...shown.map(item => nameOf(item).length));"
        "      restore.grid = {right: option.grid[0].right};"
        f"      change.grid = {{right: shown.length ? {_LEGEND_PAD} + columns * ({_LEGEND_ICON}"
        f" + {_LEGEND_CHAR} * longest) : {_LEGEND_PAD}}};"
        "    }"
        "    chart.setOption(change);"
        "  }"
        "}"
        "let url;"
        "try {"
        "  url = chart.getDataURL({type: 'png', backgroundColor: '#FFFFFF',"
        f"   pixelRatio: {int(pixel_ratio)}}});"
        "} finally { if (restore) chart.setOption(restore); }"
        f"const filename = {json.dumps(filename)};"
        "if (filename) {"
        "  const link = document.createElement('a');"
        "  link.href = url; link.download = filename; link.click(); return 'saved';"
        "}"
        "try {"
        "  if (!(window.ClipboardItem && navigator.clipboard && navigator.clipboard.write))"
        "    return 'failed: this browser cannot put a picture on the clipboard';"
        "  const blob = await (await fetch(url)).blob();"
        "  await navigator.clipboard.write([new ClipboardItem({'image/png': blob})]);"
        "  return 'image';"
        "} catch (error) { return 'failed: ' + error; }"
        "})()")
