"""ECharts options for the duration comparison plot.

Plain dicts, no nicegui - ``ui.echart`` is handed the result. That keeps the
whole chart testable, and keeps ``core/`` importable without a browser stack.

The x axis is the standard normal variate with the standard AEPs as ticks,
which is how UtilModule.plot_durations draws it (``plt.xticks(std_z, std_aeps)``).
ECharts labels a value axis at its own chosen ticks, so two things are set:
``customValues`` pins the ticks to the AEPs that were actually evaluated, and
the formatter can label *any* z, so the axis still reads as AEPs on an older
ECharts that ignores customValues.

Values are converted to plain floats on the way out. NaN cannot be represented
in JSON, so gaps become ``None``, which ECharts draws as a break in the line -
the honest rendering for an AEP a duration never reached.
"""

from __future__ import annotations

import json
import math

from .results import ABSOLUTE, format_aep, normal_variate, y_axis

# Chosen to stay distinguishable on both the light and dark themes the UI
# follows, and to survive being printed in greyscale in a report.
PALETTE = ("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1",
           "#76B7B2", "#EDC948", "#9C755F", "#FF9DA7", "#8CD17D")

ENVELOPE_COLOUR = "#888888"

# Abramowitz & Stegun 26.2.17 for the upper tail, so a tick anywhere on the
# axis can be labelled '1 in X' even when it is not one of ours.
_TAIL_JS = """
const z = Math.abs(v);
const t = 1 / (1 + 0.2316419 * z);
const d = 0.3989422804014327 * Math.exp(-z * z / 2);
let q = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
        t * (-1.821255978 + t * 1.330274429))));
if (v < 0) { q = 1 - q; }
if (!(q > 0)) { return ''; }
const a = 1 / q;
return a >= 1000 ? Math.round(a).toLocaleString()
     : (a >= 10 ? String(Math.round(a)) : a.toFixed(1));
"""


def numeric(value):
    """A JSON-safe number: NaN and infinities become None."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) or math.isinf(number) else number


def colour_for(labels) -> dict:
    return {label: PALETTE[index % len(PALETTE)]
            for index, label in enumerate(labels)}


def _axis_label_formatter(aeps) -> str:
    exact = {f"{normal_variate(aep):.4f}": format_aep(aep)
             for aep in aeps if numeric(normal_variate(aep)) is not None}
    return ("(v) => { const m = " + json.dumps(exact) + ";"
            " const k = v.toFixed(4); if (m[k] !== undefined) { return m[k]; }"
            + _TAIL_JS + "}")


def _series_data(index, values):
    return [[numeric(normal_variate(aep)), numeric(value)]
            for aep, value in zip(index, values)]


def aep_axis(aeps, *, name_gap=32) -> dict:
    """The standard normal variate axis, labelled with the AEPs it came from."""
    zs = [numeric(normal_variate(aep)) for aep in aeps]
    known = [z for z in zs if z is not None]
    return {
        "type": "value",
        "name": "AEP (1 in X)",
        "nameLocation": "middle",
        "nameGap": name_gap,
        "min": min(known) if known else None,
        "max": max(known) if known else None,
        "axisLabel": {"customValues": known, "hideOverlap": True,
                      ":formatter": _axis_label_formatter(aeps)},
        "axisTick": {"customValues": known},
    }


def duration_chart(comparison, analysis, key, *, show_envelope=True,
                   show_markup=True, title="") -> dict:
    """The frequency curves, one per duration, with the envelope over them."""
    frame = comparison.frame
    if frame.empty:
        return {"series": []}

    aeps = list(frame.index)
    label, logarithmic = y_axis(key)
    colours = colour_for(frame.columns)

    series = [{
        "name": str(column),
        "type": "line",
        "smooth": False,
        "symbolSize": 5,
        "connectNulls": False,
        "itemStyle": {"color": colours[column]},
        "lineStyle": {"width": 2},
        "emphasis": {"focus": "series"},
        "data": _series_data(aeps, frame[column]),
    } for column in frame.columns]

    if show_envelope:
        envelope = {
            "name": "envelope",
            "type": "line",
            "symbol": "none",
            "z": 1,
            "itemStyle": {"color": ENVELOPE_COLOUR},
            "lineStyle": {"width": 4, "type": "dashed", "opacity": 0.75},
            "data": _series_data(aeps, analysis.envelope),
        }
        if show_markup:
            envelope["markArea"] = _bands(analysis.bands, colours)
            points = _switch_points(analysis, frame)
            if points:
                envelope["markPoint"] = points
        series.append(envelope)

    horizontal = aep_axis(aeps)
    horizontal["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}

    return {
        "title": {"text": title, "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "cross"},
                    ":valueFormatter": "(v) => v == null ? '-' : Number(v).toPrecision(5)"},
        "legend": {"type": "scroll", "top": 24},
        "grid": {"left": 60, "right": 30, "top": 60, "bottom": 60},
        "xAxis": horizontal,
        "yAxis": {
            "type": "log" if logarithmic else "value",
            "name": label,
            "nameLocation": "middle",
            "nameGap": 45,
            "scale": True,
            "splitLine": {"lineStyle": {"opacity": 0.25}},
        },
        "series": series,
    }


def _bands(bands, colours) -> dict:
    """Shade the AEP axis by which duration owns the envelope."""
    return {
        "silent": True,
        "itemStyle": {"opacity": 0.12},
        "label": {"show": True, "position": "insideTop", "fontSize": 10,
                  "formatter": "{b}"},
        "data": [[{"name": band.label,
                   "xAxis": numeric(band.z_from),
                   "itemStyle": {"color": colours.get(band.label, ENVELOPE_COLOUR)}},
                  {"xAxis": numeric(band.z_to)}] for band in bands],
    }


def _switch_points(analysis, frame) -> dict | None:
    """Pin each real crossover. A switch inside the noise floor gets no pin.

    Marking a hop that came out of sampling noise would dress it up as a
    finding; it is reported in the warnings instead.
    """
    data = []
    for switch in analysis.switches:
        strength = numeric(switch.strength)
        if strength is not None and strength < analysis.noise_floor:
            continue
        z = numeric(normal_variate(switch.aep))
        value = numeric(analysis.envelope.get(switch.aep))
        if z is None or value is None:
            continue
        data.append({"coord": [z, value],
                     "value": f"{switch.before} → {switch.after}"})
    if not data:
        return None
    return {"symbol": "pin", "symbolSize": 48, "symbolOffset": [0, -6],
            "itemStyle": {"color": ENVELOPE_COLOUR, "opacity": 0.85},
            "label": {"fontSize": 9, "color": "#fff", "formatter": "{c}"},
            "data": data}


def critical_duration_chart(comparison, analysis) -> dict:
    """Critical duration against AEP - the transition, drawn directly.

    For lake level this should fall from left to right: long durations while the
    storage is filling, shorter ones on the rare tail as the dam behaves as a
    conveyance. A jagged line here usually means the duration curves are
    coincident and idxmax is picking up noise - check the margin column.
    """
    durations = comparison.durations
    if analysis.critical.empty or any(
            durations.get(owner) is None for owner in analysis.critical):
        return {}

    aeps = list(analysis.critical.index)
    data = [[numeric(normal_variate(aep)), numeric(durations[owner])]
            for aep, owner in analysis.critical.items()]

    return {
        "tooltip": {"trigger": "axis"},
        "grid": {"left": 60, "right": 30, "top": 20, "bottom": 50},
        "xAxis": aep_axis(aeps, name_gap=30),
        "yAxis": {"type": "log", "name": "Critical duration (h)",
                  "nameLocation": "middle", "nameGap": 40},
        "series": [{"type": "line", "step": "middle", "symbolSize": 6,
                    "itemStyle": {"color": PALETTE[0]}, "data": data}],
    }


# -- overlaying groups -------------------------------------------------------
#
# One line per group, and no mark-up: the maximum *across* groups is not a
# quantity - they are scenarios, not alternatives to be enveloped - so there is
# no envelope of envelopes, no bands and no crossover pins here. What replaces
# them is the change from a baseline, drawn underneath on the same axis.

def overlay_chart(overlay, *, title="") -> dict:
    """One group's envelope per line."""
    frame = overlay.frame
    if frame.empty:
        return {"series": []}

    aeps = list(frame.index)
    label, logarithmic = y_axis(overlay.key)
    colours = colour_for(frame.columns)

    horizontal = aep_axis(aeps)
    horizontal["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}

    return {
        "title": {"text": title, "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "cross"},
                    ":valueFormatter": "(v) => v == null ? '-' : Number(v).toPrecision(5)"},
        "legend": {"type": "scroll", "top": 24},
        "grid": {"left": 60, "right": 30, "top": 60, "bottom": 60},
        "xAxis": horizontal,
        "yAxis": {
            "type": "log" if logarithmic else "value",
            "name": label,
            "nameLocation": "middle",
            "nameGap": 45,
            "scale": True,
            "splitLine": {"lineStyle": {"opacity": 0.25}},
        },
        "series": [{
            "name": str(column),
            "type": "line",
            "smooth": False,
            "symbolSize": 5,
            "connectNulls": False,
            "itemStyle": {"color": colours[column]},
            "lineStyle": {"width": 2},
            "emphasis": {"focus": "series"},
            "data": _series_data(aeps, frame[column]),
        } for column in frame.columns],
    }


def delta_chart(overlay, deltas, *, title="") -> dict:
    """Every group as a change from the baseline, in the result type's units.

    Always a linear axis: a change is signed, and half of these are negative.
    Colours come from the overlay's full column order, so a group is the same
    colour here as it is above - taking them off this frame instead would
    recolour everything the moment the baseline changed.
    """
    frame = getattr(deltas, "frame", None)
    if frame is None or frame.empty:
        return {"series": []}

    aeps = list(frame.index)
    colours = colour_for(overlay.frame.columns)
    unit = "m" if deltas.kind == ABSOLUTE else "%"

    horizontal = aep_axis(aeps)
    horizontal["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}

    series = [{
        "name": str(column),
        "type": "line",
        "symbolSize": 5,
        "connectNulls": False,
        "itemStyle": {"color": colours.get(column, ENVELOPE_COLOUR)},
        "lineStyle": {"width": 2},
        "emphasis": {"focus": "series"},
        "data": _series_data(aeps, frame[column]),
    } for column in frame.columns]

    series[0]["markLine"] = {
        "silent": True, "symbol": "none",
        "lineStyle": {"color": ENVELOPE_COLOUR, "type": "solid", "opacity": 0.6},
        "label": {"show": False},
        "data": [{"yAxis": 0}],
    }

    return {
        "title": {"text": title, "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "cross"},
                    ":valueFormatter": "(v) => v == null ? '-' : Number(v).toPrecision(4)"},
        "legend": {"type": "scroll", "top": 24},
        "grid": {"left": 60, "right": 30, "top": 60, "bottom": 60},
        "xAxis": horizontal,
        "yAxis": {"type": "value", "name": f"Change from {deltas.baseline} ({unit})",
                  "nameLocation": "middle", "nameGap": 45, "scale": True,
                  "splitLine": {"lineStyle": {"opacity": 0.25}}},
        "series": series,
    }


def critical_overlay_chart(overlay, critical) -> dict:
    """Critical duration against AEP, one line per group.

    Whether a warmer climate or a raised dam moves the critical duration is a
    question the single-group version cannot answer.
    """
    if critical is None or critical.empty:
        return {}

    aeps = list(critical.index)
    colours = colour_for(overlay.frame.columns)

    return {
        "tooltip": {"trigger": "axis"},
        "legend": {"type": "scroll", "top": 0},
        "grid": {"left": 60, "right": 30, "top": 30, "bottom": 50},
        "xAxis": aep_axis(aeps, name_gap=30),
        "yAxis": {"type": "log", "name": "Critical duration (h)",
                  "nameLocation": "middle", "nameGap": 40},
        "series": [{
            "name": str(column),
            "type": "line",
            "step": "middle",
            "symbolSize": 6,
            "connectNulls": False,
            "itemStyle": {"color": colours.get(column, ENVELOPE_COLOUR)},
            "data": _series_data(aeps, critical[column]),
        } for column in critical.columns],
    }
