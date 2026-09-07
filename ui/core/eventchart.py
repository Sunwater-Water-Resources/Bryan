"""The AEP neutrality plot: where each candidate event sits against its loading.

Both axes are standard normal variates labelled as AEPs, so the picture is the
question being asked: rainfall rarity up, flood rarity across. The diagonal is
AEP neutrality - an event on it produced a flood exactly as rare as the rain
that caused it. Distance from the target marker is the rank key, and it is a
distance on the page because both axes are in the same units.

Points off the diagonal are the events to be careful with. Below it the flood
outran the rainfall, which means something else supplied the rarity: a full
lake, a big pre-burst, an embedded burst. Those are flagged in the table and
drawn hollow here.

Plain dicts, no nicegui - the same arrangement as core/resultchart.py, which
this borrows its AEP axis from.
"""

from __future__ import annotations

from .results import normal_variate
from .resultchart import PALETTE, aep_axis, numeric

# Enough of the standard set to label either axis over any range Bryan runs.
STANDARD_AEPS = (2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10_000,
                 20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000,
                 2_000_000, 10_000_000)

CLEAN = PALETTE[0]
FLAGGED = PALETTE[3]
CHOSEN = PALETTE[2]
NEUTRAL_LINE = "#888888"


def _ticks(values, pad=0.4) -> tuple:
    """The standard AEPs that fall inside the range of some variates."""
    known = [value for value in values if numeric(value) is not None]
    if not known:
        return STANDARD_AEPS[:6]
    low, high = min(known) - pad, max(known) + pad
    inside = [aep for aep in STANDARD_AEPS
              if low <= normal_variate(aep) <= high]
    return tuple(inside) if len(inside) >= 2 else STANDARD_AEPS[:6]


def _point(row, sim_id, colour, hollow=False) -> dict:
    flags = "; ".join(row.get("flags") or ()) or "nothing flagged"
    style = {"color": "transparent", "borderColor": colour, "borderWidth": 2} \
        if hollow else {"color": colour}
    return {
        "value": [numeric(row.get("z_result")), numeric(row.get("z_rain"))],
        "name": f"sim {int(sim_id)}",
        "itemStyle": style,
        # carried through to the tooltip
        "sim": int(sim_id),
        "rain": numeric(row.get("rain_aep")),
        "result": numeric(row.get("result_aep")),
        "delta": numeric(row.get("delta_z")),
        "flags": flags,
    }


_TOOLTIP = """
(p) => {
  const d = p.data || {};
  const aep = (v) => v == null ? '-' : (v >= 1000 ? Math.round(v).toLocaleString()
                                                  : Number(v).toPrecision(3));
  return `<b>sim ${d.sim}</b><br/>rainfall 1 in ${aep(d.rain)}`
       + `<br/>result 1 in ${aep(d.result)}`
       + `<br/>&Delta;z ${d.delta == null ? '-' : d.delta.toFixed(3)}`
       + `<br/><i>${d.flags}</i>`;
}
"""


def neutrality_chart(outcome, result_type="level", title="") -> dict:
    """Candidates against one loading, with the target and the neutral line."""
    frame = outcome.candidates
    if frame is None or frame.empty:
        return {"series": []}

    clean, flagged = [], []
    picked_id = outcome.picked_id
    for sim_id, row in frame.iterrows():
        if sim_id == picked_id:
            continue
        target_list = flagged if row.get("flags") else clean
        target_list.append(_point(row, sim_id,
                                  FLAGGED if row.get("flags") else CLEAN,
                                  hollow=bool(row.get("flags"))))

    variates = list(frame["z_result"].dropna()) + list(frame["z_rain"].dropna())
    z_target = normal_variate(outcome.aep) if outcome.aep else None
    z_rain_target = normal_variate(
        outcome.target.rain_aep or outcome.aep) if outcome.aep else None
    # Where the loading sits in this database's own realisations, which is not
    # where the design curve puts it: the curve is the envelope over the
    # durations, read as a straight line between the standard AEPs and smoothed
    # on the way. Marking only the design AEP makes that gap look like the plot
    # is wrong, so both are drawn.
    z_data = numeric(getattr(outcome, "data_z", None))
    for value in (z_target, z_rain_target, z_data):
        if numeric(value) is not None:
            variates.append(value)

    known = [value for value in variates if numeric(value) is not None]
    horizontal = aep_axis(_ticks(known))
    # aep_axis pins the range to its own ticks, and those sit *inside* the data
    # by construction here - so widen it back out, or the extreme candidates,
    # which are exactly the ones worth seeing, fall off the plot.
    if known:
        horizontal["min"] = round(min(known) - 0.25, 4)
        horizontal["max"] = round(max(known) + 0.25, 4)
    horizontal["name"] = f"AEP of the {result_type} (1 in X)"
    horizontal["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}
    vertical = dict(horizontal)
    vertical["name"] = "AEP of the rainfall (1 in X)"
    vertical["nameGap"] = 45

    series = [{
        "name": "AEP neutral",
        "type": "line",
        "symbol": "none",
        "silent": True,
        "z": 1,
        "lineStyle": {"color": NEUTRAL_LINE, "type": "dashed", "width": 1.5},
        "data": [[horizontal["min"], horizontal["min"]],
                 [horizontal["max"], horizontal["max"]]],
    }]
    if clean:
        series.append({"name": "candidates", "type": "scatter", "symbolSize": 11,
                       "data": clean, "tooltip": {":formatter": _TOOLTIP}})
    if flagged:
        series.append({"name": "flagged", "type": "scatter", "symbolSize": 11,
                       "data": flagged, "tooltip": {":formatter": _TOOLTIP}})

    picked = outcome.picked
    if picked is not None:
        series.append({
            "name": "chosen",
            "type": "scatter",
            "symbolSize": 20,
            "symbol": "diamond",
            "z": 5,
            "data": [_point(picked, picked.name, CHOSEN)],
            "tooltip": {":formatter": _TOOLTIP},
        })

    if z_data is not None:
        series.append({
            "name": f"{outcome.target.label} in this run",
            "type": "line",
            "symbol": "none",
            "silent": True,
            "z": 2,
            "lineStyle": {"color": CHOSEN, "type": "dotted", "width": 1.5},
            "data": [[z_data, vertical["min"]], [z_data, vertical["max"]]],
        })

    if numeric(z_target) is not None:
        series[0]["markLine"] = {
            "silent": True,
            "symbol": "none",
            "label": {"show": False},
            "lineStyle": {"color": NEUTRAL_LINE, "opacity": 0.6},
            "data": [{"xAxis": numeric(z_target)}, {"yAxis": numeric(z_rain_target)}],
        }

    return {
        "title": {"text": title, "left": "center", "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "item"},
        "legend": {"type": "scroll", "top": 24},
        "grid": {"left": 70, "right": 30, "top": 60, "bottom": 60},
        "xAxis": horizontal,
        "yAxis": vertical,
        "series": series,
    }


# -- the hydrograph preview --------------------------------------------------

HYDROGRAPH_COLOURS = {"inflows": PALETTE[0], "outflows": PALETTE[3],
                      "levels": PALETTE[2]}
HYDROGRAPH_LABELS = {"inflows": "Inflow", "outflows": "Outflow",
                     "levels": "Lake level"}


def hydrograph_chart(series, sim_id, title="") -> dict:
    """One realisation's stored hydrographs - the shape behind the peak.

    Flows on the left axis and the lake level on the right, because the two are
    orders of magnitude apart and the question being asked - is this one rise
    and one recession, or two - is about the shape of each, not their ratio.
    """
    from . import hydrographs                       # local: avoids a cycle

    drawn = []
    for kind in ("inflows", "outflows", "levels"):
        found = series.get(kind)
        if found is None or not len(found):
            continue
        times, values = hydrographs.for_plot(found)
        drawn.append({
            "name": HYDROGRAPH_LABELS[kind],
            "type": "line",
            "showSymbol": False,
            "smooth": False,
            "yAxisIndex": 1 if kind == "levels" else 0,
            "lineStyle": {"width": 1.8,
                          "type": "dashed" if kind == "levels" else "solid"},
            "itemStyle": {"color": HYDROGRAPH_COLOURS[kind]},
            "data": [[time, value] for time, value in zip(times, values)],
        })
    if not drawn:
        return {"series": []}

    return {
        "title": {"text": title or f"sim {int(sim_id)}", "left": "center",
                  "textStyle": {"fontSize": 13}},
        "tooltip": {"trigger": "axis"},
        "legend": {"type": "scroll", "top": 24},
        "grid": {"left": 70, "right": 70, "top": 60, "bottom": 50},
        "xAxis": {"type": "value", "name": "Time (hours)", "nameLocation": "middle",
                  "nameGap": 30, "splitLine": {"show": True,
                                               "lineStyle": {"opacity": 0.25}}},
        "yAxis": [
            {"type": "value", "name": "Flow (m\u00b3/s)", "nameGap": 45,
             "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}}},
            {"type": "value", "name": "Level (m AHD)", "nameGap": 45,
             "splitLine": {"show": False}, "scale": True},
        ],
        "series": drawn,
    }
