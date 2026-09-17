"""The house colours: Sunwater's, as Judith's window uses them.

Bryan and Judith chain in practice - Bryan produces the inflow, Judith routes it
through the reservoir while the structures fail - so they are read side by side
and are drawn as siblings. These values are copied from Judith's window
(``damfailure/ui/app.py`` in dam-failure-hydraulics) rather than shared, because
the two are separate repositories with separate environments. Change one, change
both.

Plain values, no nicegui, so the chart builders in core/ can use them.

Measured, not eyeballed (WCAG 2.1 contrast):

    INK      17.06:1 on white   headings, values
    BODY      6.44:1            body text
    MUTED     5.03:1, 4.67:1 on PAPER
                                help text. Judith's first muted grey, #74858F,
                                reached only 3.82:1 - under AA for text this size -
                                and Bryan's own text-gray-500 was 4.83:1, so taking
                                it unchanged would have made Bryan's prose harder to
                                read. Both now use this one.
    WATER     5.05:1            links, primary actions, the brand cyan darkened
    OCHRE     5.32:1            attention: problems and warnings, never a series
    READY     6.58:1, MISSING 6.54:1
                                done and not done

The window is pinned to this light palette (``dark=False``); it does not follow
the desktop theme.
"""

from __future__ import annotations

INK, BODY, MUTED = "#031E2F", "#4F616D", "#63717A"
WATER, OCHRE = "#007A8C", "#A8541B"
PAPER, SURFACE, RULE, HAIR = "#F3F7F8", "#FFFFFF", "#CBD7DC", "#E1E8EB"
WATER_SOFT = "#DCF0F3"
READY, MISSING = "#1B6B2E", "#B3261E"

# Text on the navy header. 7.85:1 on INK.
ON_INK_MUTED = "#9FB3BD"
# The brand cyan, which has no contrast requirement as a mark on navy.
BRAND_CYAN = "#00B0CA"

FONT = "Rubik, system-ui, -apple-system, sans-serif"
MONO = '"IBM Plex Mono", ui-monospace, Menlo, monospace'

# A categorical ramp built around the brand cyan, with the attention ochre kept out
# of it so a series never reads as a warning. Ten slots, because a comparison of
# storm durations routinely draws that many. Validated as a set for light mode on
# white: every slot inside the lightness band and at least 3:1 on the surface,
# worst adjacent colour-blind separation dE 9.0, worst adjacent normal-vision
# separation dE 25.1, and no two slots anywhere in the set closer than dE 11. The
# first five are the order Judith always had; the last five replace three of its
# slots that could not be told apart under protanopia.
PALETTE = ("#0089A0", "#C42E6A", "#1B6FD1", "#B8791A", "#6E4BC4",
           "#117733", "#882255", "#6B8E23", "#0F5C8C", "#B04A00")


def _axis() -> dict:
    return {
        "axisLine": {"lineStyle": {"color": RULE}},
        "axisTick": {"lineStyle": {"color": RULE}},
        "axisLabel": {"color": MUTED},
        "splitLine": {"lineStyle": {"color": HAIR}},
        "nameTextStyle": {"color": MUTED},
    }


# An ECharts theme, handed to every chart the window draws. A chart's own options
# still win wherever they name a colour.
ECHARTS_THEME = {
    "color": list(PALETTE),
    "backgroundColor": "transparent",
    "textStyle": {"fontFamily": FONT, "color": BODY},
    "title": {"textStyle": {"color": INK, "fontWeight": 500},
              "subtextStyle": {"color": MUTED}},
    "legend": {"textStyle": {"color": BODY}},
    "tooltip": {"backgroundColor": SURFACE, "borderColor": RULE,
                "textStyle": {"color": INK}},
    "valueAxis": _axis(),
    "categoryAxis": _axis(),
    "logAxis": _axis(),
    "timeAxis": _axis(),
}
