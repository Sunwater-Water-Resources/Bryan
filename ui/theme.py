"""The look of every page: fonts, the house colours, and a chart theme.

The colours live in core/palette.py, as plain values the chart builders can use.
This puts them into the page: Quasar's brand colours and a few custom ones, so
``color=muted`` on a chip and ``text-muted`` on a label mean the house grey, and
a stylesheet that names every colour Quasar would otherwise choose for itself.
"""

from __future__ import annotations

from nicegui import ui

from core.palette import (BODY, ECHARTS_THEME, FONT, INK, MISSING, MONO, MUTED, OCHRE,
                          PAPER, READY, RULE, SURFACE, WATER, WATER_SOFT)

_FONTS = ('<link rel="preconnect" href="https://fonts.googleapis.com">'
          '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
          'family=IBM+Plex+Mono:wght@400;500&family=Rubik:wght@300;400;500;600&display=swap">')

_STYLE = f"""<style>
body {{ font-family: {FONT}; background: {PAPER}; color: {BODY}; }}
.mono {{ font-family: {MONO}; font-variant-numeric: tabular-nums; }}
a, .q-link {{ color: {WATER}; }}
.q-card {{ background: {SURFACE}; color: {BODY}; }}
.q-card--bordered {{ border-color: {RULE}; }}
.text-lg.font-bold, .font-bold {{ color: {INK}; }}
.text-lg.font-bold {{ font-weight: 600; }}
.q-field__native, .q-field__prefix, .q-field__suffix, .q-field__input {{ color: {INK}; }}
.q-field__label {{ color: {MUTED}; }}
.q-field--outlined .q-field__control:before {{ border-color: {RULE}; }}
.q-table__container {{ color: {BODY}; }}
.q-table th {{ color: {MUTED}; font-weight: 500; }}
.q-table tbody td {{ color: {INK}; }}
.q-expansion-item__container .q-item__label {{ color: {INK}; }}
.q-tab {{ color: {BODY}; }}
.q-tab--active {{ color: {WATER}; }}
.q-tab__indicator {{ background: {WATER}; }}
.q-toggle__label, .q-checkbox__label {{ color: {INK}; }}
</style>"""


def apply() -> None:
    """Called at the top of every page, before anything is drawn."""
    ui.colors(primary=WATER, secondary=INK, accent=OCHRE, positive=READY,
              negative=MISSING, warning=OCHRE, info=WATER,
              ink=INK, body=BODY, muted=MUTED, attention=OCHRE, water_soft=WATER_SOFT)
    ui.add_head_html(_FONTS + _STYLE)


def house_echart(options: dict, **kwargs):
    """ui.echart in the house style: the palette, the fonts, the axis greys."""
    return ui.echart(options, theme=ECHARTS_THEME, **kwargs)
