"""The house style: Judith's palette, measured, and used everywhere.

Bryan and Judith are read side by side, so the window takes Judith's light
palette and does not follow the desktop theme. These pin the parts of that which
could drift back without anyone noticing: a Tailwind grey typed into a new page,
a chart drawn without the theme, a text colour that no longer passes AA.
"""

from __future__ import annotations

import re
from pathlib import Path

from core import palette

UI_ROOT = Path(__file__).resolve().parents[1]
PAGES = [UI_ROOT / "layout.py", *sorted((UI_ROOT / "pages").glob("*.py"))]


def contrast(one: str, other: str) -> float:
    def luminance(colour):
        red, green, blue = (int(colour.lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4))
        linear = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
                  for c in (red, green, blue)]
        return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
    high, low = sorted((luminance(one), luminance(other)), reverse=True)
    return (high + 0.05) / (low + 0.05)


def test_every_text_colour_passes_aa_on_the_ground_it_sits_on():
    text = {"INK": palette.INK, "BODY": palette.BODY, "MUTED": palette.MUTED,
            "WATER": palette.WATER, "OCHRE": palette.OCHRE,
            "READY": palette.READY, "MISSING": palette.MISSING}
    short = {name: round(contrast(colour, ground), 2)
             for name, colour in text.items()
             for ground in (palette.SURFACE, palette.PAPER)
             if contrast(colour, ground) < 4.5}
    assert not short, f"under 4.5:1: {short}"
    assert contrast(palette.ON_INK_MUTED, palette.INK) >= 4.5


def test_white_text_reads_on_every_status_chip_colour():
    for colour in (palette.MUTED, palette.BODY, palette.READY, palette.MISSING,
                   palette.OCHRE, palette.WATER):
        assert contrast("#FFFFFF", colour) >= 4.5, colour


def test_the_chart_ramp_is_ten_distinct_marks_that_stand_off_white():
    assert len(palette.PALETTE) == len(set(palette.PALETTE)) == 10
    assert all(contrast(colour, palette.SURFACE) >= 3.0 for colour in palette.PALETTE)
    # the attention colour means a problem, never a series
    assert palette.OCHRE not in palette.PALETTE


def test_no_page_names_a_tailwind_grey_or_a_quasar_grey():
    stray = {path.name: sorted(set(re.findall(r"text-gray-\d+|bg-gray-\d+|grey-\d+",
                                              path.read_text(encoding="utf-8"))))
             for path in PAGES}
    assert not {name: found for name, found in stray.items() if found}


def test_every_chart_is_drawn_with_the_theme():
    bare = [path.name for path in PAGES if "ui.echart(" in path.read_text(encoding="utf-8")]
    assert not bare, f"use theme.house_echart, not ui.echart: {bare}"


def test_the_window_does_not_follow_the_desktop_theme():
    source = (UI_ROOT / "main.py").read_text(encoding="utf-8")
    assert re.search(r"dark=False", source) and not re.search(r"dark=None", source)
