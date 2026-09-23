"""What goes on the clipboard: one table, as HTML for Word and as text."""

from __future__ import annotations

from core import wordtable as wt


def sample():
    table = wt.ReportTable(header=["Climate Horizon", "Inflow (m3/s)", "Level"],
                           align=["left", "center", "center"],
                           footnotes=["* Peak inflow for lake level critical duration."])
    table.section("Current RFSL of 215.5 m AHD")
    table.add(["Near <Term> & co", "16,301", "220.89"])
    table.add(["1,900,000\n(PMPF)", "15,247", "220.88"], shaded=True, bold=True)
    return table


def test_the_html_is_a_word_table_in_the_house_style():
    html = wt.to_html(sample())
    assert html.startswith("<table")
    assert "background:#00B0CA;color:#FFFFFF;font-weight:bold" in html
    assert "m<sup>3</sup>/s" in html
    assert 'colspan="3"' in html and "background:#C1F7FF" in html
    assert "font-family:'Rubik Light'" in html and "font-size:10pt" in html
    assert "1,900,000<br>(PMPF)" in html
    assert html.index("</table>") < html.index("* Peak inflow")


def test_cell_text_is_escaped():
    html = wt.to_html(sample())
    assert "Near &lt;Term&gt; &amp; co" in html
    assert "<Term>" not in html


def test_a_study_can_restyle_the_table():
    html = wt.to_html(sample(), {"font": "Arial", "header_fill": "#123456", "unknown": 1})
    assert "font-family:'Arial'" in html and "background:#123456" in html


def test_the_text_is_tab_separated_with_plain_units():
    lines = wt.to_text(sample()).split("\r\n")
    assert lines[0] == "Climate Horizon\tInflow (m3/s)\tLevel"
    assert lines[1] == "Current RFSL of 215.5 m AHD"
    assert lines[3] == "1,900,000 (PMPF)\t15,247\t220.88"
    assert lines[4].startswith("* Peak inflow")


def test_the_clipboard_script_carries_both_forms_and_falls_back_to_text():
    script = wt.clipboard_script(wt.clipboard_document("<table></table>"), "a\tb")
    assert "ClipboardItem" in script and "'text/html'" in script
    assert "writeText(text)" in script
    assert '"a\\tb"' in script                      # JSON-escaped, not spliced in raw
