"""The manual served by the launcher, and each page's way into it."""

from __future__ import annotations

import re

import pytest

from core import manual


def test_anchors_are_written_as_the_manual_links_them():
    assert manual.slug("Dam inputs") == "daminputs"
    assert manual.slug("What is checked before a run") == "whatischeckedbeforearun"


def test_every_link_the_manual_makes_to_itself_lands_on_a_heading():
    body, headings = manual.prepare(manual.doc_path("ui").read_text(encoding="utf-8"))
    anchors = {anchor for _, _, anchor in headings}
    for target in re.findall(r"\]\(#([^)]+)\)", body):
        assert target in anchors, f"#{target} names no heading"


def test_every_page_section_is_a_heading_of_its_doc():
    for title, (doc, heading) in manual.SECTIONS.items():
        path = manual.doc_path(doc)
        assert path is not None, f"{title}: no {doc}.md"
        if heading:
            _, headings = manual.prepare(path.read_text(encoding="utf-8"))
            assert manual.slug(heading) in {anchor for _, _, anchor in headings}, \
                f"{title}: no heading '{heading}' in {doc}.md"


def test_headings_get_anchors_and_code_is_left_alone():
    body, headings = manual.prepare("# Top\n\n```\n# not a heading\n```\n\n## Dam inputs\n"
                                    "See [the list](sim_list.md#columns).")
    assert headings == [(1, "Top", "top"), (2, "Dam inputs", "daminputs")]
    assert '<a id="daminputs"></a>' in body
    assert "# not a heading" in body and body.count("<a id=") == 2
    assert "](/manual/sim_list#columns)" in body


def test_only_the_manual_s_own_docs_are_served():
    assert manual.doc_path("ui") is not None
    for name in ("../README", "..\\ui", "", "no_such_doc", "ui.md"):
        assert manual.doc_path(name) is None


def test_a_page_is_sent_to_its_section():
    assert manual.section_for("Downstream storms") == "/manual/downstream_storms"
    assert manual.section_for("Results") == "/manual/ui#viewingresults"
    assert manual.section_for("History") is None


pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.mark.asyncio
async def test_the_manual_page_shows_the_doc(user):
    await user.open("/manual/ui")
    await user.should_see(marker="manual-contents")
    body = next(iter(user.find(marker="manual-body").elements))
    assert '<a id="thestudy"></a>' in body.content
    await user.open("/manual/nothing_here")
    await user.should_see("There is no manual page called 'nothing_here'.")


@pytest.mark.asyncio
async def test_pages_with_a_section_have_a_help_button(user):
    await user.open("/study")
    await user.should_see(marker="open-help")
    await user.open("/history")
    await user.should_not_see(marker="open-help")
