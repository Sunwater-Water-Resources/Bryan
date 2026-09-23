"""Render the Report page against a miniature study.

The core is tested without a browser; this is for what those tests cannot see:
the page raising on load, a table never reaching its preview, or a copy button
putting something other than the table on the clipboard.
"""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui.testing.user_simulation import user_simulation      # noqa: E402

from report_fixtures import GROUP, build_study                   # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        import pages
        pages.register_all()
        yield simulated


@pytest.fixture
def private_settings(monkeypatch, tmp_path):
    import settings as settings_module
    monkeypatch.setattr(settings_module, "SETTINGS_PATH", tmp_path / "ui.json")


@pytest.fixture
def opened(tmp_path, private_settings):
    from core import reporttables as rt, study as studies
    from state import STATE

    studies.forget_runs()
    rt.forget_cached()
    study = build_study(tmp_path / "study")
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Table 26 near-term", source={"run": "E099 RFSL", "group": GROUP},
                aeps=[5, 10], pmp_aep=None)
    study.put_table(spec)
    study.save()
    STATE.open_study(study.path)
    yield study
    STATE.study = None


async def preview_text(user, table_id, seconds=10.0):
    """The preview's HTML once it is drawn - tables are built off the event loop.

    ``user.find`` raises rather than returning nothing, so a miss is caught.
    """
    for _ in range(int(seconds / 0.1)):
        try:
            found = user.find(marker=f"preview-{table_id}").elements
        except AssertionError:
            found = set()
        if found:
            return next(iter(found)).content
        await asyncio.sleep(0.1)
    raise AssertionError(f"the preview of {table_id} never appeared")


@pytest.mark.asyncio
async def test_without_a_study_the_page_offers_to_open_one(user, private_settings):
    from state import STATE
    STATE.study = None
    await user.open("/report")
    await user.should_see("Study file")
    await user.should_see("New study")


@pytest.mark.asyncio
async def test_the_table_is_built_and_previewed(user, opened):
    await user.open("/report")
    await user.should_see("Table 26 near-term")
    html = await preview_text(user, "table-26-near-term")
    assert "216.00" in html and "950" in html


@pytest.mark.asyncio
async def test_copy_as_text_puts_the_table_on_the_clipboard(user, opened, monkeypatch):
    from nicegui import ui as nicegui_ui
    copied = []
    monkeypatch.setattr(nicegui_ui.clipboard, "write", copied.append)

    await user.open("/report")
    await preview_text(user, "table-26-near-term")
    user.find(marker="copy-text-table-26-near-term").click()
    await asyncio.sleep(0.2)
    assert copied and copied[0].splitlines()[2] == "10\t950\t150\t216.00\t72"


@pytest.mark.asyncio
async def test_copy_for_word_sends_the_html_and_says_what_it_managed(user, opened,
                                                                     monkeypatch):
    import pages.report as report_page
    scripts = []

    async def fake_run_javascript(code, timeout=1.0):
        scripts.append(code)
        return "html"

    monkeypatch.setattr(report_page.ui, "run_javascript", fake_run_javascript)
    await user.open("/report")
    await preview_text(user, "table-26-near-term")
    user.find(marker="copy-word-table-26-near-term").click()
    await user.should_see("Copied - paste into Word")
    assert "ClipboardItem" in scripts[0] and "216.00" in scripts[0]


@pytest.mark.asyncio
async def test_a_study_is_opened_from_the_page(user, private_settings, tmp_path):
    from state import STATE
    STATE.study = None
    study = build_study(tmp_path / "other")
    await user.open("/report")
    user.find(marker="study-path").clear().type(str(study.path))
    user.find(marker="open-study").click()
    await user.should_see(str(study.path))
    await user.should_see("1 group")
    assert STATE.study is not None and STATE.study.path == study.path
    STATE.study = None
