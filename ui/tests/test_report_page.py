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

    ``user.find`` raises rather than returning nothing, so a miss is caught. Tables
    open folded, so a folded one is opened first.
    """
    chevron = next(iter(user.find(marker=f"chevron-{table_id}").elements))
    if chevron.name == "chevron_right":
        user.find(marker=f"fold-{table_id}").click()
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
async def test_without_a_study_the_page_sends_you_to_the_study_page(user, private_settings):
    from state import STATE
    STATE.study = None
    await user.open("/report")
    await user.should_see("No study is open")
    await user.should_see("Go to Study")


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
async def test_the_study_is_named_with_a_way_to_change_it(user, opened):
    await user.open("/report")
    await user.should_see(marker="report-study")
    await user.should_see("Change study")


# -- tables that fold away -------------------------------------------------------

@pytest.fixture
def two_tables(opened):
    from core import reporttables as rt
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Table 27 long-term", source={"run": "E099 RFSL", "group": GROUP},
                aeps=[5, 10], pmp_aep=None)
    opened.put_table(spec)
    opened.save()
    from state import STATE
    STATE.open_study(opened.path)
    return [table["id"] for table in opened.tables]


@pytest.fixture
def counted_builds(monkeypatch):
    from core import reporttables as rt
    built = []
    real = rt.build

    def counting(study, spec):
        built.append(spec["id"])
        return real(study, spec)

    monkeypatch.setattr(rt, "build", counting)
    return built


def _chevron(user, table_id):
    return next(iter(user.find(marker=f"chevron-{table_id}").elements)).name


@pytest.mark.asyncio
async def test_tables_open_folded_and_named_with_nothing_built(user, two_tables,
                                                               counted_builds):
    await user.open("/report")
    await user.should_see(marker="table-contents")
    await user.should_see("1. Table 26 near-term")
    await user.should_see("2. Table 27 long-term")
    for table_id in two_tables:
        assert _chevron(user, table_id) == "chevron_right"
    await user.should_see("not built yet")
    await asyncio.sleep(0.3)
    assert counted_builds == []


@pytest.mark.asyncio
async def test_opening_one_builds_only_that_one_and_is_remembered(user, two_tables,
                                                                  counted_builds):
    from state import STATE
    first, second = two_tables
    await user.open("/report")
    await preview_text(user, second)
    assert counted_builds == [second]
    await user.should_see(marker=f"copy-word-{second}")
    await user.should_not_see(marker=f"copy-word-{first}")
    assert STATE.settings.open_tables_for(STATE.study.path) == {second}

    await user.open("/report")                  # it stays open
    assert _chevron(user, second) == "expand_more"
    assert _chevron(user, first) == "chevron_right"


@pytest.mark.asyncio
async def test_open_all_and_fold_all(user, two_tables, counted_builds):
    from state import STATE
    await user.open("/report")
    user.find(marker="open-all-tables").click()
    for table_id in two_tables:
        await preview_text(user, table_id)
    assert sorted(counted_builds) == sorted(two_tables)
    user.find(marker="fold-all-tables").click()
    assert all(_chevron(user, table_id) == "chevron_right" for table_id in two_tables)
    assert STATE.settings.open_tables_for(STATE.study.path) == set()


@pytest.mark.asyncio
async def test_the_contents_open_the_table_they_name(user, two_tables):
    second = two_tables[1]
    await user.open("/report")
    user.find(marker=f"contents-{second}").click()
    assert _chevron(user, second) == "expand_more"
    await preview_text(user, second)


@pytest.mark.asyncio
async def test_a_table_with_problems_says_so_on_its_folded_line(user, opened):
    from core import reporttables as rt
    spec = rt.new_spec(rt.DESIGN_FLOODS)
    spec.update(title="Broken", source={"run": "E099 RFSL", "group": "no such group"},
                aeps=[5], pmp_aep=None)
    stored = opened.put_table(spec)
    opened.save()
    from state import STATE
    STATE.open_study(opened.path)
    await user.open("/report")
    await preview_text(user, stored["id"])
    user.find(marker=f"fold-{stored['id']}").click()
    await user.should_see("to read")


@pytest.mark.asyncio
async def test_a_table_from_stale_results_says_so_and_again_when_copied(user, opened,
                                                                       monkeypatch):
    from nicegui import ui as nicegui_ui
    from test_staleness import make_stale
    copied = []
    monkeypatch.setattr(nicegui_ui.clipboard, "write", copied.append)
    make_stale(opened)
    table_id = "table-26-near-term"
    await user.open("/report")
    await preview_text(user, table_id)
    await user.should_see("dam.sq changed after the run")
    await user.should_see("results out of date")
    user.find(marker=f"copy-text-{table_id}").click()
    await user.should_see("Copied as text - but")
    assert copied
