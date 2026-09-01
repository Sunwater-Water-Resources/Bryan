"""Drive REAL Bryan through the whole UI path.

Reservoir routing needs no model executable, so this runs on Linux in seconds
and is the only test that proves the UI's run copies and configs are actually
acceptable to Bryan - not just to a stand-in that agrees with the UI's
assumptions.

Skipped when Bryan's own dependencies (scipy, matplotlib) are not importable,
because the UI environment deliberately does not have them.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

from conftest import write_workbook
from core import progress, runplan, runstate, runwriter
from core.bryan import BRYAN_ROOT
from core.config import load_sims_config
from core.launcher import RunManager
from core.simslist import read_sims_list
from make_rr_fixture import (SIMS_COLUMNS, build, monte_carlo_rows,
                             sims_rows)


def _bryan_python():
    """An interpreter that can import what Main.py needs, or None.

    ``BRYAN_TEST_PYTHON`` overrides the search - point it at the interpreter
    the model batch files use (``%PROJECT_ROOT%\\env\\python.exe`` on Windows)
    to run these tests against the environment that produces study results.
    """
    candidates = []
    override = os.environ.get("BRYAN_TEST_PYTHON")
    if override:
        candidates.append(Path(override))
    candidates += [BRYAN_ROOT / "env/bin/python",
                   BRYAN_ROOT / "env/Scripts/python.exe",
                   BRYAN_ROOT / ".venv/bin/python",
                   Path(sys.executable)]

    for candidate in candidates:
        if not Path(candidate).exists():
            continue
        probe = subprocess.run(
            [str(candidate), "-c",
             "import numpy, pandas, scipy, matplotlib, openpyxl"],
            capture_output=True)
        if probe.returncode == 0:
            return str(candidate)
    return None


BRYAN_PYTHON = _bryan_python()
needs_bryan = pytest.mark.skipif(
    BRYAN_PYTHON is None,
    reason="no interpreter with Bryan's dependencies (scipy, matplotlib)")


@pytest.fixture
def mini_project(tmp_path):
    """The miniature reservoir-routing model, with its sims list."""
    folder = tmp_path / "mini"
    folder.mkdir()
    config_path = build(folder)
    write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS, sims_rows())
    return config_path


@pytest.fixture
def mini_mc_project(tmp_path):
    """The same model, driven from a Monte Carlo database instead.

    Two rows over one ``Output file``, told apart only by ``Output suffix`` -
    the shape a re-routing sims list actually has.
    """
    folder = tmp_path / "mini_mc"
    folder.mkdir()
    config_path = build(folder)
    write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS, monte_carlo_rows())
    return config_path


def wait_until(predicate, timeout=180.0, interval=0.2):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(interval)
    return False


def test_the_fixture_matches_the_formats_bryan_parses(mini_project):
    """Cheap, and runs even without Bryan's dependencies."""
    sys.path.insert(0, str(BRYAN_ROOT))
    folder = mini_project.parent

    els = pd.read_csv(folder / "reservoir/dam.els")
    assert {"EL", "V"} <= set(els.columns)

    # _read_sq skips six lines then reads 'storage flow' pairs.
    lines = (folder / "reservoir/dam_base.sq").read_text().splitlines()
    assert len(lines) > 6
    for line in lines[6:]:
        if line.strip():
            storage, flow = line.split()
            float(storage), float(flow)

    database = pd.read_csv(folder / "inputs/ensemble_results.csv", index_col=0)
    assert {"duration", "tp"} <= set(database.columns), \
        "an ensemble database is detected from these columns"

    inflows = pd.read_csv(folder / "inputs/ensemble_inflows.csv", index_col=0)
    assert list(inflows.columns) == list(database.index)
    assert inflows.isna().any().any(), "an ensemble inflow file is ragged"


@needs_bryan
def test_real_bryan_runs_a_ui_written_run_folder(mini_project):
    """The whole point: plan -> write -> launch -> Bryan produces results."""
    config = load_sims_config(mini_project)
    sims = read_sims_list(config.sims_list_path)

    plan = runplan.plan_run(sims, list(sims.frame.index), n_chunks=1)
    assert not plan.blocked, [h.message for h in plan.hazards if h.blocks]

    folder = runwriter.write_run(sims, config, plan, run_id="20260805-real")
    record = runstate.from_run_folder(folder, config, max_parallel=1)

    manager = RunManager(bryan_python=BRYAN_PYTHON,
                         bryan_main=str(BRYAN_ROOT / "Main.py"))
    manager.submit(record)
    assert wait_until(lambda: (manager.poll(), not record.is_live)[1]), \
        "Bryan did not finish"

    chunk = record.chunks[0]
    console = Path(chunk.console_log).read_text(encoding="utf-8", errors="replace")
    assert chunk.returncode == 0, f"Bryan exited {chunk.returncode}\n{console[-4000:]}"

    # Bryan wrote its run log where the UI predicted.
    assert Path(chunk.run_log).is_file()
    log = pd.read_csv(chunk.run_log)
    assert list(log["Status"]) == ["completed", "completed"]

    # And the results themselves, one set per rating curve.
    results = config.project_folder / "results"
    for suffix in ("base", "raised"):
        assert (results / f"mini_enb_{suffix}.csv").is_file(), \
            f"no ensemble database for {suffix}\n{console[-4000:]}"

    # The UI's progress parser understands what Bryan actually printed.
    parsed = progress.read_chunk_progress(
        _chunk_shim(chunk), returncode=chunk.returncode, alive=False)
    assert parsed.completed == 2
    assert "All 2 simulations completed" in parsed.summary


@needs_bryan
def test_completion_sees_the_results_bryan_wrote(mini_project):
    """Closes the loop: the states on the Select page come from real outputs."""
    from core import completion

    config = load_sims_config(mini_project)
    sims = read_sims_list(config.sims_list_path)

    before = completion.assess_frame(sims.frame, config)
    assert {state.state for state in before.values()} == {completion.NOT_RUN}

    plan = runplan.plan_run(sims, list(sims.frame.index), n_chunks=1)
    folder = runwriter.write_run(sims, config, plan, run_id="20260805-real2")
    record = runstate.from_run_folder(folder, config, max_parallel=1)
    manager = RunManager(bryan_python=BRYAN_PYTHON,
                         bryan_main=str(BRYAN_ROOT / "Main.py"))
    manager.submit(record)
    assert wait_until(lambda: (manager.poll(), not record.is_live)[1])
    assert record.chunks[0].returncode == 0

    after = completion.assess_frame(sims.frame, config)
    assert {state.state for state in after.values()} == {completion.UP_TO_DATE}


class _chunk_shim:
    def __init__(self, record) -> None:
        self.index = record.index
        self.rows = tuple(record.rows)
        self.run_log = Path(record.run_log)
        self.console_log = Path(record.console_log)


# --- the Output suffix on the Monte Carlo database -------------------------

def _run(config_path, run_id):
    """Plan, write and run the whole sims list. Returns (config, sims, record)."""
    config = load_sims_config(config_path)
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, list(sims.frame.index), n_chunks=1)
    assert not plan.blocked, [h.message for h in plan.hazards if h.blocks]

    folder = runwriter.write_run(sims, config, plan, run_id=run_id)
    record = runstate.from_run_folder(folder, config, max_parallel=1)
    manager = RunManager(bryan_python=BRYAN_PYTHON,
                         bryan_main=str(BRYAN_ROOT / "Main.py"))
    manager.submit(record)
    assert wait_until(lambda: (manager.poll(), not record.is_live)[1]), \
        "Bryan did not finish"
    console = Path(record.chunks[0].console_log).read_text(encoding="utf-8",
                                                          errors="replace")
    assert record.chunks[0].returncode == 0, console[-4000:]
    return config, sims, console


@needs_bryan
def test_each_rating_curve_keeps_its_own_monte_carlo_database(mini_mc_project):
    """Both rows share an Output file, so before the suffix went on the mcdf the
    second one silently overwrote the first."""
    config, _, console = _run(mini_mc_project, "20260826-mc")
    results = config.project_folder / "results"

    databases = {suffix: results / f"mini_mc__mcdf_{suffix}.csv"
                 for suffix in ("base", "raised")}
    for suffix, path in databases.items():
        assert path.is_file(), f"no database for {suffix}\n{console[-4000:]}"
    assert not (results / "mini_mc__mcdf.csv").exists(), \
        "the unsuffixed name is the one the two rows used to fight over"

    base, raised = (pd.read_csv(path, index_col=0) for path in databases.values())
    assert len(base) == len(raised) == 20
    assert not base["level"].equals(raised["level"]), \
        "two rating curves must not produce the same routed levels"

    # And the quantile tables that were already suffixed still line up with them.
    for suffix in databases:
        for kind in ("inflow", "level", "outflow"):
            assert (results / f"mini_mc__{kind}_quantiles_{suffix}.csv").is_file()


@needs_bryan
def test_an_analysis_only_row_falls_back_to_the_pre_suffix_database(tmp_path):
    """Results routed before 26 August 2026 carry no suffix on the mcdf.

    _ensure_mcdf_loaded still finds one, so a re-analysis of an old study keeps
    working - but it says out loud that the file belongs to every suffix.
    """
    folder = tmp_path / "legacy"
    folder.mkdir()
    config_path = build(folder)

    # Route it once for a real database, then rename it as an old run left it
    # and delete the quantile tables that came with it.
    write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS,
                   monte_carlo_rows(suffixes=("base",)))
    _run(config_path, "20260826-seed")
    results = folder / "results"
    (results / "mini_mc__mcdf_base.csv").rename(results / "mini_mc__mcdf.csv")
    for kind in ("inflow", "level", "outflow"):
        (results / f"mini_mc__{kind}_quantiles_base.csv").unlink()

    rows = monte_carlo_rows(suffixes=("base",))
    rows[0]["Run models"] = "no"
    write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS, rows)
    _run(config_path, "20260826-legacy")

    # The simulation's own output is teed to its log file, not the console.
    log = (results / "mini_mc_base_log.txt").read_text(encoding="utf-8",
                                                       errors="replace")
    assert "WARNING" in log and "pre-suffix name" in log, log[-3000:]
    assert "mini_mc__mcdf.csv" in log
    for kind in ("inflow", "level", "outflow"):
        assert (results / f"mini_mc__{kind}_quantiles_base.csv").is_file(), \
            "the re-analysis still produced its own suffixed quantiles"


# --- a fixed antecedent storage over Monte Carlo input ---------------------

@needs_bryan
def test_a_fixed_adv_holds_every_realisation_at_one_storage(tmp_path):
    """The method exists to test dam operation against a fixed set of inflows,
    and the starting storage is part of the operation being tested.

    'ADV source' of 'sims list' is the only way to it: the ADV column alone is
    not read for Monte Carlo input, because the storages there come one per
    realisation from the input database.
    """
    projects = {}
    for name, rows in (
        ("default", monte_carlo_rows(suffixes=("base",))),
        ("fixed", monte_carlo_rows(suffixes=("base",), adv=2500.0,
                                   adv_source="sims list")),
    ):
        folder = tmp_path / name
        folder.mkdir()
        projects[name] = build(folder)
        write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS, rows)
        _run(projects[name], f"20260901-{name}")

    fixed = pd.read_csv(tmp_path / "fixed/results/mini_mc__mcdf_base.csv", index_col=0)
    default = pd.read_csv(tmp_path / "default/results/mini_mc__mcdf_base.csv", index_col=0)

    # Every realisation started where the sims list said, and the storage the
    # source run used is kept beside it.
    assert (fixed["ADV"] == 2500.0).all()
    assert (fixed["ADV_input"] == 5500.0).all()
    assert (default["ADV"] == 5500.0).all(), "the default still comes from the database"

    # Starting 3000 ML lower is a real difference, not a relabelled column.
    assert (fixed["level"] < default["level"]).all()

    log = (tmp_path / "fixed/results/mini_mc_base_log.txt").read_text(
        encoding="utf-8", errors="replace")
    assert "ADV source: the ADV column of the sims list" in log
    assert "conditional on that starting storage" in log, \
        "holding the lake still stops these being design flood quantiles"


@needs_bryan
def test_an_adv_the_monte_carlo_source_ignores_says_so(tmp_path):
    """An 'ADV' with no 'ADV source' does nothing at all for Monte Carlo input.

    It used to do it silently, which reads as the run having ignored the
    setting rather than never having looked at it.
    """
    folder = tmp_path / "unused"
    folder.mkdir()
    config_path = build(folder)
    write_workbook(folder / "MiniSimsList.xlsx", SIMS_COLUMNS,
                   monte_carlo_rows(suffixes=("base",), adv="fsv"))
    _run(config_path, "20260901-unused")

    mcdf = pd.read_csv(folder / "results/mini_mc__mcdf_base.csv", index_col=0)
    assert (mcdf["ADV"] == 5500.0).all(), "the database's storages, untouched"

    log = (folder / "results/mini_mc_base_log.txt").read_text(
        encoding="utf-8", errors="replace")
    assert 'which is NOT used' in log
    assert 'Set "ADV source" to "sims list"' in log
