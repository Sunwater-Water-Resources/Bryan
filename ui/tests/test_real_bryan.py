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
from make_rr_fixture import SIMS_COLUMNS, build, sims_rows


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
