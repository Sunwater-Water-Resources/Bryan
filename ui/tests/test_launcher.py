"""End to end against a stand-in Bryan - no URBS, no RORB, runs on Linux.

``fake_main.py`` takes the same command line as Main.py and imports the real
lib.RunLog and lib.LogFiles, so the run-log format cannot drift away from what
core/progress.py parses.
"""

from __future__ import annotations

import os
import sys
import time

import pytest

from conftest import TINAROO_COLUMNS, reservoir_row
from core import launcher, progress, runplan, runstate, runwriter
from core.config import load_sims_config
from core.simslist import read_sims_list

def wait_until(predicate, timeout=60.0, interval=0.05):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(interval)
    return False


def make_run(project, fake_main, rows, *, n_chunks=1, max_parallel=1,
             columns=TINAROO_COLUMNS):
    config = load_sims_config(project(columns, rows))
    sims = read_sims_list(config.sims_list_path)
    plan = runplan.plan_run(sims, list(sims.frame.index), n_chunks=n_chunks,
                            keep_groups_together=False)
    assert not plan.blocked, [h.message for h in plan.hazards if h.blocks]

    folder = runwriter.write_run(sims, config, plan, run_id="20260805-000000")
    names = {index: str(sims.frame.loc[index, "Output file"])
             for index in sims.frame.index}
    record = runstate.from_run_folder(folder, config, max_parallel=max_parallel,
                                      names=names)
    manager = launcher.RunManager(bryan_python=sys.executable,
                                  bryan_main=str(fake_main))
    return config, sims, folder, record, manager


def drain(manager, record, timeout=90.0):
    def done():
        manager.poll()
        return not record.is_live
    assert wait_until(done, timeout), (
        f"run never finished: {[(c.index, c.status) for c in record.chunks]}")


def test_a_single_chunk_runs_to_completion(project, fake_main):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}",
                                         "Log file": rf"sims_mc\log\r{d}.log"})
            for d in (18, 24)]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)

    manager.submit(record)
    drain(manager, record)

    chunk = record.chunks[0]
    assert chunk.status == runstate.COMPLETED
    assert chunk.returncode == 0

    result = progress.read_chunk_progress(
        folder.chunks[0], names=[r["Output file"] for r in rows],
        returncode=chunk.returncode, alive=False)
    assert result.completed == 2
    assert result.total == 2
    assert "All 2 simulations completed" in result.summary


def test_the_run_log_lands_in_the_run_folder_and_bryan_writes_it(project, fake_main):
    """Proves the per-chunk run log arrangement works against real RunLog."""
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)
    manager.submit(record)
    drain(manager, record)

    run_log = folder.chunks[0].run_log
    assert run_log.is_file(), f"expected a run log at {run_log}"
    assert "completed" in run_log.read_text(encoding="utf-8")


def test_a_failing_simulation_does_not_stop_the_rest(project, fake_main):
    """Main.py:97-108 catches per row; the chunk exits 1 at the end."""
    rows = [
        reservoir_row(**{"Output file": "ok_1", "Output suffix": "a"}),
        reservoir_row(**{"Output file": "bad", "Output suffix": "b",
                         "Fake result": "fail"}),
        reservoir_row(**{"Output file": "ok_2", "Output suffix": "c"}),
    ]
    columns = TINAROO_COLUMNS + ["Fake result"]
    config, sims, folder, record, manager = make_run(project, fake_main, rows,
                                                     columns=columns)
    manager.submit(record)
    drain(manager, record)

    chunk = record.chunks[0]
    assert chunk.returncode == 1
    result = progress.read_chunk_progress(folder.chunks[0], alive=False,
                                          returncode=chunk.returncode)
    assert [sim.status for sim in result.sims] == ["completed", "FAILED", "completed"]
    assert "at least one simulation failed" in progress.explain_returncode(1)


def test_a_missing_input_is_flagged_as_such(project, fake_main):
    rows = [reservoir_row(**{"Output file": "gone", "Fake result": "missing input"})]
    config, sims, folder, record, manager = make_run(
        project, fake_main, rows, columns=TINAROO_COLUMNS + ["Fake result"])
    manager.submit(record)
    drain(manager, record)

    result = progress.read_chunk_progress(folder.chunks[0], alive=False)
    assert result.sims[0].status == "FAILED - missing input"


def test_input_raises_eoferror_instead_of_hanging(project, fake_main):
    """stdin=DEVNULL. Without it this test would hang until the timeout.

    Bryan calls input() at MCScheme/EnbScheme.store_simulations and
    lib/URBSmodel.py:41.
    """
    rows = [reservoir_row(**{"Output file": "asks", "Fake result": "ask"})]
    config, sims, folder, record, manager = make_run(
        project, fake_main, rows, columns=TINAROO_COLUMNS + ["Fake result"])
    manager.submit(record)
    drain(manager, record, timeout=30)

    result = progress.read_chunk_progress(folder.chunks[0], alive=False)
    assert result.sims[0].failed
    assert "EOFError" in result.sims[0].error
    assert "open in Excel" in progress.explain_error(result.sims[0].error)


def test_parallelism_is_capped_and_the_queue_drains(project, fake_main):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}"})
            for d in (18, 24, 36)]
    config, sims, folder, record, manager = make_run(
        project, fake_main, rows, n_chunks=3, max_parallel=2)

    os.environ["FAKE_SIM_SECONDS"] = "0.6"
    try:
        manager.submit(record)
        running = sum(1 for c in record.chunks if c.status == runstate.RUNNING)
        assert running <= 2, "max_parallel must cap concurrent chunks"
        assert any(c.status == runstate.QUEUED for c in record.chunks)
        drain(manager, record)
    finally:
        os.environ["FAKE_SIM_SECONDS"] = "0.05"

    assert all(c.status == runstate.COMPLETED for c in record.chunks)


def test_cancel_kills_the_whole_process_tree(project, fake_main):
    """Killing only Bryan would orphan cmd.exe and urbs32.exe."""
    psutil = pytest.importorskip("psutil")
    rows = [reservoir_row(**{"Output file": "hangs", "Fake result": "hang"})]
    config, sims, folder, record, manager = make_run(
        project, fake_main, rows, columns=TINAROO_COLUMNS + ["Fake result"])

    manager.submit(record)
    chunk = record.chunks[0]
    assert wait_until(lambda: bool(chunk.pid) and psutil.pid_exists(chunk.pid), 30)
    parent = psutil.Process(chunk.pid)
    assert wait_until(lambda: bool(parent.children(recursive=True)), 30), \
        "the fake simulation should have spawned a grandchild"
    descendants = parent.children(recursive=True)

    manager.cancel(record.run_id)

    assert chunk.status == runstate.CANCELLED
    assert "no run-log entry" in chunk.note
    assert wait_until(lambda: not any(p.is_running() for p in descendants), 20), \
        "grandchildren must not survive the cancel"


def test_request_stop_leaves_the_running_chunk_alone(project, fake_main):
    rows = [reservoir_row(Duration=d, **{"Output file": f"out_{d}h",
                                         "Output suffix": f"o{d}"})
            for d in (18, 24)]
    config, sims, folder, record, manager = make_run(
        project, fake_main, rows, n_chunks=2, max_parallel=1)

    os.environ["FAKE_SIM_SECONDS"] = "0.5"
    try:
        manager.submit(record)
        manager.request_stop(record.run_id)
        drain(manager, record)
    finally:
        os.environ["FAKE_SIM_SECONDS"] = "0.05"

    statuses = {chunk.index: chunk.status for chunk in record.chunks}
    assert statuses[1] == runstate.COMPLETED, "the running chunk finishes"
    assert statuses[2] == runstate.CANCELLED, "the queued one never starts"


def test_reattach_finds_a_run_after_the_ui_restarts(project, fake_main):
    """Disk is authoritative - a closed browser tab must not lose the run."""
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)

    os.environ["FAKE_SIM_SECONDS"] = "1.0"
    try:
        manager.submit(record)
        assert record.chunks[0].status == runstate.RUNNING

        fresh = launcher.RunManager(bryan_python=sys.executable,
                                    bryan_main=str(fake_main))
        found = fresh.reattach(config.run_root)
        assert [r.run_id for r in found] == [record.run_id]
        assert found[0].chunks[0].status == runstate.RUNNING, \
            "a live process must not be written off as dead"

        drain(manager, record, timeout=30)
    finally:
        os.environ["FAKE_SIM_SECONDS"] = "0.05"


def test_reattach_reconciles_a_process_that_vanished(project, fake_main):
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)
    manager.submit(record)
    drain(manager, record)

    # Pretend the UI died mid-run and never reaped it.
    record.chunks[0].status = runstate.RUNNING
    record.chunks[0].pid = 999_999_999
    record.chunks[0].create_time = 1.0
    record.save()

    fresh = launcher.RunManager(bryan_python=sys.executable,
                                bryan_main=str(fake_main))
    found = fresh.reattach(config.run_root)[0]
    assert found.chunks[0].status == runstate.COMPLETED, \
        "the console log's closing summary says it finished"


def test_a_stale_pid_is_not_mistaken_for_this_run(project, fake_main):
    """The create-time token - a bare pid check would reattach to a stranger."""
    psutil = pytest.importorskip("psutil")
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)

    chunk = runstate.ChunkRecord(
        index=9, config="", sims_list="", console_log="", run_log="",
        status=runstate.RUNNING, pid=os.getpid(),
        create_time=psutil.Process(os.getpid()).create_time() - 5000,
    )
    assert not manager._external_alive(chunk)

    chunk.create_time = psutil.Process(os.getpid()).create_time()
    assert manager._external_alive(chunk)


def test_launch_error_when_bryan_is_not_configured(project, fake_main):
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)
    manager.bryan_main = ""
    with pytest.raises(launcher.LaunchError, match="has not been configured"):
        manager.submit(record)


def test_prune_removes_old_finished_runs(project, fake_main):
    rows = [reservoir_row(**{"Output file": "out_18h"})]
    config, sims, folder, record, manager = make_run(project, fake_main, rows)
    manager.submit(record)
    drain(manager, record)

    record.created = time.time() - 400 * 86400
    record.save()
    removed = runstate.prune(config.run_root, older_than_days=30)
    assert removed == [record.run_id]
    assert not folder.folder.exists()
