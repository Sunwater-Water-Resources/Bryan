"""Long steps run to the end, or until cancelled, and say for how long."""

from __future__ import annotations

import sys
import threading
import time

import pytest

from core import processes

SLEEPER = "import time\nprint('started', flush=True)\ntime.sleep(60)\n"


def test_a_finished_process_gives_its_code_and_output(tmp_path):
    done = processes.run([sys.executable, "-c", "import sys; print('hello'); sys.exit(3)"])
    assert done.returncode == 3 and "hello" in done.output
    assert not done.cancelled and not done.timed_out


def test_cancel_stops_it_promptly(tmp_path):
    cancel = threading.Event()
    threading.Timer(1.0, cancel.set).start()
    started = time.monotonic()
    done = processes.run([sys.executable, "-c", SLEEPER], cancel=cancel)
    assert time.monotonic() - started < 15
    assert done.cancelled and done.returncode != 0
    assert "started" in done.output and done.output.endswith("Cancelled.")


def test_a_process_past_its_time_is_stopped():
    done = processes.run([sys.executable, "-c", SLEEPER], timeout=1.0)
    assert done.timed_out and not done.cancelled
    assert "ran past 1 s" in done.output


def test_a_program_that_is_not_there_says_so(tmp_path):
    done = processes.run([tmp_path / "nothing.exe"])
    assert done.returncode == 1 and "could not start" in done.output


def test_a_long_log_does_not_stall_it():
    code = "import sys\nfor i in range(200000): sys.stdout.write('x' * 50 + '\\n')\n"
    done = processes.run([sys.executable, "-c", code], timeout=60)
    assert done.returncode == 0 and done.output.count("\n") >= 200000


def test_the_lake_record_steps_can_be_cancelled(tmp_path):
    from core import lakerecord
    script = tmp_path / "slow.py"
    script.write_text(SLEEPER, encoding="utf-8")
    cancel = threading.Event()
    threading.Timer(1.0, cancel.set).start()
    result = lakerecord.run_script(script, {"out": str(tmp_path)}, tmp_path / "job.json",
                                   sys.executable, cancel)
    assert result.cancelled and not result.ok


pytest.importorskip("nicegui")
pytest_asyncio = pytest.importorskip("pytest_asyncio")

from nicegui import ui                                            # noqa: E402
from nicegui.testing.user_simulation import user_simulation      # noqa: E402


@pytest_asyncio.fixture
async def user():
    async with user_simulation() as simulated:
        yield simulated


@pytest.mark.asyncio
async def test_the_running_line_counts_up_and_cancels(user):
    from widgets import Running
    made = {}

    @ui.page("/running")
    def page() -> None:
        holder = ui.row()
        made["running"] = Running(holder, "homogenising")
        made["holder"] = holder

    await user.open("/running")
    await user.should_see("homogenising")
    running = made["running"]
    running.started -= 75                         # as if a minute and a quarter had gone
    running._tick()
    await user.should_see("homogenising - 1:15")
    user.find(marker="running-cancel").click()
    assert running.cancel.is_set()
    await user.should_see("stopping...")
    running.done()
    await user.should_not_see(marker="running")
