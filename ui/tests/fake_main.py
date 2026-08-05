"""A stand-in for Bryan's Main.py, for testing the launcher without URBS.

Takes the same command line, reads the same config and sims list, and imports
the **real** ``lib.RunLog`` and ``lib.LogFiles`` so the run-log format cannot
drift away from what the UI parses. The banners it prints come from the same
regexes core/progress.py matches against, and test_progress.py separately greps
the real Main.py for that literal text - so if Bryan's wording changes, the
suite says so rather than the progress display quietly going blank.

An extra ``Fake result`` column drives the behaviour of each row:

    ok             completes (the default)
    fail           raises, so the run log records FAILED
    missing input  raises FileNotFoundError -> 'FAILED - missing input'
    hang           spawns a long-lived grandchild, so tree-kill is really tested
    ask            calls input(), so the stdin=DEVNULL -> EOFError path is pinned
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from lib import LogFiles, RunLog  # noqa: E402

SLEEP_SECONDS = float(os.environ.get("FAKE_SIM_SECONDS", "0.05"))


def main() -> int:
    simulation_file = sys.argv[1]
    print("Opening simulation configuration file: ", simulation_file)
    with open(simulation_file) as stream:
        sim_config = json.load(stream)

    project_folder = os.path.dirname(simulation_file)
    if "project_folder" in sim_config:
        if str(sim_config["project_folder"]).lower() != "default":
            project_folder = sim_config["project_folder"]

    sim_list_file = os.path.join(project_folder, sim_config["simulation_list"])
    sim_df = pd.read_excel(sim_list_file, sheet_name=0)

    filepaths = sim_config.get("filepaths", {})
    folder = os.path.dirname(simulation_file) or os.getcwd()
    for key, relpath in list(filepaths.items()):
        filepaths[key] = os.path.normpath(os.path.join(folder, relpath))

    log_filepath = RunLog.log_filepath(sim_config)
    included = sim_df[sim_df["Include"] == "yes"]
    included = LogFiles.resolve_duplicates(included)

    failures = []
    total = len(included)
    for position, (_, sim) in enumerate(included.iterrows(), start=1):
        entry = RunLog.start_entry(sim, filepaths, executable=sys.argv[0],
                                   config_file=sys.argv[1])
        # The exact banner Main.py:83-86 prints.
        print(f'\n{"=" * 90}')
        print(f'Simulation {position} of {total}: {sim["Output file"]}')
        print(f'  method: {sim["Method"]}    method config: {sim["Config file"]}')
        print(f'{"=" * 90}')

        try:
            _simulate(sim, project_folder)
            status, error = RunLog.COMPLETED, ""
        except Exception as exception:      # noqa: BLE001 - mirrors Main.py:97
            status = RunLog.status_for(exception)
            error = f"{type(exception).__name__}: {exception}"
            failures.append((sim["Output file"], status, error))
            print(f'\n{"!" * 60}')
            print(f'{status}: {sim["Output file"]}')
            traceback.print_exc()
            print("Moving on to the next simulation in the list.")
            print(f'{"!" * 60}')

        RunLog.write_entry(RunLog.close_entry(entry, status, error), log_filepath)

    RunLog.print_summary(failures, total)
    return 1 if failures else 0


def _simulate(sim, project_folder) -> None:
    behaviour = str(sim.get("Fake result", "ok") or "ok").strip().lower()

    log_file = sim.get("Log file")
    if isinstance(log_file, str) and log_file.strip():
        path = Path(project_folder) / log_file.replace("\\", "/")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as stream:
            stream.write(f"log for {sim['Output file']}\n")
            stream.flush()
            time.sleep(SLEEP_SECONDS)
            stream.write("done\n")
    else:
        time.sleep(SLEEP_SECONDS)

    if behaviour == "fail":
        raise RuntimeError("the fake simulation was told to fail")
    if behaviour == "missing input":
        raise FileNotFoundError("no such file: pretend.csv")
    if behaviour == "ask":
        # MCScheme.store_simulations does this when a CSV is open in Excel. With
        # stdin=DEVNULL it raises EOFError instead of hanging.
        input("Could not save the simulation file. Please close it and press enter.")
    if behaviour == "hang":
        # A grandchild that outlives a naive kill of the parent only.
        subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
        time.sleep(600)

    # Write something that looks like a results database, so completion
    # detection has a file to find.
    output_file = str(sim.get("Output file") or "")
    results_folder = sim.get("Results folder")
    if output_file:
        if isinstance(results_folder, str) and results_folder.strip():
            base = (Path(project_folder) / results_folder.replace("\\", "/")
                    / f"{Path(output_file).name}__mcdf.csv")
        else:
            base = Path(project_folder) / f"{output_file}__mcdf.csv".replace("\\", "/")
        base.parent.mkdir(parents=True, exist_ok=True)
        base.write_text("index,inflow,level,outflow\n0,1,2,3\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
