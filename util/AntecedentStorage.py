"""Antecedent storage from the homogenised lake, and the lake configs Bryan samples.

    python util/AntecedentStorage.py <job>.json

The job names the homogenisation job, the catchment rainfall, the IFD and the
method settings (see lib/antecedent/job.py). For each target rating it writes
``<out>/<target>/antecedent_storage.csv``, ``scurve_params.csv`` and one
``lake_config_<basis>_<target>_01.json`` per basis; and ``<out>/summary.json``
for the launcher's Lake record page, which writes the job and runs this with
Bryan's interpreter.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lib.antecedent import job as jobs                     # noqa: E402
from lib.homogenise.job import JobError as HomogeniseError   # noqa: E402


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    try:
        summary = jobs.load_and_run(argv[0])
    except (jobs.JobError, HomogeniseError, ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}")
        return 1
    for note in summary["notes"]:
        print(f"NOTE: {note}")
    for target in summary["targets"]:
        print(f"{target['name']}: {target['qualified']} of {target['years']} water years "
              f"give a sample; FSV {target['fsv_ML']:,.0f} ML")
        for basis, fit in target["fits"].items():
            print(f"  {basis}: k={fit['k']} z0={fit['z0']} Vf={fit['Vf']:,.0f} "
                  f"Vc={fit['Vc']:,.0f} -> {target['configs'][basis]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
