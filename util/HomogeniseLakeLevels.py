"""Homogenise a dam's recorded lake levels onto one or more target ratings.

    python util/HomogeniseLakeLevels.py <job>.json

The job names the gauge exports, storage table, rating register, evaporation
and target ratings (see lib/homogenise/job.py). For each target it writes
``<out>/<target>/``: ``trace.csv.gz`` (every routing step), ``daily.csv``,
``ams.csv`` (the annual maxima by water year) and ``peaks.csv``; and
``<out>/summary.json`` for the launcher. The launcher's Lake record page writes
the job and runs this with Bryan's interpreter.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lib.homogenise import job as jobs                     # noqa: E402


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    try:
        job = jobs.load(argv[0])
        summary = jobs.run(job)
    except (jobs.JobError, ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}")
        return 1
    for note in summary["notes"]:
        print(f"NOTE: {note}")
    for target in summary["targets"]:
        print(f"wrote {target['name']}: {target['years']} water years, highest "
              f"{target['max_level']:.3f} m AHD -> {Path(target['files']['ams']).parent}")
    print(json.dumps(summary["record"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
