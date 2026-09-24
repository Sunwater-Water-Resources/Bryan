"""The inflow record: annual maximum inflow and event hydrographs from the level record.

    python util/InflowRecord.py <job>.json

The inflow is stage 1 of the homogenisation - the net inflow derived against the
rating in force at each step (lib/homogenise/inflow.py) - so the job names the
homogenisation job and how to read the inflow:

    {
      "homogenise": "homogenise_job.json",       # or the job itself, inline
      "evaporation": false,                       # add lake evaporation back in
      "recession_correction": true,               # see inflow.py: peaks unaffected
      "smoothing": "1h",                          # blank for the native intervals
      "durations_h": [24, 36, 48, 72],
      "catchment_km2": 519,                       # for runoff depths
      "rainfall": "lake_record/rainfall.csv",     # optional: rain beside each volume
      "top_events": 10, "before_days": 3, "after_days": 7,
      "events": [{"name": "Jan 2013", "start": "2013-01-20", "end": "2013-02-05"}],
      "out": "lake_record/inflow"
    }

Writes ``inflow_ams.csv`` (peak inflow, burst volumes, runoff depths, rainfall by
water year), ``inflow_intervals.csv.gz`` (every interval), ``hydrographs/*.csv``
(one per event, both the corrected and the uncorrected inflow, and the intervals
where the release is uncertain), and ``summary.json`` for the launcher.
"""

from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd                                          # noqa: E402

from lib.homogenise import inflow, job as jobs               # noqa: E402

UNSAFE = re.compile(r"[^A-Za-z0-9_.-]+")


def _path(folder: Path, value) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else folder / path


def run(job: dict, folder: Path) -> dict:
    source = job.get("homogenise")
    hjob = (jobs.from_dict(source, folder, routing=False) if isinstance(source, dict)
            else jobs.load(_path(folder, source), routing=False))
    with_evaporation = bool(job.get("evaporation", False))
    inputs = jobs.load_inputs(hjob, with_evaporation=with_evaporation)
    derived = inputs.dam().derive_inflow()
    frame = inflow.intervals(derived, evaporation=with_evaporation,
                             recession_correction=bool(job.get("recession_correction", True)))
    window = job.get("smoothing", "1h") or ""
    frame = inflow.with_smoothing(frame, window)

    rainfall = None
    if job.get("rainfall"):
        path = _path(folder, job["rainfall"])
        if path.is_file():
            rain = pd.read_csv(path, comment="#", parse_dates=["date"])
            rainfall = rain.set_index("date")["rain_mm"].astype(float).sort_index()
    durations = [int(d) for d in job.get("durations_h") or inflow.DURATIONS_H]
    start_month = int(hjob["water_year_start"])
    ams = inflow.annual_maxima(frame, durations, start_month,
                               catchment_km2=job.get("catchment_km2"), rainfall=rainfall)

    out = _path(folder, job.get("out") or "inflow")
    (out / "hydrographs").mkdir(parents=True, exist_ok=True)
    ams.to_csv(out / "inflow_ams.csv", float_format="%.4f")
    frame.to_csv(out / "inflow_intervals.csv.gz", float_format="%.4f")

    events = inflow.event_windows(ams[ams["Complete"]], job.get("top_events") or 0,
                                  job.get("before_days", 3.0), job.get("after_days", 7.0))
    events += [(item["name"], item["start"], item["end"]) for item in job.get("events") or []]
    written = []
    for name, start, end in events:
        part = inflow.hydrograph(frame, start, end, window)
        if part.empty:
            continue
        path = out / "hydrographs" / f"{UNSAFE.sub('_', str(name)).strip('_')}.csv"
        part.to_csv(path, float_format="%.4f")
        written.append({"name": str(name), "file": str(path),
                        "start": f"{pd.Timestamp(start):%Y-%m-%d %H:%M}",
                        "end": f"{pd.Timestamp(end):%Y-%m-%d %H:%M}",
                        "peak_m3s": float(part.get("Inflow_smoothed_m3s",
                                                   part["Inflow_m3s"]).max()),
                        "uncertain_share": float(part["Release_uncertain"].mean())})

    top = ams.sort_values("Peak_inflow_m3s", ascending=False).head(10)
    summary = {
        "record": {"start": f"{frame.index.min():%Y-%m-%d %H:%M}",
                   "end": f"{frame.index.max():%Y-%m-%d %H:%M}", "intervals": int(len(frame))},
        "notes": inputs.notes,
        "settings": {"evaporation": with_evaporation,
                     "recession_correction": bool(job.get("recession_correction", True)),
                     "smoothing": window, "durations_h": durations,
                     "catchment_km2": job.get("catchment_km2")},
        "years": int(len(ams)), "complete_years": int(ams["Complete"].sum()),
        "largest": [{"period": row["Period"], "peak_m3s": float(row["Peak_inflow_m3s"]),
                     "at": f"{pd.Timestamp(row['Peak_time']):%Y-%m-%d %H:%M}"}
                    for _, row in top.iterrows()],
        "ams": str(out / "inflow_ams.csv"),
        "hydrographs": written,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__)
        return 2
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    path = Path(argv[0])
    try:
        summary = run(json.loads(path.read_text(encoding="utf-8")), path.parent)
    except (jobs.JobError, ValueError, FileNotFoundError, KeyError) as exc:
        print(f"ERROR: {exc}")
        return 1
    for note in summary["notes"]:
        print(f"NOTE: {note}")
    print(f"{summary['years']} water years ({summary['complete_years']} complete); "
          f"{len(summary['hydrographs'])} hydrograph(s) -> {Path(summary['ams']).parent}")
    for item in summary["largest"][:5]:
        print(f"  {item['period']}: {item['peak_m3s']:,.0f} m3/s at {item['at']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
