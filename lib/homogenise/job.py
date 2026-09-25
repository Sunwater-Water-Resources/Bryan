"""A homogenisation job: every input and setting, and the run.

The job is a JSON object; relative paths in it resolve against the job file's
folder. What Callide held as constants and a repository layout is here:

    {
      "gauges": ["gauge/130314A.csv", "gauge/130314C.csv"],   # operating order
      "overlay": {"file": "gauge/130314B_hydstra.csv",       # optional
                  "below": 200.90, "reconnect_margin": 0.10},
      "storage": "storage/CALLIDE_STORAGE.els",               # EL,A,V
      "register": "ratings/RatingCurves.xlsx",                # Register + a sheet per rating
      "register_fsl": null,                                   # one rating only, see below
      "evaporation": "climate/silo_-24.35_150.65.txt",        # SILO Data Drill
      "pan_factors": {"1": 0.82, ...},                        # optional, month -> factor
      "step": "1h",
      "recession_correction": true,
      "water_year_start": 10,
      "separation_days": 5, "drop_m": 0.5,
      "targets": [{"name": "rfsl", "rating": "ratings/CALLIDE_RFSL.sq", "fsl": 215.5}],
      "out": "homogenised"
    }

``overlay`` is Callide's intake splice made optional: a second gauge that
replaces the chain wherever it reads below ``below`` (a pool that partitions
below a sediment bar), its last value held until the chain climbs
``reconnect_margin`` clear.

``register`` is a workbook - a ``Register`` sheet of ``Rating, from, to, FSL``
and one ``level, flow`` sheet per rating - or, for a dam whose spillway has not
changed, one rating for the whole record: a URBS ``.rat`` (level, flow) or
``.sq``, or a ``level,flow`` csv. ``register_fsl`` is its full supply level where
the file does not state one; a ``.sq`` without one in its header needs it.

A target rating is a URBS ``.sq`` (storage above full supply against outflow,
tied to the FSL its header declares - never re-based) or a ``level,flow`` csv
starting at the FSL; ``curves.load_rating`` refuses a mismatch with ``fsl``.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import curves, evaporation as evaporation_module, gauges, model, peaks

log = logging.getLogger("bryan.homogenise.job")

DEFAULTS = {
    "gauges": [], "overlay": None, "storage": "", "register": "", "register_fsl": None,
    "evaporation": "",
    "pan_factors": None, "step": "1h", "recession_correction": True,
    "water_year_start": peaks.DEFAULT_WATER_YEAR_START,
    "separation_days": peaks.DEFAULT_SEPARATION_DAYS, "drop_m": peaks.DEFAULT_DROP_M,
    "targets": [], "out": "homogenised",
}


class JobError(Exception):
    """A job that cannot be run, with what to change."""


@dataclass
class Job:
    settings: dict
    folder: Path

    def path(self, key: str) -> Path:
        """The file a setting names, resolved against the job's folder."""
        return self.resolve(self.settings.get(key), what=key)

    def resolve(self, value, what="file") -> Path:
        if not value:
            raise JobError(f"the job names no {what}")
        path = Path(str(value))
        return path if path.is_absolute() else (self.folder / path)

    @property
    def out(self) -> Path:
        return self.path("out")

    def __getitem__(self, key):
        return self.settings[key]


def load(path, routing: bool = True) -> Job:
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise JobError(f"{path} could not be read ({exc})") from exc
    return from_dict(data, path.parent, routing)


def from_dict(data: dict, folder, routing: bool = True) -> Job:
    """``routing=False`` for stage 1 alone (the inflow record): no target rating is
    needed, and the evaporation only where the inflow keeps it."""
    settings = dict(DEFAULTS)
    settings.update(data or {})
    if not settings["gauges"]:
        raise JobError("the job names no gauge exports")
    for key in ("storage", "register") + (("evaporation",) if routing else ()):
        if not settings[key]:
            raise JobError(f"the job names no {key} file")
    if routing and not settings["targets"]:
        raise JobError("the job names no target rating to route through")
    return Job(settings=settings, folder=Path(folder))


# -- inputs --------------------------------------------------------------------

@dataclass
class Inputs:
    volume_of_level: object
    area_of_level: object
    register: object
    ratings: dict
    evaporation: object
    record: pd.DataFrame
    recession_correction: bool = True
    notes: list = field(default_factory=list)

    def dam(self, start=None, end=None):
        return model.LakeModel(self.record, self.evaporation, self.register, self.ratings,
                               self.volume_of_level, self.area_of_level, start=start,
                               end=end, recession_correction=self.recession_correction)


def water_year_start(analysis_job: dict, homogenise_job) -> int:
    """The month an analysis built on a homogenisation starts its water year.

    Its own when its job gives one (``"water_year_start"``), else the
    homogenisation job's. The launcher keeps one water year for the study and
    lets each analysis override it, so the antecedent storage and the inflow
    record need not label their years as the homogenisation does.
    """
    return int(analysis_job.get("water_year_start") or homogenise_job["water_year_start"])


def load_inputs(job: Job, with_evaporation: bool = True) -> Inputs:
    """Read every shared input; clip the record to where they all cover it.

    ``with_evaporation=False`` is for the inflow record without losses: the
    evaporation is then zero over the whole level record, so the record is
    clipped only where the rating register ends, not where SILO does.
    """
    notes = []
    for key in ("storage", "register") + (("evaporation",) if with_evaporation else ()):
        if not job.path(key).is_file():
            raise JobError(f"{key} file not found: {job.path(key)}")
    volume_of_level, area_of_level = curves.read_storage(job.path("storage"))
    fsl = job["register_fsl"]
    register, ratings = curves.read_rating_source(
        job.path("register"), volume_of_level, float(fsl) if fsl not in (None, "") else None)

    record = gauges.read_chain([job.resolve(item, "gauge export") for item in job["gauges"]])
    if with_evaporation:
        evaporation = evaporation_module.read_evaporation(job.path("evaporation"),
                                                          job["pan_factors"])
    else:
        days = pd.date_range(record.index.min().normalize() - pd.Timedelta(days=1),
                             record.index.max().normalize() + pd.Timedelta(days=1), freq="D")
        evaporation = evaporation_module.Evaporation(pd.Series(0.0, index=days))
    overlay = job["overlay"]
    if overlay:
        overlay_path = job.resolve(overlay.get("file"), "overlay gauge")
        if not overlay_path.is_file():
            raise JobError(f"overlay gauge not found: {overlay_path}")
        record = gauges.splice_intake(
            record, gauges.read_wmip_export(overlay_path),
            dead_storage_level=float(overlay["below"]),
            reconnect_margin=float(overlay.get("reconnect_margin",
                                               gauges.INTAKE_RECONNECT_MARGIN)),
            intake_name=overlay_path.stem)
    record = gauges.cap_timestep(record, job["step"])

    usable_end = min(evaporation.end, register["to"].max())
    if record.index.max() > usable_end:
        limit = ("evaporation" if evaporation.end <= register["to"].max()
                 else "rating register")
        notes.append(f"the level record runs to {record.index.max():%Y-%m-%d}; routed to "
                     f"{usable_end:%Y-%m-%d}, where the {limit} ends")
        record = record.loc[:usable_end]
    usable_start = max(evaporation.start, register["from"].min())
    if record.index.min() < usable_start:
        notes.append(f"the level record starts {record.index.min():%Y-%m-%d}; routed from "
                     f"{usable_start:%Y-%m-%d}, where the evaporation or register starts")
        record = record.loc[usable_start:]
    return Inputs(volume_of_level=volume_of_level, area_of_level=area_of_level,
                  register=register, ratings=ratings, evaporation=evaporation,
                  record=record, recession_correction=bool(job["recession_correction"]),
                  notes=notes)


# -- the run -------------------------------------------------------------------

def run_target(job: Job, inputs: Inputs, target: dict, derived=None):
    """Route one target rating. Returns (trace, rating, derived inflow record)."""
    rating_path = job.resolve(target.get("rating"), "target rating")
    if not rating_path.is_file():
        raise JobError(f"target rating not found: {rating_path}")
    rating = curves.load_rating(rating_path, inputs.volume_of_level,
                                fsl=target.get("fsl"), name=target.get("name"))
    dam = inputs.dam()
    if derived is None:
        derived = dam.derive_inflow()
    return dam.simulate(rating, record=derived), rating, derived


def write_outputs(job: Job, name: str, result, rating) -> dict:
    """The run_scenario outputs, one folder per target."""
    folder = job.out / name
    folder.mkdir(parents=True, exist_ok=True)
    written = {}
    trace = folder / "trace.csv.gz"
    result[model.TRACE_COLUMNS].to_csv(trace, float_format="%.4f")
    written["trace"] = trace
    written["daily"] = folder / "daily.csv"
    model.daily_summary(result).to_csv(written["daily"], float_format="%.4f")
    ams = peaks.annual_maxima(result, start_month=int(job["water_year_start"]))
    written["ams"] = folder / "ams.csv"
    ams.to_csv(written["ams"], float_format="%.4f")
    independent = peaks.independent_peaks(
        result, separation_days=float(job["separation_days"]), drop=float(job["drop_m"]),
        threshold=rating.fsl, start_month=int(job["water_year_start"]))
    written["peaks"] = folder / "peaks.csv"
    independent.to_csv(written["peaks"], float_format="%.4f")
    return written


def run(job: Job) -> dict:
    """Every target, sharing one derived inflow. Returns a summary for the launcher."""
    inputs = load_inputs(job)
    derived, targets = None, []
    for target in job["targets"]:
        name = str(target.get("name") or Path(str(target.get("rating"))).stem)
        result, rating, derived = run_target(job, inputs, target, derived)
        files = write_outputs(job, name, result, rating)
        ams = pd.read_csv(files["ams"], index_col=0)
        targets.append({"name": name, "fsl": float(rating.fsl),
                        "files": {key: str(path) for key, path in files.items()},
                        "years": int(len(ams)),
                        "max_level": float(ams["New_Level_max"].max())})
    summary = {"record": {"start": f"{inputs.record.index.min():%Y-%m-%d %H:%M}",
                          "end": f"{inputs.record.index.max():%Y-%m-%d %H:%M}",
                          "steps": int(len(inputs.record)),
                          "gauges": sorted(map(str, inputs.record["Gauge"].unique()))},
               "notes": inputs.notes, "targets": targets}
    job.out.mkdir(parents=True, exist_ok=True)
    (job.out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
