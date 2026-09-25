"""The Lake levels page's logic: settings, the analysis job, the chart and exports.

The page reads a headwater level record (Hydstra or WMIP exports) or an annual
maximum CSV, shows the annual maxima, fits frequency curves through them, lays
the Monte Carlo design floods over them, and exports the report figure and the
series. Its settings live in ``lake_frequency.json`` **beside the project's
sims_config.json**, not in the user's settings: they are inputs to an analysis,
so they belong with the project and can travel with it to a colleague.

Two halves, split on what each needs:

- the record and its annual maxima come from ``lib.LakeLevelRecord``, which is
  pandas-only and imported here directly, so the points, the CSV and the water
  year scores need no subprocess;
- the curve fits, their resampled bands and the report figure need scipy and
  matplotlib, so they run as ``util/LakeLevelFrequency.py`` under Bryan's
  interpreter - the pattern ``core/critexport`` set. It writes a results file
  named by the job's fingerprint into ``_lake_frequency/`` beside the settings,
  and the page draws from that.

Both halves are handed the same job (``LakeLevelRecord.JOB_DEFAULTS``), so the
page and the export cannot disagree about which maxima went in.
"""

from __future__ import annotations

import copy
import json
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import dam as dams, processes
from . import study as studies
from .bryan import BRYAN_ROOT, lake_level_record
from .palette import BODY, MUTED, PALETTE, SURFACE
from .paths import atomic_write_json, clean_path_text, normalise_sep, read_json

RECORD = lake_level_record()

SETTINGS_NAME = "lake_frequency.json"

# With a study open, what of the settings the study keeps: they are about the
# dam's record. The level record itself is the study's dam inputs, and the water
# year the study's unless overridden; the design comparison and exports stay here.
STUDY_KEY = "lake_levels"
STUDY_FIELDS = ("min_coverage", "include_incomplete", "carryover_days", "tie_tolerance",
                "fsl", "fsl_label", "reference_levels", "fit", "axes")
CACHE_FOLDER = "_lake_frequency"
SCRIPT = BRYAN_ROOT / "util" / "LakeLevelFrequency.py"

# The shouldered form resamples in about half a minute; a long record with every
# duration's mcdf to read can take a few.
TIMEOUT_SECONDS = 900

UNSAFE_NAME = re.compile(r"[\\/:*?\"<>|]")

FORM_LABELS = {"shouldered": "Shouldered plateau", "logistic": "Logistic",
               "none": "No curve"}

UPPER_JOIN_LABELS = {"free": "Free - a step up from the plateau is allowed",
                     "fsl": "At full supply - continuous"}

CURVE_HELP = (
    "A gated dam steps up where the gates stop holding the lake, so leave the curve "
    "above FSL free; an uncontrolled spillway is continuous. A plateau resting on one "
    "or two years is fragile - widen the tolerance, or set it to 0. A degree above 1 "
    "above FSL needs 4 maxima above the plateau per coefficient, which a dam that "
    "spills most years has and a dry-belt dam does not. Choose on how the dam works, "
    "not on the rmse: more terms always lower it.")

# What the page adds to a job to make it a settings file: which design floods by
# group and duration rather than by path, since the paths follow from the sims
# list, and where exports go.
PAGE_DEFAULTS = {
    "design": {"include": False, "group": "", "durations": None},
    "export": {"folder": "", "name": ""},
}

# The record's resolution is worth a warning past this: a maximum read off
# readings this far apart can miss the peak.
COARSE_READING_MINUTES = 360


# -- the settings file ------------------------------------------------------------

def settings_path(config_path) -> Path:
    return Path(config_path).resolve().parent / SETTINGS_NAME


def default_settings() -> dict:
    settings = RECORD.complete_job({})
    settings["design"] = dict(PAGE_DEFAULTS["design"])
    settings["export"] = dict(PAGE_DEFAULTS["export"])
    return settings


def complete(settings: dict | None) -> dict:
    out = RECORD.complete_job(settings or {})
    for key, default in PAGE_DEFAULTS.items():
        merged = dict(default)
        merged.update((settings or {}).get(key) or {})
        out[key] = merged
    return out


def load_settings(config_path, study=None) -> dict:
    """The page's settings: this sims list's own, or, with a study open, the study's.

    With a study, the level record is the study's gauge exports (``core/dam.py``)
    or this analysis's own - a homogenised series, say (``record_source``) -
    the water year is the study's unless the analysis overrides it
    (``"water_year"``), and the fit and annual-maximum settings are kept in the
    study (``STUDY_KEY``) - they are about the dam's record, not a run. The design
    flood comparison and the exports stay with the sims list: they name its runs.
    A study that has no Lake levels settings yet starts from this sims list's.
    """
    settings = complete(read_json(settings_path(config_path), default={}) or {})
    if study is None:
        return settings
    base = settings_path(config_path).parent
    part = study.extra.get(STUDY_KEY)
    if not isinstance(part, dict):
        part = _study_part(settings, study, base)
    for key in STUDY_FIELDS:
        if key in part:
            settings[key] = copy.deepcopy(part[key])
    dam = dams.settings(study)
    source = part.get("record_source") or DAM_RECORD
    files = dams.gauges(dam) if source == DAM_RECORD else part.get("files") or []
    settings["record"] = {
        "files": [str(studies.resolve(study.folder, item)) for item in files],
        "ams_csv": (str(studies.resolve(study.folder, part["ams_csv"]))
                    if part.get("ams_csv") else ""),
        "level_column": part.get("level_column", ""),
    }
    settings["record_source"] = source
    settings["water_year"] = part.get("water_year")
    settings["water_year_start"] = int(part.get("water_year") or dam["water_year_start"])
    return settings


# Where the Lake levels record comes from, with a study open. The study's gauge
# exports are the recorded levels; a dam whose spillway changed is compared on a
# homogenised series instead (Callide's is _FFA/CLD_FSL_homogenised.csv), which is
# this analysis's own. A homogenised series is written in the gauges' export
# layout, so the two cannot be told apart from the file - the source is said.
DAM_RECORD = "dam"
OWN_RECORD = "own"


def _study_part(settings: dict, study, base: Path) -> dict:
    """The study's share of the settings, as a sims list's settings give it.

    A list that names a record keeps it as its own - it may be a homogenised
    series, and swapping it for the raw gauges would change the analysis.
    """
    part = {key: copy.deepcopy(settings[key]) for key in STUDY_FIELDS}
    files = [resolve(base, text) for text in settings["record"]["files"] if str(text).strip()]
    part["record_source"] = OWN_RECORD if files else DAM_RECORD
    part["files"] = [studies.portable(study.folder, str(path)) for path in files]
    ams_csv = resolve(base, settings["record"]["ams_csv"])
    part["ams_csv"] = studies.portable(study.folder, str(ams_csv)) if ams_csv else ""
    part["level_column"] = settings["record"]["level_column"]
    shared = int(dams.settings(study)["water_year_start"])
    own = int(settings["water_year_start"])
    part["water_year"] = None if own == shared else own
    return part


def make_study_gauges(settings: dict, study) -> str:
    """Make this analysis's own record the study's gauge exports, and read those.

    For a record that is the recorded levels (Kroombit's 130360A), not a
    homogenised series - which is why it is asked for, never done unasked.
    Returns what was done, to be said.
    """
    files = [path for path in settings["record"]["files"] if str(path).strip()]
    dam = dams.settings(study)
    dam["gauges"] = [studies.portable(study.folder, path) for path in files]
    dams.store(study, dam)
    settings["record_source"] = DAM_RECORD
    return (f"{', '.join(Path(path).name for path in files)} are now the study's gauge "
            f"exports, which every analysis of the lake level record reads.")


def save_settings(config_path, settings: dict, study=None) -> Path:
    """Write the settings, keeping paths inside the project relative to it.

    With a study, its share goes to the study as well. ``lake_frequency.json`` is
    still written whole, so a launcher that has not been updated finds what it
    always did.
    """
    base = settings_path(config_path).parent
    stored = complete(settings)
    stored["record"] = dict(stored["record"])
    if study is not None:
        # With a study, a path in the record is relative to the study; this file's
        # are relative to the sims list, so they go in as full paths first.
        stored["record"]["files"] = [str(studies.resolve(study.folder, text))
                                     for text in stored["record"]["files"]
                                     if str(text).strip()]
        if stored["record"]["ams_csv"]:
            stored["record"]["ams_csv"] = str(studies.resolve(study.folder,
                                                              stored["record"]["ams_csv"]))
    stored["record"]["files"] = [portable(base, text)
                                 for text in stored["record"]["files"] if str(text).strip()]
    stored["record"]["ams_csv"] = portable(base, stored["record"]["ams_csv"])
    stored["export"]["folder"] = portable(base, stored["export"]["folder"])
    path = settings_path(config_path)
    atomic_write_json(path, stored)
    if study is not None:
        part = {key: copy.deepcopy(settings[key]) for key in STUDY_FIELDS}
        part["ams_csv"] = (studies.portable(study.folder, settings["record"]["ams_csv"])
                           if settings["record"]["ams_csv"] else "")
        part["level_column"] = settings["record"]["level_column"]
        part["water_year"] = settings.get("water_year")
        part["record_source"] = settings.get("record_source") or DAM_RECORD
        part["files"] = ([studies.portable(study.folder, text)
                          for text in settings["record"]["files"] if str(text).strip()]
                         if part["record_source"] == OWN_RECORD else [])
        study.extra[STUDY_KEY] = part
        study.save()
    return path


def portable(base: Path, text) -> str:
    """A path as it should be stored: relative when it sits under ``base``."""
    text = clean_path_text(text or "")
    if not text:
        return ""
    path = Path(normalise_sep(text))
    if not path.is_absolute():
        return text
    try:
        return Path(os.path.relpath(path, base)).as_posix() \
            if path.resolve().is_relative_to(base.resolve()) else text
    except (OSError, ValueError):
        return text


def resolve(base: Path, text) -> Path | None:
    text = clean_path_text(text or "")
    if not text:
        return None
    path = Path(normalise_sep(text))
    return path if path.is_absolute() else (Path(base) / path).resolve()


# -- design floods ------------------------------------------------------------------

@dataclass(frozen=True)
class DesignChoice:
    duration: float
    path: Path
    label: str


def design_choices(groups: dict, group: str) -> list:
    """The Monte Carlo databases a group offers, one per duration.

    ``groups`` is ``core.events.sources_by_group(project)``. A row with no
    duration cannot sit on a duration envelope and is left out.
    """
    return [DesignChoice(float(source.duration), Path(source.path), source.label)
            for source in groups.get(group, []) if source.duration is not None]


# -- the job ------------------------------------------------------------------------

@dataclass
class JobPlan:
    job: dict = field(default_factory=dict)
    problems: list = field(default_factory=list)       # stop the record being read
    fit_problems: list = field(default_factory=list)   # stop only the curves

    @property
    def can_read(self) -> bool:
        return not self.problems

    @property
    def can_fit(self) -> bool:
        return not self.problems and not self.fit_problems


def build_job(config_path, settings: dict, groups: dict | None = None) -> JobPlan:
    """The job ``util/LakeLevelFrequency.py`` and the record module both read.

    Paths are made absolute here and nowhere else. Problems are reported rather
    than raised, so the page can show all of them at once.
    """
    settings = complete(settings)
    base = settings_path(config_path).parent
    plan = JobPlan()
    job = RECORD.complete_job({key: value for key, value in settings.items()
                               if key in RECORD.JOB_DEFAULTS})
    files = [resolve(base, text) for text in settings["record"]["files"]]
    files = [path for path in files if path is not None]
    ams_csv = resolve(base, settings["record"]["ams_csv"])
    job["record"] = {"files": [str(path) for path in files],
                     "ams_csv": "" if files or ams_csv is None else str(ams_csv),
                     "level_column": settings["record"]["level_column"]}

    if not files and ams_csv is None:
        plan.problems.append("Name a level record export, or an annual maximum CSV.")
    for path in files if files else ([ams_csv] if ams_csv else []):
        if not path.is_file():
            plan.problems.append(f"Not found: {path}")

    fit = job["fit"]
    if fit["form"] == "shouldered" and job["fsl"] is None:
        plan.fit_problems.append("The shouldered form needs the full supply level.")

    design = settings["design"]
    sources = []
    if design["include"]:
        choices = design_choices(groups or {}, design["group"])
        wanted = design["durations"]
        if wanted is not None:
            wanted = {float(value) for value in wanted}
            choices = [choice for choice in choices if choice.duration in wanted]
        if not choices:
            plan.fit_problems.append("Design floods are switched on, but no Monte Carlo "
                                 "database is chosen.")
        sources = [{"duration": choice.duration, "path": str(choice.path)}
                   for choice in choices]
    job["design"] = {"include": bool(design["include"] and sources),
                     "label": design["group"], "sources": sources}
    plan.job = job
    return plan


def cache_folder(config_path) -> Path:
    return settings_path(config_path).parent / CACHE_FOLDER


def results_path(config_path, job) -> Path:
    return cache_folder(config_path) / f"{RECORD.fingerprint(job)}.json"


def cached_results(config_path, job) -> dict | None:
    """The fitted results for exactly this job and these inputs, if there are any."""
    path = results_path(config_path, job)
    saved = read_json(path)
    if saved and saved.get("fingerprint") == RECORD.fingerprint(job):
        return saved
    return None


def quick_job(job) -> dict:
    """The same job with no resampling: the curves in seconds, without their bands.

    Trying settings is a matter of the fit, which is a fraction of a second; the
    resampling is most of a minute. So **Fit curves** runs this, and the bands are
    resampled once the settings are chosen. Its own fingerprint keeps the two
    results files apart.
    """
    quick = json.loads(json.dumps(job))
    quick["fit"]["draws"] = 0
    return quick


def shown_results(config_path, job) -> tuple:
    """(results, banded): the resampled results for these settings if they exist,
    else the quick fit, else None."""
    banded = cached_results(config_path, job)
    if banded is not None:
        return banded, True
    return cached_results(config_path, quick_job(job)), False


# -- reading the record, cached ------------------------------------------------------

_RECORDS: dict = {}
_RECORD_LIMIT = 2


def _record_key(job) -> tuple:
    stamps = []
    for path in job["record"]["files"] or [job["record"]["ams_csv"]]:
        try:
            info = os.stat(path)
            stamps.append((path, info.st_size, info.st_mtime_ns))
        except OSError:
            stamps.append((path, None, None))
    return tuple(stamps) + (job["record"]["level_column"],)


@dataclass
class RecordView:
    ams: object                  # every water year, LakeLevelRecord.annual_maxima
    positions: object            # the maxima in the curve, with plotting positions
    level: object = None         # the level series, when read from exports
    sites: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def read_view(job) -> RecordView:
    """The annual maxima for a job. Reading the export is the slow part, so the
    level series is kept per file set; the maxima are re-derived every call,
    because the water year and carry-over settings change them."""
    key = _record_key(job)
    level_record = None
    if job["record"]["files"]:
        level_record = _RECORDS.get(key)
        if level_record is None:
            level_record = RECORD.read_record(job["record"]["files"])
            _RECORDS[key] = level_record
            while len(_RECORDS) > _RECORD_LIMIT:
                _RECORDS.pop(next(iter(_RECORDS)))
        ams = RECORD.annual_maxima(level_record.level, int(job["water_year_start"]),
                                   float(job["min_coverage"]),
                                   float(job["carryover_days"]),
                                   float(job["tie_tolerance"]))
    else:
        ams, _ = RECORD.ams_for_job(job)
    positions = RECORD.with_positions(ams, bool(job["include_incomplete"]))
    notes = list(level_record.notes) if level_record is not None else []
    coarse = positions[positions["reading_interval_min"] > COARSE_READING_MINUTES]
    if len(coarse):
        notes.append(f"{len(coarse)} annual maxim{'a were' if len(coarse) > 1 else 'um was'} "
                     f"read off a record logged "
                     f"less often than every {COARSE_READING_MINUTES // 60} hours "
                     f"({', '.join(coarse['period'].head(6))}"
                     f"{'...' if len(coarse) > 6 else ''}); a peak between readings "
                     f"would be missed")
    incomplete = ams[~ams["complete"]]
    if len(incomplete):
        verb = "included" if job["include_incomplete"] else "left out"
        many = len(incomplete) > 1
        notes.append(f"{len(incomplete)} water year{'s are' if many else ' is'} not fully "
                     f"covered ({', '.join(incomplete['period'])}) and "
                     f"{'are' if many else 'is'} {verb}")
    return RecordView(ams=ams, positions=positions,
                      level=None if level_record is None else level_record.level,
                      sites=[] if level_record is None else level_record.sites,
                      notes=notes)


def is_quick(job) -> bool:
    """Whether reading the view is quick enough to do on the event loop: an
    annual maximum CSV, or exports already read."""
    return not job["record"]["files"] or _record_key(job) in _RECORDS


def forget_records() -> None:
    _RECORDS.clear()


# -- running the util script ---------------------------------------------------------

@dataclass
class ScriptResult:
    returncode: int
    output: str
    written: tuple = ()
    cancelled: bool = False

    @property
    def ok(self) -> bool:
        return self.returncode == 0


def write_job(config_path, job) -> Path:
    path = cache_folder(config_path) / f"{RECORD.fingerprint(job)}_job.json"
    atomic_write_json(path, job)
    return path


def command(python, job_path, *, results=None, png=None, without_design=False,
            ams_csv=None) -> list:
    argv = [str(python), "-u", str(SCRIPT), "--job", str(job_path)]
    if results:
        argv += ["--results", str(results)]
    if png:
        argv += ["--png", str(png)]
    if without_design:
        argv.append("--without-design")
    if ams_csv:
        argv += ["--ams-csv", str(ams_csv)]
    return argv


def interpreter_problem(python) -> str:
    if not python:
        return "Bryan's Python interpreter has not been set - set it on the Project page."
    if not Path(python).exists():
        return f"Interpreter not found: {python}"
    if not SCRIPT.is_file():
        return f"Script not found: {SCRIPT}"
    return ""


def run(argv, outputs=(), cancel=None) -> ScriptResult:
    """Run the util script. Blocking - call it off the UI thread; ``cancel`` (a
    threading.Event) stops it."""
    environment = dict(os.environ)
    environment["PYTHONUNBUFFERED"] = "1"
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    finished = processes.run(argv, env=environment, timeout=TIMEOUT_SECONDS, cancel=cancel)
    return ScriptResult(finished.returncode, finished.output,
                        tuple(Path(path) for path in outputs if Path(path).is_file()),
                        finished.cancelled)


# -- exports ----------------------------------------------------------------------------

@dataclass
class ExportPaths:
    png: Path | None = None
    csv: Path | None = None
    problems: list = field(default_factory=list)


def default_export(config_path, settings) -> tuple:
    """(folder, name) - the saved ones, or the project folder and the site."""
    base = settings_path(config_path).parent
    export = complete(settings)["export"]
    folder = resolve(base, export["folder"]) or base
    return folder, export["name"] or "lake_level_frequency"


def export_paths(folder, name, *, with_design) -> ExportPaths:
    paths = ExportPaths()
    name = (name or "").strip()
    if not name or UNSAFE_NAME.search(name):
        paths.problems.append("The name is empty or holds a path separator - it becomes "
                              "a filename, so it cannot contain \\ / : * ? \" < > |.")
    if not folder:
        paths.problems.append("No output folder.")
    if paths.problems:
        return paths
    folder = Path(folder)
    paths.png = folder / f"{name}_{'validation' if with_design else 'record'}.png"
    paths.csv = folder / f"{name}_ams.csv"
    return paths


def write_ams_csv(view: RecordView, job, path) -> Path:
    table = view.positions.sort_values("water_year").reset_index(drop=True)
    level_record = None
    if view.level is not None:
        level_record = RECORD.Record(level=view.level, sites=view.sites)
    text = RECORD.ams_csv_text(table, RECORD.ams_metadata(job, level_record))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="")
    return path


# -- the chart ------------------------------------------------------------------------

RECORD_COLOUR = PALETTE[2]      # the house blue
DESIGN_COLOUR = PALETTE[9]      # and its burnt orange: the same pairing as before

EY_TICKS = (4, 3, 2, 1)
AEP_TICKS = (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 5e-3, 2e-3, 1e-3, 5e-4, 1e-4, 1e-5, 1e-6)


def _num(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) or math.isinf(number) else number


def frequency_label(aep) -> str:
    """EY inside 50% AEP, "1 in X" past it - the report figure's convention."""
    if aep > 0.5:
        return f"{-math.log(1.0 - aep):.3g}EY"
    return f"1 in {1.0 / aep:,.0f}"


def frequency_axis(z_min, z_max) -> dict:
    marks = [1.0 - math.exp(-ey) for ey in EY_TICKS] + list(AEP_TICKS)
    ticks = {}
    for aep in marks:
        z = RECORD.normal_variate_of_probability(aep)
        if z_min - 1e-9 <= z <= z_max + 1e-9:
            ticks[f"{z:.4f}"] = frequency_label(aep)
    return {
        "type": "value", "min": z_min, "max": z_max,
        "name": "Exceedance frequency (EY, then 1 in X AEP)",
        "nameLocation": "middle", "nameGap": 34,
        # z = 0 is 50% AEP, and the level axis would otherwise be drawn through it.
        "axisLine": {"onZero": False},
        "splitLine": {"show": True, "lineStyle": {"opacity": 0.25}},
        "axisLabel": {"customValues": [float(k) for k in ticks], "hideOverlap": True,
                      ":formatter": "(v) => { const m = " + json.dumps(ticks)
                                    + "; return m[v.toFixed(4)] ?? ''; }"},
        "axisTick": {"customValues": [float(k) for k in ticks]},
    }


def _line(name, xs, ys, keep=None, **style) -> dict:
    data = [[_num(x), _num(y)] for i, (x, y) in enumerate(zip(xs, ys))
            if keep is None or keep[i]]
    series = {"name": name, "type": "line", "symbol": "none", "smooth": False,
              "connectNulls": False, "data": data}
    series.update(style)
    return series


def chart_options(positions, results=None, settings=None, *, show_design=True,
                  whole_record=False) -> dict:
    """The page's chart: the maxima, and whatever of the fit and design floods exist.

    ``positions`` is ``RecordView.positions``. ``results`` is the util script's
    output for the same job, or None when the curves have not been fitted.
    ``whole_record`` widens the frame from the report figure's (1EY onwards) to
    every maximum.
    """
    settings = complete(settings)
    rows = positions.to_dict("records")
    storm = [row for row in rows if not row["carried_over"]]
    carried = [row for row in rows if row["carried_over"]]

    def point(row, z_key="z"):
        when = row.get("level_at")
        return {"value": [_num(row[z_key]), _num(row["level"])],
                "name": f"{row['period']}: {row['level']:.3f} m, "
                        f"{frequency_label(row['aep'])}"
                        + (f", {str(when)[:16]}" if when is not None and when == when else "")
                        + (" (carried over)" if row["carried_over"] else "")}

    series = []
    grid_z = (results or {}).get("grid", {}).get("z") or []
    fits = (results or {}).get("fits") or {}
    whole, storm_fit = fits.get("all") or {}, fits.get("storm") or {}

    def extent(block):
        z_max = block.get("z_max") if block.get("form") == "shouldered" else None
        return [z_max is None or z <= z_max for z in grid_z]

    if whole.get("band_lo"):
        keep = extent(whole)
        lo = [value for value, k in zip(whole["band_lo"], keep) if k]
        width = [h - l for h, l, k in zip(whole["band_hi"], whole["band_lo"], keep) if k]
        zs = [z for z, k in zip(grid_z, keep) if k]
        series.append(_line("band", zs, lo, stack="band", silent=True,
                            lineStyle={"opacity": 0}, tooltip={"show": False}))
        series.append(_line(f"90% band, all maxima ({whole['draws_used']:,} resamples)",
                            zs, width, stack="band", silent=True,
                            lineStyle={"opacity": 0}, tooltip={"show": False},
                            itemStyle={"color": RECORD_COLOUR},
                            areaStyle={"color": RECORD_COLOUR, "opacity": 0.13}))

    design = (results or {}).get("design") if show_design else None
    if design:
        for hours, curve in design["durations"].items():
            keep = [aep <= 0.5 for aep in curve["aep"]]
            series.append(_line("Monte Carlo durations", curve["z"], curve["level"], keep,
                                silent=True, itemStyle={"color": MUTED},
                                lineStyle={"width": 1, "opacity": 0.45},
                                tooltip={"show": False}))
        series.append(_line("Design flood envelope", grid_z, design["envelope"],
                            itemStyle={"color": DESIGN_COLOUR}, lineStyle={"width": 3}))

    series.append({"name": "Storm-driven maxima", "type": "scatter", "symbolSize": 8,
                   "itemStyle": {"color": RECORD_COLOUR},
                   "data": [point(row) for row in storm], "z": 5})
    series.append({"name": "Carried over the water year", "type": "scatter",
                   "symbolSize": 8,
                   "itemStyle": {"color": SURFACE, "borderColor": RECORD_COLOUR,
                                 "borderWidth": 1.5},
                   "data": [point(row) for row in carried], "z": 5})

    if whole.get("curve"):
        series.append(_line(f"Fit to all maxima (RMSE {whole['rmse']:.2f} m)", grid_z,
                            whole["curve"], extent(whole),
                            itemStyle={"color": RECORD_COLOUR}, lineStyle={"width": 2.5}))
    if storm_fit.get("curve") and carried:
        censored = [row for row in storm if _num(row.get("storm_z")) is not None]
        z_floor = min(row["storm_z"] for row in censored) if censored else -99
        keep = [k and z >= z_floor for k, z in zip(extent(storm_fit), grid_z)]
        series.append({"name": "…at censored positions", "type": "scatter",
                       "symbol": "path://M0,0L10,10M10,0L0,10", "symbolSize": 7,
                       "itemStyle": {"color": "rgba(0,0,0,0)", "borderColor": RECORD_COLOUR,
                                     "borderWidth": 1.2, "opacity": 0.6},
                       "data": [point(row, "storm_z") for row in censored], "z": 4})
        if storm_fit.get("band_lo"):
            for position, edge in enumerate((storm_fit["band_lo"], storm_fit["band_hi"])):
                series.append(_line("90% band, storm-driven", grid_z, edge, keep,
                                    silent=True, tooltip={"show": False},
                                    itemStyle={"color": RECORD_COLOUR},
                                    lineStyle={"width": 1, "type": "dotted",
                                               "opacity": 0.8}))
        series.append(_line(f"Fit to storm-driven (RMSE {storm_fit['rmse']:.2f} m)",
                            grid_z, storm_fit["curve"], keep,
                            itemStyle={"color": RECORD_COLOUR},
                            lineStyle={"width": 2, "type": "dashed"}))

    references = [(item.get("label", ""), _num(item.get("level")))
                  for item in settings.get("reference_levels") or []]
    if _num(settings.get("fsl")) is not None:
        references.append((settings.get("fsl_label") or "FSL", _num(settings["fsl"])))
    references = [(label, level) for label, level in references if level is not None]
    if references:
        anchor = next(item for item in series if item["name"] == "Storm-driven maxima")
        anchor["markLine"] = {
                "silent": True, "symbol": "none",
                "lineStyle": {"color": BODY, "type": "dotted", "width": 1},
                "label": {"position": "insideStartTop", "formatter": "{b}",
                          "color": "#333", "fontSize": 11},
                "data": [{"name": f"{label} {level:g} m", "yAxis": level}
                         for label, level in references]}

    zs = [row["z"] for row in rows if _num(row["z"]) is not None]
    ey1 = RECORD.normal_variate_of_probability(1.0 - math.exp(-1.0))
    z_max = max(zs + [grid_z[-1]] if design and grid_z else zs) + 0.1
    axes = settings.get("axes") or {}
    y_min, y_max = _num(axes.get("level_min")), _num(axes.get("level_max"))
    if whole_record:
        z_min = min([ey1] + zs) - 0.1
    else:
        # The report figure's frame: from 1EY, and levels from what is in it -
        # the maxima, the all-maxima band and the envelope - so the frequent
        # maxima and the storm-driven band's tail cannot squash the comparison.
        z_min = ey1
        seen = [row["level"] for row in rows if row["z"] >= z_min]
        if whole.get("band_lo"):
            keep = extent(whole)
            seen += [v for v, z, k in zip(whole["band_lo"], grid_z, keep) if k and z >= z_min]
            seen += [v for v, z, k in zip(whole["band_hi"], grid_z, keep) if k and z >= z_min]
        if design:
            seen += design["envelope"]
        seen += [level for _, level in references]
        seen = [value for value in seen if _num(value) is not None]
        if seen and y_min is None:
            y_min = math.floor(min(seen) - 0.3)
        if seen and y_max is None:
            y_max = math.ceil((max(seen) + 0.8) * 5) / 5
    legend_names = []
    for item in series:
        if item["name"] != "band" and item["name"] not in legend_names:
            legend_names.append(item["name"])
    return {
        # Thousands of design flood points redrawn on every setting change; the
        # entry animation only delays the picture.
        "animation": False,
        "tooltip": {"trigger": "item", ":formatter": "(p) => p.name || ''"},
        "legend": {"type": "scroll", "top": 0, "data": legend_names},
        "grid": {"left": 64, "right": 24, "top": 56, "bottom": 64},
        "xAxis": frequency_axis(z_min, z_max),
        "yAxis": {"type": "value", "scale": True, "axisLine": {"show": True, "onZero": False},
                  "name": "Annual maximum lake level (m AHD)", "nameLocation": "middle",
                  "nameGap": 48, "min": y_min, "max": y_max,
                  "splitLine": {"lineStyle": {"opacity": 0.25}}},
        # Scroll to zoom, drag to pan; nothing is filtered out as it leaves view.
        "dataZoom": [{"type": "inside", "xAxisIndex": 0, "filterMode": "none"},
                     {"type": "inside", "yAxisIndex": 0, "filterMode": "none"}],
        "series": series,
    }


# -- water year options ---------------------------------------------------------------

SCORE_LABELS = {
    "start_month": "Water year starts",
    "years": "Years",
    "carried_over": "Carried over",
    "carried_over_high": "Carried over, high",
    "closest_storm_max_d": "Closest storm maximum (d)",
    "cuts_rising": "Cuts a rise",
    "cuts_above_fsl": "Cuts above FSL",
    "promoted": "Promoted",
    "level_50pct_aep": "Level at 50% AEP",
    "level_20pct_aep": "Level at 20% AEP",
    "level_10pct_aep": "Level at 10% AEP",
}

SCORE_HELP = (
    "Carried over: maxima inherited from the year before (high: above the level "
    "given). Closest storm maximum: how near any storm-driven maximum comes to a "
    "year boundary - the margin against a flood being split. Cuts a rise: "
    "boundaries where the lake rose over the day either side. Promoted: storm-"
    "driven maxima below a level reached within six months, so a maximum only "
    "because of where the year was cut. Levels: the empirical curve, read "
    "between plotting positions.")


def score_rows(frame, adopted_month: int) -> list:
    rows = []
    adopted = RECORD.MONTH_NAMES[adopted_month - 1]
    for record in frame.to_dict("records"):
        row = {}
        for key, value in record.items():
            if isinstance(value, float):
                value = None if math.isnan(value) else round(value, 2)
            row[key] = value
        row["adopted"] = record["start_month"] == adopted
        rows.append(row)
    return rows


def seasonal_chart(monthly, counts) -> dict:
    """Median level by calendar month, with its 10-90% range, and the month the
    storm-driven maxima fall in."""
    months = list(monthly["month"])
    return {
        "tooltip": {"trigger": "axis"},
        "legend": {"top": 0},
        "grid": {"left": 64, "right": 56, "top": 36, "bottom": 36},
        "xAxis": {"type": "category", "data": months},
        "yAxis": [{"type": "value", "scale": True, "name": "Daily level (m AHD)",
                   "nameLocation": "middle", "nameGap": 48,
                   "splitLine": {"lineStyle": {"opacity": 0.25}}},
                  {"type": "value", "name": "Maxima", "nameLocation": "middle",
                   "nameGap": 32, "minInterval": 1, "splitLine": {"show": False}}],
        "series": [
            {"name": "10%", "type": "line", "symbol": "none", "stack": "range",
             "lineStyle": {"opacity": 0}, "data": [_num(v) for v in monthly["p10"]],
             "tooltip": {"show": False}},
            {"name": "10-90% of days", "type": "line", "symbol": "none", "stack": "range",
             "lineStyle": {"opacity": 0},
             "areaStyle": {"color": RECORD_COLOUR, "opacity": 0.12},
             "data": [_num(h - l) for h, l in zip(monthly["p90"], monthly["p10"])]},
            {"name": "Median daily level", "type": "line",
             "itemStyle": {"color": RECORD_COLOUR},
             "data": [_num(v) for v in monthly["median"]]},
            {"name": "Storm-driven maxima", "type": "bar", "yAxisIndex": 1,
             "itemStyle": {"color": DESIGN_COLOUR, "opacity": 0.6},
             "data": [int(v) for v in counts.reindex(months, fill_value=0)]},
        ],
    }
