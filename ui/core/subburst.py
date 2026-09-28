"""The Results page's Sub-bursts tab: are the embedded bursts in the design storms neutral?

For every storm duration of a group, Bryan records the wettest window of each shorter
duration in each storm's main burst, and the IFD depth for that window at the same
AEP (``subburst_<d>h`` and ``ifd_<d>h`` in the Monte Carlo database). The total
probability theorem applied to those windows, divided by the IFD, is the
*neutrality margin* of Manual/SubDocs/sub_burst_check.md. Above 1, embedded bursts of
that duration occur more often in the design storms than the rainfall statistics
allow - evidence for filtering them, or for down-weighting the patterns that carry
them. The tab tests one group, optionally beside another (unfiltered against
filtered, say), and can calibrate the pattern weights to neutrality and show what
they do to the design level curve.

The arithmetic needs scipy, so it runs as ``util/SubBurstCheck.py`` under Bryan's
interpreter - the pattern ``core/lakefreq`` set - using Bryan's own
``Simulator.analyse_sub_bursts`` arithmetic and ``CalibrateTpWeights.calibrate``. The
results file is named by a fingerprint of the job and of the databases it reads, in
``_subburst/`` beside the sims_config.json, so a result is shown again until a run
changes it. Nothing here imports nicegui; ``pages/subbursts.py`` is the view.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import events, overlay, processes
from .bryan import BRYAN_ROOT
from .outputs import resolve_value
from .palette import MUTED, OCHRE, PALETTE
from .paths import atomic_write_json, read_json
from .results import format_aep, normal_variate
from .resultchart import aep_axis, numeric

SCRIPT = BRYAN_ROOT / "util" / "SubBurstCheck.py"
CACHE_FOLDER = "_subburst"

# Margins only take seconds per duration; calibrating nine durations takes minutes.
TIMEOUT_SECONDS = 3600

# A margin above this is a breach; between 1 and this, close to neutral. The TPT
# curve of a finite sample wanders a few per cent either side of the truth.
BREACH = 1.05

TESTED, COMPARE, CALIBRATED = "tested", "compare", "calibrated"


# -- what there is to test -----------------------------------------------------------

@dataclass
class Candidate:
    group: str
    databases: dict              # duration label ('72h') -> Path
    recorded: bool               # the databases hold sub-burst depths
    config_file: Path | None = None


_HEADERS: dict = {}


def _has_subbursts(path: Path) -> bool:
    """Does the database record sub-burst depths? Read from its header, cached by mtime."""
    try:
        stat = path.stat()
    except OSError:
        return False
    key = (str(path), stat.st_mtime_ns)
    if key not in _HEADERS:
        try:
            if path.suffix.lower() == ".parquet":
                columns = []                      # the UI cannot read parquet; say unrecorded
            else:
                columns = list(pd.read_csv(path, nrows=0).columns)
        except (OSError, ValueError):
            columns = []
        _HEADERS[key] = any(str(c).startswith("subburst_") for c in columns)
    return _HEADERS[key]


def _by_duration(sources) -> dict:
    """{'24h': path} for one group. The source labels are made unique across every group
    ('24h (CLD_mc_24h_GWL1p3)' once two groups both have a 24h run), so the plain duration
    is taken back where it is unique within the group."""
    plain = [f"{s.duration:g}h" if s.duration is not None else s.label for s in sources]
    return {(name if plain.count(name) == 1 else source.label): source.path
            for name, source in zip(plain, sources)}


def _hours(label) -> float:
    """Sort key for a duration label: '24h' and '24h (the output name)' both sort as 24."""
    try:
        return float(str(label).split("h")[0])
    except ValueError:
        return float("inf")


def candidates(project) -> dict:
    """Every group with a Monte Carlo database, and whether it recorded sub-bursts."""
    frame = project.frame
    folder = project.config.project_folder
    out = {}
    for group, sources in events.sources_by_group(project).items():
        databases = _by_duration(sources)
        config = None
        for source in sources:
            config = resolve_value(folder, frame.loc[source.row_index].get("Config file"))
            if config is not None and config.is_file():
                break
        out[group] = Candidate(group, databases,
                               any(_has_subbursts(path) for path in databases.values()), config)
    return out


def short_label(group) -> str:
    """The last part of a group key, which names the weights files: 'sims_mc\\results\\CLD_mc_GWL1p3'
    is 'CLD_mc_GWL1p3'."""
    return re.split(r"[\\/]", str(group).rstrip("\\/"))[-1] or str(group)


def display_names(groups) -> dict:
    """{group key: the part that tells it from the others}, as the Groups tab's legend has it."""
    return overlay.distinguishing_labels(list(groups))


def scheme(config_file) -> tuple[dict | None, str]:
    """The sampling scheme from a Monte Carlo config: ``scheme_config`` first, then the
    top level, as lib/ReservoirRouting.py reads it (see core/completion.py)."""
    if config_file is None or not Path(config_file).is_file():
        return None, "the group's Monte Carlo config file was not found ('Config file')"
    try:
        data = json.loads(Path(config_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return None, f"{Path(config_file).name} could not be read ({exc})"
    source = data.get("scheme_config", {})
    out = {}
    for key in ("lower_aep", "upper_aep", "number_of_main_divisions", "number_of_sub_divisions",
                "number_of_temporal_patterns", "aep_of_pmp"):
        value = source.get(key, data.get(key))
        if value is not None:
            out[key] = value
    missing = [k for k in ("lower_aep", "upper_aep", "number_of_main_divisions",
                           "number_of_sub_divisions") if k not in out]
    if missing:
        return None, f"{Path(config_file).name} has no {', '.join(missing)}"
    out.setdefault("number_of_temporal_patterns", 10)
    return out, ""


# -- the job -------------------------------------------------------------------------

@dataclass
class Plan:
    job: dict = field(default_factory=dict)
    problems: list = field(default_factory=list)

    @property
    def fingerprint(self) -> str:
        return self.job.get("fingerprint", "")


def cache_folder(config_path) -> Path:
    return Path(config_path).parent / CACHE_FOLDER


def build_job(project, tested: Candidate, compare: Candidate | None = None, *,
              calibrate: bool = False, calibration: dict | None = None) -> Plan:
    plan = Plan()
    if tested is None:
        plan.problems.append("Choose a group to test.")
        return plan
    if not tested.recorded:
        plan.problems.append(f"{tested.group}: its databases hold no sub-burst depths - they were "
                             "made by a Bryan version from before sub-burst tracking.")
    scheme_values, problem = scheme(tested.config_file)
    if problem:
        plan.problems.append(f"{tested.group}: {problem}")
    job = {
        "scheme": scheme_values or {},
        TESTED: {"label": short_label(tested.group),
                 "databases": {d: str(p) for d, p in tested.databases.items()}},
        COMPARE: None if compare is None else {
            "label": short_label(compare.group), "databases": {d: str(p) for d, p in compare.databases.items()}},
        "calibrate": bool(calibrate),
        "calibration": dict(calibration or {}),
    }
    stamp = []
    for role in (TESTED, COMPARE):
        for path in (job[role] or {}).get("databases", {}).values():
            try:
                stat = Path(path).stat()
                stamp.append([path, stat.st_mtime_ns, stat.st_size])
            except OSError:
                stamp.append([path, None, None])
    text = json.dumps({"job": job, "files": stamp}, sort_keys=True)
    job["fingerprint"] = hashlib.sha256(text.encode("utf-8")).hexdigest()[:20]
    job["weights_folder"] = str(cache_folder(project.config.config_path) / "weights"
                                / job["fingerprint"])
    plan.job = job
    return plan


def results_path(config_path, fingerprint) -> Path:
    return cache_folder(config_path) / f"{fingerprint}.json"


def cached(config_path, plan: Plan) -> dict | None:
    saved = read_json(results_path(config_path, plan.fingerprint))
    if saved and saved.get("fingerprint") == plan.fingerprint:
        return saved
    return None


def shown(project, tested, compare, calibration=None) -> dict | None:
    """The calibrated results for this pair if they exist, else the margins alone."""
    config_path = project.config.config_path
    for calibrate in (True, False):
        saved = cached(config_path, build_job(project, tested, compare, calibrate=calibrate,
                                              calibration=calibration))
        if saved is not None:
            return saved
    return None


def write_job(config_path, plan: Plan) -> Path:
    path = cache_folder(config_path) / f"{plan.fingerprint}_job.json"
    atomic_write_json(path, plan.job)
    return path


def command(python, job_path, results) -> list:
    return [str(python), "-u", str(SCRIPT), "--job", str(job_path), "--results", str(results)]


def interpreter_problem(python) -> str:
    if not python:
        return "Bryan's Python interpreter has not been set - set it in the settings panel."
    if not Path(python).exists():
        return f"Interpreter not found: {python}"
    if not SCRIPT.is_file():
        return f"Script not found: {SCRIPT}"
    return ""


def run(argv, cancel=None) -> processes.Finished:
    """Run the check. Blocking - call it off the UI thread."""
    environment = dict(os.environ)
    environment["PYTHONUNBUFFERED"] = "1"
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    return processes.run(argv, env=environment, timeout=TIMEOUT_SECONDS, cancel=cancel)


# -- reading the results ---------------------------------------------------------------

def _aep(key) -> float:
    return float(key)


def standard_aeps(results) -> list:
    return [float(a) for a in (results or {}).get("standard_aeps", [])]


def window(results, low=None, high=None) -> list:
    """The standard AEPs from ``low`` to ``high``; by default all but the first and the
    last, which the TPT's edges pin (the calibration's own default)."""
    aeps = standard_aeps(results)
    if low is None and high is None:
        return aeps[1:-1]
    return [a for a in aeps if (low is None or a >= low) and (high is None or a <= high)]


def _durations(results, role) -> dict:
    return (((results or {}).get("groups") or {}).get(role) or {}).get("durations") or {}


def parents(results, role=TESTED) -> list:
    return sorted(_durations(results, role), key=_hours)


def _subs(margins) -> list:
    return sorted(margins, key=_hours)


def worst(results, role, aeps) -> dict:
    """{parent: {sub: worst margin in the window}}, and the single worst overall."""
    wanted = {float(a) for a in aeps}
    table = {}
    for parent, entry in _durations(results, role).items():
        row = {}
        for sub, series in entry.get("margins", {}).items():
            values = [v for k, v in series.items() if _aep(k) in wanted and v is not None]
            row[sub] = max(values) if values else None
        table[parent] = row
    return table


def worst_overall(results, role, aeps):
    """(margin, parent, sub, aep) of the largest margin in the window, or None."""
    wanted = {float(a) for a in aeps}
    best = None
    for parent, entry in _durations(results, role).items():
        for sub, series in entry.get("margins", {}).items():
            for key, value in series.items():
                if value is not None and _aep(key) in wanted and (best is None or value > best[0]):
                    best = (value, parent, sub, _aep(key))
    return best


def verdict(results, aeps) -> tuple[str, str, str]:
    """(severity, message, hint) for the tested group, read one-sided as the manual reads it."""
    found = worst_overall(results, TESTED, aeps)
    if found is None:
        return "info", "No sub-burst margins in this AEP range.", ""
    value, parent, sub, aep = found
    where = f"the {sub} windows of the {parent} storms reach {value:.2f} times the IFD at 1 in {format_aep(aep)}"
    if value > BREACH:
        return ("warn", f"Embedded bursts occur more often than the IFD implies: {where}.",
                "The patterns over-generate embedded bursts of that duration. Filtering them, or "
                "calibrating the pattern weights, brings the ensemble back to neutral.")
    if value > 1.0:
        return ("info", f"Close to neutral: {where}.",
                f"Within {BREACH - 1:.0%} of the IFD, which is within the sampling wander of the curves.")
    return ("info", f"Neutral: no sub-burst curve reaches the IFD ({where.replace('reach', 'come to')}).",
            "Leaving the embedded bursts unfiltered is not contradicted - though not proven either, "
            "since the check runs one storm duration at a time.")


def matrix_rows(results, aeps) -> tuple[list, list]:
    """Rows for the worst-margin table: one per parent duration, a column per sub-duration,
    each cell 'tested' or 'tested -> compare'."""
    tested = worst(results, TESTED, aeps)
    compare = worst(results, COMPARE, aeps) if _durations(results, COMPARE) else {}
    subs = _subs({s for row in tested.values() for s in row})
    rows = []
    for parent in parents(results):
        row = {"parent": parent}
        for sub in subs:
            a = tested.get(parent, {}).get(sub)
            b = compare.get(parent, {}).get(sub) if compare else None
            if a is None:
                row[sub] = ""
            elif compare:
                row[sub] = f"{a:.2f} → {b:.2f}" if b is not None else f"{a:.2f} → -"
            else:
                row[sub] = f"{a:.2f}"
            row[f"{sub}_breach"] = a is not None and a > BREACH
        rows.append(row)
    return subs, rows


def margin_chart(results, parent) -> dict:
    """Margin against AEP for one parent duration: a line per sub-duration, the
    comparison group dashed in the same colour, and the IFD at 1."""
    tested = _durations(results, TESTED).get(parent, {}).get("margins", {})
    compare = _durations(results, COMPARE).get(parent, {}).get("margins", {})
    aeps = standard_aeps(results)
    series = []
    for index, sub in enumerate(_subs(tested)):
        colour = PALETTE[index % len(PALETTE)]
        for role, margins, dashed in ((TESTED, tested, False), (COMPARE, compare, True)):
            if sub not in margins:
                continue
            data = [[numeric(normal_variate(_aep(k))), numeric(v)]
                    for k, v in sorted(margins[sub].items(), key=lambda kv: _aep(kv[0]))]
            series.append({
                "name": f"{sub}" + (" (compared)" if dashed else ""),
                "type": "line", "symbolSize": 4, "connectNulls": False,
                "itemStyle": {"color": colour},
                "lineStyle": {"width": 2, "type": "dashed" if dashed else "solid"},
                "data": data,
            })
    if series:
        series[0]["markLine"] = {
            "silent": True, "symbol": "none",
            "lineStyle": {"color": OCHRE, "width": 1.5},
            "label": {"formatter": "IFD", "color": OCHRE},
            "data": [{"yAxis": 1.0}],
        }
    horizontal = aep_axis(aeps)
    horizontal["splitLine"] = {"show": True, "lineStyle": {"opacity": 0.25}}
    return {
        "tooltip": {"trigger": "axis",
                    ":valueFormatter": "(v) => v == null ? '-' : Number(v).toFixed(3)"},
        "legend": {"type": "scroll", "top": 0, "textStyle": {"color": MUTED}},
        "grid": {"left": 60, "right": 30, "top": 40, "bottom": 60},
        "xAxis": horizontal,
        # always take in the IFD line at 1, which a scaled axis would leave off the chart
        "yAxis": {"type": "value", "name": "Sub-burst depth / IFD", "nameLocation": "middle",
                  "nameGap": 45, "scale": True, "splitLine": {"lineStyle": {"opacity": 0.25}},
                  ":min": "(v) => Math.floor(Math.min(v.min, 1) * 20 - 1) / 20",
                  ":max": "(v) => Math.ceil(Math.max(v.max, 1) * 20 + 1) / 20"},
        "series": series,
    }


# -- calibration ----------------------------------------------------------------------

def calibration_rows(results) -> list:
    rows = []
    durations = ((results or {}).get("calibration") or {}).get("durations") or {}
    for duration in sorted(durations, key=_hours):
        entry = durations[duration]
        floored = ", ".join(f"{g}: {n}/{len(entry['weights'].get(g, []))}"
                            for g, n in entry.get("at_floor", {}).items() if n)
        rows.append({
            "duration": duration,
            "converged": "yes" if entry.get("converged") else "no",
            "iterations": entry.get("iterations"),
            "margin": (f"{entry['worst_before']:.2f} → {entry['worst_after']:.2f}"
                       if entry.get("worst_before") is not None else ""),
            "floored": floored or "none",
            "weights_file": entry.get("weights_file", ""),
        })
    return rows


def has_curves(results) -> bool:
    return bool(((results or {}).get("calibration") or {}).get("curves"))


def curve_aep(results, case, level) -> float | None:
    """AEP (1 in X) at which the case's design curve reaches ``level``."""
    calibration = (results or {}).get("calibration") or {}
    grid, curve = calibration.get("grid"), (calibration.get("curves") or {}).get(case)
    if not grid or not curve or level is None:
        return None
    grid, p = np.asarray(grid, dtype=float), np.asarray(curve["p"], dtype=float)
    if level < grid[0] or level > grid[-1]:
        return None
    # in log probability, which is near linear in level on the tail; linear in p is not
    positive = p > 0
    if not positive.any():
        return None
    probability = float(np.exp(np.interp(level, grid[positive], np.log(p[positive]))))
    return 1.0 / probability


def level_rows(results, level=None) -> tuple[list, list]:
    """(cases present, rows): the design level at each standard AEP per case, and a last
    row with the AEP of ``level`` on each."""
    curves = ((results or {}).get("calibration") or {}).get("curves") or {}
    cases = [c for c in (COMPARE, TESTED, CALIBRATED) if c in curves]
    rows = []
    for aep in standard_aeps(results):
        key = str(int(aep)) if float(aep).is_integer() else str(aep)
        row = {"aep": format_aep(aep)}
        for case in cases:
            value = curves[case]["levels"].get(key)
            row[case] = "" if value is None else f"{value:.2f}"
        rows.append(row)
    if level is not None:
        row = {"aep": f"AEP of {level:.2f}"}
        for case in cases:
            value = curve_aep(results, case, level)
            row[case] = f"1 in {format_aep(value)}" if value else "off the curve"
        rows.append(row)
    return cases, rows


CASE_LABELS = {COMPARE: "Compared group", TESTED: "Tested group", CALIBRATED: "Tested, calibrated weights"}
