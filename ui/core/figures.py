"""Report figures: frequency curves drawn from a study's runs.

The successor to ``PlotFrequencyCurves_v03.py`` and its plot-list workbook. A
figure lives in the study file (``"figures"``) and names its curves by **run and
group**, never by path, so moving the figures from E012 to E013 is one change to
a run rather than a new plot list with every path retyped; and each curve
carries its **own label for that figure**, because the same curve is 'URBS'
against the FFA and 'GWL 1.3 °C' among the horizons.

Three kinds of curve:

``group``  A group's design curve - the envelope over its durations - read the
           way the Results page reads it (``results.compare``/``analyse``).
``file``   A table with AEP ('1 in X') in its first column and the result as a
           column: a previous study's adopted curve (Sunwater 2020).
``ffa``    An RMC Bestfit export: the posterior mode and/or mean, the 90%
           credible interval, the annual maxima, and any paleoflood lower
           bounds - what the v03 script drew from its 'FFA' sheet.

The UI builds the data here (pandas only) and draws the preview; **Export**
hands exactly that data to ``util/ReportFigure.py``, which draws the PNG with
matplotlib under Bryan's interpreter in the v03 script's style. The preview and
the PNG are therefore the same numbers and the same labels.

Nothing here imports nicegui.
"""

from __future__ import annotations

import copy
import json
import math
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import results
from .bryan import BRYAN_ROOT
from .paths import atomic_write_json
from .study import Study, StudyError, resolve

SCRIPT = BRYAN_ROOT / "util" / "ReportFigure.py"
TIMEOUT_SECONDS = 300
FIGURES_KEY = "figures"
UNSAFE_NAME = re.compile(r"[\\/:*?\"<>|]")

GROUP, FILE, FFA = "group", "file", "ffa"
POSTERIORS = ("Mode", "Mean", "Both", "AMS only")

RESULT_LABELS = {"level": "Lake level (m AHD)", "inflow": "Flow (m³/s)",
                 "outflow": "Flow (m³/s)"}

TEMPLATE = {
    "id": "", "title": "", "filename": "", "type": "level",
    "curves": [],                   # [{kind, label, run, group | file, column | ...}]
    "reference_levels": [],         # [{label, level, colour}]
    "min_aep": 2, "max_aep": 2_000_000, "aep_of_pmp": None,
    "y_min": None, "y_max": None, "show_title": False,
}


def result_label(key: str) -> str:
    if key in RESULT_LABELS:
        return RESULT_LABELS[key]
    match = re.match(r"^(?:inflow)?Vol(\d+(?:[._]\d+)?)h$", key or "", re.IGNORECASE)
    if match:
        return f"{match.group(1).replace('_', '.')} hr volume (ML)"
    return key


def logarithmic(key: str) -> bool:
    return key != "level"


def standard_aeps(lower, upper) -> list:
    """2, 5, 10, 20, 50 ... from ``lower`` to ``upper``, as v03's get_standard_aeps."""
    out, aep = [], float(lower)
    while aep <= float(upper):
        out.append(aep)
        aep = aep * 5 // 2 if math.log10(aep / 2) % 1 == 0 else aep * 2
    return out


# -- the study keeps them --------------------------------------------------------

def figures(study: Study) -> list:
    return [spec for spec in study.extra.get(FIGURES_KEY) or [] if isinstance(spec, dict)]


def new_spec(**overrides) -> dict:
    spec = copy.deepcopy(TEMPLATE)
    spec.update(overrides)
    return spec


def put(study: Study, spec: dict) -> dict:
    spec = copy.deepcopy(spec)
    stored = figures(study)
    if not spec.get("id"):
        taken = {each.get("id") for each in stored}
        base = re.sub(r"[^a-z0-9]+", "-", (spec.get("filename") or spec.get("title")
                                           or "figure").lower()).strip("-") or "figure"
        spec["id"], n = base, 2
        while spec["id"] in taken:
            spec["id"], n = f"{base}-{n}", n + 1
    for position, each in enumerate(stored):
        if each.get("id") == spec["id"]:
            stored[position] = spec
            break
    else:
        stored.append(spec)
    study.extra[FIGURES_KEY] = stored
    return spec


def remove(study: Study, figure_id: str) -> None:
    study.extra[FIGURES_KEY] = [each for each in figures(study) if each.get("id") != figure_id]


def figure(study: Study, figure_id: str) -> dict | None:
    return next((each for each in figures(study) if each.get("id") == figure_id), None)


# -- the data ----------------------------------------------------------------------

@dataclass
class Series:
    """One line or set of markers, as (AEP '1 in X', value) pairs."""

    label: str
    points: list                     # [(aep, value)], ascending AEP
    style: str = "line"              # line, ffa-mode, ffa-mean, ffa-ci, ams, paleo
    legend: bool = True


@dataclass
class FigureData:
    spec: dict
    series: list = field(default_factory=list)
    problems: list = field(default_factory=list)

    @property
    def key(self) -> str:
        return self.spec.get("type") or "level"

    def to_job(self, png: Path) -> dict:
        spec = self.spec
        return {"png": str(png), "type": self.key, "y_label": result_label(self.key),
                "log_y": logarithmic(self.key), "title": spec.get("title") or "",
                "show_title": bool(spec.get("show_title")),
                "min_aep": spec.get("min_aep") or 2, "max_aep": spec.get("max_aep") or 2e6,
                "aep_of_pmp": spec.get("aep_of_pmp"),
                "y_min": spec.get("y_min"), "y_max": spec.get("y_max"),
                "reference_levels": spec.get("reference_levels") or [],
                "series": [{"label": s.label, "style": s.style, "legend": s.legend,
                            "points": [[float(a), float(v)] for a, v in s.points]}
                           for s in self.series]}


def _trimmed(series: pd.Series, spec) -> list:
    low, high = float(spec.get("min_aep") or 2), float(spec.get("max_aep") or 2e6)
    out = []
    for aep, value in series.items():
        try:
            aep, value = float(aep), float(value)
        except (TypeError, ValueError):
            continue
        if low <= aep <= high and math.isfinite(value):
            out.append((aep, value))
    return sorted(out)


def group_envelope(study: Study, run_name: str, group: str, key: str) -> pd.Series:
    run = study.open_run(run_name)
    rows = run.rows_in(group)
    if not rows:
        raise StudyError(f"{run_name} has no group {group!r}")
    sources = results.sources_for_rows(run.frame, run.project_folder, rows).get(key, [])
    if not sources:
        raise StudyError(f"{group}: no {key} results on disk")
    comparison = results.compare(sources)
    if comparison.is_empty:
        raise StudyError(f"{group}: the {key} results could not be read")
    return results.analyse(comparison).envelope


def file_curve(path: Path, column: str) -> pd.Series:
    """A table with AEP ('1 in X') in its first column."""
    try:
        frame = pd.read_csv(path)
    except (OSError, ValueError) as exc:
        raise StudyError(f"{Path(path).name} could not be read ({exc})") from exc
    if column not in frame.columns:
        raise StudyError(f"{Path(path).name} has no {column!r} column "
                         f"(it has {', '.join(map(str, frame.columns[1:]))})")
    return pd.Series(pd.to_numeric(frame[column], errors="coerce").to_numpy(),
                     index=pd.to_numeric(frame.iloc[:, 0], errors="coerce"))


def read_ffa(path: Path) -> dict:
    """An RMC Bestfit export, each series as AEP '1 in X' -> value.

    Columns come in ``<name>_x``/``<name>_y`` pairs, x the exceedance
    probability. The paleoflood rows are censored - the flood reached *at least*
    ``Interval Data_yLower`` - so only that lower bound is kept, as v03 drew it.
    """
    try:
        frame = pd.read_csv(path)
    except (OSError, ValueError) as exc:
        raise StudyError(f"{Path(path).name} could not be read ({exc})") from exc

    def pair(x_column, y_column):
        if x_column not in frame.columns or y_column not in frame.columns:
            return pd.Series(dtype=float)
        part = frame[[x_column, y_column]].apply(pd.to_numeric, errors="coerce").dropna()
        part = part[part[x_column] > 0]
        return pd.Series(part[y_column].to_numpy(), index=1.0 / part[x_column].to_numpy())

    out = {"Mode": pair("Posterior Mode_x", "Posterior Mode_y"),
           "Mean": pair("Posterior Mean_x", "Posterior Mean_y"),
           "upper": pair("90% Credible Intervals_x", "90% Credible Intervals_y"),
           "lower": pair("90% Credible Intervals_x2", "90% Credible Intervals_y2"),
           "AMS": pair("Exact Data_x", "Exact Data_y"),
           "paleo": pair("Interval Data_x", "Interval Data_yLower")}
    if out["Mode"].empty and out["Mean"].empty and out["AMS"].empty:
        raise StudyError(f"{Path(path).name} does not look like an RMC Bestfit export "
                         f"(no 'Posterior Mode_x' or 'Exact Data_x' column)")
    return out


def build(study: Study, spec: dict) -> FigureData:
    data = FigureData(spec=spec)
    key = spec.get("type") or "level"
    ffa_count = sum(1 for curve in spec.get("curves") or [] if curve.get("kind") == FFA)
    for curve in spec.get("curves") or []:
        kind = curve.get("kind", GROUP)
        label = curve.get("label") or curve.get("group") or Path(str(curve.get("file", ""))).stem
        try:
            if kind == GROUP:
                envelope = group_envelope(study, curve.get("run", ""), curve.get("group", ""),
                                          key)
                data.series.append(Series(label, _trimmed(envelope, spec)))
            elif kind == FILE:
                path = resolve(study.folder, curve.get("file"))
                series = file_curve(path, curve.get("column") or key)
                data.series.append(Series(label, _trimmed(series, spec)))
            elif kind == FFA:
                _add_ffa(data, study, curve, label, spec, single=ffa_count == 1)
            else:
                data.problems.append(f"{label}: unknown curve kind {kind!r}")
        except StudyError as exc:
            data.problems.append(f"{label}: {exc}")
    return data


def _add_ffa(data, study, curve, label, spec, *, single) -> None:
    path = resolve(study.folder, curve.get("file"))
    if path is None:
        raise StudyError("no FFA file given")
    ffa = read_ffa(path)
    posterior = curve.get("posterior") or "Both"
    named = label or "FFA"
    data.series.append(Series(curve.get("ams_label") or "AMS",
                              _trimmed(ffa["AMS"], spec), style="ams"))
    if posterior != "AMS only":
        # the credible interval is drawn once, and only when there is one FFA:
        # two sets of grey bands on one plot cannot be told apart
        if single:
            data.series.append(Series("90% credible interval", _trimmed(ffa["upper"], spec),
                                      style="ffa-ci", legend=False))
            data.series.append(Series("90% credible interval", _trimmed(ffa["lower"], spec),
                                      style="ffa-ci", legend=False))
        for kind in ("Mode", "Mean"):
            if posterior in (kind, "Both") and not ffa[kind].empty:
                data.series.append(Series(f"{named}: Posterior {kind}",
                                          _trimmed(ffa[kind], spec),
                                          style=f"ffa-{kind.lower()}"))
    if not ffa["paleo"].empty:
        data.series.append(Series("Paleoflood (lower bound)", _trimmed(ffa["paleo"], spec),
                                  style="paleo"))


# -- the export ----------------------------------------------------------------------

def output_path(study: Study, spec: dict) -> Path:
    folder = resolve(study.folder, spec.get("folder")) or (study.folder / "figures")
    name = UNSAFE_NAME.sub("_", spec.get("filename") or spec.get("id") or "figure")
    return folder / f"{name}.png"


@dataclass
class ExportResult:
    png: Path
    returncode: int
    output: str = ""

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and self.png.is_file()


def export(data: FigureData, png: Path, python) -> ExportResult:
    """Draw the PNG with util/ReportFigure.py. Blocking - call it off the UI thread.

    The job is written beside the PNG as ``<name>.json``, so a figure in the
    report can always be traced back to the numbers and labels it was drawn from.
    """
    png = Path(png)
    png.parent.mkdir(parents=True, exist_ok=True)
    job_path = png.with_suffix(".json")
    atomic_write_json(job_path, data.to_job(png))
    environment = dict(os.environ)
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    try:
        finished = subprocess.run(
            [str(python), "-u", str(SCRIPT), str(job_path)], capture_output=True,
            text=True, timeout=TIMEOUT_SECONDS, env=environment,
            stdin=subprocess.DEVNULL, cwd=str(png.parent))
    except subprocess.TimeoutExpired:
        return ExportResult(png, 1, f"timed out after {TIMEOUT_SECONDS} s")
    except OSError as exc:
        return ExportResult(png, 1, f"could not start {SCRIPT.name}: {exc}")
    return ExportResult(png, finished.returncode,
                        (finished.stdout or "") + (finished.stderr or ""))


def job_json(data: FigureData, png: Path) -> str:
    return json.dumps(data.to_job(png), indent=2)


# -- what the page's text boxes hold ----------------------------------------------

def parse_reference_levels(text) -> list:
    """One 'label = level' or 'label = level = colour' per line."""
    out = []
    for line in str(text or "").splitlines():
        parts = [part.strip() for part in line.split("=")]
        if len(parts) < 2:
            continue
        try:
            level = float(parts[1])
        except ValueError:
            continue
        item = {"label": parts[0], "level": level}
        if len(parts) > 2 and parts[2]:
            item["colour"] = parts[2]
        out.append(item)
    return out


def reference_levels_text(levels) -> str:
    lines = []
    for item in levels or []:
        parts = (item.get("label", ""), item.get("level"), item.get("colour"))
        lines.append(" = ".join(str(part) for part in parts if part not in (None, "")))
    return "\n".join(lines)
