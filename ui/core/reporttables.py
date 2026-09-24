"""The design flood report's result tables, built from a study's runs.

Each table the report carries is a *kind* here, and a study file holds one spec
per table naming the run and group it is drawn from. Re-routing a revised rating
then means re-running the group and pressing Copy, not re-typing four hundred
cells - which is what the Callide report's tables were, and why half of them had
drifted from the results on disk by the time this was written (see below).

The kinds, and what each reproduces:

``design_floods``  Tables 1 and 26-31. One group. Per standard AEP: the level
                   envelope over the durations and its critical duration, with
                   the inflow and outflow **read at the level's critical
                   duration**, not at their own - that is what the report's
                   "* Peak inflow for lake level critical duration" says. A
                   duration that does not spill has no outflow quantile at that
                   AEP, and the table says 0, as the report does. Table 1 adds
                   two rows in AEP order: the dam crest flood (``dcf``: its
                   AEP by ``aep_of_level``, no flows, shaded) and the PMF
                   (``pmf``: the highest ensemble event at the notional AEP
                   adopted on the PMF page).
``flood_levels``   Table 32. The AEP at which the level envelope reaches each of
                   a list of levels (the dam crest, each embankment crest),
                   rounded to 10. By default from the critical duration's own
                   realisations, as the Callide DCF was estimated; see
                   ``aep_of_level``.
``peak_at_aep``    Table 33. Level, inflow and outflow at one AEP (the AEP of
                   the PMP) for several groups.
``ensemble_peak``  Table 34. The PMF: the largest level in an ensemble database,
                   with the inflow and outflow of the event that produced it
                   (``pick: highest``). ``pick: median`` takes Bryan's median
                   pattern instead, for any other ensemble result - see
                   core/ensemble.py.

The multi-group kinds hold ``sections`` of rows, each row either a {run, group}
or fixed ``values`` - the Sunwater 2020 baseline rows are a previous study's
numbers, not a run in this one.

**Where flows are read.** Tables 33 and 34 have one duration column headed
"Level critical duration", and the scripts that filled them took each result's
*own* maximum instead - the inflow's critical duration is shorter than the
level's, so the PMPF inflow in Table 33 was a different storm from the PMPF
level beside it, and not the value in Table 26's PMPF row. ``flows_at`` says
which: ``level`` (the default, consistent with the header and with Tables 1 and
26-31) or ``own`` (what the old scripts did).

Validated against E012 on 23 September 2026: Tables 26 and 28 reproduce the
report to every digit. Everywhere the report differs, it is the report that is
behind - Tables 29-31 predate the FSL re-run of 21 September, and Table 32 and 33
predate E012 altogether; E012's own ``_00h_`` critical-duration files agree with
this module, not with the report.

pandas and the standard library only; nothing here imports nicegui.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import pandas as pd

from . import ensemble, results
from .bryan import representative_events
from .study import Study, StudyError
from .wordtable import ReportTable

EVENTS = representative_events()

DESIGN_FLOODS = "design_floods"
FLOOD_LEVELS = "flood_levels"
PEAK_AT_AEP = "peak_at_aep"
ENSEMBLE_PEAK = "ensemble_peak"
REPRESENTATIVE = "representative_events"

AT_LEVEL = "level"
AT_OWN = "own"

# How a level becomes an AEP - see aep_of_level.
MCDF = "mcdf"
ENVELOPE = "envelope"

DASH = "–"
NO_VALUE = DASH

# The standard AEPs a design flood table carries, 1 in 5 up to the AEP of the
# PMP. Rows the run did not produce are reported, not invented.
DEFAULT_AEPS = [5, 10, 20, 50, 100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000,
                50_000, 100_000, 200_000, 500_000, 1_000_000]
DEFAULT_PMP_AEP = 1_900_000

INFLOW_FOOTNOTE = "* Peak inflow for lake level critical duration."


@dataclass(frozen=True)
class Kind:
    key: str
    label: str
    help: str
    multi: bool                  # sections of rows, or one source
    template: dict


KINDS = {
    DESIGN_FLOODS: Kind(
        DESIGN_FLOODS, "Design flood estimates",
        "One group: peak inflow, outflow and level by AEP, flows at the level's "
        "critical duration (report Tables 1, 26-31).",
        False,
        {"kind": DESIGN_FLOODS, "title": "", "source": {"run": "", "group": ""},
         "aeps": list(DEFAULT_AEPS), "pmp_aep": DEFAULT_PMP_AEP, "pmp_label": "PMPF",
         "thousands": True, "footnote": INFLOW_FOOTNOTE, "dcf": None, "pmf": None}),
    FLOOD_LEVELS: Kind(
        FLOOD_LEVELS, "AEP of given lake levels",
        "The AEP at which each group's level curve reaches the dam crest, an "
        "embankment crest or any other level (report Table 32).",
        True,
        {"kind": FLOOD_LEVELS, "title": "", "first_column": "Climate Horizon",
         "levels": [{"label": "DCF", "level": 219.13}], "method": MCDF,
         "round_to": 10, "duration": True, "sections": []}),
    PEAK_AT_AEP: Kind(
        PEAK_AT_AEP, "Peaks at one AEP",
        "Lake level, inflow and outflow at one AEP - the AEP of the PMP for the "
        "PMPF - for several groups (report Table 33).",
        True,
        {"kind": PEAK_AT_AEP, "title": "", "first_column": "Climate Horizon",
         "aep": DEFAULT_PMP_AEP, "flows_at": AT_LEVEL, "thousands": True,
         "sections": []}),
    ENSEMBLE_PEAK: Kind(
        ENSEMBLE_PEAK, "Ensemble peak (PMF)",
        "The largest lake level in an ensemble run's database, with that event's "
        "inflow and outflow (report Table 34).",
        True,
        {"kind": ENSEMBLE_PEAK, "title": "", "first_column": "Climate Horizon",
         "pick": "highest", "flows_at": AT_LEVEL, "thousands": True, "sections": []}),
    REPRESENTATIVE: Kind(
        REPRESENTATIVE, "Representative events",
        "The events chosen on the Events page, one section per group: each loading's "
        "AEP, lake level, trigger, simulation and duration, with the PMF's event "
        "(report Tables 35-36).",
        True,
        {"kind": REPRESENTATIVE, "title": "", "method": MCDF, "round_to": 10,
         "pmp_aep": DEFAULT_PMP_AEP, "pmp_label": "PMPF",
         "triggers": [{"label": "DCF", "level": 219.13}],
         "sections": []}),     # [{heading, run, group, pmf: {run, group} | None}]
}


# -- formatting ----------------------------------------------------------------

def _finite(value) -> bool:
    try:
        return not math.isnan(float(value))
    except (TypeError, ValueError):
        return False


def fmt_aep(aep) -> str:
    return results.format_aep(aep)


def fmt_flow(value, thousands=True) -> str:
    if not _finite(value):
        return NO_VALUE
    return f"{float(value):,.0f}" if thousands else f"{float(value):.0f}"


def fmt_level(value) -> str:
    return f"{float(value):.2f}" if _finite(value) else NO_VALUE


def fmt_duration(label) -> str:
    """'72h' -> '72', '4.5h' -> '4.5'. A label with no duration is kept."""
    if label is None or (isinstance(label, float) and math.isnan(label)):
        return NO_VALUE
    text = str(label).strip()
    number = text[:-1] if text.lower().endswith("h") else text
    try:
        return f"{float(number):g}"
    except ValueError:
        return text


def fmt_rounded_aep(aep, round_to=10) -> str:
    if not _finite(aep):
        return NO_VALUE
    step = float(round_to or 1)
    return f"{round(float(aep) / step) * step:,.0f}"


# -- reading a group -----------------------------------------------------------

@dataclass
class GroupCurves:
    """One group's quantile curves, compared across its durations."""

    run: str
    group: str
    comparisons: dict = field(default_factory=dict)   # result type -> Comparison
    analyses: dict = field(default_factory=dict)      # result type -> Analysis
    problems: list = field(default_factory=list)

    def envelope(self, key) -> pd.Series:
        analysis = self.analyses.get(key)
        return analysis.envelope if analysis is not None else pd.Series(dtype=float)

    def critical(self, key) -> pd.Series:
        analysis = self.analyses.get(key)
        return analysis.critical if analysis is not None else pd.Series(dtype=object)

    def at(self, key, aep, label) -> float:
        """One duration's value at one AEP; NaN when the duration has none there."""
        comparison = self.comparisons.get(key)
        if comparison is None or comparison.is_empty:
            return math.nan
        frame = comparison.frame
        if aep not in frame.index or label not in frame.columns:
            return math.nan
        value = frame.loc[aep, label]
        return float(value) if _finite(value) else math.nan

    def has(self, key) -> bool:
        comparison = self.comparisons.get(key)
        return comparison is not None and not comparison.is_empty


_CURVE_CACHE: dict = {}


def group_curves(study: Study, run_name: str, group: str) -> GroupCurves:
    """Read a group's level, inflow and outflow curves, once per change on disk."""
    curves = GroupCurves(run=run_name, group=group)
    if not run_name or not group:
        curves.problems.append("no run and group chosen")
        return curves
    try:
        run = study.open_run(run_name)
    except StudyError as exc:
        curves.problems.append(str(exc))
        return curves
    rows = run.rows_in(group)
    if not rows:
        curves.problems.append(f"{run_name} has no group {group!r}")
        return curves

    sources = results.sources_for_rows(run.frame, run.project_folder, rows)
    stamp = []
    for key in ("level", "inflow", "outflow"):
        for source in sources.get(key, []):
            try:
                stamp.append((str(source.path), source.path.stat().st_mtime_ns))
            except OSError:
                stamp.append((str(source.path), None))
    cache_key = (str(run.config.config_path), group, tuple(stamp))
    if cache_key in _CURVE_CACHE:
        return _CURVE_CACHE[cache_key]

    for key in ("level", "inflow", "outflow"):
        found = sources.get(key, [])
        if not found:
            curves.problems.append(f"{group}: no {key} results on disk")
            continue
        comparison = results.compare(found)
        for label, why in comparison.problems:
            curves.problems.append(f"{group} {key} {label}: {why}")
        curves.comparisons[key] = comparison
        curves.analyses[key] = results.analyse(comparison)
    _CURVE_CACHE[cache_key] = curves
    return curves


def forget_cached() -> None:
    _CURVE_CACHE.clear()


# -- the kinds -----------------------------------------------------------------

def build(study: Study, spec: dict) -> ReportTable:
    kind = spec.get("kind")
    builder = BUILDERS.get(kind)
    if builder is None:
        table = ReportTable(header=["-"], title=spec.get("title", ""))
        table.problems.append(f"unknown table kind {kind!r}")
        return table
    table = builder(study, spec)
    table.title = spec.get("title", "")
    # A group with no results is met once per row and once per level; say it once.
    table.problems = list(dict.fromkeys(table.problems))
    footnote = str(spec.get("footnote") or "").strip()
    if footnote:
        table.footnotes.append(footnote)
    return table


def _design_floods(study: Study, spec: dict) -> ReportTable:
    table = ReportTable(
        header=["AEP (1 in x)", "Peak inflow*\n(m3/s)", "Peak outflow (m3/s)",
                "Peak level (m AHD)", "Level critical duration (h)"],
        align=["right"] * 5)
    source = spec.get("source") or {}
    curves = group_curves(study, source.get("run", ""), source.get("group", ""))
    table.problems.extend(curves.problems)
    if not curves.has("level"):
        return table

    thousands = spec.get("thousands", True)
    pmp_aep = spec.get("pmp_aep")
    aeps = [float(aep) for aep in (spec.get("aeps") or DEFAULT_AEPS)]
    if _finite(pmp_aep) and float(pmp_aep) not in aeps:
        aeps.append(float(pmp_aep))

    rows = []                      # (AEP, cells, row options), sorted at the end
    envelope, critical = curves.envelope("level"), curves.critical("level")
    for aep in aeps:
        label = fmt_aep(aep)
        if _finite(pmp_aep) and aep == float(pmp_aep) and spec.get("pmp_label"):
            label = f"{label}\n({spec['pmp_label']})"
        if aep not in envelope.index or not _finite(envelope.get(aep)):
            rows.append((aep, [label, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE], {}))
            table.problems.append(f"no level quantile at 1 in {fmt_aep(aep)}")
            continue
        duration = critical.get(aep)
        inflow = curves.at("inflow", aep, duration)
        outflow = curves.at("outflow", aep, duration)
        if not _finite(outflow) and curves.has("outflow"):
            outflow = 0.0          # that duration did not spill at this AEP
        rows.append((aep, [label, fmt_flow(inflow, thousands), fmt_flow(outflow, thousands),
                           fmt_level(envelope[aep]), fmt_duration(duration)], {}))

    rows += _dcf_row(study, curves, spec.get("dcf"), table)
    rows += _pmf_row(study, spec.get("pmf"), thousands, table)
    for _, cells, options in sorted(rows, key=lambda item: item[0]):
        table.add(cells, **options)
    return table


def _dcf_row(study, curves, dcf, table) -> list:
    """Table 1's dam crest flood row: its AEP and level, no flows, shaded."""
    if not dcf or not _finite(dcf.get("level")):
        return []
    level = float(dcf["level"])
    aep, duration, problem = aep_of_level(study, curves, level, dcf.get("method", MCDF))
    if problem:
        table.problems.append(f"DCF ({level:.2f} m) {problem}")
    if not _finite(aep):
        return []
    label = f"{fmt_rounded_aep(aep, dcf.get('round_to', 10))} ({dcf.get('label') or 'DCF'})"
    return [(aep, [label, "", "", fmt_level(level), fmt_duration(duration)],
             {"bold": True, "shaded": True})]


def _pmf_row(study, pmf, thousands, table) -> list:
    """Table 1's PMF row: the highest ensemble event at the adopted notional AEP."""
    if not pmf or not pmf.get("group"):
        return []
    aep = pmf.get("aep") or ensemble.settings(study).get("adopted_aep")
    if not _finite(aep):
        table.problems.append("PMF row: no notional AEP - adopt one on the PMF page")
        return []
    try:
        frame = ensemble.load(study, pmf.get("run", ""), pmf["group"])
    except StudyError as exc:
        table.problems.append(f"PMF row: {exc}")
        return []
    chosen = ensemble.pick(frame, pmf.get("pick", ensemble.HIGHEST))
    if not chosen.found:
        table.problems.append("PMF row: no levels in the database")
        return []
    label = f"{fmt_aep(float(aep))} ({pmf.get('label') or 'PMF'})"
    return [(float(aep), [label, fmt_flow(chosen.inflow, thousands),
                          fmt_flow(chosen.outflow, thousands), fmt_level(chosen.level),
                          fmt_duration(chosen.duration)], {})]


def _rows(spec):
    """(section heading or None, row) in order, sections flattened."""
    for section in spec.get("sections") or []:
        yield section.get("heading") or "", None
        for row in section.get("rows") or []:
            yield None, row


def _fixed(row, width) -> list:
    values = [str(value) for value in (row.get("values") or [])]
    return (values + [""] * width)[:width]


def _flood_levels(study: Study, spec: dict) -> ReportTable:
    levels = [item for item in spec.get("levels") or []
              if _finite(item.get("level"))]
    header = [spec.get("first_column") or "Climate Horizon"]
    header += [item.get("label") or f"{float(item['level']):.2f} m" for item in levels]
    if spec.get("duration", True):
        header.append("Critical duration (h)")
    table = ReportTable(header=header, align=["left"] + ["center"] * (len(header) - 1))
    round_to = spec.get("round_to", 10)

    for heading, row in _rows(spec):
        if row is None:
            if heading:
                table.section(heading)
            continue
        label = row.get("label", "")
        if "values" in row:
            table.add([label] + _fixed(row, len(header) - 1))
            continue
        curves = group_curves(study, row.get("run", ""), row.get("group", ""))
        table.problems.extend(curves.problems)
        cells = [label]
        first_duration = None
        for item in levels:
            aep, duration, problem = aep_of_level(study, curves, item["level"],
                                                  spec.get("method", MCDF))
            cells.append(fmt_rounded_aep(aep, round_to))
            if problem:
                table.problems.append(f"{label}: {item.get('label') or item['level']} "
                                      f"({float(item['level']):.2f} m) {problem}")
            if first_duration is None:
                first_duration = duration
        if spec.get("duration", True):
            cells.append(fmt_duration(first_duration))
        table.add(cells)
    return table


def aep_of_level(study: Study, curves: GroupCurves, level: float, method: str = MCDF):
    """(AEP, critical duration, problem) of one lake level in one group.

    ``envelope`` reads the level off the design curve - the envelope over the
    durations, linear in (log level, z) between the standard AEPs
    (util/DesignFloodInterpolation.py). ``mcdf`` - the default, and how the
    Callide DCF was estimated - takes the **critical duration** from that curve
    and then reads the level off that duration's own realisations
    (``AEPofDCF_v2.py``): thousands of events rather than a straight line between
    two standard AEPs a decade apart. On E012 near-term RFSL the two give 1 in
    25,350 and 1 in 27,230 for the dam crest. Where the duration has no database
    on disk the envelope answer is given, and the problem says so.
    """
    lookup = EVENTS.aep_for_level(curves.envelope("level"), level)
    if not lookup.found:
        where = ("is above the top of the level curve - rarer than anything the run "
                 "produced" if lookup.above_curve
                 else "is below the bottom of the level curve" if "below" in lookup.note
                 else f"- {lookup.note or 'not on the level curve'}")
        return math.nan, None, where
    duration = _critical_near(curves, lookup.aep)
    if method == ENVELOPE:
        return lookup.aep, duration, ""
    try:
        databases = ensemble.mc_databases(study, curves.run, curves.group)
    except StudyError as exc:
        return lookup.aep, duration, f"read off the curve: {exc}"
    path = databases.get(str(duration))
    if path is None:
        return lookup.aep, duration, (f"read off the curve - no {duration} results "
                                      f"database to interpolate")
    try:
        real = ensemble.read_realisations(path)
    except (StudyError, OSError) as exc:
        return lookup.aep, duration, f"read off the curve: {exc}"
    aep = ensemble.aep_at_value(real, level)
    if not _finite(aep):
        return lookup.aep, duration, (f"read off the curve - the {duration} "
                                      f"realisations do not reach it")
    return aep, duration, ""


def _critical_near(curves: GroupCurves, aep):
    """The level's critical duration at the standard AEP nearest ``aep``."""
    critical = curves.critical("level").dropna()
    if critical.empty or not _finite(aep):
        return None
    z = EVENTS.normal_variate(aep)
    nearest = min(critical.index, key=lambda each: abs(EVENTS.normal_variate(each) - z))
    return critical.loc[nearest]


def _peak_header(spec, level_title="Lake level\n(m AHD)", inflow="Inflow (m3/s)",
                 outflow="Outflow (m3/s)"):
    duration = ("Level critical duration (h)" if spec.get("flows_at", AT_LEVEL) == AT_LEVEL
                else "Critical duration (h)")
    return [spec.get("first_column") or "Climate Horizon", level_title, inflow,
            outflow, duration]


def _peak_at_aep(study: Study, spec: dict) -> ReportTable:
    header = _peak_header(spec)
    table = ReportTable(header=header, align=["left"] + ["center"] * 4)
    thousands = spec.get("thousands", True)
    try:
        aep = float(spec.get("aep") or DEFAULT_PMP_AEP)
    except (TypeError, ValueError):
        table.problems.append(f"AEP {spec.get('aep')!r} is not a number")
        return table
    own = spec.get("flows_at", AT_LEVEL) == AT_OWN

    for heading, row in _rows(spec):
        if row is None:
            if heading:
                table.section(heading)
            continue
        label = row.get("label", "")
        if "values" in row:
            table.add([label] + _fixed(row, 4))
            continue
        curves = group_curves(study, row.get("run", ""), row.get("group", ""))
        table.problems.extend(curves.problems)
        envelope = curves.envelope("level")
        if aep not in envelope.index:
            table.add([label, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE])
            table.problems.append(f"{label}: no level quantile at 1 in {fmt_aep(aep)}")
            continue
        duration = curves.critical("level").get(aep)
        if own:
            inflow = curves.envelope("inflow").get(aep, math.nan)
            outflow = curves.envelope("outflow").get(aep, math.nan)
        else:
            inflow = curves.at("inflow", aep, duration)
            outflow = curves.at("outflow", aep, duration)
        if not _finite(outflow) and curves.has("outflow"):
            outflow = 0.0
        table.add([label, fmt_level(envelope[aep]), fmt_flow(inflow, thousands),
                   fmt_flow(outflow, thousands), fmt_duration(duration)])
    return table


def _ensemble_peak(study: Study, spec: dict) -> ReportTable:
    header = _peak_header(spec, "Peak lake level\n(m AHD)", "Peak dam inflow (m3/s)",
                          "Peak dam outflow (m3/s)")
    table = ReportTable(header=header, align=["left"] + ["center"] * 4)
    thousands = spec.get("thousands", True)
    own = spec.get("flows_at", AT_LEVEL) == AT_OWN
    convention = spec.get("pick", ensemble.HIGHEST)

    for heading, row in _rows(spec):
        if row is None:
            if heading:
                table.section(heading)
            continue
        label = row.get("label", "")
        if "values" in row:
            table.add([label] + _fixed(row, 4))
            continue
        try:
            frame = ensemble.load(study, row.get("run", ""), row.get("group", ""))
        except StudyError as exc:
            table.add([label, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE])
            table.problems.append(f"{label}: {exc}")
            continue
        chosen = ensemble.pick(frame, convention)
        if not chosen.found:
            table.add([label, NO_VALUE, NO_VALUE, NO_VALUE, NO_VALUE])
            table.problems.append(f"{label}: no levels in the database")
            continue
        inflow = frame["inflow"].max() if own else chosen.inflow
        table.add([label, fmt_level(chosen.level), fmt_flow(inflow, thousands),
                   fmt_flow(chosen.outflow, thousands), fmt_duration(chosen.duration)])
    return table


def selection_for(study: Study, run_name: str, group: str):
    """The Events page's saved selection for a group: (path, targets)."""
    from .events import load_targets, selection_path      # local: events is heavy
    from .outputs import find_database

    run = study.open_run(run_name)
    rows = run.rows_in(group)
    if not rows:
        raise StudyError(f"{run_name} has no group {group!r}")
    for index in rows:
        database = find_database(run.frame.loc[index], run.project_folder)
        if database is not None:
            path = selection_path(database.parent, group)
            if not path.is_file():
                raise StudyError(f"{group}: no events chosen yet - pick them on the "
                                 f"Events page and press Save ({path.name})")
            targets, _ = load_targets(path)
            return path, targets
    raise StudyError(f"{group}: no results database on disk")


def _event_duration(run, target) -> str:
    """The duration of the run the event was taken from, '72h'."""
    wanted = str(target.output_file or "").replace("\\", "/").rsplit("/", 1)[-1]
    for index in run.frame.index:
        row = run.frame.loc[index]
        name = str(row.get("Output file") or "").replace("\\", "/").rsplit("/", 1)[-1]
        if name and name == wanted:
            duration = results.duration_of(row, name)
            return f"{duration:g}h" if duration is not None else ""
    return ""


def _representative(study: Study, spec: dict) -> ReportTable:
    table = ReportTable(header=["AEP\n(1 in Y)", "Lake Level (m AHD)", "Trigger",
                                "Representative Event", "Critical duration"],
                        align=["center"] * 5)
    pmp_aep = spec.get("pmp_aep")
    triggers = {round(float(item["level"]), 3): item.get("label", "")
                for item in spec.get("triggers") or [] if _finite(item.get("level"))}
    for section in spec.get("sections") or []:
        if section.get("heading"):
            table.section(section["heading"])
        run_name, group = section.get("run", ""), section.get("group", "")
        try:
            _, targets = selection_for(study, run_name, group)
            run = study.open_run(run_name)
        except StudyError as exc:
            table.problems.append(str(exc))
            targets, run = [], None
        curves = group_curves(study, run_name, group) if run is not None else None
        rows = []
        for target in targets:
            if target.picked is None:
                table.problems.append(f"{group}: no event picked for {target.label}")
                continue
            duration = _event_duration(run, target)
            if target.kind == "level":
                level = float(target.value)
                aep, _, problem = aep_of_level(study, curves, level,
                                               spec.get("method", MCDF))
                if problem:
                    table.problems.append(f"{group}: {level:.2f} m {problem}")
                trigger = (target.comment or "").strip() or triggers.get(round(level, 3), "")
                rows.append((aep if _finite(aep) else math.inf,
                             [fmt_rounded_aep(aep, spec.get("round_to", 10)),
                              fmt_level(level), trigger, str(target.picked), duration]))
            else:
                aep = float(target.value)
                level = EVENTS.value_for_aep(curves.envelope("level"), aep)
                label = (spec.get("pmp_label") or "PMPF"
                         if _finite(pmp_aep) and aep == float(pmp_aep) else fmt_aep(aep))
                rows.append((aep, [label, fmt_level(level), "", str(target.picked),
                                   duration]))
        pmf = section.get("pmf")
        if pmf and pmf.get("group"):
            try:
                chosen = ensemble.pick(ensemble.load(study, pmf.get("run", ""),
                                                     pmf["group"]), ensemble.HIGHEST)
            except StudyError as exc:
                table.problems.append(f"PMF: {exc}")
                chosen = None
            if chosen is not None and chosen.found:
                rows.append((math.inf, [pmf.get("label") or "PMF", fmt_level(chosen.level),
                                        "", str(chosen.event),
                                        f"{chosen.duration:g}h"]))
        for _, cells in sorted(rows, key=lambda item: item[0]):
            table.add(cells)
    return table


BUILDERS = {
    DESIGN_FLOODS: _design_floods,
    FLOOD_LEVELS: _flood_levels,
    PEAK_AT_AEP: _peak_at_aep,
    ENSEMBLE_PEAK: _ensemble_peak,
    REPRESENTATIVE: _representative,
}


def new_spec(kind: str) -> dict:
    import copy
    return copy.deepcopy(KINDS[kind].template)


# -- what the page's text boxes hold -------------------------------------------

def parse_aeps(text) -> list:
    """'5, 10, 1,000' is ambiguous, so AEPs are one per line or space-separated;
    thousands separators inside a number are allowed ('1,000')."""
    out = []
    for token in str(text or "").replace(";", " ").split():
        try:
            value = float(token.replace(",", "").replace("_", ""))
        except ValueError:
            continue
        if value > 1:
            out.append(int(value) if value.is_integer() else value)
    return out


def aeps_text(aeps) -> str:
    return " ".join(f"{aep:,.0f}" if float(aep).is_integer() else f"{aep:g}"
                    for aep in aeps or [])


def parse_levels(text) -> list:
    """One 'label = level' per line, as the Lake levels page takes them."""
    out = []
    for line in str(text or "").splitlines():
        if "=" not in line:
            continue
        label, _, value = line.rpartition("=")
        try:
            out.append({"label": label.strip(), "level": float(value.strip())})
        except ValueError:
            continue
    return out


def levels_text(levels) -> str:
    return "\n".join(f"{item.get('label', '')} = {item.get('level')}" for item in levels or [])


def parse_values(text) -> list:
    """A fixed row's cells, separated by '|' so a value may hold a comma."""
    return [cell.strip() for cell in str(text or "").split("|")]
