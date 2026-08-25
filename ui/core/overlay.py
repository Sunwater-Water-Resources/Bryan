"""Overlaying one group's envelope against another's.

The Durations tab asks a question *inside* one group: which storm duration is
critical. This asks one *between* groups: what does a warmer climate, or a
raised full supply level, or a different antecedent storage distribution do to
the design flood. The curve that carries that answer is the **envelope** - the
maximum over the durations, which is the design quantile
(``util/MaxQuantiles.py:96``, and ``results.analyse`` does the same thing) - so
a group contributes exactly one line here however many durations it ran.

Two decisions shape the module:

- **The overlay reuses ``results.Comparison``'s shape**, columns being groups
  instead of durations. The AEP index, the standard normal variate axis, the
  range trim and the table renderer then work unaltered.
- **What it must *not* reuse is the mark-up.** The maximum across groups is
  meaningless: groups are scenarios, not alternatives to be enveloped, so there
  is no envelope-of-envelopes, no critical-group bands and no crossover pins.
  The comparison people actually want is against a **baseline** group, and that
  is ``deltas``.

Deltas carry the same unit convention as the critical-duration margins, from
the same ``results.margin_scale``: percent for flows and volumes, **metres for
level**. A climate uplift of 0.4 m on a lake level is 0.06% of a number on an
arbitrary datum, which is not a quantity anyone can act on.

Overlaying envelopes hides how each one was built, so ``notes`` says out loud
what the picture cannot: envelopes over different numbers of durations, a group
whose envelope is pinned to the end of its own duration range (a lower bound,
so a difference measured against it is biased), a group mixing two methods, and
groups that do not reach the same AEP.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

from . import results
from .columns import normalise_method

# Where a group key may be cut. Keys are built by ``grouping.group_key`` out of
# name parts joined with '|', on top of whatever the study's own naming
# convention uses, so the label is trimmed on those boundaries and never
# mid-token: 'GWL1p3' and 'GWL1p7' must not come out as '3' and '7'.
SEPARATOR = re.compile(r"([_|\-\s/\\]+)")
TRIM = "_|-/\\ "


@dataclass(frozen=True)
class GroupCurve:
    """One group's envelope, and enough provenance to say what it is."""

    key: str                 # the group key, as grouping derived it
    label: str               # the trimmed label, distinguishing it from the rest
    envelope: pd.Series
    critical: pd.Series      # which duration owns the envelope, per AEP
    durations: dict          # duration label -> hours
    sources: tuple = ()      # CurveSource, for 'Files read'
    pinned: tuple = ()       # this group's pinned-end warnings
    methods: tuple = ()      # the distinct Method values that contributed

    @property
    def duration_count(self) -> int:
        return len(self.durations)


@dataclass
class Overlay:
    """Every selected group's envelope, on one AEP index."""

    frame: pd.DataFrame = field(default_factory=pd.DataFrame)
    curves: tuple = ()
    key: str = ""
    problems: list = field(default_factory=list)   # (label, why) - unreadable
    notes: list = field(default_factory=list)      # about the comparison itself

    @property
    def is_empty(self) -> bool:
        return self.frame.empty

    @property
    def labels(self) -> list:
        return list(self.frame.columns)

    def curve(self, label):
        for curve in self.curves:
            if curve.label == label:
                return curve
        return None


@dataclass
class Deltas:
    """Every group against the baseline, in the units of the result type."""

    frame: pd.DataFrame = field(default_factory=pd.DataFrame)
    baseline: str = ""
    kind: str = results.PERCENT
    label: str = "change %"
    undefined: tuple = ()          # AEPs where the baseline is zero or missing

    @property
    def is_empty(self) -> bool:
        return self.frame.empty


# -- labelling ---------------------------------------------------------------

def distinguishing_labels(keys) -> dict:
    """Group keys trimmed to the part that tells them apart.

    ``TFD_mc_C030_ebf_L20-4_GWL1p3|exg`` is not a legend entry. What
    distinguishes it from its neighbour is ``GWL1p3``, so the common leading
    and trailing parts come off - on separator boundaries, so a shared
    ``GWL1p`` can never be trimmed down to ``3``.

    Falls back to the full keys whenever trimming would leave nothing, which is
    what happens when one key is a prefix of another.
    """
    keys = list(keys)
    if len(keys) < 2:
        return {key: key for key in keys}

    parts = [SEPARATOR.split(key) for key in keys]
    head = _common_run(parts)
    tail = _common_run([list(reversed(items)) for items in parts])
    if not head and not tail:
        return {key: key for key in keys}

    labels = {}
    for key, items in zip(keys, parts):
        if head + tail >= len(items):
            return {key: key for key in keys}
        labels[key] = "".join(items[head:len(items) - tail]).strip(TRIM)
    if any(not label for label in labels.values()):
        return {key: key for key in keys}
    return labels


def _common_run(parts) -> int:
    """How many leading items every list shares, cut to a token boundary.

    ``re.split`` with a capturing group alternates token, separator, token, so
    an even count always ends on a boundary - which is the whole point.
    """
    shortest = min(len(items) for items in parts)
    count = 0
    while count < shortest and len({items[count] for items in parts}) == 1:
        count += 1
    return count - (count % 2)


# -- building ----------------------------------------------------------------

def result_keys(available, group_keys=None) -> list:
    """Every result type any of these groups offers, in the page's order."""
    keys = []
    for group, found in available.items():
        if group_keys is not None and group not in group_keys:
            continue
        for key in found:
            if key not in keys:
                keys.append(key)
    return keys


def groups_with(available, result_key) -> list:
    """The groups that have this result type, in the order they were scanned."""
    return [group for group, found in available.items() if found.get(result_key)]


def build(available, group_keys, result_key, sims_frame=None,
          labels=None) -> Overlay:
    """One envelope per group, read off the files that are already on disk.

    ``labels`` overrides the trimming. The page passes the labels it computed
    over *every* group offering this result type, so that ticking a group on
    and off does not rename the rest of them mid-comparison.
    """
    group_keys = [key for key in group_keys if key in available]
    labels = dict(labels) if labels else distinguishing_labels(group_keys)

    columns, curves, problems = [], [], []
    for group in group_keys:
        label = labels.get(group, group)
        sources = available.get(group, {}).get(result_key, [])
        if not sources:
            problems.append((label, f"has no {result_key} results"))
            continue
        comparison = results.compare(sources)
        problems.extend((f"{label} - {who}", why)
                        for who, why in comparison.problems)
        if comparison.is_empty:
            continue
        analysis = results.analyse(comparison)
        columns.append(analysis.envelope.rename(label))
        curves.append(GroupCurve(
            key=group, label=label, envelope=analysis.envelope,
            critical=analysis.critical, durations=dict(comparison.durations),
            sources=tuple(sources),
            pinned=tuple(results.pinned_end_warnings(comparison,
                                                     analysis.critical)),
            methods=_methods(sources, sims_frame)))

    if not columns:
        return Overlay(key=result_key, problems=problems)

    frame = pd.concat(columns, axis=1).sort_index()
    # An AEP no group reached is not plottable, and the outer join is
    # deliberate: groups genuinely stop at different AEPs, because the AEP of
    # the PMP comes from the storm config.
    frame = frame.dropna(axis=0, how="all")
    frame.index.name = results.AEP_COLUMN
    return Overlay(frame=frame, curves=tuple(curves), key=result_key,
                   problems=problems, notes=_notes(curves))


def _methods(sources, sims_frame) -> tuple:
    """The distinct Method values behind a group's curves."""
    if sims_frame is None or "Method" not in getattr(sims_frame, "columns", []):
        return ()
    seen = []
    for source in sources:
        if source.row_index is None or source.row_index not in sims_frame.index:
            continue
        method = normalise_method(sims_frame.loc[source.row_index].get("Method"))
        if method and method not in seen:
            seen.append(method)
    return tuple(seen)


def _notes(curves) -> list:
    """What the overlay cannot show, said out loud."""
    notes = []
    if len(curves) > 1:
        counts = {curve.label: curve.duration_count for curve in curves}
        if len(set(counts.values())) > 1:
            notes.append(
                "These envelopes are not built from the same number of "
                "durations ("
                + ", ".join(f"{label}: {count}" for label, count in counts.items())
                + "). An envelope over fewer durations is a lower bound, so the "
                  "comparison is biased towards the better covered group.")

        reach = {curve.label: _last_aep(curve.envelope) for curve in curves}
        known = {value for value in reach.values() if value == value}
        if len(known) > 1:
            notes.append(
                "The groups do not reach the same AEP ("
                + ", ".join(f"{label}: 1 in {results.format_aep(aep)}"
                            for label, aep in reach.items())
                + "), so the curves end at different places.")

    for curve in curves:
        if len(curve.methods) > 1:
            notes.append(
                f"{curve.label} draws on more than one method "
                f"({', '.join(curve.methods)}), so its envelope mixes them.")
    for curve in curves:
        notes.extend(f"{curve.label}: {warning}" for warning in curve.pinned)
    return notes


def _last_aep(envelope: pd.Series):
    finite = envelope.dropna()
    return float(finite.index[-1]) if len(finite) else float("nan")


# -- comparing against a baseline --------------------------------------------

def deltas(overlay: Overlay, baseline: str) -> Deltas:
    """Every other group's envelope as a change from the baseline's.

    Percent for flows and volumes, metres for level - ``results.margin_scale``,
    the same convention the critical-duration margins use.
    """
    kind, _floor, _heading = results.margin_scale(overlay.key)
    heading = "change (m)" if kind == results.ABSOLUTE else "change %"
    if overlay.is_empty or baseline not in overlay.frame.columns:
        return Deltas(baseline=baseline, kind=kind, label=heading)

    base = overlay.frame[baseline]
    others = overlay.frame.drop(columns=[baseline])
    if others.empty:
        return Deltas(baseline=baseline, kind=kind, label=heading)

    difference = others.sub(base, axis=0)
    if kind == results.PERCENT:
        # A dam that does not spill has an outflow quantile of zero at the
        # frequent end, and a percentage change from zero is not a number. The
        # AEPs where that happens are reported rather than drawn as infinity.
        difference = difference.div(base.abs(), axis=0) * 100.0
        difference = difference.replace([float("inf"), float("-inf")], float("nan"))
        undefined = tuple(aep for aep, value in base.items()
                          if not value == value or value == 0)
    else:
        undefined = tuple(aep for aep, value in base.items() if not value == value)

    return Deltas(frame=difference, baseline=baseline, kind=kind, label=heading,
                  undefined=undefined)


# -- the critical duration, per group ----------------------------------------

def critical_frame(overlay: Overlay) -> tuple:
    """(AEP x group frame of critical durations in hours, groups left out).

    A group whose durations are not all known contributes nothing rather than a
    line that cannot be plotted - the same guard
    ``resultchart.critical_duration_chart`` already applies to one group.
    """
    columns, skipped = [], []
    for curve in overlay.curves:
        if curve.critical.empty or any(curve.durations.get(owner) is None
                                       for owner in curve.critical):
            skipped.append(curve.label)
            continue
        hours = curve.critical.map(
            lambda owner, table=curve.durations: table[owner])
        columns.append(pd.to_numeric(hours).rename(curve.label))
    if not columns:
        return pd.DataFrame(), skipped
    frame = pd.concat(columns, axis=1).sort_index()
    frame.index.name = results.AEP_COLUMN
    return frame, skipped


def table(overlay: Overlay, deltas_result: Deltas = None) -> pd.DataFrame:
    """The on-screen table: one column per group, then the changes."""
    if overlay.is_empty:
        return pd.DataFrame()
    frame = overlay.frame.copy()
    if deltas_result is not None and not deltas_result.is_empty:
        for column in deltas_result.frame.columns:
            frame[f"{column} {deltas_result.label}"] = deltas_result.frame[column]
    return frame
