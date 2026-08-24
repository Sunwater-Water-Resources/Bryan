"""Comparing frequency curves across storm durations.

The question this answers is whether the durations that were run bracket the
critical one, so the plot on ``/results`` can say 'run longer storms' rather
than leaving it to be read off a picture.

Two things about critical durations are worth stating, because they shape what
counts as a finding here:

- **The direction is not the same for every result type.** For lake level the
  critical duration tends to be *long* at frequent AEPs, because it takes
  rainfall volume to fill and charge the storage, and *shorter* on the rare
  tail as the dam starts behaving as a conveyance rather than a volume system.
  Peak inflow follows catchment response instead and does no such thing. So
  nothing here assumes a direction: every AEP is checked independently for
  whether its critical duration is pinned to either end of the range that was
  actually run.
- **``idxmax`` hops on noise.** Where the duration curves are nearly coincident
  - the rare tail of a level curve especially, once everything is spilling -
  the winning duration changes on Monte Carlo sampling noise alone, and that
  looks exactly like a real crossover. Every critical duration is reported with
  the margin by which it beat the runner-up, and a switch inside
  ``MARGIN_NOISE_PERCENT`` is called out as noise rather than marked up.

No scipy: ``statistics.NormalDist`` covers the standard normal variate over the
dozen or so standard AEPs, which keeps this importable in the UI environment.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist

import pandas as pd

from .grouping import DURATION_TOKEN_RE
from .outputs import quantile_files
from .paths import cell_text, normalise_sep

AEP_COLUMN = "aep (1 in x)"

# Below this, the 'critical' duration beat the runner-up by less than the Monte
# Carlo noise floor, so the switch is not evidence of anything.
#
# The UNIT matters, and it is not the same for every result type. Flows and
# volumes are ratio scales - zero means no flow - so a percentage of the
# runner-up is the natural measure. **Lake level is not**: it is an interval
# scale on an arbitrary datum, so a percentage of a number like 217 m AHD says
# nothing. At Callide the durations separate by 0.01-0.10 m across the whole
# frequency range, which is 0.005-0.05% - every crossover would be dismissed as
# noise, on the one result type where the critical duration is most interesting.
# So level is measured in metres, against a floor in metres.
MARGIN_NOISE_PERCENT = 1.0

# A rule of thumb, not a physical constant - roughly the level difference worth
# acting on, and well inside the freeboard decisions these curves feed. The
# page lets it be changed per plot.
MARGIN_NOISE_METRES = 0.05

PERCENT = "percent"
ABSOLUTE = "absolute"

# Level is plotted linear, everything else logarithmic - as
# UtilModule.plot_durations does it.
LINEAR_TYPES = ("level",)

Y_LABELS = {
    "inflow": "Peak inflow (m³/s)",
    "outflow": "Peak outflow (m³/s)",
    "level": "Peak lake level (m AHD)",
}


def margin_scale(key) -> tuple:
    """(kind, default floor, column heading) for a result type's margin."""
    if key in LINEAR_TYPES:
        return ABSOLUTE, MARGIN_NOISE_METRES, "margin (m)"
    return PERCENT, MARGIN_NOISE_PERCENT, "margin %"


def format_margin(value, kind) -> str:
    if value != value:                      # NaN - only one duration reaches here
        return "-"
    return f"{value:.3f} m" if kind == ABSOLUTE else f"{value:.2f}%"


@dataclass(frozen=True)
class CurveSource:
    """One duration's quantile file, and how to label it."""

    label: str
    duration: float | None
    path: Path
    column: str
    key: str
    output_name: str = ""
    row_index: object = None

    @property
    def sort_key(self) -> tuple:
        return (0, self.duration) if self.duration is not None else (1, 0.0)


@dataclass
class Comparison:
    """Every selected duration's curve, on one AEP index."""

    frame: pd.DataFrame = field(default_factory=pd.DataFrame)
    durations: dict = field(default_factory=dict)
    problems: list = field(default_factory=list)   # (label, why)
    key: str = ""                                  # which result type these are

    @property
    def labels(self) -> list:
        return list(self.frame.columns)

    @property
    def is_empty(self) -> bool:
        return self.frame.empty


@dataclass(frozen=True)
class Switch:
    aep: float
    before: str
    after: str
    strength: float        # see below - NOT the margin at the switch itself


@dataclass(frozen=True)
class Band:
    """A run of AEPs owned by one duration - the shading under the envelope.

    ``peak_margin`` is the most the owner beat the runner-up by anywhere in the
    band, and it is what says whether the band is real. The margin *at* a
    crossover is near zero by definition - the curves are equal there - so it
    tells you nothing about whether the crossover means anything. What
    separates a real crossover from an idxmax hop is whether the incoming
    duration then goes on to win by something worth having.
    """

    label: str
    z_from: float
    z_to: float
    aep_from: float
    aep_to: float
    peak_margin: float = float("nan")


@dataclass
class Analysis:
    envelope: pd.Series = field(default_factory=pd.Series)
    critical: pd.Series = field(default_factory=pd.Series)
    margin: pd.Series = field(default_factory=pd.Series)
    switches: list = field(default_factory=list)
    bands: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    never_critical: list = field(default_factory=list)
    margin_kind: str = PERCENT
    margin_label: str = "margin %"
    noise_floor: float = MARGIN_NOISE_PERCENT


# -- the standard normal variate axis ----------------------------------------

def normal_variate(aep) -> float:
    """z for an AEP expressed as '1 in X'. ndtri(1 - 1/aep), without scipy."""
    try:
        aep = float(aep)
    except (TypeError, ValueError):
        return math.nan
    if not aep > 1:
        return math.nan
    return NormalDist().inv_cdf(1.0 - 1.0 / aep)


def format_aep(aep) -> str:
    try:
        aep = float(aep)
    except (TypeError, ValueError):
        return str(aep)
    if aep >= 1000:
        return f"{aep:,.0f}"
    if abs(aep - round(aep)) < 1e-9:
        return f"{aep:.0f}"
    return f"{aep:.4g}"


def y_axis(key: str) -> tuple:
    """(label, log?) for a result key."""
    if key in Y_LABELS:
        return Y_LABELS[key], key not in LINEAR_TYPES
    match = re.match(r"^(inflow|outflow)Vol(\d+(?:[._]\d+)?)h$", key, re.IGNORECASE)
    if match:
        window = match.group(2).replace("_", ".")
        return f"{match.group(1).capitalize()} volume, {window} h window (ML)", True
    return key, True


# -- finding the files -------------------------------------------------------

def _duration_of(row, output_name: str):
    value = row.get("Duration")
    try:
        number = float(value)
        if not pd.isna(number):
            return number
    except (TypeError, ValueError):
        pass
    # A row whose Duration is blank may still carry it in the output name, the
    # way grouping.strip_duration_token finds it.
    match = DURATION_TOKEN_RE.search(output_name)
    if match:
        try:
            return float(match.group(1).replace("_", "."))
        except ValueError:
            return None
    return None


def sources_for_rows(frame, project_folder, rows=None) -> dict:
    """Every plottable curve, keyed by result key ('inflow', 'inflowVol24h').

    Rows that have not been analysed simply contribute nothing - the files are
    probed, never assumed.
    """
    if frame is None or len(frame) == 0:
        return {}
    indices = list(frame.index) if rows is None else [i for i in rows if i in frame.index]

    found: dict = {}
    for index in indices:
        row = frame.loc[index]
        output_name = cell_text(row.get("Output file"))
        # Sims lists hold Windows paths, and 'sims_mc\\results\\X'.name is the
        # whole string on POSIX - normalise before taking the basename, or the
        # labels and the exported filename carry the folders with them.
        display_name = Path(normalise_sep(output_name)).name if output_name else ""
        duration = _duration_of(row, output_name)
        for key, quantile in quantile_files(row, project_folder).items():
            found.setdefault(key, []).append(CurveSource(
                label=f"{duration:g}h" if duration is not None else (
                    display_name or f"row {index + 2}"),
                duration=duration,
                path=quantile.path,
                column=quantile.column,
                key=key,
                output_name=display_name,
                row_index=index,
            ))
    return {key: _uniquely_labelled(sorted(sources, key=lambda s: s.sort_key))
            for key, sources in found.items()}


def _uniquely_labelled(sources: list) -> list:
    """Two rows of the same duration must not both be called '24h'.

    It happens whenever one duration is run twice under different options, and
    a duplicate label would silently drop a curve when the frame is built.
    """
    counts: dict = {}
    for source in sources:
        counts[source.label] = counts.get(source.label, 0) + 1

    resolved, used = [], set()
    for source in sources:
        label = source.label
        if counts[label] > 1:
            extra = source.output_name or str(source.path.stem)
            label = f"{label} ({extra})"
        while label in used:
            label = f"{label}'"
        used.add(label)
        resolved.append(CurveSource(
            label=label, duration=source.duration, path=source.path,
            column=source.column, key=source.key,
            output_name=source.output_name, row_index=source.row_index))
    return resolved


# -- reading and comparing ---------------------------------------------------

def read_curve(source: CurveSource) -> pd.Series:
    """One quantile file as a series indexed by AEP.

    Monte carlo (lib/MCScheme.py:353) and reservoir routing
    (lib/ReservoirRouting.py:310) write the same three columns, so one reader
    does both. Raises ValueError with something worth showing the user.
    """
    try:
        frame = pd.read_csv(source.path)
    except OSError as exc:
        raise ValueError(f"could not be opened ({exc.strerror or exc})") from exc
    except Exception as exc:                     # noqa: BLE001 - report, not crash
        raise ValueError(f"could not be read ({exc})") from exc

    if AEP_COLUMN not in frame.columns:
        raise ValueError(f"has no '{AEP_COLUMN}' column")
    if source.column not in frame.columns:
        raise ValueError(f"has no '{source.column}' column")

    series = pd.to_numeric(frame[source.column], errors="coerce")
    series.index = pd.to_numeric(frame[AEP_COLUMN], errors="coerce")
    series = series[series.index.notna()]
    series = series[~series.index.duplicated(keep="first")]
    series.name = source.label
    return series.sort_index()


def compare(sources) -> Comparison:
    """Read the selected curves onto one AEP index.

    The AEP sets can genuinely differ - ``get_standard_aeps`` extends to the
    AEP of the PMP, which comes from the storm config - so this is an outer
    join and the gaps stay NaN rather than being interpolated over.
    """
    columns, durations, problems = [], {}, []
    for source in sources:
        try:
            columns.append(read_curve(source))
        except ValueError as exc:
            problems.append((source.label, f"{source.path.name} {exc}"))
            continue
        durations[source.label] = source.duration

    if not columns:
        return Comparison(problems=problems, key=sources[0].key if sources else "")
    frame = pd.concat(columns, axis=1).sort_index()
    # An AEP no duration reached is not plottable, and pandas 2.x raises on an
    # all-NA row in idxmax rather than returning NaN.
    frame = frame.dropna(axis=0, how="all")
    frame.index.name = AEP_COLUMN
    return Comparison(frame=frame, durations=durations, problems=problems,
                      key=sources[0].key if sources else "")


def analyse(comparison: Comparison, noise_floor=None) -> Analysis:
    """The envelope, who owns it, and by how much.

    ``noise_floor`` overrides the default for the result type - in percent for
    flows and volumes, in metres for level.
    """
    kind, default_floor, label = margin_scale(comparison.key)
    floor = default_floor if noise_floor is None else float(noise_floor)
    frame = comparison.frame
    if frame.empty:
        return Analysis(margin_kind=kind, margin_label=label, noise_floor=floor)

    envelope = frame.max(axis=1, skipna=True)
    critical = frame.idxmax(axis=1, skipna=True)
    margin = _margins(frame, envelope, kind)

    bands = _bands(critical, envelope.index, margin)
    switches = _switches(critical, bands)
    warnings = _warnings(comparison, critical, switches, kind, floor)
    never = _never_critical(comparison, critical)

    return Analysis(envelope=envelope, critical=critical, margin=margin,
                    switches=switches, bands=bands, warnings=warnings,
                    never_critical=never, margin_kind=kind,
                    margin_label=label, noise_floor=floor)


def _margins(frame: pd.DataFrame, envelope: pd.Series, kind: str) -> pd.Series:
    """How far the winner beat the runner-up, in the units of the result type."""
    if frame.shape[1] < 2:
        return pd.Series(float("nan"), index=frame.index)
    runner_up = frame.apply(
        lambda values: values.dropna().nlargest(2).iloc[-1]
        if values.notna().sum() >= 2 else float("nan"), axis=1)
    difference = envelope - runner_up
    if kind == ABSOLUTE:
        return difference
    margin = difference / runner_up.abs() * 100.0
    return margin.replace([float("inf"), float("-inf")], float("nan"))


def _bands(critical: pd.Series, aeps, margin: pd.Series) -> list:
    """Contiguous runs of one owner, extended to the midpoints between them.

    Midpoints, so the shading is continuous: the changeover happens somewhere
    between the two AEPs that were evaluated, not at either of them.
    """
    labels = [value for value in critical.tolist()]
    zs = [normal_variate(aep) for aep in aeps]
    bands, start = [], 0
    for position in range(1, len(labels) + 1):
        if position < len(labels) and labels[position] == labels[start]:
            continue
        if isinstance(labels[start], str):
            z_from = zs[start] if start == 0 else (zs[start - 1] + zs[start]) / 2
            z_to = (zs[position - 1] if position == len(labels)
                    else (zs[position - 1] + zs[position]) / 2)
            within = margin.iloc[start:position].dropna()
            bands.append(Band(label=labels[start], z_from=z_from, z_to=z_to,
                              aep_from=float(aeps[start]),
                              aep_to=float(aeps[position - 1]),
                              peak_margin=float(within.max()) if len(within)
                              else float("nan")))
        start = position
    return bands


def _switches(critical: pd.Series, bands: list) -> list:
    """Where the envelope changes hands, and how convincing the takeover was.

    ``strength`` is the incoming band's peak margin, not the margin at the
    switch - see Band. A switch whose new owner never gets clear of the
    runner-up is an idxmax hop, and is reported as noise rather than marked up.
    """
    return [Switch(aep=later.aep_from, before=earlier.label,
                   after=later.label, strength=later.peak_margin)
            for earlier, later in zip(bands, bands[1:])]


def _extreme_labels(comparison: Comparison) -> tuple:
    """The labels of the shortest and longest durations actually plotted."""
    known = {label: duration for label, duration in comparison.durations.items()
             if duration is not None and label in comparison.frame.columns}
    if len(known) < 2:
        return None, None
    shortest = min(known, key=known.get)
    longest = max(known, key=known.get)
    return shortest, longest


def _warnings(comparison, critical, switches, kind, floor) -> list:
    """What the plot should say out loud."""
    warnings = []
    shortest, longest = _extreme_labels(comparison)
    if shortest is None:
        return warnings

    for label, direction, advice in (
        (longest, "longest", "Longer storms would show whether the curve has "
                             "turned over beyond it."),
        (shortest, "shortest", "Shorter storms would show whether the curve has "
                               "turned over below it."),
    ):
        pinned = [aep for aep, owner in critical.items() if owner == label]
        if not pinned:
            continue
        warnings.append(
            f"The {direction} duration you ran ({label}) is critical at "
            f"1 in {', '.join(format_aep(aep) for aep in pinned)}. {advice}")

    for switch in switches:
        if switch.strength == switch.strength and switch.strength < floor:
            warnings.append(
                f"{switch.after} takes over from {switch.before} at 1 in "
                f"{format_aep(switch.aep)} but never gets more than "
                f"{format_margin(switch.strength, kind)} clear of the next "
                f"duration - inside sampling noise, so not necessarily a real "
                f"crossover.")
    return warnings


def _never_critical(comparison, critical) -> list:
    if len(comparison.frame.columns) < 3:
        return []
    owners = {owner for owner in critical.tolist() if isinstance(owner, str)}
    return [label for label in comparison.frame.columns if label not in owners]


def table(comparison: Comparison, analysis: Analysis) -> pd.DataFrame:
    """The on-screen critical duration table: durations, max, who, by how much."""
    if comparison.is_empty:
        return pd.DataFrame()
    frame = comparison.frame.copy()
    frame["max"] = analysis.envelope
    frame["critical duration"] = analysis.critical
    frame[analysis.margin_label] = analysis.margin
    return frame
