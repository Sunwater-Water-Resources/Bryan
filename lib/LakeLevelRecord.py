"""The recorded lake level: reading it, and its annual maximum series.

The first half of the lake level frequency analysis - everything that can be
done with pandas alone. ``lib/LakeLevelFrequency.py`` holds the other half, the
curve fits and their resampled bands, which need scipy.

**Dependencies are deliberately pandas and the standard library only.** The run
launcher imports this module directly (``ui/core/bryan.py``) so the page can show
the annual maxima, write them out and score the water year options without
starting a subprocess, and the util script that fits and plots them reads the
record through the same functions. The annual maxima on the page and in the
report figure are therefore one series, not two that ought to agree.

What it covers:

- **Reading a headwater level export.** Both layouts the Queensland exports come
  in: Hydstra's own (four header rows, point values, the site and variable
  descriptions running down the third column beside the first rows of data, and
  no quality code against any value) and the WMIP web export (the same header,
  plus a quality column). Several files are joined in the order given, each
  owning the record from its first reading to the next file's first reading, so a
  gauge replaced mid-record reads as one series.
- **Water years**, labelled by the calendar year they end in: with a 1 October
  start, water year 2011 is October 2010 to September 2011, printed "2010-11".
  The start month is a parameter throughout, never a constant.
- **The annual maximum series**, with the two things about each maximum that a
  frequency curve needs to know: whether the water year was covered well enough
  for its maximum to be the year's, and whether the maximum was *carried over* -
  the level the year opened at, left from the wet season before, rather than a
  rise a storm produced inside the year.
- **Plotting positions** (Cunnane), including the censored positions a subset of
  the maxima takes on the full record's probability scale.
- **Scoring the candidate water year starts** against the level record alone.

Levels are metres; AEP is carried both as a probability (``aep``) and as
``1 in X`` (``aep_1_in_x``), and the column name always says which.
"""

from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

_NORMAL = NormalDist()

# Cunnane's plotting position, (rank - a) / (n + 1 - 2a).
CUNNANE_A = 0.4

DEFAULT_WATER_YEAR_START = 10

# A water year covered for less than this share of its length is kept but
# flagged: its maximum is the maximum of what was recorded, which cannot be
# assumed to be the year's.
DEFAULT_MIN_COVERAGE = 0.9

# A maximum falling within this many days of the start of its water year is
# taken to be carried over from the year before. Read off the Callide record,
# where the maxima inside the first 1.5 days are followed by nothing until day 9;
# check the days-into-year column of another dam before relying on it.
DEFAULT_CARRYOVER_DAYS = 3.0

# A maximum inside that window is still storm-driven if the lake falls away from
# it and later comes back to within this of it, after the window. The case is a
# gated lake held at full supply: it opens the year on the plateau and is put
# back on it by the next flood, and a record read to the millimetre makes the two
# equal, so which one "is" the maximum is decided by rounding. Both conditions
# matter - a lake held flat on the plateau for a fortnight from 1 October has not
# been put back there by anything, and a millimetre of ripple in the first hours
# is not a flood. One reading's resolution: this reproduces the Callide carry-over
# classification from its record read to 0.1 mm and to 1 mm alike, where 5 mm
# starts reclassifying years that genuinely fell away.
DEFAULT_TIE_TOLERANCE_M = 0.001

# The Hydstra header is four rows: site, variable code, variable name, and the
# kind of value (Point, Mean, Max...). WMIP's fourth row names its columns.
HEADER_ROWS = 4

# The date layouts the exports use. Hydstra writes the day and hour unpadded
# ("1/10/1989 8:00"), which strptime accepts under the padded directives.
TIMESTAMP_FORMATS = ("%d/%m/%Y %H:%M", "%H:%M:%S %d/%m/%Y", "%d/%m/%Y %H:%M:%S")

# WMIP quality codes that carry no usable value. A Hydstra export has no quality
# column at all, so nothing can be screened out of one.
UNUSABLE_QUALITY = (151, 255)

# The export's own footer is a couple of lines; more than this failing to parse
# means the layout has changed, and a record silently half-read looks exactly
# like a short record.
MAX_UNPARSED_ROWS = 5


# -- reading an export --------------------------------------------------------

@dataclass
class Record:
    """A headwater level record, and what reading it found worth saying."""

    level: pd.Series                         # metres, on a sorted DatetimeIndex
    sites: list = field(default_factory=list)       # (site id, description)
    value_kinds: list = field(default_factory=list)  # 'Point', per file
    notes: list = field(default_factory=list)        # things the user should see
    screened: bool = False                   # quality codes were available

    @property
    def start(self):
        return self.level.index.min()

    @property
    def end(self):
        return self.level.index.max()

    @property
    def site_label(self) -> str:
        return " + ".join(site for site, _ in self.sites) or "record"


def _header(path: Path) -> list:
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as stream:
        reader = csv.reader(stream)
        return [next(reader, []) for _ in range(HEADER_ROWS)]


def _cell(row, position) -> str:
    return row[position].strip().strip('"').strip() if len(row) > position else ""


def _parse_stamps(text: pd.Series) -> pd.Series:
    best = None
    for layout in TIMESTAMP_FORMATS:
        parsed = pd.to_datetime(text, format=layout, errors="coerce")
        if best is None or parsed.notna().sum() > best.notna().sum():
            best = parsed
    return best


def read_export(path) -> Record:
    """Read one Hydstra or WMIP level export.

    Raises ValueError when the file is not recognisably an export, or when more
    than a footer's worth of rows fail to parse.
    """
    path = Path(path)
    header = _header(path)
    if len(header) < HEADER_ROWS or not header[0]:
        raise ValueError(f"{path.name}: too short to be a level export")

    site = _cell(header[0], 1)
    kind = _cell(header[3], 1)
    labels = [_cell(header[3], i).lower() for i in range(len(header[3]))]
    quality_at = labels.index("quality") if "quality" in labels else None
    # The descriptions run down the first free column after the values.
    meta_at = 3 if quality_at is not None else 2

    frame = pd.read_csv(path, skiprows=HEADER_ROWS, header=None, dtype=str,
                        names=range(8), keep_default_na=False,
                        encoding="utf-8-sig", encoding_errors="replace")
    stamps = _parse_stamps(frame[0].str.strip())
    unparsed = int(stamps.isna().sum() - (frame[0].str.strip() == "").sum())
    if stamps.notna().sum() == 0:
        raise ValueError(f"{path.name}: no row has a recognisable timestamp")
    if unparsed > MAX_UNPARSED_ROWS:
        raise ValueError(f"{path.name}: {unparsed} rows have no recognisable "
                         f"timestamp; the export layout may have changed")

    notes = []
    description = ""
    meta = frame[meta_at].str.strip().tolist()
    if "Sites:" in meta:
        following = meta.index("Sites:") + 1
        if following < len(meta):
            description = meta[following]

    level = pd.to_numeric(frame[1].str.strip(), errors="coerce")
    data = pd.DataFrame({"level": level.to_numpy()}, index=stamps)
    if quality_at is not None:
        data["quality"] = pd.to_numeric(frame[quality_at].str.strip(),
                                        errors="coerce").to_numpy()
    data = data[data.index.notna() & data["level"].notna()].sort_index()

    if quality_at is not None:
        unusable = data["quality"].isin(UNUSABLE_QUALITY)
        if unusable.any():
            notes.append(f"{path.name}: dropped {int(unusable.sum())} values with "
                         f"unusable quality codes")
            data = data[~unusable]
    else:
        notes.append(f"{path.name}: no quality codes in this export, so no value "
                     f"could be screened out")

    duplicated = data.index.duplicated(keep="last")
    if duplicated.any():
        notes.append(f"{path.name}: {int(duplicated.sum())} duplicate timestamps, "
                     f"the last of each kept")
        data = data[~duplicated]

    if kind and kind.lower() != "point":
        notes.append(f"{path.name}: values are '{kind}', not point readings - an "
                     f"averaged series understates every peak")
    series = data["level"].astype(float)
    series.index.name = "time"
    series.name = "level"
    return Record(level=series, sites=[(site or path.stem, description)],
                  value_kinds=[kind], notes=notes, screened=quality_at is not None)


def read_record(paths) -> Record:
    """Join one or more exports into a single record, in the order given.

    Each file owns the record from its own first reading up to the first reading
    of the next, which is how a replaced gauge hands over. Files have to be given
    in the order the gauges operated.
    """
    paths = [Path(p) for p in ([paths] if isinstance(paths, (str, Path)) else paths)]
    if not paths:
        raise ValueError("no level export given")
    records = [read_export(path) for path in paths]
    starts = [record.start for record in records]
    if starts != sorted(starts):
        raise ValueError("the level exports are not in the order the gauges "
                         "operated - give the earliest first")

    pieces, notes = [], []
    for position, record in enumerate(records):
        piece = record.level.loc[starts[position]:]
        if position + 1 < len(records):
            piece = piece.loc[piece.index < starts[position + 1]]
        pieces.append(piece)
        notes += record.notes
        if position:
            before = pieces[position - 1]
            if len(before) and len(piece):
                step = float(piece.iloc[0] - before.iloc[-1])
                if abs(step) > 0.05:
                    notes.append(f"handover to {paths[position].name} is a "
                                 f"{step:+.3f} m jump; check the two share a datum")

    level = pd.concat(pieces).sort_index()
    level.name = "level"
    return Record(level=level,
                  sites=[site for record in records for site in record.sites],
                  value_kinds=[k for record in records for k in record.value_kinds],
                  notes=notes, screened=all(r.screened for r in records))


# -- water years ---------------------------------------------------------------

MONTH_NAMES = ("January", "February", "March", "April", "May", "June", "July",
               "August", "September", "October", "November", "December")


def water_year(index, start_month=DEFAULT_WATER_YEAR_START) -> np.ndarray:
    """The calendar year each timestamp's water year ends in."""
    index = pd.DatetimeIndex(index)
    if start_month == 1:
        return index.year.to_numpy()
    return (index.year + (index.month >= start_month)).to_numpy()


def period_label(year, start_month=DEFAULT_WATER_YEAR_START) -> str:
    if start_month == 1:
        return str(int(year))
    return f"{int(year) - 1}-{str(int(year))[-2:]}"


def water_year_start(year, start_month=DEFAULT_WATER_YEAR_START) -> pd.Timestamp:
    return pd.Timestamp(year=int(year) - (start_month > 1), month=start_month, day=1)


# -- the annual maximum series -------------------------------------------------

AMS_COLUMNS = ("water_year", "period", "level", "level_at", "days_into_year",
               "carried_over", "coverage", "complete", "reading_interval_min")


def annual_maxima(level: pd.Series, start_month=DEFAULT_WATER_YEAR_START,
                  min_coverage=DEFAULT_MIN_COVERAGE,
                  carryover_days=DEFAULT_CARRYOVER_DAYS,
                  tie_tolerance=DEFAULT_TIE_TOLERANCE_M) -> pd.DataFrame:
    """One row per water year: the maximum, when it fell, and how far to trust it.

    ``level_at`` is the first time the maximum was reached, except for a year
    whose lake leaves that maximum and returns to it after the carry-over window
    (see :data:`DEFAULT_TIE_TOLERANCE_M`), where it is the return - the flood's.

    ``reading_interval_min`` is the typical gap between readings in the day
    either side of the maximum. A record read once a day misses peaks that fall
    between readings, and the early parts of most records were.
    """
    level = level.dropna().sort_index()
    if level.empty:
        return pd.DataFrame(columns=list(AMS_COLUMNS))
    years = water_year(level.index, start_month)
    stamps = level.index.to_numpy()
    rows = []
    for year, block in level.groupby(years):
        at = block.idxmax()
        opened = water_year_start(year, start_month)
        closed = opened + pd.DateOffset(years=1)
        covered = (block.index.max() - block.index.min()).total_seconds()
        length = (closed - opened).total_seconds()

        lo = np.searchsorted(stamps, np.datetime64(at - pd.Timedelta("1D")))
        hi = np.searchsorted(stamps, np.datetime64(at + pd.Timedelta("1D")),
                             side="right")
        gaps = np.diff(stamps[lo:hi]).astype("timedelta64[s]").astype(float) / 60.0
        interval = float(np.median(gaps)) if len(gaps) else math.nan

        days_in = (at - opened).total_seconds() / 86400.0
        if days_in < carryover_days:
            peak = block.max()
            after = block[block.index > at]
            left = after.index[after < peak - tie_tolerance]
            if len(left):
                again = after[(after.index > left[0])
                              & (after.index >= opened + pd.Timedelta(days=carryover_days))
                              & (after >= peak - tie_tolerance)]
                if len(again):
                    at = again.index[0]
                    days_in = (at - opened).total_seconds() / 86400.0
        rows.append({
            "water_year": int(year),
            "period": period_label(year, start_month),
            "level": float(block.max()),
            "level_at": at,
            "days_into_year": days_in,
            "carried_over": bool(days_in < carryover_days),
            "coverage": covered / length,
            "complete": bool(covered / length >= min_coverage),
            "reading_interval_min": interval,
        })
    return pd.DataFrame(rows, columns=list(AMS_COLUMNS))


def read_ams_csv(path, level_column=None, start_month=DEFAULT_WATER_YEAR_START,
                 carryover_days=DEFAULT_CARRYOVER_DAYS,
                 min_coverage=DEFAULT_MIN_COVERAGE) -> pd.DataFrame:
    """An annual maximum series someone has already derived.

    For a record that cannot be used as it was measured - a dam whose full supply
    level or operation changed part way through, where the maxima have to be
    re-routed to one condition first. Reads the file this module writes, and the
    common alternatives: ``WaterYear``/``Level_max``/``Level_max_at`` and the
    ``New_`` homogenised columns. ``level_column`` picks one explicitly.

    Without a timestamp column nothing can be said about carry-over, so every
    maximum is taken as storm-driven, and the caller is told so by an all-NaN
    ``days_into_year``.
    """
    frame = pd.read_csv(path, comment="#")
    lower = {column.lower(): column for column in frame.columns}

    def pick(*names):
        for name in names:
            if name.lower() in lower:
                return lower[name.lower()]
        return None

    level_name = level_column or pick("level", "new_level_max", "level_max")
    if level_name is None or level_name not in frame.columns:
        raise ValueError(f"{Path(path).name}: no level column found - name it "
                         f"(columns are {list(frame.columns)})")
    at_name = pick(f"{level_name}_at", "level_at", "new_level_max_at", "level_max_at")
    year_name = pick("water_year", "wateryear", "year")

    out = pd.DataFrame({"level": pd.to_numeric(frame[level_name], errors="coerce")})
    out["level_at"] = (pd.to_datetime(frame[at_name], format="ISO8601", errors="coerce")
                       if at_name else pd.NaT)
    if year_name is not None:
        out["water_year"] = pd.to_numeric(frame[year_name], errors="coerce")
    elif at_name:
        out["water_year"] = water_year(out["level_at"], start_month)
    else:
        raise ValueError(f"{Path(path).name}: needs a water year or a timestamp column")
    out = out.dropna(subset=["level", "water_year"]).copy()
    out["water_year"] = out["water_year"].astype(int)
    out["period"] = [period_label(y, start_month) for y in out["water_year"]]
    opened = pd.to_datetime([water_year_start(y, start_month) for y in out["water_year"]])
    out["days_into_year"] = ((out["level_at"] - opened).dt.total_seconds() / 86400.0
                             if at_name else math.nan)
    out["carried_over"] = (out["days_into_year"] < carryover_days).fillna(False) \
        .astype(bool)
    coverage_name = pick("coverage")
    out["coverage"] = (pd.to_numeric(frame.loc[out.index, coverage_name], errors="coerce")
                       if coverage_name else 1.0)
    out["complete"] = out["coverage"].fillna(1.0) >= min_coverage
    out["reading_interval_min"] = math.nan
    return out[list(AMS_COLUMNS)].reset_index(drop=True)


# -- plotting positions --------------------------------------------------------

def normal_variate_of_probability(p) -> float:
    """z for an exceedance probability: the inverse normal of 1 - p."""
    p = float(p)
    if not 0.0 < p < 1.0:
        return math.nan
    return _NORMAL.inv_cdf(1.0 - p)


def exceedance_of_variate(z) -> float:
    return 1.0 - _NORMAL.cdf(float(z))


def plotting_positions(values, n_years=None, variate=None) -> tuple:
    """Cunnane exceedance probability and z for each value, in input order.

    ``n_years`` defaults to the number of values. Giving it separately places a
    subset of the annual maxima on the full record's probability scale - see
    :func:`censored_positions`.

    ``variate`` maps an array of exceedance probabilities to z. The default is
    ``statistics.NormalDist``; ``lib/LakeLevelFrequency.py`` passes scipy's, which
    differs only in the sixteenth figure - but the constrained shoulder fit is
    sensitive enough to that for a resample or two to converge under one and not
    the other, so the resampling uses the implementation the curves were
    developed with.
    """
    values = np.asarray(values, dtype=float)
    n_years = len(values) if n_years is None else int(n_years)
    order = np.argsort(-values, kind="stable")
    rank = np.empty(len(values))
    rank[order] = np.arange(1, len(values) + 1)
    p = (rank - CUNNANE_A) / (n_years + 1 - 2 * CUNNANE_A)
    if variate is not None:
        return p, np.asarray(variate(p), dtype=float)
    z = np.array([normal_variate_of_probability(x) for x in p])
    return p, z


def censored_positions(values, n_years, variate=None) -> tuple:
    """(z, level) for a subset of the maxima, ascending, on the record's scale.

    Ranking the subset among itself would slide every point to a more frequent
    AEP and flatter any curve fitted through it. Keeping the whole record's
    length as the denominator leaves the frequent end of the scale empty, which
    is where the excluded years belong.
    """
    values = np.sort(np.asarray(values, dtype=float))[::-1]
    _, z = plotting_positions(values, n_years, variate)
    order = np.argsort(z, kind="stable")
    return z[order], values[order]


def with_positions(ams: pd.DataFrame, include_incomplete=True) -> pd.DataFrame:
    """The maxima that go into the curve, with their plotting positions.

    Sorted by ascending z. ``storm_z`` is each storm-driven maximum's censored
    position (NaN for the carried-over ones); see :func:`censored_positions`.
    """
    used = ams if include_incomplete else ams[ams["complete"]]
    used = used.dropna(subset=["level"]).copy()
    p, z = plotting_positions(used["level"])
    used["aep"] = p
    used["aep_1_in_x"] = 1.0 / p
    used["z"] = z
    used = used.sort_values("z", kind="stable").reset_index(drop=True)

    storm = ~used["carried_over"].astype(bool)
    storm_values = used.loc[storm, "level"].to_numpy()
    _, storm_z = plotting_positions(storm_values, len(used))
    used["storm_z"] = math.nan
    used.loc[storm, "storm_z"] = storm_z
    return used


# -- writing the series out ----------------------------------------------------

def ams_csv_text(ams: pd.DataFrame, metadata: dict | None = None) -> str:
    """The annual maxima as CSV, with what produced them in '#' lines above.

    The water year start, the carry-over cut and the source files are part of
    the series - change any of them and the numbers change - so they travel
    with it. :func:`read_ams_csv` skips the '#' lines.
    """
    buffer = io.StringIO()
    for key, value in (metadata or {}).items():
        buffer.write(f"# {key}: {value}\n")
    table = ams.copy()
    if "level_at" in table:
        table["level_at"] = pd.to_datetime(table["level_at"]).dt.strftime("%Y-%m-%d %H:%M")
    for column, places in (("level", 3), ("days_into_year", 2), ("coverage", 4),
                           ("reading_interval_min", 0), ("aep", 5),
                           ("aep_1_in_x", 2), ("z", 4), ("storm_z", 4)):
        if column in table:
            table[column] = table[column].astype(float).round(places)
    table.to_csv(buffer, index=False, lineterminator="\n")
    return buffer.getvalue()


# -- choosing the water year ----------------------------------------------------

# Empirical levels at these AEPs are reported per candidate start, to show
# whether the choice moves the curve over the range it is read.
SCORE_AEPS = (0.5, 0.2, 0.1)

# Two maxima within this window of each other are the same wet season.
PROMOTION_WINDOW = pd.Timedelta("183D")


def _empirical_level(levels, aep) -> float:
    levels = np.asarray(levels, dtype=float)
    if len(levels) < 3:
        return math.nan
    _, z = plotting_positions(levels)
    order = np.argsort(z)
    target = normal_variate_of_probability(aep)
    if not z[order][0] <= target <= z[order][-1]:
        return math.nan
    return float(np.interp(target, z[order], levels[order]))


def water_year_scores(level: pd.Series, start_months=range(1, 13), *,
                      carryover_days=DEFAULT_CARRYOVER_DAYS,
                      min_coverage=DEFAULT_MIN_COVERAGE,
                      high_level=None, full_supply=None, rise_m=0.10,
                      promotion_tolerance_m=0.05,
                      include_incomplete=True) -> pd.DataFrame:
    """Score each candidate water year start against the level record.

    None of these asks which start has been adopted, so they can be used to
    choose one:

    ``carried_over``, ``carried_over_high``
        maxima inherited from the year before, in all and above ``high_level``.
        An inherited maximum is real exposure, but it means one flood supplies
        two years' maxima; what matters most is how many land high enough to
        sit in the range the curve is read over.
    ``closest_storm_max_d``
        how near any storm-driven maximum comes to a boundary. The margin
        against a future flood being split between two years.
    ``cuts_rising``
        boundaries where the lake rose more than ``rise_m`` over the day either
        side - a flood in progress at the cut.
    ``cuts_above_fsl``
        boundaries with the lake above ``full_supply`` (when given): spilling.
    ``promoted``
        storm-driven maxima more than ``promotion_tolerance_m`` below a level
        reached within six months of them, and so a maximum only because of
        where the year was cut.
    ``level_<aep>``
        the empirical level at 50, 20 and 10% AEP.

    Water years at either end of the record are only counted where covered to
    ``min_coverage`` unless ``include_incomplete``.
    """
    level = level.dropna().sort_index()
    stamps = level.index
    values = level.to_numpy()
    start, end = stamps.min(), stamps.max()
    rows = []
    for month in start_months:
        ams = annual_maxima(level, month, min_coverage, carryover_days)
        if not include_incomplete:
            ams = ams[ams["complete"]]
        storm = ams[~ams["carried_over"]]

        cuts = [water_year_start(year, month) for year in
                range(start.year - 1, end.year + 2)]
        cuts = [cut for cut in cuts if start < cut < end]

        def span(cut, before, after):
            lo = stamps.searchsorted(cut - before)
            hi = stamps.searchsorted(cut + after, side="right")
            return values[lo:hi]

        rising = above = 0
        for cut in cuts:
            window = span(cut, pd.Timedelta("1D"), pd.Timedelta("1D"))
            if len(window) >= 2 and window[-1] - window[0] > rise_m:
                rising += 1
            if full_supply is not None:
                at_cut = span(cut, pd.Timedelta("1D"), pd.Timedelta(0))
                if len(at_cut) and at_cut[-1] > full_supply:
                    above += 1

        closest = math.nan
        if len(storm) and cuts:
            closest = min(min(abs((at - cut).total_seconds()) for cut in cuts)
                          for at in storm["level_at"]) / 86400.0

        promoted = 0
        for at, peak in zip(storm["level_at"], storm["level"]):
            nearby = span(at, PROMOTION_WINDOW, PROMOTION_WINDOW)
            if len(nearby) and nearby.max() > peak + promotion_tolerance_m:
                promoted += 1

        row = {
            "start_month": MONTH_NAMES[month - 1],
            "years": int(len(ams)),
            "carried_over": int(ams["carried_over"].sum()),
            "carried_over_high": (int((ams["carried_over"]
                                       & (ams["level"] > high_level)).sum())
                                  if high_level is not None else None),
            "closest_storm_max_d": None if math.isnan(closest) else round(closest, 1),
            "cuts_rising": rising,
            "cuts_above_fsl": above if full_supply is not None else None,
            "promoted": promoted,
        }
        for aep in SCORE_AEPS:
            row[f"level_{int(round(aep * 100))}pct_aep"] = _empirical_level(
                ams["level"], aep)
        rows.append(row)
    return pd.DataFrame(rows)


def monthly_levels(level: pd.Series) -> pd.DataFrame:
    """Median, 10th and 90th percentile daily level by calendar month.

    Free of any water year choice, so it can justify one: where the drawdown
    bottoms out and when the lake starts to rise.
    """
    daily = level.dropna().resample("1D").mean().dropna()
    grouped = daily.groupby(daily.index.month)
    frame = pd.DataFrame({
        "month": [MONTH_NAMES[m - 1][:3] for m in grouped.median().index],
        "median": grouped.median().to_numpy(),
        "p10": grouped.quantile(0.10).to_numpy(),
        "p90": grouped.quantile(0.90).to_numpy(),
    })
    return frame


def maxima_by_month(ams: pd.DataFrame) -> pd.Series:
    """How many storm-driven maxima fall in each calendar month.

    Depends on the water year it was cut with, so it describes a choice rather
    than justifying one; :func:`monthly_levels` is the independent evidence.
    """
    storm = ams[~ams["carried_over"]]
    months = pd.to_datetime(storm["level_at"]).dt.month
    counts = months.value_counts().reindex(range(1, 13), fill_value=0)
    counts.index = [name[:3] for name in MONTH_NAMES]
    return counts


# -- the analysis job ------------------------------------------------------------
#
# The launcher page and util/LakeLevelFrequency.py read one description of the
# analysis - the job - so the maxima on screen, in the exported CSV and in the
# report figure are derived by the same call. Paths in a job are absolute; the
# page resolves the project-relative ones before writing it.

JOB_DEFAULTS = {
    "record": {"files": [], "ams_csv": "", "level_column": ""},
    "water_year_start": DEFAULT_WATER_YEAR_START,
    "min_coverage": DEFAULT_MIN_COVERAGE,
    "include_incomplete": True,
    "carryover_days": DEFAULT_CARRYOVER_DAYS,
    "tie_tolerance": DEFAULT_TIE_TOLERANCE_M,
    "fsl": None,
    "fsl_label": "FSL",
    "reference_levels": [],          # [{"label": ..., "level": ...}]
    "fit": {"form": "shouldered", "degree": 4, "plateau_tolerance": 0.025,
            "plateau_gap": 0.100, "storm_driven": True, "draws": 400,
            "seed": 20260826},
    "design": {"include": False, "label": "", "sources": []},   # [{"duration", "path"}]
    "axes": {"rare_aep_1_in_x": 2000, "level_min": None, "level_max": None},
}


def complete_job(job: dict) -> dict:
    """The job with every missing key given its default. Never mutates ``job``."""
    out = {}
    for key, default in JOB_DEFAULTS.items():
        value = (job or {}).get(key, default)
        if isinstance(default, dict):
            merged = dict(default)
            merged.update(value or {})
            value = merged
        out[key] = value
    return out


def job_inputs(job: dict) -> list:
    """Every file the job's results depend on."""
    job = complete_job(job)
    paths = list(job["record"]["files"] or [])
    if job["record"]["ams_csv"]:
        paths.append(job["record"]["ams_csv"])
    paths += [source["path"] for source in job["design"]["sources"] or []]
    return [str(path) for path in paths]


def fingerprint(job: dict) -> str:
    """A hash of the job and of the size and time of every input it reads.

    What lets a results file be reused: the same fingerprint means the same
    numbers, and a re-exported record or a re-run duration changes it.
    """
    import hashlib
    import json
    import os

    job = complete_job(job)
    stamps = []
    for path in job_inputs(job):
        try:
            info = os.stat(path)
            stamps.append([path, info.st_size, info.st_mtime_ns])
        except OSError:
            stamps.append([path, None, None])
    text = json.dumps({"job": job, "inputs": stamps}, sort_keys=True, default=str)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def ams_for_job(job: dict) -> tuple:
    """``(ams, record)`` for a job. ``record`` is None when read from an AMS CSV."""
    job = complete_job(job)
    source = job["record"]
    start = int(job["water_year_start"])
    if source["files"]:
        level_record = read_record(source["files"])
        ams = annual_maxima(level_record.level, start, float(job["min_coverage"]),
                            float(job["carryover_days"]), float(job["tie_tolerance"]))
        return ams, level_record
    if source["ams_csv"]:
        ams = read_ams_csv(source["ams_csv"], source["level_column"] or None, start,
                           float(job["carryover_days"]), float(job["min_coverage"]))
        return ams, None
    raise ValueError("the job names neither a level record nor an annual maximum CSV")


def ams_metadata(job: dict, level_record=None) -> dict:
    """What produced a series, for the '#' lines above it in the CSV."""
    job = complete_job(job)
    source = job["record"]
    meta = {
        "source": "; ".join(source["files"]) if source["files"] else source["ams_csv"],
        "water_year_start": MONTH_NAMES[int(job["water_year_start"]) - 1],
        "water_year_label": "the calendar year the water year ends in",
        "carried_over": f"maximum within {job['carryover_days']:g} days of the start "
                        f"of the water year, and not reached again later in it",
        "complete": f"water year covered for at least {job['min_coverage']:.0%}",
        "plotting_position": f"Cunnane (a = {CUNNANE_A}), over "
                             f"{'all' if job['include_incomplete'] else 'complete'} "
                             f"water years",
    }
    if level_record is not None:
        meta["sites"] = "; ".join(f"{site} {text}".strip()
                                  for site, text in level_record.sites)
        meta["record"] = f"{level_record.start} to {level_record.end}"
    return meta
