"""Monotone lookups: storage curves and spillway ratings.

Two kinds of rating appear in this analysis and they are not interchangeable:

*   A **stage-discharge** table (``RatingCurves.xlsx``, ``rating_6_215.5.csv``) gives
    release as a function of lake level.
*   A **storage-discharge** table (``CALLIDE_RFSL.sq``, ``CALLIDE_FSL.sq``) gives
    release as a function of storage above full supply.  This is the form URBS
    uses, and the form the hydrologic model's gate operation is defined in.

The two sq files are the 2026 rating review (rating 10.3) at the two operating
datums, and each has a level-indexed twin (``CLD_OUTFLOW_RFSL.rat``,
``CLD_OUTFLOW_FSL.rat``) carrying the same flows.  Both were generated from
``CALLIDE_STORAGE.els``, so every pair lands exactly on a row of the storage
table: converting a pair back through ``V(EL)`` reproduces its storage to the
last ML.  That is what makes the sq and the rat interchangeable descriptions of
one curve rather than two curves that happen to be close.

The two datums differ only in where flow starts being counted.  From 216.204 m
up the tables are identical, because it is the same structure and only the gate
operation's trigger level moves.  Below that the RFSL file carries the gates
already working (147 m3/s at 216.1 m) while the FSL file is still at zero.

The sq curve releases **142 m3/s at 120 ML above full supply** -- 10 mm of
level.  That is the auto gate operation ("no pre-release"), and it makes the
rating very nearly vertical at the origin.  Nothing here smooths that: a
level-indexed table would have to resolve millimetres to represent it, so the sq
rating is evaluated in storage space, where the segment is an ordinary one.
The sq curve also extends to 223.0 m and 30,552 m3/s, against the stage-discharge
218.0 m and 5,265 m3/s.
"""

from __future__ import annotations

import logging
import re
from bisect import bisect_right
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

log = logging.getLogger("bryan.homogenise.curves")

# Resolution of the dense tables the routing loop reads.  The storage curve has
# V' ~ 12,000 ML/m and V'' ~ 700 ML/m2, so linear interpolation at this spacing
# costs about 1e-11 m of level -- far below anything that matters, and it buys
# an O(1) scalar lookup in the inner loop.
DENSE_STEP = 0.0005


def monotone_interpolator(x, y, name, strict=False):
    """Shape-preserving interpolator over (x, y), with extrapolation disabled.

    Rating curves and storage tables are monotonic by construction.  A cubic
    spline through them is not: fitting ``scipy.interpolate.CubicSpline`` to
    these rating tables overshoots to 10,602 m3/s inside the flat FSL-to-lip
    section of rating 4 -- twice the maximum flow anywhere in that table.  PCHIP
    is monotone by construction and cannot do this.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if not np.all(np.diff(x) > 0):
        raise ValueError(f"{name}: x values must be strictly increasing")
    dy = np.diff(y)
    if strict and not np.all(dy > 0):
        raise ValueError(f"{name}: y values must be strictly increasing")
    if not np.all(dy >= 0):
        raise ValueError(f"{name}: y values must be non-decreasing")
    return PchipInterpolator(x, y, extrapolate=False)


class Curve:
    """A monotone PCHIP lookup with explicit, counted end-clamping.

    ``PchipInterpolator(extrapolate=False)`` returns NaN outside the fitted
    range, which would silently poison the balance.  Clamping is the defensible
    choice here (flow below FSL is zero; above the top of the table it is at
    least the top-of-table value), but it is counted and reported rather than
    hidden.
    """

    def __init__(self, x, y, name, strict=False):
        self.name = name
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self._interp = monotone_interpolator(self.x, self.y, name, strict=strict)
        self.clamps = {}

    @property
    def lo(self):
        return float(self.x[0])

    @property
    def hi(self):
        return float(self.x[-1])

    def _count(self, below, above):
        if below.any() or above.any():
            self.clamps[self.name] = (
                self.clamps.get(self.name, 0) + int(below.sum() + above.sum())
            )

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        below, above = x < self.x[0], x > self.x[-1]
        self._count(below, above)
        out = self._interp(np.clip(x, self.x[0], self.x[-1]))
        out = np.where(below, self.y[0], out)
        out = np.where(above, self.y[-1], out)
        return out if out.ndim else float(out)


class DenseCurve(Curve):
    """A ``Curve`` resampled onto a uniform grid for O(1) scalar lookup.

    The routing loop solves an implicit equation at every timestep and there are
    of the order of half a million of them, so the storage and area curves are
    evaluated hundreds of thousands of times.  A ``PchipInterpolator`` call on a
    Python float costs microseconds; a uniform-grid lookup costs tens of
    nanoseconds.

    Both the vectorised path (stage 1) and the scalar path (stage 2) read the
    *same* table, so the two stages are exactly self-consistent.  That is what
    makes the no-op test -- new rating equal to the rating in force reproduces
    the record -- hold to solver tolerance rather than to interpolation
    tolerance.
    """

    def __init__(self, x, y, name, strict=False, step=DENSE_STEP):
        super().__init__(x, y, name, strict=strict)
        count = int(round((self.hi - self.lo) / step)) + 1
        self.grid = np.linspace(self.lo, self.hi, count)
        self.values = self._interp(self.grid)
        self.values[0], self.values[-1] = self.y[0], self.y[-1]
        self._step = (self.hi - self.lo) / (count - 1)
        self._inv_step = 1.0 / self._step
        self._table = self.values.tolist()
        self._last = count - 1

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        below, above = x < self.lo, x > self.hi
        self._count(below, above)
        out = np.interp(np.clip(x, self.lo, self.hi), self.grid, self.values)
        return out if out.ndim else float(out)

    def at(self, x):
        """Scalar evaluation.  Hot path -- deliberately free of numpy."""
        if x <= self.lo:
            return self._table[0]
        if x >= self.hi:
            return self._table[-1]
        position = (x - self.lo) * self._inv_step
        index = int(position)
        if index >= self._last:
            return self._table[-1]
        low = self._table[index]
        return low + (self._table[index + 1] - low) * (position - index)

    def slope_at(self, x):
        """Scalar first derivative, from the same table the value comes from."""
        if x <= self.lo or x >= self.hi:
            return 0.0
        index = int((x - self.lo) * self._inv_step)
        if index >= self._last:
            index = self._last - 1
        return (self._table[index + 1] - self._table[index]) * self._inv_step


# --------------------------------------------------------------------------
# Ratings
# --------------------------------------------------------------------------

class StageDischargeRating:
    """Release as a tabulated function of lake level."""

    kind = "stage-discharge"

    def __init__(self, level, flow, name, fsl=None):
        self.name = name
        self.curve = DenseCurve(level, flow, name)
        self.fsl = float(np.min(level)) if fsl is None else float(fsl)
        self.top_level = self.curve.hi
        self.top_flow = float(self.curve.y[-1])

    @property
    def clamps(self):
        return self.curve.clamps

    def __call__(self, level):
        level = np.asarray(level, dtype=float)
        flow = self.curve(level)
        return np.where(level > self.fsl, flow, 0.0)

    def at(self, level):
        return self.curve.at(level) if level > self.fsl else 0.0

    def slope_at(self, level):
        return self.curve.slope_at(level) if level > self.fsl else 0.0

    def describe(self):
        return (f"{self.name}: stage-discharge, FSL {self.fsl:.3f} m AHD, "
                f"to {self.top_level:.3f} m / {self.top_flow:,.0f} m3/s")


class StorageDischargeRating:
    """Release as a tabulated function of storage above full supply (URBS sq).

    Evaluated in storage space via the storage curve, so the near-vertical
    segment at the origin -- 0 to 142 m3/s within 120 ML, which is 10 mm of
    level -- is represented exactly instead of being flattened into a
    level-indexed table that cannot resolve it.

    The pairs are joined linearly rather than by PCHIP.  URBS interpolates its
    own tabulated pairs linearly, so this is what the hydrologic model does; it
    also cannot overshoot, which on a curve this steep is worth having.
    """

    kind = "storage-discharge"

    def __init__(self, storage, flow, volume_of_level, fsl, name):
        self.name = name
        self.storage = np.asarray(storage, dtype=float)
        self.flow = np.asarray(flow, dtype=float)
        if not np.all(np.diff(self.storage) > 0):
            raise ValueError(f"{name}: storage must be strictly increasing")
        if not np.all(np.diff(self.flow) >= 0):
            raise ValueError(f"{name}: flow must be non-decreasing")
        self.volume_of_level = volume_of_level
        self.fsl = float(fsl)
        self.fsv = float(volume_of_level(self.fsl))
        self._storage = self.storage.tolist()
        self._flow = self.flow.tolist()
        self._last = len(self._storage) - 1
        self.top_storage = float(self.storage[-1])
        self.top_flow = float(self.flow[-1])
        self.clamps = {}

    @property
    def top_level(self):
        grid = self.volume_of_level.grid
        return float(np.interp(self.fsv + self.top_storage,
                               self.volume_of_level.values, grid))

    def _interp(self, storage):
        # Deliberately uncounted: the routing loop probes this at the top of the
        # storage table (230 m, far above the rating) to bracket the root, so
        # counting clamps here would report one on every single step.  Whether
        # the *solution* left the table is checked afterwards, against the
        # levels that were actually reached.
        if storage <= 0.0:
            return 0.0
        if storage >= self.top_storage:
            return self._flow[-1]
        index = bisect_right(self._storage, storage) - 1
        if index >= self._last:
            return self._flow[-1]
        x0, x1 = self._storage[index], self._storage[index + 1]
        y0, y1 = self._flow[index], self._flow[index + 1]
        return y0 + (y1 - y0) * (storage - x0) / (x1 - x0)

    def __call__(self, level):
        level = np.asarray(level, dtype=float)
        storage = self.volume_of_level(level) - self.fsv
        above = storage > self.top_storage
        if np.any(above):
            self.clamps[self.name] = self.clamps.get(self.name, 0) + int(above.sum())
        flow = np.interp(np.clip(storage, 0.0, self.top_storage),
                         self.storage, self.flow)
        return np.where(storage > 0.0, flow, 0.0)

    def at(self, level):
        return self._interp(self.volume_of_level.at(level) - self.fsv)

    def slope_at(self, level):
        """dQ/dLevel, by the chain rule through the storage curve."""
        storage = self.volume_of_level.at(level) - self.fsv
        if storage <= 0.0 or storage >= self.top_storage:
            return 0.0
        index = bisect_right(self._storage, storage) - 1
        if index >= self._last:
            return 0.0
        x0, x1 = self._storage[index], self._storage[index + 1]
        y0, y1 = self._flow[index], self._flow[index + 1]
        return (y1 - y0) / (x1 - x0) * self.volume_of_level.slope_at(level)

    def describe(self):
        return (f"{self.name}: storage-discharge, FSL {self.fsl:.3f} m AHD "
                f"(FSV {self.fsv:,.0f} ML), to {self.top_storage:,.0f} ML above "
                f"FSV / {self.top_flow:,.0f} m3/s")


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------

def read_storage(els_file, step=DENSE_STEP):
    """Read the elevation / area / storage table.

    Only the forward curves are returned.  Fitting a separate ``EL(V)`` curve
    and using it as the inverse of ``V(EL)`` is not safe: two independent PCHIP
    fits are not exact inverses of each other, and routing hundreds of thousands
    of steps through the round trip accumulates drift.  The simulation therefore
    solves in level space and never needs an inverse.

    A table that opens with several rows of no storage - below the lowest outlet,
    or where the survey starts; Kroombit's reads 0 ML from 248.6 to 249.0 m - has
    no volume curve until the last of them, so the rows before it are dropped. A
    flat stretch anywhere else is still refused.
    """
    els = pd.read_csv(els_file).sort_values("EL").reset_index(drop=True)
    volumes = els["V"].to_numpy(dtype=float)
    first = 0
    while first + 1 < len(volumes) and volumes[first + 1] <= volumes[0]:
        first += 1
    if first:
        log.info("Storage: %d rows of %g ML below %.2f m AHD left out; the volume "
                 "curve starts there", first, volumes[0], els["EL"].iloc[first])
        els = els.iloc[first:]
    volume_of_level = DenseCurve(els["EL"], els["V"], "storage V(EL)",
                                 strict=True, step=step)
    area_of_level = DenseCurve(els["EL"], els["A"], "storage A(EL)", step=step)
    log.info("Storage: %d points, %.2f to %.2f m AHD",
             len(els), els["EL"].min(), els["EL"].max())
    return volume_of_level, area_of_level


def read_sq(path):
    """Read a URBS storage-discharge file.

    The format is a free-text title, ``*``-prefixed comments, a line reading
    ``<n> PAIRS:`` and then n whitespace-separated pairs of storage above full
    supply (ML) and release (m3/s).  The declared count is checked, not
    trusted: a truncated table would otherwise cap the release at whatever the
    last surviving pair happened to be.
    """
    pairs, comments = _read_pairs(path)
    table = pd.DataFrame(pairs, columns=["storage_ML", "flow_m3s"])
    log.info("URBS sq: %s, %d pairs, 0 to %s ML above FSV, 0 to %s m3/s",
             Path(path).name, len(table), f"{table['storage_ML'].max():,.0f}",
             f"{table['flow_m3s'].max():,.0f}")
    return table, comments


def read_rat(path):
    """Read a URBS level-discharge file (``.rat``): the ``.sq`` layout, with the
    first column the lake level (m AHD) rather than the storage above full
    supply. Returns the table as ``level``, ``flow`` and the header's comments."""
    pairs, comments = _read_pairs(path)
    table = pd.DataFrame(pairs, columns=["level", "flow"])
    log.info("URBS rat: %s, %d pairs, %.3f to %.3f m AHD, 0 to %s m3/s",
             Path(path).name, len(table), table["level"].min(), table["level"].max(),
             f"{table['flow'].max():,.0f}")
    return table, comments


def _read_pairs(path):
    """The pairs and the header comments of a URBS ``.sq`` or ``.rat`` file."""
    lines = [line.strip() for line in Path(path).read_text().splitlines()]
    header = [i for i, line in enumerate(lines) if "PAIRS" in line.upper()]
    if not header:
        raise ValueError(f"{path}: no 'PAIRS:' line found")
    start = header[0]
    declared = int(lines[start].split()[0])

    pairs = []
    for line in lines[start + 1:]:
        if not line or line.startswith("*"):
            continue
        fields = line.split()
        if len(fields) < 2:
            break
        pairs.append((float(fields[0]), float(fields[1])))
        if len(pairs) == declared:
            break
    if len(pairs) != declared:
        raise ValueError(f"{path}: declares {declared} pairs but {len(pairs)} "
                         f"could be read")
    comments = [line for line in lines[:start] if line]
    return pairs, comments


# Full supply written into the sq header.  Both spellings have to be accepted:
# the reduced-FSL file says "RFSL is 215.5 mAHD" and the design-FSL file says
# "FSL is 216.1 mAHD".  Matching only the first silently returns None for the
# second, which switches the re-basing guard in ``load_rating`` off for exactly
# the file that most needs it -- the one at the other datum.
_SQ_FSL_PATTERN = r"\bR?FSL\s+is\s+([0-9]+(?:\.[0-9]+)?)"

FSL_TOLERANCE_M = 1e-6


def declared_sq_fsl(comments):
    """The full supply level a ``.sq`` header states, if it states one."""
    for line in comments:
        match = re.search(_SQ_FSL_PATTERN, line, re.IGNORECASE)
        if match:
            return float(match.group(1))
    return None


def load_rating(path, volume_of_level, fsl=None, name=None):
    """Build a rating from a URBS ``.sq`` or ``.rat`` file, or a level/flow CSV.

    ``fsl`` is checked against the file rather than applied to it.  That
    distinction matters: a ``.sq`` is indexed by storage *above full supply*, so
    handing it a different full supply level does not move the gate operation,
    it slides the entire structure rating -- ogee included.  Re-basing
    ``CALLIDE_RFSL.sq`` from 215.5 to 216.1 m reports 573 m3/s at 217.0 m
    against a true capacity of 2,885.  Each full supply level needs its own
    ``.sq``, generated from the gate operation rules for that level -- which is
    what ``CALLIDE_FSL.sq`` is.

    The check is only as good as the header it reads, so ``_SQ_FSL_PATTERN``
    has to match both spellings; see the note there.

    A stage-discharge table is checked the same way, against its first row,
    because there the full supply level *is* the origin of the table.
    """
    path = Path(path)
    name = name or path.name

    if path.suffix.lower() == ".sq":
        table, comments = read_sq(path)
        declared = declared_sq_fsl(comments)
        if declared is None and fsl is None:
            raise ValueError(
                f"{path}: the full supply level is not stated in the header "
                f"and none was supplied"
            )
        if declared is not None and fsl is not None and abs(declared - fsl) > FSL_TOLERANCE_M:
            raise ValueError(
                f"{path} declares a full supply level of {declared:.3f} m AHD "
                f"but {fsl:.3f} m was requested. A storage-discharge table "
                f"cannot be re-based: its storage axis is measured above full "
                f"supply, so shifting the datum moves the whole structure "
                f"rating rather than only the gate operation. Supply a .sq "
                f"generated for {fsl:.3f} m."
            )
        if declared is not None:
            fsl = declared
        rating = StorageDischargeRating(table["storage_ML"], table["flow_m3s"],
                                        volume_of_level, fsl, name)
    elif path.suffix.lower() == ".rat":
        table, comments = read_rat(path)
        declared = declared_sq_fsl(comments)
        if declared is not None and fsl is not None and abs(declared - fsl) > FSL_TOLERANCE_M:
            raise ValueError(
                f"{path} declares a full supply level of {declared:.3f} m AHD but "
                f"{fsl:.3f} m was given. Leave the full supply level blank to use "
                f"the file's, or supply a rating for {fsl:.3f} m.")
        claimed = declared if declared is not None else fsl
        origin = float(table["level"].min())
        if claimed is not None and abs(origin - claimed) > FSL_TOLERANCE_M:
            raise ValueError(
                f"{path} starts at {origin:.3f} m AHD but its full supply level is "
                f"{claimed:.3f} m. In a level-discharge table the full supply level "
                f"is the first row, so these have to agree.")
        rating = StageDischargeRating(table["level"], table["flow"], name, fsl=origin)
    else:
        table = pd.read_csv(path)
        origin = float(table["level"].min())
        if fsl is not None and abs(origin - fsl) > FSL_TOLERANCE_M:
            raise ValueError(
                f"{path} starts at {origin:.3f} m AHD but {fsl:.3f} m was "
                f"requested. In a stage-discharge table the full supply level "
                f"is the first row, so these have to agree."
            )
        rating = StageDischargeRating(table["level"], table["flow"], name,
                                      fsl=fsl if fsl is not None else origin)
    log.info("New rating -- %s", rating.describe())
    return rating


def read_ratings(workbook):
    """Read the rating register and one stage-discharge curve per rating."""
    register = pd.read_excel(workbook, sheet_name="Register", index_col=0)
    register["from"] = pd.to_datetime(register["from"], dayfirst=True)
    register["to"] = pd.to_datetime(register["to"], dayfirst=True)
    register = register.sort_values("from")

    ratings = {}
    for rating_id in register.index:
        frame = pd.read_excel(workbook, sheet_name=str(rating_id))
        register_fsl = float(register.loc[rating_id, "FSL"])
        ratings[rating_id] = StageDischargeRating(
            frame["level"], frame["flow"], f"rating {rating_id}", fsl=register_fsl,
        )
        table_fsl = float(frame["level"].min())
        if not np.isclose(table_fsl, register_fsl, atol=1e-6):
            log.warning("Rating %s: table starts at %.3f m but register FSL is "
                        "%.3f m", rating_id, table_fsl, register_fsl)
    log.info("Ratings: %d curves, register %s to %s", len(ratings),
             register["from"].min().date(), register["to"].max().date())
    return register, ratings


REGISTER_SUFFIXES = (".xlsx", ".xlsm", ".xls")

# A single rating is in force over any record: from long before a dam could have
# been gauged to long after. The register's span also clips the record, so it has
# to be wide enough never to.
SINGLE_RATING_SPAN = (pd.Timestamp("1800-01-01"), pd.Timestamp("2200-01-01"))


def single_rating(path, volume_of_level, fsl=None):
    """One rating for the whole record, as a one-row register.

    For a dam whose spillway has not changed: a URBS ``.rat`` or ``.sq``, or a
    ``level,flow`` csv, read by ``load_rating`` with the same full supply checks
    as a target rating. Everything downstream reads the register as it would a
    workbook's, so the homogenisation is unchanged.
    """
    rating = load_rating(path, volume_of_level, fsl=fsl, name=Path(path).name)
    register = pd.DataFrame({"from": [SINGLE_RATING_SPAN[0]], "to": [SINGLE_RATING_SPAN[1]],
                             "FSL": [float(rating.fsl)]},
                            index=pd.Index([1], name="Rating"))
    log.info("Ratings: one rating for the whole record, %s", rating.describe())
    return register, {1: rating}


def read_rating_source(path, volume_of_level, fsl=None):
    """The ratings in force over the record: a register workbook, or one rating.

    ``fsl`` is used only for a single rating whose file does not state one.
    """
    if Path(path).suffix.lower() in REGISTER_SUFFIXES:
        return read_ratings(path)
    return single_rating(path, volume_of_level, fsl)
