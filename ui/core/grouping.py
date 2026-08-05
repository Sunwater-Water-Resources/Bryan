"""Grouping simulations that differ only by storm duration.

A 'group' is a set of rows that form one analysis - the same case run across
several storm durations, so a critical duration can be picked. Selecting and
running them together is the natural unit of work.

Bryan has no notion of this and no sims list in the Tinaroo or Callide projects
carries a ``Group`` column, so one is derived. The rule is deliberately narrow:
strip a duration token from the row's name **only when its value equals that
row's own Duration**. That stops ``GWL1p7`` or a fixed ``_24h`` baked into a
model name from being mistaken for the varying part.

Worked against all three real conventions:

    TFD_SimsList_01        TFD_mc_C030_ebf_18h_L20-4_GWL1p3, Duration 18
                           -> '18h' matches -> TFD_mc_C030_ebf_L20-4_GWL1p3|exg
    CLD_RFSL_mc_sims_01    E010_no-pbp_L105-1.7_GWL0p3_RGN, Duration 120
                           -> no token -> unchanged (24 durations, one group)
    TFD_SimsList_LongList_02   Output file blank (uncached formulas)
                           -> FR-4C-no_raise|670.42|reservoir routing
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

from .paths import cell_text, is_blank

GROUP_COLUMN = "Group"

# 18h / 24 hr / 1.5hrs / 120hour. The lookarounds stop it matching inside a
# longer token, so 'GWL1p7' and 'L20-4' are safe.
DURATION_TOKEN_RE = re.compile(
    r"(?<![A-Za-z0-9])(\d+(?:[._]\d+)?)\s*(?:h|hr|hrs|hour|hours)(?![A-Za-z0-9])",
    re.IGNORECASE,
)

# Provenance values, shown in the UI so a derived key is never mistaken for one
# the sims list actually declared.
FROM_COLUMN = "Group column"
FROM_NAME = "derived from the output name"
FROM_NAME_STRIPPED = "derived from the output name (duration removed)"
FROM_FALLBACK = "derived from Output suffix (no output name)"


@dataclass(frozen=True)
class GroupKey:
    key: str
    provenance: str

    def __str__(self) -> str:
        return self.key


@dataclass(frozen=True)
class GroupReport:
    """Whether the derivation produced something worth showing."""

    total_rows: int
    group_count: int
    from_column: bool
    degenerate: str = ""     # '' when the grouping looks useful

    @property
    def is_useful(self) -> bool:
        return not self.degenerate


def _as_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if pd.isna(number) else number


def strip_duration_token(name: str, duration) -> tuple[str, bool]:
    """Remove duration tokens from ``name`` that match ``duration``.

    Returns (name, stripped). A token is only removed when it is numerically
    the row's own duration - that equality check is what makes this safe.
    """
    target = _as_number(duration)
    if target is None or not name:
        return name, False

    stripped = False

    def replace(match):
        nonlocal stripped
        token = _as_number(match.group(1).replace("_", "."))
        if token is not None and token == target:
            stripped = True
            return "\x00"        # a placeholder, tidied below
        return match.group(0)

    result = DURATION_TOKEN_RE.sub(replace, name)
    if not stripped:
        return name, False
    return _tidy(result.replace("\x00", "")), True


def _tidy(text: str) -> str:
    """Close the gap a removed token leaves, without mangling the rest."""
    text = re.sub(r"[_\-]{2,}", lambda m: m.group(0)[0], text)
    text = re.sub(r"\s{2,}", " ", text)
    return text.strip(" _-")


def group_key(row) -> GroupKey:
    """The group this row belongs to."""
    declared = row.get(GROUP_COLUMN) if hasattr(row, "get") else None
    if not is_blank(declared):
        return GroupKey(cell_text(declared), FROM_COLUMN)

    # 'Basename' is the spreadsheet's own idea of the name without the run
    # suffix, so prefer it where a list has one (Callide does).
    base = cell_text(row.get("Basename")) or cell_text(row.get("Output file"))
    suffix = cell_text(row.get("Output suffix"))

    if not base:
        # Every row of TFD_SimsList_LongList_02 lands here: its Output file is
        # an uncached formula. The option name plus the dam being routed is
        # still enough to tell the groups apart.
        parts = [suffix, cell_text(row.get("FSL")), cell_text(row.get("Method"))]
        return GroupKey("|".join(p for p in parts if p) or "(ungrouped)",
                        FROM_FALLBACK)

    base, stripped = strip_duration_token(base, row.get("Duration"))

    # Two options over one base name (the Tinaroo reservoir-routing case, where
    # every option re-routes the same source run) must not collapse together.
    parts = [base]
    if suffix and suffix != base:
        parts.append(suffix)
    results_folder = cell_text(row.get("Results folder"))
    if results_folder:
        parts.append(results_folder)

    return GroupKey("|".join(parts),
                    FROM_NAME_STRIPPED if stripped else FROM_NAME)


def add_group_keys(frame: pd.DataFrame) -> pd.Series:
    """A group key per row, aligned to ``frame``'s index."""
    if frame.empty:
        return pd.Series(dtype=object)
    return pd.Series(
        [group_key(frame.loc[index]).key for index in frame.index],
        index=frame.index, dtype=object,
    )


def group_provenance(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        return pd.Series(dtype=object)
    return pd.Series(
        [group_key(frame.loc[index]).provenance for index in frame.index],
        index=frame.index, dtype=object,
    )


def derived_group_report(frame: pd.DataFrame) -> GroupReport:
    """Say plainly when the derivation has not produced a useful grouping."""
    total = len(frame)
    if total == 0:
        return GroupReport(0, 0, False, "the sims list has no rows")

    keys = add_group_keys(frame)
    count = keys.nunique()
    from_column = GROUP_COLUMN in frame.columns and not frame[GROUP_COLUMN].isna().all()

    degenerate = ""
    if not from_column:
        if count == total and total > 1:
            degenerate = (
                f"every row came out in its own group ({count} of {total}). "
                f"The output names probably differ by more than the duration - "
                f"add a '{GROUP_COLUMN}' column to say what belongs together."
            )
        elif count == 1 and total > 1:
            degenerate = (
                f"all {total} rows came out in one group. The output names "
                f"probably do not distinguish the cases - add a "
                f"'{GROUP_COLUMN}' column."
            )
    return GroupReport(total, count, from_column, degenerate)


def groups_in_order(frame: pd.DataFrame) -> list[str]:
    """Distinct group keys, in the order they first appear in the sheet."""
    keys = add_group_keys(frame)
    seen, order = set(), []
    for key in keys:
        if key not in seen:
            seen.add(key)
            order.append(key)
    return order


def add_group_column(frame: pd.DataFrame) -> pd.DataFrame:
    """A copy of ``frame`` with a real ``Group`` column written in.

    Used by the 'write a Group column' action, which then goes through the
    ordinary save-as flow - the master workbook is never touched.
    """
    out = frame.copy()
    out[GROUP_COLUMN] = add_group_keys(frame)
    return out
