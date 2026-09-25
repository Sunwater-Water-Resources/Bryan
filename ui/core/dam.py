"""The dam's inputs: what several analyses read, kept once in the study file.

The lake level record (gauge exports and an optional overlay gauge), the storage
table, the ratings, the SILO evaporation and its pan factors, the catchment
polygon and area, and the month the water year starts. The Lake record page's
four steps all read some of them, and so will Lake levels; they are entered once,
on the Study page, and kept under ``"dam"`` in the study file with paths relative
to it.

**Where they were before.** Each lived in the Lake record's own section - the
record under ``homogenise``, the polygon under ``rainfall``, the area under
``inflow``. Two things follow from that:

- a study with no ``"dam"`` section yet takes its values from there (``settings``),
  so nothing is retyped;
- every save writes them back there too (``write_through``), so a launcher that
  predates this, on a machine that has not pulled, still finds them. That is for
  one release; after it, only ``"dam"`` is read.

pandas-free and nicegui-free: plain dictionaries.
"""

from __future__ import annotations

import copy
from pathlib import Path

from .study import Study

KEY = "dam"
CALLIDE_PAN_FACTORS = [0.82, 0.82, 0.83, 0.79, 0.75, 0.71, 0.76, 0.81, 0.79, 0.81, 0.82, 0.83]

DEFAULTS = {
    "gauges": [], "overlay": None, "storage": "", "register": "", "register_fsl": None,
    "evaporation": "",
    "pan_factors": list(CALLIDE_PAN_FACTORS),
    "shapefile": "", "field": "", "value": "", "catchment_km2": None,
    "water_year_start": 10,
}

# Where each dam input lived in the Lake record's section, as (step, key). Read to
# fill a new "dam" section, and written on every save for older launchers.
LEGACY = {
    "gauges": ("homogenise", "gauges"),
    "overlay": ("homogenise", "overlay"),
    "storage": ("homogenise", "storage"),
    "register": ("homogenise", "register"),
    "register_fsl": ("homogenise", "register_fsl"),
    "evaporation": ("homogenise", "evaporation"),
    "pan_factors": ("homogenise", "pan_factors"),
    "water_year_start": ("homogenise", "water_year_start"),
    "shapefile": ("rainfall", "shapefile"),
    "field": ("rainfall", "field"),
    "value": ("rainfall", "value"),
    "catchment_km2": ("inflow", "catchment_km2"),
}

LAKE_RECORD = "lake_record"


def settings(study: Study) -> dict:
    """The study's dam inputs, completed with the defaults.

    A study saved before the dam section existed takes them from its Lake record
    settings, where they used to be kept.
    """
    out = copy.deepcopy(DEFAULTS)
    stored = study.extra.get(KEY)
    if isinstance(stored, dict):
        out.update(copy.deepcopy(stored))
        return out
    legacy = study.extra.get(LAKE_RECORD) or {}
    for key, (step, name) in LEGACY.items():
        value = (legacy.get(step) or {}).get(name)
        if value not in (None, "", []):
            out[key] = copy.deepcopy(value)
    return out


def write_through(section: dict, dam: dict) -> dict:
    """Lay the dam inputs over a Lake record section, where its jobs read them."""
    for key, (step, name) in LEGACY.items():
        section.setdefault(step, {})[name] = copy.deepcopy(dam[key])
    return section


def store(study: Study, dam: dict) -> None:
    """Keep the dam inputs in the study - and, for older launchers, where they were."""
    study.extra[KEY] = copy.deepcopy(dam)
    write_through(study.extra.setdefault(LAKE_RECORD, {}), dam)


def gauges(dam: dict) -> list:
    return [item for item in dam.get("gauges") or [] if str(item).strip()]


# A register is a workbook; anything else in its place is one rating for the whole
# record (lib/homogenise/curves.read_rating_source).
REGISTER_SUFFIXES = (".xlsx", ".xlsm", ".xls")


def single_rating(dam: dict) -> bool:
    """Whether the ratings are one rating (.rat, .csv, .sq) rather than a register."""
    text = str(dam.get("register") or "").strip()
    return bool(text) and Path(text).suffix.lower() not in REGISTER_SUFFIXES
