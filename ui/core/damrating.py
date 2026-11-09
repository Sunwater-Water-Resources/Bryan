"""The dam's rating, for the outflow at a level no run reports as a quantile.

The design flood table's dam crest flood row has a level (the crest) and an AEP
read off the realisations, but the run's outflow quantiles sit at the standard
AEPs, not at the crest. The outflow at the crest is a property of the dam, not
of the storms: it is read off the rating the group was routed through.

**Which rating.** The one the group's results were made with:

- a ``reservoir routing`` row names its own: ``ELS file``, ``SQ file`` and
  ``FSL`` in the sims list, relative to the project folder (Bryan runs there);
- a ``monte carlo`` or ``ensemble`` row was routed by URBS, through the
  ``DAM ROUTE`` line of the model config's ``vec_file`` whose ``location=`` is
  ``max_keys.level`` (as ``Simulator.initialise_model`` finds it): ``FSL=``,
  ``datafile=`` (the ELS) and ``file=`` (the SQ), both in ``ratings_folder``.

**How it is read**, as ``lib/ReservoirRouting.py`` routes: the volume at the
level off the ELS, linearly; the SQ's storages are ML above the full supply
volume; the outflow linearly off the SQ; nothing below the FSL. A level above
the top of either curve is refused rather than extrapolated.

The files are read as they are on disk now, so a rating edited after the run
gives the new rating's outflow; ``core/staleness.py`` is what says a run is
older than its inputs.

pandas, numpy and the standard library only; nothing here imports nicegui.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .bryan import config_paths
from .paths import normalise_sep

ROUTED_BY_URBS = ("monte carlo", "ensemble")
RESERVOIR_ROUTING = "reservoir routing"


class RatingError(Exception):
    """A rating that cannot be found or read."""


@dataclass(frozen=True)
class DamRating:
    els: Path
    sq: Path
    fsl: float
    el: np.ndarray            # ELS levels, ascending (m AHD)
    volume: np.ndarray        # ELS volumes (ML)
    storage: np.ndarray       # SQ storage above the FSV (ML), ascending
    flow: np.ndarray          # SQ outflow (m3/s)
    where: str                # where the rating was named, for a message

    @property
    def fsv(self) -> float:
        return float(np.interp(self.fsl, self.el, self.volume))

    def outflow_at(self, level: float) -> float:
        """The outflow (m3/s) at a lake level; RatingError off the top of a curve."""
        level = float(level)
        if level <= self.fsl:
            return 0.0
        if level > self.el[-1]:
            raise RatingError(f"{level:.2f} m is above the top of the storage curve "
                              f"{self.els.name} ({self.el[-1]:.2f} m)")
        above = float(np.interp(level, self.el, self.volume)) - self.fsv
        if above > self.storage[-1]:
            raise RatingError(f"{level:.2f} m is above the top of the rating {self.sq.name} "
                              f"({self.storage[-1]:,.0f} ML above full supply)")
        return float(np.interp(above, self.storage, self.flow))


# -- reading the files ---------------------------------------------------------

def read_els(path: Path):
    """(levels, volumes) of an ELS file (CSV with EL and V), ascending in level."""
    try:
        frame = pd.read_csv(path)
    except (OSError, ValueError) as exc:
        raise RatingError(f"the storage curve {path.name} could not be read ({exc})") from exc
    missing = [column for column in ("EL", "V") if column not in frame.columns]
    if missing:
        raise RatingError(f"the storage curve {path.name} has no {' or '.join(missing)} column")
    frame = frame[["EL", "V"]].astype(float).sort_values("EL")
    if len(frame) < 2:
        raise RatingError(f"the storage curve {path.name} needs at least two rows")
    return frame["EL"].to_numpy(), frame["V"].to_numpy()


def read_sq(path: Path):
    """(storage above FSV, outflow) of an SQ file: five header lines, then pairs -
    as ``lib/ReservoirRouting._read_sq`` reads it."""
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise RatingError(f"the rating {Path(path).name} could not be read ({exc})") from exc
    storage, flow = [], []
    for number, line in enumerate(lines[5:], start=6):
        if not line.strip():
            continue
        parts = re.split(r"\s+", line.strip())
        try:
            s, q = (float(part) for part in parts[:2])
        except ValueError:
            raise RatingError(f"the rating {Path(path).name}: line {number} is not a "
                              f"storage and a flow ({line.strip()!r})") from None
        storage.append(s)
        flow.append(q)
    if len(storage) < 2:
        raise RatingError(f"the rating {Path(path).name} needs at least two pairs")
    storage, flow = np.array(storage), np.array(flow)
    if np.any(np.diff(storage) <= 0):
        raise RatingError(f"the rating {Path(path).name}: storage does not rise down the file")
    return storage, flow


def _in_folder(folder: Path, name: str) -> Path:
    """A rating file named on a DAM ROUTE line. URBS runs on Windows, where
    ``CALLIDE.SQ`` is ``callide.sq``; match the case the same way anywhere."""
    path = folder / name
    if path.is_file():
        return path
    try:
        for each in folder.iterdir():
            if each.name.lower() == name.lower():
                return each
    except OSError:
        pass
    return path


# -- finding the rating --------------------------------------------------------

def dam_route(vec_path: Path, location: str) -> dict:
    """FSL, ELS and SQ off the DAM ROUTE line naming ``location``, found as
    ``UrbsModel.find_dam_routing_line`` finds it: a commented line (``{...}``)
    does not start with DAM ROUTE, so it is never taken."""
    try:
        text = Path(vec_path).read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise RatingError(f"the vec file {Path(vec_path).name} could not be read "
                          f"({exc})") from exc
    for line in text.splitlines():
        if line.strip().upper().startswith("DAM ROUTE") and location in line:
            found = {"fsl": None, "els": None, "sq": None}
            for part in line.split():
                key, _, value = part.partition("=")
                key = key.lower()
                if key == "fsl":
                    found["fsl"] = float(value)
                elif key == "datafile":
                    found["els"] = value
                elif key == "file":
                    found["sq"] = value
                elif key == "vbf":
                    raise RatingError(f"the dam {location} is routed by volume (VBF=) in "
                                      f"{Path(vec_path).name}, which gives no full supply "
                                      f"level to read the rating from")
            for key, token in (("fsl", "FSL="), ("els", "datafile="), ("sq", "file=")):
                if found[key] is None:
                    raise RatingError(f"the DAM ROUTE line for {location} in "
                                      f"{Path(vec_path).name} has no {token}")
            return found
    raise RatingError(f"no DAM ROUTE line for {location!r} in {Path(vec_path).name}")


def _urbs_rating(model_config: Path) -> tuple:
    """(els, sq, fsl, where) for a run URBS routed."""
    try:
        data = json.loads(Path(model_config).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RatingError(f"the model config {Path(model_config).name} could not be read "
                          f"({exc})") from exc
    folder = str(Path(model_config).parent)
    try:
        location = data["max_keys"]["level"]
        vec = Path(config_paths().resolve(config_paths().resolve(folder, data["model_folder"]),
                                          data["vec_file"]))
        ratings = Path(config_paths().resolve(folder, data["ratings_folder"]))
    except (KeyError, TypeError) as exc:
        raise RatingError(f"the model config {Path(model_config).name} has no "
                          f"{exc.args[0] if exc.args else 'max_keys.level'}") from exc
    route = dam_route(vec, location)
    return (_in_folder(ratings, route["els"]), _in_folder(ratings, route["sq"]),
            route["fsl"], f"the DAM ROUTE line for {location} in {vec.name}")


def _sims_list_rating(row, project_folder: Path) -> tuple:
    """(els, sq, fsl, where) for a reservoir routing row."""
    def path_of(column):
        value = row.get(column)
        if value is None or (isinstance(value, float) and math.isnan(value)) \
                or not str(value).strip():
            raise RatingError(f"the reservoir routing row gives no {column}")
        local = Path(normalise_sep(str(value).strip()))
        return local if local.is_absolute() else Path(project_folder) / local
    try:
        fsl = float(row.get("FSL"))
    except (TypeError, ValueError):
        raise RatingError("the reservoir routing row gives no FSL") from None
    if math.isnan(fsl):
        raise RatingError("the reservoir routing row gives no FSL")
    return path_of("ELS file"), path_of("SQ file"), fsl, "the sims list"


_CACHE: dict = {}


def _load(els: Path, sq: Path, fsl: float, where: str) -> DamRating:
    for path, what in ((els, "storage curve"), (sq, "rating")):
        if not path.is_file():
            raise RatingError(f"the {what} named by {where} is not on disk ({path})")
    key = (str(els), els.stat().st_mtime_ns, str(sq), sq.stat().st_mtime_ns, fsl)
    if key not in _CACHE:
        el, volume = read_els(els)
        storage, flow = read_sq(sq)
        _CACHE[key] = DamRating(els=els, sq=sq, fsl=fsl, el=el, volume=volume,
                                storage=storage, flow=flow, where=where)
    return _CACHE[key]


def rating_for(run, rows) -> DamRating:
    """The rating a group's rows were routed through (``run`` is a study RunData).

    Every row of a group should share one; rows that name different ratings are
    refused rather than one of them picked."""
    if not rows:
        raise RatingError("the group has no rows")
    found = set()
    for index in rows:
        row = run.frame.loc[index]
        method = str(row.get("Method") or "").strip().lower()
        if method == RESERVOIR_ROUTING:
            found.add(_sims_list_rating(row, run.project_folder))
        elif method in ROUTED_BY_URBS:
            model_config = run.config.filepaths.get("model_config")
            if model_config is None:
                raise RatingError(f"{run.config.config_path.name} names no model_config")
            found.add(_urbs_rating(model_config))
        else:
            raise RatingError(f"a {method or 'blank'} method row names no rating")
    if len({(str(els), str(sq), fsl) for els, sq, fsl, _ in found}) > 1:
        names = sorted({f"{Path(sq).name} (FSL {fsl:g})" for _, sq, fsl, _ in found})
        raise RatingError(f"the group's rows route through different ratings: "
                          f"{', '.join(names)}")
    return _load(*found.pop())


def forget_cached() -> None:
    _CACHE.clear()
