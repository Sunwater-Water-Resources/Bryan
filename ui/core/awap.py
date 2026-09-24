"""Catchment-average daily rainfall from gridded AWAP / AWRA-L netCDF files.

The antecedent storage analysis needs one daily rainfall series for the dam's
catchment. It used to arrive pre-made from another project (the JPA areal series
for Callide); this makes it from the grids themselves: point it at the catchment
polygon (a shapefile) and the folder of yearly netCDF files as downloaded -
``rain_day_1990.nc``, one per year, the AWRA-L layout of 0.05 degree cells over
Australia - and it writes ``date, rain_mm``.

**Weighting.** Each grid cell counts in proportion to the area of it the
catchment covers, found by testing a mesh of points in the cell (``SUBDIVISIONS``
a side) against the polygon, and to the cell's own area, which shrinks with
cos(latitude). A catchment of a few hundred cells is weighted to well under a
per cent either way. ``centre`` weighting instead counts every cell whose centre
is in the catchment equally, which is what the ``Average_<year>_Callide_mask.csv``
files in callide-fsl-reinstate are - a plain mask average - and is here so those
can be reproduced.

**Dates are the file's own.** An AWAP day is the 24 hours to 9 am, and which
date the file stamps on it is the file's convention, not something this changes;
the series carries the dates as read, and the analysis that uses it decides.

Reads only the catchment's window of each file, so a century of daily grids is
a century of small reads. Needs ``netCDF4`` and ``pyshp`` (and ``pyproj`` for a
polygon not in latitude and longitude), all in ui/requirements-ui.txt; imported
where used, so the rest of the UI runs without them. No nicegui.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

SUBDIVISIONS = 10
AREA, CENTRE = "area", "centre"
LAT_NAMES = ("latitude", "lat", "y")
LON_NAMES = ("longitude", "lon", "x")
YEAR_RE = re.compile(r"(18|19|20)\d{2}")


class AwapError(Exception):
    """Inputs that cannot give a rainfall series."""


def _need(module: str, purpose: str):
    try:
        return __import__(module)
    except ImportError as exc:
        raise AwapError(f"{purpose} needs the '{module}' package in the launcher's "
                        f"environment: pip install -r ui/requirements-ui.txt") from exc


# -- the catchment ----------------------------------------------------------------

@dataclass
class Catchment:
    """The polygon as rings of (longitude, latitude), outer and holes alike.

    Even-odd filling treats a hole as a hole whichever way it winds, so the
    rings need not be told apart.
    """

    rings: list
    name: str = ""

    @property
    def bounds(self) -> tuple:
        points = np.concatenate(self.rings)
        return (float(points[:, 0].min()), float(points[:, 1].min()),
                float(points[:, 0].max()), float(points[:, 1].max()))


def read_catchment(path, field_name: str = "", value: str = "") -> Catchment:
    """A polygon out of a shapefile, in longitude and latitude.

    With ``field_name``/``value`` only the matching records are taken (a regions
    file holding several catchments, as ``CLD_JPA_Regions_v3.shp`` does);
    without, every polygon in the file.
    """
    shapefile = _need("shapefile", "Reading a shapefile")
    path = Path(path)
    if not path.is_file():
        raise AwapError(f"shapefile not found: {path}")
    try:
        reader = shapefile.Reader(str(path))
    except Exception as exc:                          # noqa: BLE001 - report, never crash
        raise AwapError(f"{path.name} could not be read ({exc})") from exc
    fields = [f[0] for f in reader.fields[1:]]
    if field_name and field_name not in fields:
        raise AwapError(f"{path.name} has no field {field_name!r} (it has "
                        f"{', '.join(fields)})")
    rings = []
    for record in reader.iterShapeRecords():
        if field_name and str(record.record[field_name]).strip() != str(value).strip():
            continue
        shape = record.shape
        if shape.shapeType not in (5, 15, 25):       # polygon, polygonZ, polygonM
            continue
        parts = list(shape.parts) + [len(shape.points)]
        for start, end in zip(parts, parts[1:]):
            rings.append(np.asarray(shape.points[start:end], dtype=float)[:, :2])
    if not rings:
        raise AwapError(f"no polygon in {path.name}"
                        + (f" with {field_name} = {value!r}" if field_name else ""))
    rings = _to_geographic(path, rings)
    return Catchment(rings=rings, name=str(value or path.stem))


def _to_geographic(path: Path, rings: list) -> list:
    prj = path.with_suffix(".prj")
    wkt = prj.read_text(encoding="utf-8", errors="ignore") if prj.is_file() else ""
    if not wkt or wkt.lstrip().upper().startswith("GEOGCS"):
        x = np.concatenate(rings)[:, 0]
        if wkt or (np.all(np.abs(x) <= 180)):
            return rings
        raise AwapError(f"{path.name} has no .prj and its coordinates are not "
                        f"longitudes - give the projection file with it")
    pyproj = _need("pyproj", "A projected catchment polygon")
    try:
        transformer = pyproj.Transformer.from_crs(pyproj.CRS.from_wkt(wkt), "EPSG:4326",
                                                  always_xy=True)
    except Exception as exc:                          # noqa: BLE001
        raise AwapError(f"the projection of {path.name} could not be read ({exc})") from exc
    out = []
    for ring in rings:
        lon, lat = transformer.transform(ring[:, 0], ring[:, 1])
        out.append(np.column_stack([lon, lat]))
    return out


def inside(catchment: Catchment, lon, lat) -> np.ndarray:
    """Even-odd point-in-polygon over every ring, vectorised over the points."""
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    result = np.zeros(lon.shape, dtype=bool)
    for ring in catchment.rings:
        x0, y0 = ring[:-1, 0], ring[:-1, 1]
        x1, y1 = ring[1:, 0], ring[1:, 1]
        for ax, ay, bx, by in zip(x0, y0, x1, y1):
            crosses = (ay > lat) != (by > lat)
            if not crosses.any():
                continue
            x_at = ax + (lat - ay) * (bx - ax) / ((by - ay) if by != ay else 1e-300)
            result ^= crosses & (lon < x_at)
    return result


# -- the grid ------------------------------------------------------------------------

@dataclass
class Grid:
    lats: np.ndarray
    lons: np.ndarray
    lat_name: str
    lon_name: str
    variable: str


def describe(dataset, variable: str = "") -> Grid:
    names = {name.lower(): name for name in dataset.variables}
    lat = next((names[n] for n in LAT_NAMES if n in names), None)
    lon = next((names[n] for n in LON_NAMES if n in names), None)
    if lat is None or lon is None:
        raise AwapError(f"no latitude/longitude coordinates (it has "
                        f"{', '.join(dataset.variables)})")
    if variable:
        if variable not in dataset.variables:
            raise AwapError(f"no variable {variable!r} (it has {', '.join(dataset.variables)})")
    else:
        gridded = [name for name, var in dataset.variables.items()
                   if len(var.dimensions) == 3]
        if len(gridded) != 1:
            raise AwapError(f"name the rainfall variable - the file holds "
                            f"{', '.join(gridded) or 'no'} daily grids")
        variable = gridded[0]
    return Grid(lats=np.asarray(dataset.variables[lat][:], dtype=float),
                lons=np.asarray(dataset.variables[lon][:], dtype=float),
                lat_name=lat, lon_name=lon, variable=variable)


@dataclass
class Weights:
    """The catchment's cells in one grid: index windows and the weight of each."""

    lat_slice: slice
    lon_slice: slice
    weights: np.ndarray          # (lat, lon) over the window, summing to 1
    cells: int                   # cells with any weight

    def signature(self) -> tuple:
        return (self.lat_slice.start, self.lat_slice.stop, self.lon_slice.start,
                self.lon_slice.stop)


def _window(values: np.ndarray, low: float, high: float, step: float) -> slice:
    inside_range = np.where((values >= low - step) & (values <= high + step))[0]
    if not len(inside_range):
        raise AwapError("the catchment lies outside the grid")
    return slice(int(inside_range.min()), int(inside_range.max()) + 1)


def cell_weights(catchment: Catchment, grid: Grid, weighting: str = AREA,
                 subdivisions: int = SUBDIVISIONS) -> Weights:
    west, south, east, north = catchment.bounds
    d_lat = float(np.median(np.abs(np.diff(grid.lats))))
    d_lon = float(np.median(np.abs(np.diff(grid.lons))))
    lat_slice = _window(grid.lats, south, north, d_lat)
    lon_slice = _window(grid.lons, west, east, d_lon)
    lats, lons = grid.lats[lat_slice], grid.lons[lon_slice]
    centre_lon, centre_lat = np.meshgrid(lons, lats)

    if weighting == CENTRE:
        weights = inside(catchment, centre_lon, centre_lat).astype(float)
    else:
        offsets = (np.arange(subdivisions) + 0.5) / subdivisions - 0.5
        count = np.zeros(centre_lon.shape)
        for dy in offsets:
            for dx in offsets:
                count += inside(catchment, centre_lon + dx * d_lon, centre_lat + dy * d_lat)
        weights = count / subdivisions ** 2 * np.cos(np.radians(centre_lat))
    total = weights.sum()
    if total <= 0:
        raise AwapError("no grid cell falls in the catchment - is the polygon the "
                        "right one, and in the same place as the grid?")
    return Weights(lat_slice=lat_slice, lon_slice=lon_slice, weights=weights / total,
                   cells=int((weights > 0).sum()))


# -- the series ------------------------------------------------------------------------

def rainfall_files(folder, pattern: str = "*.nc") -> list:
    """The yearly files in a folder, in year order where the names carry one."""
    files = sorted(Path(folder).glob(pattern))
    if not files:
        raise AwapError(f"no {pattern} files in {folder}")

    def year(path):
        match = YEAR_RE.search(path.stem)
        return (int(match.group(0)) if match else 9999, path.name)
    return sorted(files, key=year)


@dataclass
class Series:
    frame: pd.DataFrame                              # date, rain_mm
    files: list = field(default_factory=list)
    cells: int = 0
    problems: list = field(default_factory=list)


def _dates(dataset) -> pd.DatetimeIndex:
    netCDF4 = _need("netCDF4", "Reading netCDF")
    time = dataset.variables.get("time")
    if time is None:
        raise AwapError("no 'time' variable")
    values = netCDF4.num2date(time[:], time.units,
                              getattr(time, "calendar", "standard"),
                              only_use_cftime_datetimes=False,
                              only_use_python_datetimes=True)
    return pd.DatetimeIndex([pd.Timestamp(value.year, value.month, value.day)
                             for value in values])


def catchment_rainfall(files, catchment: Catchment, variable: str = "",
                       weighting: str = AREA, progress=None) -> Series:
    """The weighted daily mean over the catchment, one file after another.

    Cells with no value on a day (the grid's fill) drop out of that day and the
    rest are reweighted; a day with no cell at all is left out and reported.
    """
    netCDF4 = _need("netCDF4", "Reading netCDF")
    parts, problems, weights, used = [], [], None, []
    for number, path in enumerate(files):
        if progress:
            progress(number, len(files), Path(path).name)
        try:
            with netCDF4.Dataset(str(path)) as dataset:
                grid = describe(dataset, variable)
                if weights is None or not _same_grid(weights, grid):
                    weights = cell_weights(catchment, grid, weighting)
                    weights.grid = (grid.lats.copy(), grid.lons.copy())
                var = dataset.variables[grid.variable]
                var.set_auto_mask(True)
                block = var[:, weights.lat_slice, weights.lon_slice]
                dates = _dates(dataset)
        except AwapError as exc:
            problems.append(f"{Path(path).name}: {exc}")
            continue
        except Exception as exc:                      # noqa: BLE001 - report, never crash
            problems.append(f"{Path(path).name} could not be read ({exc})")
            continue
        data = np.ma.filled(np.ma.asarray(block, dtype=float), np.nan)
        valid = np.isfinite(data) & (weights.weights > 0)
        weight = np.where(valid, weights.weights, 0.0)
        totals = weight.sum(axis=(1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            means = np.nansum(np.where(valid, data, 0.0) * weight, axis=(1, 2)) / totals
        empty = totals <= 0
        if empty.any():
            problems.append(f"{Path(path).name}: {int(empty.sum())} day(s) with no value "
                            f"in the catchment, left out")
        parts.append(pd.DataFrame({"date": dates[~empty], "rain_mm": means[~empty]}))
        used.append(Path(path))
    if not parts:
        raise AwapError("no rainfall could be read: " + "; ".join(problems[:3]))
    frame = pd.concat(parts, ignore_index=True).sort_values("date")
    duplicated = frame["date"].duplicated(keep="first")
    if duplicated.any():
        problems.append(f"{int(duplicated.sum())} date(s) in more than one file - the "
                        f"first kept")
        frame = frame[~duplicated]
    return Series(frame=frame.reset_index(drop=True), files=used,
                  cells=weights.cells if weights else 0, problems=problems)


def _same_grid(weights: Weights, grid: Grid) -> bool:
    lats, lons = getattr(weights, "grid", (None, None))
    return (lats is not None and lats.shape == grid.lats.shape
            and lons.shape == grid.lons.shape and np.allclose(lats, grid.lats)
            and np.allclose(lons, grid.lons))


def write_series(series: Series, path, *, source: str = "") -> Path:
    """``date,rain_mm`` with what produced it in '#' lines above, as the AMS CSV does."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = series.frame.copy()
    frame["date"] = frame["date"].dt.strftime("%Y-%m-%d")
    frame["rain_mm"] = frame["rain_mm"].round(3)
    header = [f"# catchment-average daily rainfall, {len(frame)} days, "
              f"{series.cells} grid cells",
              f"# from {len(series.files)} file(s): {series.files[0].name} .. "
              f"{series.files[-1].name}" if series.files else "# from no files"]
    if source:
        header.append(f"# {source}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("\n".join(header) + "\n")
        frame.to_csv(stream, index=False)
    return path


def read_series(path) -> pd.Series:
    """A series written by ``write_series`` (or any date,rain_mm CSV)."""
    frame = pd.read_csv(path, comment="#", parse_dates=["date"])
    return frame.set_index("date")["rain_mm"]
