"""Catchment-average rainfall from gridded netCDF (core/awap.py).

The grids here are written in the AWRA-L layout - yearly files, a 3-D daily
variable on latitude and longitude - with rain chosen so the right average is
known: every cell rains its own column number, so a catchment covering columns
2 and 3 equally averages 2.5 whatever the weighting.
"""

from __future__ import annotations

from datetime import date

import numpy as np
import pandas as pd
import pytest

netCDF4 = pytest.importorskip("netCDF4")
shapefile = pytest.importorskip("shapefile")

from core import awap                                               # noqa: E402

STEP = 0.05
LONS = 150.0 + STEP * np.arange(8)          # cell centres 150.00 .. 150.35
LATS = -24.0 - STEP * np.arange(6)          # descending, as AWRA-L stores them


def write_year(path, year, days=3, fill_cell=None, value=None):
    with netCDF4.Dataset(str(path), "w") as dataset:
        dataset.createDimension("time", None)
        dataset.createDimension("latitude", len(LATS))
        dataset.createDimension("longitude", len(LONS))
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = f"days since {year}-01-01"
        time.calendar = "standard"
        time[:] = np.arange(days)
        dataset.createVariable("latitude", "f8", ("latitude",))[:] = LATS
        dataset.createVariable("longitude", "f8", ("longitude",))[:] = LONS
        rain = dataset.createVariable("rain_day", "f4", ("time", "latitude", "longitude"),
                                      fill_value=-999.0)
        grid = np.tile(np.arange(len(LONS), dtype=float), (days, len(LATS), 1))
        if value is not None:
            grid[:] = value
        if fill_cell is not None:
            grid[:, fill_cell[0], fill_cell[1]] = -999.0
        rain[:] = grid
    return path


def write_polygon(path, rings, prj=True):
    with shapefile.Writer(str(path), shapeType=shapefile.POLYGON) as writer:
        writer.field("NAME", "C")
        writer.poly(rings)
        writer.record("CLD")
        writer.poly([[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]])
        writer.record("OTHER")
    if prj:
        path.with_suffix(".prj").write_text(
            'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",'
            '6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925]]')
    return path


def square(west, east, south, north):
    return [[(west, south), (west, north), (east, north), (east, south), (west, south)]]


@pytest.fixture
def box(tmp_path):
    # columns 2 and 3 (centres 150.10, 150.15), rows 1-3, cell edges exactly
    return write_polygon(tmp_path / "cld.shp",
                         square(150.075, 150.175, -24.175, -24.025))


def test_a_catchment_over_two_columns_averages_them(box, tmp_path):
    files = [write_year(tmp_path / "rain_day_2001.nc", 2001),
             write_year(tmp_path / "rain_day_2000.nc", 2000)]
    catchment = awap.read_catchment(box, "NAME", "CLD")
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path), catchment)
    assert list(series.frame["date"].dt.year.unique()) == [2000, 2001]   # year order
    assert np.allclose(series.frame["rain_mm"], 2.5)
    assert series.cells == 6 and not series.problems
    assert [path.name for path in series.files] == ["rain_day_2000.nc", "rain_day_2001.nc"]


def test_a_cell_half_in_the_catchment_counts_half(tmp_path):
    # columns 2 and 3 fully, column 4 half: (2 + 3 + 0.5 x 4) / 2.5 = 2.8
    shp = write_polygon(tmp_path / "half.shp", square(150.075, 150.2, -24.175, -24.025))
    write_year(tmp_path / "rain_day_2000.nc", 2000)
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path, "rain_*.nc"),
                                     awap.read_catchment(shp, "NAME", "CLD"))
    assert series.frame["rain_mm"].iloc[0] == pytest.approx(2.8, abs=0.01)


def test_centre_weighting_counts_only_the_cells_whose_centre_is_in(tmp_path):
    shp = write_polygon(tmp_path / "half.shp", square(150.075, 150.2, -24.175, -24.025))
    write_year(tmp_path / "rain_day_2000.nc", 2000)
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path, "rain_*.nc"),
                                     awap.read_catchment(shp, "NAME", "CLD"),
                                     weighting=awap.CENTRE)
    assert series.frame["rain_mm"].iloc[0] == pytest.approx(2.5)     # column 4 left out


def test_a_missing_cell_is_left_out_and_the_rest_reweighted(box, tmp_path):
    write_year(tmp_path / "rain_day_2000.nc", 2000, fill_cell=(1, 2))
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path),
                                     awap.read_catchment(box, "NAME", "CLD"))
    # three cells of column 3 at 3 mm, two of column 2 at 2 mm
    assert series.frame["rain_mm"].iloc[0] == pytest.approx((2 * 2 + 3 * 3) / 5, abs=0.01)


def test_overlapping_files_keep_the_first_and_say_so(box, tmp_path):
    write_year(tmp_path / "rain_day_2000.nc", 2000, value=1.0)
    write_year(tmp_path / "rain_day_2000_v2.nc", 2000, value=9.0)
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path),
                                     awap.read_catchment(box, "NAME", "CLD"))
    assert len(series.frame) == 3 and np.allclose(series.frame["rain_mm"], 1.0)
    assert any("more than one file" in problem for problem in series.problems)


def test_a_catchment_off_the_grid_is_refused(tmp_path):
    shp = write_polygon(tmp_path / "far.shp", square(120.0, 120.1, -30.1, -30.0))
    write_year(tmp_path / "rain_day_2000.nc", 2000)
    with pytest.raises(awap.AwapError, match="outside the grid|no rainfall"):
        awap.catchment_rainfall(awap.rainfall_files(tmp_path),
                                awap.read_catchment(shp, "NAME", "CLD"))


def test_a_field_that_is_not_there_names_the_ones_that_are(box):
    with pytest.raises(awap.AwapError, match="has no field 'ID' \\(it has NAME\\)"):
        awap.read_catchment(box, "ID", "1")


def test_a_projected_polygon_is_brought_back_to_latitude_and_longitude(tmp_path):
    pyproj = pytest.importorskip("pyproj")
    to_mga = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:28356", always_xy=True)
    ring = [to_mga.transform(lon, lat) for lon, lat in square(150.075, 150.175,
                                                               -24.175, -24.025)[0]]
    shp = write_polygon(tmp_path / "mga.shp", [ring], prj=False)
    shp.with_suffix(".prj").write_text(pyproj.CRS("EPSG:28356").to_wkt("WKT1_ESRI"))
    catchment = awap.read_catchment(shp, "NAME", "CLD")
    west, south, east, north = catchment.bounds
    assert west == pytest.approx(150.075, abs=1e-6) and north == pytest.approx(-24.025, abs=1e-6)


def test_the_series_round_trips_through_its_csv(box, tmp_path):
    write_year(tmp_path / "rain_day_2000.nc", 2000)
    series = awap.catchment_rainfall(awap.rainfall_files(tmp_path),
                                     awap.read_catchment(box, "NAME", "CLD"))
    path = awap.write_series(series, tmp_path / "out" / "cld_rain.csv", source="test")
    text = path.read_text(encoding="utf-8")
    assert text.startswith("# catchment-average daily rainfall, 3 days, 6 grid cells")
    back = awap.read_series(path)
    assert back.index[0] == pd.Timestamp(date(2000, 1, 1)) and back.iloc[0] == 2.5
