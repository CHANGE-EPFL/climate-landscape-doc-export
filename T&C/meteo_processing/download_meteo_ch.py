"""Download and prepare gridded meteorological inputs for the Swiss sites.

Projected Swiss catchment boundaries are converted to WGS84 before requesting
Earth Engine data. Outputs are written as <site>/meteo/<site>/ relative to the
repository root, matching the model launch files.
"""

import argparse
import os
from pathlib import Path

import ee
import geemap
import numpy as np
import pandas as pd
import shapefile
import xarray as xr
from pyproj import CRS as PyprojCRS
from pyproj import Transformer
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
SHAPE_ROOT = SCRIPT_DIR / "shapefiles" / "CH"

SITE_PERIODS = {
    "2203": (2015, 2020),
    "2420": (2015, 2020),
    "2469": (2015, 2020),
    "2481": (2015, 2020),
}

# Fixed 250-m forcing grids used for the study simulations. Explicit grid
# origins and dimensions prevent software-dependent spatial realignment.
SITE_GRIDS = {
    # site: (EPSG code, first x centre, first y centre, nx, ny)
    "2203": (32632, 343811.6770772488, 5129638.422444288, 78, 41),
    "2420": (32632, 504006.15659088286, 5112770.951430169, 78, 154),
    "2469": (32632, 383961.84730189416, 5140919.38460363, 111, 113),
    "2481": (32632, 448275.63010847435, 5179251.521634782, 72, 93),
}

OUTPUT_VARIABLES = ["Pr", "Pre", "S", "U", "Ws", "Ta"]
ERA5_OUTPUT_VARIABLES = [
    "total_precipitation",
    "surface_pressure",
    "surface_solar_radiation_downwards_hourly",
    "relative_humidity",
    "wind_10m",
    "temperature_2m",
]


def initialize_earth_engine(project):
    """Initialize Earth Engine, requesting authentication when required."""
    try:
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception:
        ee.Authenticate()
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()


def shapefile_bbox_wgs84(shapefile_path):
    """Return a shapefile bounding box as WGS84 longitude/latitude."""
    reader = shapefile.Reader(str(shapefile_path))
    xmin, ymin, xmax, ymax = reader.bbox
    projection_path = Path(shapefile_path).with_suffix(".prj")

    if not projection_path.exists():
        if (
            -180 <= xmin <= 180
            and -180 <= xmax <= 180
            and -90 <= ymin <= 90
            and -90 <= ymax <= 90
        ):
            return [xmin, ymin, xmax, ymax]
        raise ValueError(
            f"{shapefile_path} has projected coordinates but no .prj file"
        )

    source_crs = PyprojCRS.from_wkt(projection_path.read_text())
    if source_crs.is_geographic:
        return [xmin, ymin, xmax, ymax]

    transformer = Transformer.from_crs(
        source_crs, "EPSG:4326", always_xy=True
    )
    corners = [
        transformer.transform(xmin, ymin),
        transformer.transform(xmin, ymax),
        transformer.transform(xmax, ymin),
        transformer.transform(xmax, ymax),
    ]
    longitudes, latitudes = zip(*corners)
    return [
        min(longitudes),
        min(latitudes),
        max(longitudes),
        max(latitudes),
    ]


def utm_crs_from_bbox(roi_bounds):
    utm_crs_list = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=roi_bounds[0],
            south_lat_degree=roi_bounds[1],
            east_lon_degree=roi_bounds[2],
            north_lat_degree=roi_bounds[3],
        ),
    )
    if not utm_crs_list:
        raise ValueError(f"Could not determine a UTM CRS for {roi_bounds}")
    return "EPSG:" + utm_crs_list[0].code, utm_crs_list[0].code


def production_grid(site_id, scale, crs):
    """Return the explicit xee grid used by the production simulation."""
    if scale != 250:
        raise ValueError(
            "The archived production grids are defined at 250 m; "
            f"received --scale {scale}."
        )

    epsg_code, x_first, y_first, nx, ny = SITE_GRIDS[site_id]
    expected_crs = f"EPSG:{epsg_code}"
    if crs != expected_crs:
        raise ValueError(
            f"UTM CRS mismatch for {site_id}: derived {crs}, "
            f"expected {expected_crs}."
        )

    x_origin = x_first - scale / 2
    y_top_centre = y_first + (ny - 1) * scale
    y_origin = y_top_centre + scale / 2
    return {
        "crs_transform": (scale, 0, x_origin, 0, -scale, y_origin),
        "shape_2d": (nx, ny),
    }


def wc_get(crs, grid):
    """Retrieve WorldClim v1 monthly temperature and precipitation."""
    worldclim = ee.ImageCollection("WORLDCLIM/V1/MONTHLY")
    worldclim_data = geemap.ee_to_xarray(
        worldclim, crs=crs, **grid
    ).load()
    return worldclim_data


def process_single_month(year, month, crs, grid, worldclim_data):
    """Download and process one calendar month."""
    if month != 12:
        start = ee.Date(f"{year:04d}-{month:02d}-01")
        end = ee.Date(f"{year:04d}-{month + 1:02d}-01")
    else:
        start = ee.Date(f"{year:04d}-{month:02d}-01")
        end = ee.Date(f"{year + 1:04d}-01-01")

    era5 = (
        ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
        .filter(ee.Filter.date(start, end))
        .select(
            [
                "temperature_2m",
                "dewpoint_temperature_2m",
                "total_precipitation",
                "surface_solar_radiation_downwards_hourly",
                "surface_pressure",
                "u_component_of_wind_10m",
                "v_component_of_wind_10m",
                "total_precipitation_hourly",
            ]
        )
    )
    era5_data = geemap.ee_to_xarray(
        era5, crs=crs, **grid
    ).load()

    era5_data["surface_solar_radiation_downwards_hourly"] /= 3600
    era5_data["surface_pressure"] /= 100
    era5_data["total_precipitation_hourly"] *= 1000

    wind_speed = np.sqrt(
        np.power(era5_data["u_component_of_wind_10m"].values, 2)
        + np.power(era5_data["v_component_of_wind_10m"].values, 2)
    )
    era5_data["u_component_of_wind_10m"].values = wind_speed
    era5_data = era5_data.rename({"u_component_of_wind_10m": "wind_10m"})
    era5_data = era5_data.drop_vars("v_component_of_wind_10m")

    dewpoint_c = era5_data["dewpoint_temperature_2m"].values - 273.15
    temperature_c = era5_data["temperature_2m"].values - 273.15
    relative_humidity = np.exp(17.625 * dewpoint_c / (243.04 + dewpoint_c)) / np.exp(
        17.625 * temperature_c / (243.04 + temperature_c)
    )
    era5_data["dewpoint_temperature_2m"].values = relative_humidity
    era5_data = era5_data.rename(
        {"dewpoint_temperature_2m": "relative_humidity"}
    )

    era5_hourly = era5_data["total_precipitation_hourly"]
    era_time = pd.to_datetime(era5_data["time"].values)
    era5_data = era5_data.assign_coords(time=era_time)
    era5_hourly = era5_data["total_precipitation_hourly"]
    era5_daily = era5_hourly.groupby("time.date").sum(dim="time")

    chirps = (
        ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
        .filterDate(start, end)
        .select("precipitation")
    )
    chirps_data = geemap.ee_to_xarray(
        chirps, crs=crs, **grid
    ).load()
    chirps_daily = chirps_data["precipitation"].assign_coords(
        time=pd.to_datetime(chirps_data["time"].values).date
    )

    era5_daily = era5_daily.transpose("date", "y", "x")
    chirps_daily = chirps_daily.transpose("time", "y", "x")
    chirps_factor = np.divide(
        chirps_daily.values,
        era5_daily.values,
        out=np.ones_like(chirps_daily.values),
        where=era5_daily.values != 0,
    )

    scaled_hourly = era5_hourly.copy()
    unique_dates = np.unique(era5_data["time"].dt.date.values)
    for day_index, day in enumerate(unique_dates):
        day_mask = era5_data["time"].dt.date == day
        scaled_hourly[day_mask, :, :] *= chirps_factor[day_index, :, :]
    era5_data["total_precipitation"] = scaled_hourly

    monthly_total = np.asarray(era5_data["total_precipitation"].sum(axis=0))
    worldclim_precipitation = np.asarray(
        worldclim_data["prec"][month - 1, :, :].values
    )
    worldclim_factor = np.divide(
        worldclim_precipitation,
        monthly_total,
        out=np.ones_like(worldclim_precipitation),
        where=monthly_total != 0,
    )
    era5_data["total_precipitation"] *= worldclim_factor

    worldclim_temperature = (
        0.1 * worldclim_data["tavg"][month - 1, :, :].values + 273.15
    )
    era5_monthly_temperature = np.mean(era5_data["temperature_2m"], axis=0)
    unique_temperature, inverse_index = np.unique(
        era5_monthly_temperature, return_inverse=True
    )
    inverse_index = inverse_index.reshape(worldclim_temperature.shape)
    aggregated_worldclim_temperature = np.zeros_like(worldclim_temperature)
    for temperature_index in range(len(unique_temperature)):
        temperature_mask = inverse_index == temperature_index
        aggregated_worldclim_temperature[temperature_mask] = np.mean(
            worldclim_temperature[temperature_mask]
        )
    temperature_adjustment = (
        worldclim_temperature - aggregated_worldclim_temperature
    )
    era5_data["temperature_2m"] = (
        era5_data["temperature_2m"] + temperature_adjustment - 273.15
    )

    output = {}
    for output_name, source_name in zip(
        OUTPUT_VARIABLES, ERA5_OUTPUT_VARIABLES
    ):
        output[output_name] = (
            era5_data[source_name]
            .sortby("y")
            .transpose("time", "x", "y")
            .rename(output_name)
            .to_dataset()
        )
    return output


def download_site(site_id, start_year, end_year, scale, overwrite=False):
    """Process one six-year Swiss site."""
    shapefile_path = SHAPE_ROOT / f"{site_id}.shp"
    if not shapefile_path.is_file():
        raise FileNotFoundError(f"Missing catchment shapefile: {shapefile_path}")

    output_dir = REPO_ROOT / site_id / "meteo" / site_id
    output_dir.mkdir(parents=True, exist_ok=True)
    period = f"{start_year}{end_year}yr"
    expected = [output_dir / f"{name}_{period}.nc" for name in OUTPUT_VARIABLES]
    expected.append(output_dir / "DEM.nc")
    if not overwrite and all(path.is_file() for path in expected):
        print(f"Skipping {site_id}: all expected download outputs already exist")
        return

    roi_bounds = shapefile_bbox_wgs84(shapefile_path)
    crs, epsg_code = utm_crs_from_bbox(roi_bounds)
    grid = production_grid(site_id, scale, crs)

    print(f"\n=== Processing Swiss site {site_id}: {start_year}-{end_year} ===")
    print(f"SHP: {shapefile_path}")
    print(f"OUT: {output_dir}")

    worldclim_data = wc_get(crs, grid)

    hydro_dem = ee.Image("WWF/HydroSHEDS/15CONDEM")
    dem_data = geemap.ee_to_xarray(
        hydro_dem, crs=crs, **grid
    ).load()
    dem_output = (
        dem_data["b1"]
        .isel(time=0)
        .sortby("y")
        .transpose("x", "y")
        .rename("DEM")
        .to_dataset()
    )
    dem_output.attrs["crs"] = epsg_code
    dem_output.to_netcdf(output_dir / "DEM.nc")

    all_data = {name: [] for name in OUTPUT_VARIABLES}
    dates = pd.date_range(
        f"{start_year}-01-01", f"{end_year}-12-01", freq="MS"
    )
    for date in dates:
        print(f"Processing {date.year}-{date.month:02d}")
        monthly_data = process_single_month(
            date.year, date.month, crs, grid, worldclim_data
        )
        for variable in OUTPUT_VARIABLES:
            all_data[variable].append(monthly_data[variable])

    for variable in OUTPUT_VARIABLES:
        print(f"Saving {variable}_{period}.nc")
        combined = xr.concat(all_data[variable], dim="time")
        combined.attrs["crs"] = epsg_code
        combined.to_netcdf(output_dir / f"{variable}_{period}.nc")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sites",
        nargs="+",
        choices=sorted(SITE_PERIODS),
        default=sorted(SITE_PERIODS),
        help="Site IDs to process; the default is all Swiss study sites.",
    )
    parser.add_argument(
        "--project",
        default=os.environ.get("GOOGLE_CLOUD_PROJECT"),
        help="Google Cloud project used for Earth Engine.",
    )
    parser.add_argument("--scale", type=int, default=250, help="Output grid size in m.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate outputs even when all expected files exist.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    initialize_earth_engine(args.project)
    for site_id in args.sites:
        start_year, end_year = SITE_PERIODS[site_id]
        download_site(
            site_id,
            start_year,
            end_year,
            args.scale,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
