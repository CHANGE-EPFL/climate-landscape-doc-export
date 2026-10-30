# T&C catchment simulations

This repository contains the model configuration, spatial inputs and initial
states required to reproduce the contemporary and late-century T&C simulations
for eight catchments in the contiguous United States (CONUS) and four
catchments in Switzerland. Large model-output files are not included.

## Repository structure

```text
T&C/
|-- meteo_processing/        Meteorological download and preprocessing scripts
|-- 01054200/                CONUS site configuration
|-- ...
|-- 12447390/                CONUS site configuration
|-- 2203/                    Swiss site configuration
|-- ...
`-- 2481/                    Swiss site configuration
```

Each site directory contains the model code, spatial inputs, parameter files,
initial states, climate Factors of Change and MATLAB launch scripts. The
provided initial states allow the simulations reported in the study to be run
without repeating the computationally intensive spin-up procedure.

## Requirements

The meteorological download scripts require Python, a Google Earth Engine
account and a Google Cloud project with Earth Engine enabled. The principal
Python packages are `earthengine-api`, `geemap`, `xee`, `xarray`, `numpy`,
`pandas`, `pyshp`, `pyproj` and `netCDF4`.

The preprocessing and simulations require MATLAB. The Mapping Toolbox and
Statistics and Machine Learning Toolbox are used for coordinate transformations
and nearest-neighbour mapping.

## Reproduction workflow

### 1. Download the gridded meteorological data

Enter `meteo_processing` and run the CONUS and Swiss download scripts:

```powershell
cd meteo_processing
python download_meteo_us.py --project YOUR_EARTH_ENGINE_PROJECT
python download_meteo_ch.py --project YOUR_EARTH_ENGINE_PROJECT
```

The first run may start the Earth Engine authentication procedure. By default,
each script processes all sites in its region. Individual sites can be selected,
for example:

```powershell
python download_meteo_us.py --sites 01434025 --project YOUR_EARTH_ENGINE_PROJECT
python download_meteo_ch.py --sites 2203 --project YOUR_EARTH_ENGINE_PROJECT
```

Use `--overwrite` to replace existing download outputs. Catchment shapefiles
needed for the requests are provided in `meteo_processing/shapefiles`.

The scripts retrieve the following products through Google Earth Engine:

- [ERA5-Land hourly](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY) for the hourly meteorological fields;
- [CHIRPS daily precipitation](https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY) and [WorldClim v1 monthly climatology](https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_MONTHLY) for precipitation and temperature adjustment; and
- the [HydroSHEDS hydrologically conditioned DEM](https://developers.google.com/earth-engine/datasets/catalog/WWF_HydroSHEDS_15CONDEM).

The resulting files are written directly to
`<site>/meteo/<site>/`, using the periods and grids defined for each site in the
download scripts.

### 2. Create the model-ready meteorological inputs

From MATLAB, with `meteo_processing` as the current directory, run:

```matlab
meteo_data_postprocessing
```

This creates the derived radiation, cloud and contemporary atmospheric CO2
inputs required by T&C. A selected site can be processed with:

```matlab
meteo_data_postprocessing("01434025")
```

To regenerate an existing set of derived files, use:

```matlab
meteo_data_postprocessing("01434025", true)
```

### 3. Create the future atmospheric CO2 input

After `Ca.nc` has been generated for the sites, run:

```matlab
Create_future_Ca_file
```

This creates `<site>/meteo/<site>/Ca_585_88_93.nc`, containing the annual
SSP5-8.5 atmospheric CO2 concentrations prescribed for 2088-2093.

### 4. Run the simulations

Enter the directory of the required catchment and run its MATLAB launch script.
For example:

```matlab
cd(fullfile('..', '01434025'))
TC_2D_hpc_BG_01434025
```

`TC_2D_hpc_BG_<site>.m` launches the contemporary simulation. The four future
simulations are launched with:

```text
TC_2D_hpc_BG_<site>_ACCESS.m    ACCESS-CM2
TC_2D_hpc_BG_<site>_FGOALS.m   FGOALS-g3
TC_2D_hpc_BG_<site>_KACE.m     KACE-1-0-G
TC_2D_hpc_BG_<site>_NorESM2.m  NorESM2-MM
```

The associated `Factor_<site>_*.mat` files contain the monthly climate Factors
of Change used to perturb precipitation and air temperature. Additional
`_5var` files at site `01434025` correspond to the meteorological-forcing
sensitivity experiment.

## Additional forcing data supplied with the repository

For the CONUS sites, `<site>/meteo/<site>/<site>.nc` contains the `DateTime` and
`Rainf` variables extracted for the corresponding six-year period from
[CAMELSH](https://doi.org/10.5281/zenodo.15066778). These hourly catchment-mean
precipitation data are used by the CONUS run templates for the final
precipitation correction.

For the Swiss sites, `<site>_pr_corr_intent.mat` contains the correction factors
and fallback precipitation derived from daily catchment-mean precipitation in
[CAMELS-CH](https://doi.org/10.5281/zenodo.7784632). These files are read by the
Swiss run templates.

The climate Factors of Change supplied as `Factor_<site>_*.mat` were derived
from NASA NEX-GDDP-CMIP6 projections under SSP5-8.5. The raw NEX-GDDP-CMIP6
archive is not redistributed here.

## Citation

If you use this code or the accompanying data, please cite our manuscript and
the associated datasets referenced therein.
