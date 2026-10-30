# Extended Data Figure 1

Reproduce the catchment climate/vegetation context and climate-projection
comparison from two local CSV files. The workflow is **cleaned catchment-level
data -> display analysis -> plotting**. Median markers, covariance ellipses
and axis limits are calculated in the script. No intermediate CSV files are
written.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_extended_data_fig1.py
```

Paths are relative to the script, so it can also be launched from another
directory. No simulation files, other figure folders or downloads are needed.
PNG, PDF, SVG and TIFF files are saved in `output/`; rerunning replaces them.
The figure retains its original 183 x 80 mm layout and 600 dpi raster resolution.
Use Arial for the same text appearance.

## Data

### `catchment_climate.csv`

One row for each of the 12 study catchments, used in panel a:

| Column | Meaning |
|---|---|
| `gauge_id` | Catchment identifier; read as a string to preserve leading zeros |
| `region` | `CONUS` or `Alps_CH`, determining the marker shape |
| `biome_group` | Vegetation group, determining the point colour |
| `ta_mean_degC` | Contemporary mean air temperature, degrees Celsius |
| `pr_total_mm_yr` | Contemporary annual precipitation, mm yr^-1 |

These are the cleaned contemporary climate summaries used in the original
figure. The table excludes source-file paths, quality-control bookkeeping and
variables not used in the plot.

### `climate_changes.csv`

372 rows: 31 climate-model configurations x 12 catchments, used in panel b:

| Column | Meaning |
|---|---|
| `model_member` | Climate-model configuration identifier |
| `site_id` | Catchment identifier matching `gauge_id` in the first table |
| `dta_K` | Projected change in mean air temperature, K |
| `dpr_percent` | Projected precipitation change, percent of the reference value |

Changes compare 2081-2100 with 1995-2014 under SSP5-8.5. They retain the original
catchment aggregation: monthly climatologies were weighted by month length,
then climate grid cells were combined using their fractions of catchment model
points. Temperature change is future minus reference temperature; precipitation
change is `(future / reference - 1) * 100`. This package starts from those
cleaned basin-level changes and does not require the upstream monthly climate
grids or change-factor files.

All 12 catchments and all 372 climate records are retained. The script validates
complete catchment coverage and finite values before plotting.

## Analysis and interpretation

Panel a uses colour for vegetation and shape for region. Panel b highlights
the four projections used in the study:

| Climate model | Figure legend |
|---|---|
| FGOALS-g3 | Low warming |
| NorESM2-MM | Moderate warming |
| ACCESS-CM2 | High warming |
| KACE-1-0-G | Very high warming |

The highlighted 48 records are selected directly from the complete table;
there is no duplicate four-model CSV. The remaining 324 records are grey.
These four configurations are fixed study choices, not the output of a new
ranking or selection calculation during plotting.

The script calculates an 80% bivariate-normal covariance contour for each
highlighted model and one for all 372 model-catchment records. Covariance is
normalized by n (`bias=True`), matching the original figure, and the contour
radius is `sqrt(-2 * log(1 - 0.80))`. Contours are centred on arithmetic means;
they describe the spread of climate responses and are not confidence intervals
for the ensemble mean. Each model-catchment record receives equal weight.

The coloured diamonds mark coordinate-wise medians across the 12 catchments
for each highlighted model. They need not coincide with the ellipse centres.
Ellipse coordinates, median values and axis limits are calculated in memory.
There are no separate curve, median, selected-model or ranking tables.

See the manuscript Methods for the climate dataset, catchment aggregation
and rationale for the four study projections. If you use these data or scripts,
please cite our manuscript and the original datasets referenced therein.
