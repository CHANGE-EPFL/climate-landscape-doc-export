# Figure 5

Reproduce the CONUS DOC-export vulnerability maps from cleaned basin predictors
and a small regression training table. The script refits the selected regression,
calculates basin predictions, and draws the maps and component scatter plot.
Intermediate prediction and summary tables remain in memory.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_fig5.py
```

Paths are relative to the script. No model outputs, other figure folders,
Earth Engine connection or downloads are needed. PNG, PDF, SVG and TIFF files
are saved in `output/`; rerunning replaces them. The original 230 x 198 mm
layout is retained, with 600 dpi raster output. Use Arial for matching text.

## Basin inputs

`data/basin_inputs.csv` contains the 38,639 basins included in the figure and
only seven columns:

| Column | Meaning |
|---|---|
| `coarse_basin_id` | Basin identifier linking the table to the spatial grids |
| `mean_dTa_degC_4models` | Mean projected temperature change across four climate models, degrees Celsius |
| `mean_dPr_mm_yr_4models` | Mean projected annual precipitation change across four climate models, mm yr^-1 |
| `mean_elevation_m` | Mean basin elevation, m |
| `mean_distance_to_channel_m` | Mean flow-path distance to the channel, m |
| `slope_p90_p10_deg` | Difference between the 90th and 10th percentiles of slope angle, degrees |
| `fenneman_division_id` | Physiographic division used to group basins in panel d |

The table retains the basin screening used in the final figure, including the
large-lake and no-stream exclusions. It contains the 38,639 basins with finite
mapped components out of the original 39,850 basins. Excluded basins are not
included as empty rows. Basin area, baseflow index, unused climate variables,
individual climate-model outputs and precomputed prediction components are not
needed by this figure and are omitted.

The four highlighted Fenneman divisions are Appalachian Highlands (1), Interior
Highlands (3), Interior Plains (4) and Rocky Mountain System (8). Other division
IDs remain in the basin table so their basins contribute to the full scatter
and overall percentages, but they do not receive an ellipse or coloured outline.

## Regression and prediction

`data/regression_training.csv` contains 48 records: 12 study catchments x four
climate models. Gauge IDs are strings with leading zeros preserved. The table
contains temperature change (`delta_t_degC`), annual precipitation change
(`delta_pr_mm_yr`), current DOC-export rate (`current_doc_flux_mgC_m2_h`), its
absolute change (`delta_doc_flux_mgC_m2_h`), and the same three landscape
attributes listed above. `future_model` identifies each climate projection.
The DOC-export rates are in mg C m^-2 h^-1 and retain the six-year analysis
summaries used for Figure 4.

The response is calculated as absolute DOC-export change divided by current
DOC-export rate, multiplied by 100. Ordinary least squares is then fitted with
an intercept and the five predictors. This refits the relative-change model
selected in Figure 4; predictor selection is not repeated here. The package
includes its own training data and does not load Figure 4 files.

The projection separates the fitted response into:

- **Climate pressure:** temperature coefficient x temperature change, plus
  precipitation coefficient x precipitation change.
- **Landscape component:** intercept plus the three landscape terms.
- **Predicted relative change:** the sum of those two components.

Because the regression is additive and linear, applying it to the four-model
mean climate changes gives the mean of the four model-specific predictions.
The script recalculates all three components; none are stored as redundant
columns in the input CSV.

## Map files

The CSV values alone do not define basin locations. Three small spatial files
in `data/maps/` supply exactly the geographic information needed to draw them:

- `basin_ids_1.tif` and `basin_ids_2.tif`: basin-ID grids in EPSG:5070. These
  retain the exact nearest-neighbour samples used for the original figure,
  with each tile's longest dimension limited to 1,200 cells. Zero denotes
  excluded basins or background. They are display grids, not the native
  250 m terrain data, and contain no pre-rendered DOC values or colours.
- `boundaries.gpkg`: the CONUS outline (`conus_boundary` layer) and the four
  highlighted Fenneman outlines (`fenneman_divisions` layer), also in EPSG:5070.

Basin IDs derive from the study's aggregated MERIT-Basins map. The CONUS outline
derives from TIGER/2018/States, and the regional outlines from the Fenneman
physiographic divisions. See the manuscript Methods for the original sources
and spatial processing. The supplied geometry and screening reproduce the
figure without requiring the full upstream geospatial datasets.

## Plotting summaries

Boxplots separate positive and negative basin values, with median labels and
1.5-IQR whiskers. All three maps share a symmetric colour scale whose bound
is the largest 98th percentile of absolute values among the three components.
Extreme values saturate at the colour limits; the underlying values are unchanged.

Panel d plots all 38,639 basins. Percentages and quadrant counts are based on
basin counts, not area weighting. Regional centres are arithmetic means, and
ellipses show the 50% bivariate-normal covariance contours calculated from
basin values; they are not confidence intervals for the regional mean. These
summaries use all basins, including points outside the displayed axis limits.

No separate ellipse, quadrant, boxplot, prediction or pixel-coordinate CSV
files are written. If you use these data or scripts, please cite our manuscript
and the original datasets referenced therein.
