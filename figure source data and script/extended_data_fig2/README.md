# Extended Data Figure 2

Recalculate the robustness tests from cleaned catchment-model data, then draw
the two-panel figure. The script selects the regression predictors and runs
10,000 site-block bootstrap samples and 10,000 site-level permutations for
each response. Replicates and summaries remain in memory; no intermediate
CSV files are saved.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_extended_data_fig2.py
```

This figure takes longer than a simple redraw because all four sets of 10,000
replicates are recalculated. Progress messages appear after each set finishes.
Paths are relative to the script. No simulation outputs, other figure folders
or network access are needed.

PNG, PDF, SVG and TIFF versions are saved in `output/`; rerunning replaces them.
The original layout is 219 x 91 mm, with 600 dpi raster output. Use Arial for
the same text appearance.

## Basic data

`data/regression_inputs.csv` is the only input table. It contains 48 records:
12 catchments x four climate models. Each record summarizes a six-year simulation
comparison, rather than a single year.

| Column | Meaning and units |
|---|---|
| `gauge_id` | Catchment identifier, read as a string to preserve leading zeros |
| `future_model` | Climate projection identifier |
| `delta_t_degC` | Change in mean air temperature, degrees Celsius |
| `delta_pr_mm_yr` | Change in mean annual precipitation, mm yr^-1 |
| `current_doc_flux_mgC_m2_h` | Contemporary mean DOC-export rate, mg C m^-2 h^-1 |
| `delta_doc_flux_mgC_m2_h` | Future minus contemporary mean DOC-export rate, mg C m^-2 h^-1 |
| `mean_elevation_m` | Catchment mean elevation, m |
| `mean_distance_to_channel_m` | Mean flow-path distance to the channel, m |
| `slope_p90_p10_deg` | Difference between the 90th and 10th percentiles of slope angle, degrees |

These are the cleaned analysis inputs used for the manuscript regressions.
Relative DOC-export change is calculated as absolute change divided by the
contemporary rate, multiplied by 100. All 48 records are retained. The script
checks for missing values, duplicate records and invalid contemporary fluxes.

## Regression selection

For each response, ordinary least squares includes an intercept and both
climate predictors. All eight subsets of the three landscape predictors are
compared by AICc, with BIC and predictor count used to break ties. The likelihood
parameter count includes the intercept, slopes and residual variance.

The selected models are:

- Absolute change: temperature change, precipitation change, mean distance
  to channel and slope range.
- Relative change: the same predictors plus mean elevation.

Selection is performed once using the complete input table. These predictor
sets remain fixed during bootstrap and permutation analysis, matching the
manuscript workflow. The package does not read saved selection results or any
files from Figure 4.

## Panel a: site-block bootstrap

Each replicate draws 12 catchments with replacement. All four climate-model
records of every sampled catchment are retained together, including repeated
copies when a catchment is drawn more than once. This preserves the within-site
grouping of climate scenarios.

The response and predictors are standardized within each replicate using sample
standard deviations (`ddof=1`), then OLS is fitted with an intercept. The raw
design is also checked for full rank before accepting a replicate. Rank-deficient
or unstandardizable samples are redrawn until 10,000 valid replicates are obtained.

Random seeds are 20260623 for absolute change and 20260624 for relative change.
Catchments are sorted by their identifiers before sampling. The boxes span the
25th-75th percentiles, white lines and dots indicate medians, and whiskers use
the 2.5th and 97.5th percentile settings. Outliers are not drawn.

## Panel b: site-topography pairing test

Selected landscape attributes are shuffled together among the 12 catchments,
preserving their joint structure. The four climate/response records associated
with each catchment stay together; only their pairing with the landscape vector
changes. The full regression is refitted for each of 10,000 permutations.
Each response uses its own seed listed above, reset independently of the
bootstrap calculation.

Grey boxes show the distribution of fitted R-squared under random pairings,
with quartile boxes, median lines and 2.5th/97.5th percentile whisker settings.
Black diamonds show the fitted R-squared with the observed pairing; dashed
lines show the corresponding climate-only fit. These are in-sample fit
statistics, not cross-validation scores.

The script recalculates coefficients, R-squared values and all plotting
summaries. No bootstrap-replicate, permutation-replicate, coefficient-summary
or boxplot-statistics CSV files are included or generated.

See the manuscript Methods for the simulation summaries and landscape
attributes. If you use these data or scripts, please cite our manuscript and
the original datasets referenced therein.
