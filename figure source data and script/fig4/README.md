# Figure 4

Reproduce the climate-only and climate-plus-landscape regression comparison
from two local CSV files. The workflow is **cleaned six-year summaries ->
analysis -> plotting**. The script calculates changes, selects predictors,
fits regressions and runs leave-one-site-out validation. Intermediate tables
remain in memory; only figure files are saved.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_fig4.py
```

Paths are relative to the script, so it can also be launched from another
directory. No model outputs, external project folders or network access are
needed. PNG, PDF, SVG and TIFF versions are saved in `output/`; rerunning
replaces them. The original 230 x 168 mm layout and 600 dpi raster resolution
are retained. Use Arial for the same text appearance.

## Data

| File | Contents |
|---|---|
| `scenario_summaries.csv` | 48 records: 12 catchments x four climate models, with contemporary and future temperature, precipitation and DOC export summaries |
| `sites.csv` | 12 catchments, their US1-US8 / CH1-CH4 labels and three landscape attributes |

Read `gauge_id` as a string to preserve leading zeros. `future_model` identifies
the climate projection. Contemporary values repeat across the four projections
for a given catchment. Unlike the annual process-attribution analysis in
Figure 3, each Figure 4 record summarizes a full six-year simulation.

| Column | Meaning and units |
|---|---|
| `current_t_degC`, `future_t_degC` | Mean air temperature, degrees Celsius |
| `current_pr_mm_yr`, `future_pr_mm_yr` | Mean annual precipitation, mm yr^-1 |
| `current_doc_flux_mgC_m2_h`, `future_doc_flux_mgC_m2_h` | Mean DOC-export rate, mg C m^-2 h^-1 |
| `mean_elevation_m` | Catchment mean elevation, m |
| `mean_distance_to_channel_m` | Mean flow-path distance to the channel, m |
| `slope_p90_p10_deg` | Difference between the 90th and 10th percentiles of slope angle, degrees |

Temperature and DOC-export summaries retain the original analysis convention:
means of the 72 monthly means. Precipitation is the six-year accumulated total
divided by six. DOC export here is an area-normalized rate, not an annual total.
The script calculates future minus contemporary changes and relative DOC-export
change as `(future - current) / current * 100`.

## Analysis

Both temperature and precipitation change are always included with an intercept.
For each DOC response (absolute and relative change), the script compares all
eight subsets of the three landscape predictors. It selects the minimum-AICc
model, breaking ties by BIC and then predictor count. The Gaussian-likelihood
parameter count includes the intercept, slopes and residual variance.

Ordinary least squares is fitted in the original variable units. Standardized
coefficients are calculated as `coefficient * SD(predictor) / SD(response)`,
using sample standard deviations (`ddof=1`). All 48 records are retained.

Leave-one-site-out (LOSO) validation withholds all four climate projections of
one catchment together, refits coefficients on the remaining 11 catchments and
predicts the withheld catchment. Predictor sets are selected on the full dataset
and remain fixed across folds, matching the manuscript analysis; selection is
not repeated within each fold. LOSO R-squared is calculated from all 48 held-out
predictions using the pooled response mean.

## Panels and interpretation

- **a, b:** climate-only predictions of absolute and relative changes in DOC export.
- **c, d:** predictions from the corresponding AICc-selected landscape models.
- **e:** standardized coefficients for the four models. Hatched bars indicate
  climate-only models; filled bars indicate climate-plus-landscape models.
  Predictors absent from a selected model have no bar.

Scatter points show full-data fitted predictions, with one point per catchment
and climate projection. Their x coordinates are process-model outputs, not field
observations. The annotated R-squared and overall F-test P value describe the
full-data regression; the separately labelled LOSO R-squared describes held-out
predictions.

**The P values beside the filled bars in panel e come from separate univariate
regressions of each predictor against the response.** They are not the partial
coefficient P values from the multivariable models. This preserves the original
figure's analysis and annotations.

The script prints the selected predictors and LOSO scores to the terminal.
No separate CSV files of fitted values, coefficients, model rankings or plotting
coordinates are produced. See the manuscript Methods for the simulation and
landscape-attribute calculations.

If you use these data or scripts, please cite our manuscript and the original
datasets referenced therein.
