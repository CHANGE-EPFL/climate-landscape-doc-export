# Figure 3

Reproduce the six-panel process-attribution figure from the CSV files in
`data/`. No simulation outputs, external project folders or downloads are
needed. The script fits the regressions and calculates Shapley contributions
from the supplied annual anomalies before drawing the figure. Intermediate
results are kept in memory and are not saved as additional CSV files.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_fig3.py
```

The script also works when launched by its full path from another directory.
It saves PNG, PDF, SVG and TIFF versions in `output/`; rerunning replaces them.
The figure is 183 x 136 mm, with 600 dpi raster output. Use Arial, the original
figure font, for the same text appearance.

## Data

| CSV file | Contents |
|---|---|
| `annual_anomalies.csv` | 288 annual records: 12 catchments x four climate models x six years. Each record contains both DOC responses and all three process predictors. |
| `sites.csv` | Gauge IDs, US1-US8 / CH1-CH4 labels and vegetation groups for the 12 catchments |

Gauge IDs must be read as strings to preserve leading zeros. `Scenario` and
`YearIndex` identify the climate projection and paired simulation year (1-6)
for each record. `Scenario` retains the original analysis identifiers; the
climate-model name identifies the projection. The five numerical columns are:

| Column | Meaning |
|---|---|
| `dDOC_concentration` | Change in annual mean DOC concentration |
| `dDOC_export` | Change in annual DOC export |
| `dPlant_input` | Change in the catchment-mean plant carbon-input indicator |
| `dQlat_out` | Change in the catchment-mean lateral water-outflow indicator |
| `dDOC_uptake` | Change in the catchment-mean microbial DOC-uptake indicator |

The predictor values retain their original analysis scale and are standardized
within each catchment before fitting. They have not been converted to annual
totals in this package.

## Analysis

For each catchment and response, the script standardizes the 24 annual records
using sample standard deviations (`ddof=1`) and fits ordinary least squares
with an intercept and all three predictors. Predictions are converted back to
the response's original units for panels a and b. Shapley decomposition uses
all predictor subsets to partition the fitted R-squared for panels c and d.
Separate standardized univariate regressions supply the comparison coefficients
in panels e and f. All 288 annual records are retained.

## Units and interpretation

- DOC-export anomalies are in mg C m^-2 yr^-1; concentration anomalies are
  in mg C L^-1. Anomalies are future minus the corresponding contemporary year.
- In panels a and b, **Observed** means the response supplied to the regression:
  an anomaly calculated from the process-based simulations. It does not mean
  an independent field measurement. **Predicted** is the fitted regression value.
  These panels show in-sample fits of catchment-specific regressions.
- Pooled R-squared is 1 minus the summed squared prediction residuals divided
  by the summed squared deviations of response anomalies from their respective
  catchment means. This is not squared Pearson correlation or a cross-validation
  score. It is recalculated from the fitted annual predictions.
- Shapley shares are dimensionless fractions of explained variance, normalized
  to sum to one within each catchment and response. The ternary vertices represent
  plant input (left), lateral flow (right) and DOC uptake (top).
- Coefficients in panels e and f are dimensionless standardized regression
  coefficients. Open symbols show univariate estimates; filled symbols show
  three-variable estimates. Red lines join pairs with opposite, nonzero signs.
  The script counts these sign changes from the newly fitted coefficients.
- Colours identify catchments in panels a-d. Marker shapes identify vegetation
  groups across all panels. Fixed vertical offsets separate coefficient pairs;
  no jitter-coordinate tables are needed.

The manuscript Methods describe the simulations, regression fitting and
Shapley decomposition. If you use these data or scripts, please cite our
manuscript and the original datasets referenced therein.
