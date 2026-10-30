# Figure 2

Reproduce Figure 2 from the two CSV files in `data/`. No model outputs,
external project folders, or network access are needed for plotting.

## Run

With Python 3.12 installed, open a terminal in this folder and run:

```powershell
python -m pip install -r requirements.txt
python plot_fig2.py
```

The script also works when launched by its full path from another directory.
PNG, PDF, SVG and TIFF files are saved in `output/`; rerunning replaces them.
Figure dimensions are 183 x 140 mm, with 600 dpi raster output. Use Arial,
the original figure font, to retain the same text appearance.

## Data

`ensemble_responses.csv` contains 300 rows: 12 catchments x 5 metrics x
5 simulations (Current plus four future climate models). The ratios and C-Q
slopes are already calculated; the plotting script reads these values directly.

| Panel | `metric` | Meaning of `value` |
|---|---|---|
| a | `p_ratio` | Future/current precipitation ratio |
| b | `Qlat_out` | Future/current lateral water outflow ratio |
| c | `Q_sim_mmh` | Future/current outlet discharge ratio |
| d | `DOC_flux` | Future/current DOC export ratio |
| e | `cq_slope_b` | Fitted slope of log10(DOC concentration) against log10(discharge) |

All five metrics are dimensionless. Current ratios equal 1; C-Q slopes are
absolute fitted values for both current and future simulations. `scenario`
identifies Current or the climate model; `scenario_group` distinguishes Current
from Future. `normalization` records whether a value is a ratio or a log-log
slope. Gauge IDs must be read as strings to preserve leading zeros.

`sites.csv` supplies the row order, gauge IDs, US1-US8 / CH1-CH4 labels, and
the precomputed `aridity_modis` values plotted as AI in panel f. No aridity
estimation is performed by this script.

The model-to-legend mapping is:

| Climate model | Figure legend |
|---|---|
| FGOALS-g3 | Low warming |
| NorESM2-MM | Moderate warming |
| ACCESS-CM2 | High warming |
| KACE-1-0-G | Very high warming |

For each catchment, coloured points show the four future simulations; the solid
line is their median and the dashed line is the current value. The All boxes
pool the 48 future catchment-model values, with default 1.5-IQR whiskers and
all values overlaid as points. Their dashed lines show the median current value
across catchments. The script calculates these ordinary plotting summaries
and uses fixed seeds for reproducible point jitter. No separate boxplot,
median, or jitter-coordinate files are needed.

The values were taken from the existing manuscript Figure 2 analysis. See the
manuscript Methods for the simulation periods and scientific definitions.
If you use these data or scripts, please cite our manuscript and the original
datasets referenced therein.
