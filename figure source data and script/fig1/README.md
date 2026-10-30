# Figure 1

Reproduce the manuscript's six-panel model-validation figure using the CSV
files in `data/`. No simulation outputs, external data folders or downloads
are needed.

## Environment

The figure was reproduced with Python 3.12 and the package versions listed in `requirements.txt`.
Install these dependencies with:

```powershell
python -m pip install -r requirements.txt
```

## Run

Open a terminal in this folder and run:

```powershell
python plot_fig1.py
```

The script saves PNG, PDF, SVG and TIFF versions in `output/` (183 x 146 mm,
600 dpi for raster images). It can also be launched from another directory:
paths are relative to the script. Keep `plot_helpers.py` beside `plot_fig1.py`.
Use Arial, the original figure font, for the same text appearance.

## Data

| Panel | CSV file | Contents |
|---|---|---|
| a | `gpp_distribution_us8.csv` | Temporally averaged pixel GPP values for the eight CONUS catchments, grouped by MODIS, SLOPE and Simulation |
| b | `gpp_distribution_ch4.csv` | Corresponding MODIS and simulated GPP values for the four Swiss catchments |
| a, b | `gpp_spatial_correlations_us8.csv`, `gpp_spatial_correlations_ch4.csv` | Precomputed spatial correlations displayed in the legends |
| c | `validation_relative_difference_kge_boxplot.csv` | Catchment-level relative differences (%) and KGE values |
| d, e, f | `doc_1to1_12sites.csv` | Matched observed and simulated DOC concentrations and fluxes, and observed discharge |
| d, e | `scatter_statistics.csv` | Precomputed R-squared and sample counts for the scatter plots |
| f | `site_labels.csv` | Catchment identifiers and their US1-US8 / CH1-CH4 labels |
| Map insets | `map_boundaries.csv`, `map_stations.csv` | Country outlines and station positions already projected for plotting |

The plotting script calculates GPP density curves and DOC half-violin curves
from the supplied values. It also calculates the means and medians shown on
the violins and the quartiles/whiskers for the boxplots. These are ordinary
plotting operations; no separate curve, box-statistics or jitter tables are
required. Point jitter uses a fixed seed to preserve the original appearance.

KGE, relative differences, spatial correlations and scatter-plot R-squared
are read directly from their CSV tables and are not recalculated.

## Units and interpretation

- GPP `value`: g C m^-2 d^-1. The legend correlations use spatially matched
  pixels; do not infer pixel pairing from the order of the pooled GPP tables.
- `doc_obs`, `doc_sim`: mg C L^-1.
- `doc_flux_obs`, `doc_flux_sim`: mg C m^-2 h^-1.
- `q_obs_mmh`: mm h^-1. Scatter-point colours use log10(Q_sim/Q_obs), with
  Q_sim obtained from simulated DOC flux divided by concentration. Ratios
  are clipped to 0.1-10 for the colour scale.
- Scatter `r2`: squared Pearson correlation, using linear concentration
  values in panel d and log10 flux values in panel e.
- Gauge IDs are strings; preserve their leading zeros. Missing numeric values
  are blank. Nonfinite pairs are excluded as in the original figure, and flux
  scatter plots require positive values for logarithmic axes.
- Map `x`, `y`: metres in EPSG:5070 (US) or EPSG:2056 (Switzerland).
  Boundaries use WKT geometry and include the coordinate-system identifier.

Original observational products and analysis methods are described in the
manuscript. Country outlines are derived from Natural Earth administrative
boundaries (public domain: https://www.naturalearthdata.com/about/terms-of-use/).
Station positions come from the study's CAMELS-US / CAMELS-CH-Chem inputs.

If you use these data or scripts, please cite our manuscript and the original
datasets referenced therein.
