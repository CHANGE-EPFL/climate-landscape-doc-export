"""Analyse six-year catchment summaries and reproduce Figure 4.

Read the two local CSV files, select landscape predictors using AICc, fit the
regressions and run leave-one-site-out validation before plotting. Intermediate
analysis tables stay in memory; only the figure files are written to disk.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
import re

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Resolve all input and output paths relative to this script.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
SITE_LABELS_CSV = DATA_DIR / "sites.csv"

SITE_SCATTER_EDGE_COLOR = "#3F3F3F"
SITE_SCATTER_ALPHA = 0.73
SITE_SCATTER_EDGE_WIDTH = 0.22
PREDICTORS = ["delta_t_degC", "delta_pr_mm_yr"]
TOPO_PREDICTORS = ["mean_elevation_m", "mean_distance_to_channel_m", "slope_p90_p10_deg"]
ABSOLUTE_DOC_FLUX_COL = "delta_doc_flux_mgC_m2_h"
RELATIVE_DOC_FLUX_COL = "relative_delta_doc_flux_percent"
COMPARISON_PREDICTOR_ORDER = ["delta_pr_mm_yr", "delta_t_degC", *TOPO_PREDICTORS]
PREDICTOR_LABELS = {
    "delta_t_degC": r"$\Delta$T",
    "delta_pr_mm_yr": r"$\Delta$P",
    "mean_elevation_m": "Mean elevation",
    "mean_distance_to_channel_m": "Mean distance to channel",
    "slope_p90_p10_deg": "Slope range",
}


def load_analysis_data() -> pd.DataFrame:
    """Calculate climate and DOC changes from the six-year scenario summaries."""
    data = pd.read_csv(DATA_DIR / "scenario_summaries.csv", dtype={"gauge_id": str},
                       float_precision="round_trip")
    sites = pd.read_csv(SITE_LABELS_CSV, dtype={"gauge_id": str}, float_precision="round_trip")
    data = data.merge(sites[["gauge_id", *TOPO_PREDICTORS]], on="gauge_id",
                      how="left", validate="many_to_one", sort=False)
    for key in ("t_degC", "pr_mm_yr", "doc_flux_mgC_m2_h"):
        data["delta_" + key] = data["future_" + key] - data["current_" + key]
    baseline = data["current_doc_flux_mgC_m2_h"]
    if not np.isfinite(baseline).all() or (baseline <= 0).any():
        raise ValueError("Relative DOC change requires positive finite current DOC fluxes.")
    data[RELATIVE_DOC_FLUX_COL] = data[ABSOLUTE_DOC_FLUX_COL] / baseline * 100.0
    required = [*PREDICTORS, *TOPO_PREDICTORS, ABSOLUTE_DOC_FLUX_COL, RELATIVE_DOC_FLUX_COL]
    if not np.isfinite(data[required].to_numpy(dtype=float)).all():
        raise ValueError("The analysis inputs contain missing or nonfinite values.")
    return data


def site_display_labels(gauge_ids: list[str]) -> dict[str, str]:
    """Read US1-US8 and CH1-CH4 labels, retaining leading zeros in gauge IDs."""
    sites = pd.read_csv(SITE_LABELS_CSV, dtype={"gauge_id": str}).set_index("gauge_id")
    return sites.loc[gauge_ids, "plot_label"].to_dict()


def setup_style() -> None:
    """Set manuscript fonts and editable vector-export settings."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 9.2,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 8.8,
            "ytick.labelsize": 8.8,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def fit_ols(data: pd.DataFrame, response_col: str, predictors: list[str], model_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit OLS with an intercept and calculate coefficients, predictions and fit statistics."""
    use_cols = [response_col, *predictors]
    fit_data = data.dropna(subset=use_cols).copy()
    y = fit_data[response_col].to_numpy(dtype=float)
    x = fit_data[predictors].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(fit_data)), x])
    terms = ["Intercept", *predictors]

    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    yhat = design @ beta
    residual = y - yhat
    n = len(y)
    p = design.shape[1]
    df_resid = n - p
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - sse / sst
    adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / df_resid
    sigma2 = sse / df_resid
    cov_beta = sigma2 * np.linalg.inv(design.T @ design)
    se = np.sqrt(np.diag(cov_beta))
    t_values = beta / se
    p_values = 2.0 * stats.t.sf(np.abs(t_values), df_resid)
    df_model = p - 1
    if df_model > 0 and df_resid > 0 and sse > 0 and sst > sse:
        ms_model = ((sst - sse) / df_model)
        ms_error = sse / df_resid
        f_value = ms_model / ms_error
        overall_p = float(stats.f.sf(f_value, df_model, df_resid))
    else:
        f_value = np.nan
        overall_p = np.nan

    x_std = x.std(axis=0, ddof=1)
    y_std = y.std(ddof=1)
    # Convert fitted slopes to changes per one standard deviation of each predictor.
    standardized = [np.nan, *(beta[1:] * x_std / y_std)]

    result_rows = []
    for term, coef, se_value, t_value, p_value, std_coef in zip(terms, beta, se, t_values, p_values, standardized):
        result_rows.append(
            {
                "response": response_col,
                "model": model_name,
                "term": term,
                "coef": coef,
                "std_error": se_value,
                "t_value": t_value,
                "p_value": p_value,
                "standardized_coef": std_coef,
                "n": n,
                "df_resid": df_resid,
                "df_model": df_model,
                "R2": r2,
                "adj_R2": adj_r2,
                "RMSE": float(np.sqrt(sse / n)),
                "F_value": f_value,
                "overall_p": overall_p,
            }
        )

    predictions = fit_data[["gauge_id", "future_model", *predictors, response_col]].copy()
    predictions["response"] = response_col
    predictions["model"] = model_name
    predictions["observed"] = y
    predictions["predicted"] = yhat
    predictions["residual"] = residual
    return pd.DataFrame.from_records(result_rows), predictions


def predict_with_ols(train_data: pd.DataFrame, test_data: pd.DataFrame, response_col: str, predictors: list[str]) -> np.ndarray:
    """Fit on the training catchments and predict the withheld catchment."""
    train = train_data.dropna(subset=[response_col, *predictors]).copy()
    test = test_data.dropna(subset=[response_col, *predictors]).copy()
    if train.empty or test.empty:
        return np.array([], dtype=float)
    x_train = train[predictors].to_numpy(dtype=float)
    y_train = train[response_col].to_numpy(dtype=float)
    design_train = np.column_stack([np.ones(len(train)), x_train])
    beta = np.linalg.lstsq(design_train, y_train, rcond=None)[0]
    design_test = np.column_stack([np.ones(len(test)), test[predictors].to_numpy(dtype=float)])
    return design_test @ beta


def ols_model_selection_row(data: pd.DataFrame, response_col: str, predictors: list[str], model_id: str) -> dict[str, float | int | str | bool]:
    """Calculate Gaussian-likelihood AICc for one candidate predictor set."""
    use_cols = [response_col, *predictors]
    fit_data = data.dropna(subset=use_cols).copy()
    y = fit_data[response_col].to_numpy(dtype=float)
    x = fit_data[predictors].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(fit_data)), x])

    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    yhat = design @ beta
    residual = y - yhat
    n = len(y)
    p_design = design.shape[1]
    df_resid = n - p_design
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / df_resid if df_resid > 0 and np.isfinite(r2) else np.nan
    rmse = float(np.sqrt(sse / n))
    df_model = p_design - 1
    if df_model > 0 and df_resid > 0 and sse > 0 and sst > sse:
        f_value = ((sst - sse) / df_model) / (sse / df_resid)
        overall_p = float(stats.f.sf(f_value, df_model, df_resid))
    else:
        f_value = np.nan
        overall_p = np.nan

    # Gaussian OLS likelihood; parameter count includes regression parameters
    # plus residual variance.
    k_likelihood = p_design + 1
    sigma2_mle = max(sse / n, np.finfo(float).tiny)
    log_likelihood = float(-0.5 * n * (np.log(2.0 * np.pi) + 1.0 + np.log(sigma2_mle)))
    aic = float(-2.0 * log_likelihood + 2.0 * k_likelihood)
    bic = float(-2.0 * log_likelihood + np.log(n) * k_likelihood)
    if n - k_likelihood - 1 > 0:
        aicc = float(aic + (2.0 * k_likelihood * (k_likelihood + 1.0)) / (n - k_likelihood - 1.0))
    else:
        aicc = np.nan

    optional = [p for p in predictors if p not in PREDICTORS]
    return {
        "response": response_col,
        "model_id": model_id,
        "predictors": " + ".join(predictors),
        "topographic_predictors": "none" if not optional else " + ".join(optional),
        "n": n,
        "n_predictors": len(predictors),
        "n_topographic_predictors": len(optional),
        "rank_deficient": bool(np.linalg.matrix_rank(design) < p_design),
        "df_resid": df_resid,
        "df_model": df_model,
        "R2": r2,
        "adj_R2": adj_r2,
        "RMSE": rmse,
        "F_value": f_value,
        "overall_p": overall_p,
        "log_likelihood": log_likelihood,
        "AIC": aic,
        "AICc": aicc,
        "BIC": bic,
    }


def site_colors(gauge_ids: list[str]) -> dict[str, tuple[float, float, float, float]]:
    """Assign the same catchment colours used in the manuscript."""
    cmap = plt.get_cmap("tab20")
    colors = {gauge_id: cmap(ii % 20) for ii, gauge_id in enumerate(site_order_by_label(gauge_ids))}
    return deepen_light_site_colors(colors)


def deepen_light_site_colors(
    colors: dict[str, tuple[float, float, float, float]],
    luminance_threshold: float = 0.55,
    darken_factor: float = 0.78,
) -> dict[str, tuple[float, float, float, float]]:
    """Darken pale colours to keep points visible on a white background."""
    adjusted: dict[str, tuple[float, float, float, float]] = {}
    for gauge_id, color in colors.items():
        rgba = mpl.colors.to_rgba(color)
        rgb = np.asarray(rgba[:3], dtype=float)
        luminance = float(0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2])
        if luminance >= luminance_threshold:
            rgb = np.clip(rgb * darken_factor, 0.0, 1.0)
        adjusted[gauge_id] = (float(rgb[0]), float(rgb[1]), float(rgb[2]), float(rgba[3]))
    return adjusted


def site_label_sort_key(label: str) -> tuple[int, int, str]:
    """Sort US labels first, then Swiss labels, in numeric order."""
    match = re.fullmatch(r"(US|CH)(\d+)", str(label))
    if match:
        region_rank = 0 if match.group(1) == "US" else 1
        return region_rank, int(match.group(2)), str(label)
    return 2, 999, str(label)


def site_order_by_label(gauge_ids: list[str]) -> list[str]:
    """Return catchments in the shared legend order."""
    labels = site_display_labels(gauge_ids)
    return sorted(gauge_ids, key=lambda gauge_id: site_label_sort_key(labels[gauge_id]))


def format_panel_p_value(p_value: float, broad_threshold: bool = False) -> str:
    """Format the significance thresholds displayed in the figure."""
    if not np.isfinite(p_value):
        return r"$P$ = n.a."
    if broad_threshold and p_value < 0.05:
        return r"$P$ < 0.05"
    if p_value < 0.01:
        return r"$P$ < 0.01"
    if p_value < 0.05:
        return r"$P$ < 0.05"
    return r"$P$ > 0.05"


def univariate_slope_p_values(data: pd.DataFrame, response_col: str, predictors: list[str]) -> dict[str, float]:
    """Test each predictor separately against a DOC response."""
    p_values: dict[str, float] = {}
    for predictor in predictors:
        fit_data = data[[response_col, predictor]].dropna()
        if len(fit_data) < 3 or fit_data[predictor].std(ddof=1) == 0 or fit_data[response_col].std(ddof=1) == 0:
            p_values[predictor] = np.nan
            continue
        result = stats.linregress(
            fit_data[predictor].to_numpy(dtype=float),
            fit_data[response_col].to_numpy(dtype=float),
        )
        p_values[predictor] = float(result.pvalue)
    return p_values


@mpl.rc_context({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7.7,
    "ytick.labelsize": 7.7,
    "legend.fontsize": 7.6,
    "axes.linewidth": 0.72,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.68,
    "ytick.major.width": 0.68,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
def plot_model_comparison_figure(
    predictions: pd.DataFrame,
    regression: pd.DataFrame,
    loso_summary: pd.DataFrame,
    relative_flux_climate_predictions: pd.DataFrame,
    relative_flux_climate_regression: pd.DataFrame,
    relative_flux_predictions: pd.DataFrame,
    relative_flux_regression: pd.DataFrame,
    relative_flux_loso_summary: pd.DataFrame,
    doc_flux_univariate_p: dict[str, float],
    relative_flux_univariate_p: dict[str, float],
    n_sites: int,
) -> None:
    """Draw the four prediction panels and standardized-coefficient comparison."""
    predictor_order = COMPARISON_PREDICTOR_ORDER
    response_colors = {
        "delta_doc_flux_mgC_m2_h": "#C94E3C",
        RELATIVE_DOC_FLUX_COL: "#7A6BAE",
    }
    response_labels = {
        "delta_doc_flux_mgC_m2_h": r"absolute change in $\boldsymbol{F}_{\mathbf{DOC}}$",
        RELATIVE_DOC_FLUX_COL: r"relative change in $\boldsymbol{F}_{\mathbf{DOC}}$",
    }
    model_labels = {
        "climate_only": "Climate-only",
        "relative_flux_climate_only": "Climate-only",
        "selected_model": "Climate-plus-landscape",
        "relative_flux_selected_model": "Climate-plus-landscape",
    }

    all_predictions = pd.concat(
        [predictions, relative_flux_climate_predictions, relative_flux_predictions],
        ignore_index=True,
    )
    all_regression = pd.concat(
        [regression, relative_flux_climate_regression, relative_flux_regression],
        ignore_index=True,
    )
    all_loso_summary = pd.concat([loso_summary, relative_flux_loso_summary], ignore_index=True)
    site_order = site_order_by_label(list(pd.unique(all_predictions["gauge_id"].astype(str))))
    colors = site_colors(site_order)
    display_labels = site_display_labels(site_order)

    scatter_specs = [
        {
            "panel": "a",
            "model": "climate_only",
            "response": "delta_doc_flux_mgC_m2_h",
            "unit": r"mg C m$^{-2}$ h$^{-1}$",
            "title": response_labels["delta_doc_flux_mgC_m2_h"],
            "data": predictions,
        },
        {
            "panel": "b",
            "model": "relative_flux_climate_only",
            "response": RELATIVE_DOC_FLUX_COL,
            "unit": "%",
            "title": response_labels[RELATIVE_DOC_FLUX_COL],
            "data": relative_flux_climate_predictions,
        },
        {
            "panel": "c",
            "model": "selected_model",
            "response": "delta_doc_flux_mgC_m2_h",
            "unit": r"mg C m$^{-2}$ h$^{-1}$",
            "title": response_labels["delta_doc_flux_mgC_m2_h"],
            "data": predictions,
        },
        {
            "panel": "d",
            "model": "relative_flux_selected_model",
            "response": RELATIVE_DOC_FLUX_COL,
            "unit": "%",
            "title": response_labels[RELATIVE_DOC_FLUX_COL],
            "data": relative_flux_predictions,
        },
    ]

    response_lims: dict[str, tuple[float, float]] = {}
    for response_col in ["delta_doc_flux_mgC_m2_h", RELATIVE_DOC_FLUX_COL]:
        subset = all_predictions[all_predictions["response"].eq(response_col)]
        lim_min = float(np.nanmin([subset["observed"].min(), subset["predicted"].min()]))
        lim_max = float(np.nanmax([subset["observed"].max(), subset["predicted"].max()]))
        pad = (lim_max - lim_min) * 0.10 if lim_max > lim_min else 1.0
        response_lims[response_col] = (lim_min - 2.2 * pad, lim_max + 2.3 * pad)

    # Keep panel geometry unchanged, with headroom for the enlarged panel letters.
    width_mm, height_mm = 230.0, 168.0
    panel_scale = width_mm / 183.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))

    def axes_mm(left: float, bottom: float, width: float, height: float):
        return fig.add_axes([left / width_mm, bottom / height_mm,
                             width / width_mm, height / height_mm])

    scatter_axes = {
        "a": axes_mm(19, 105, 52, 52),
        "b": axes_mm(86, 105, 52, 52),
        "c": axes_mm(19, 32, 52, 52),
        "d": axes_mm(86, 32, 52, 52),
    }
    ax_coef = axes_mm(163, 32, 62, 125)

    for spec in scatter_specs:
        model_name = spec["model"]
        response_col = spec["response"]
        ax = scatter_axes[spec["panel"]]
        subset = spec["data"][(spec["data"]["model"] == model_name) & (spec["data"]["response"] == response_col)]
        for gauge_id, site in subset.groupby("gauge_id", sort=True):
            ax.scatter(
                site["observed"],
                site["predicted"],
                s=12,
                color=colors[gauge_id],
                edgecolor=SITE_SCATTER_EDGE_COLOR,
                linewidth=SITE_SCATTER_EDGE_WIDTH,
                alpha=SITE_SCATTER_ALPHA,
                label=gauge_id,
                zorder=3,
            )
        lims = response_lims[response_col]
        ax.plot(lims, lims, color="#333333", lw=0.75, ls=(0, (4, 3)), zorder=2)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#DADADA", linewidth=0.42, alpha=0.70)
        ax.set_xlabel(f"Process-model simulated ({spec['unit']})", labelpad=3)
        ax.set_ylabel(f"Regression predicted ({spec['unit']})", labelpad=3)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4))
        if spec["title"]:
            ax.set_title(f"{model_labels[model_name]}\n{spec['title']}",
                         pad=4, fontsize=8.5,
                         fontweight="bold", linespacing=1.05)

        r2 = all_regression[
            (all_regression["response"] == response_col)
            & (all_regression["model"] == model_name)
            & (all_regression["term"] == "Intercept")
        ]["R2"].iloc[0]
        overall_p = all_regression[
            (all_regression["response"] == response_col)
            & (all_regression["model"] == model_name)
            & (all_regression["term"] == "Intercept")
        ]["overall_p"].iloc[0]
        loso_r2 = all_loso_summary[
            (all_loso_summary["response"] == response_col)
            & (all_loso_summary["model"] == model_name)
        ]["LOSO_R2"].iloc[0]
        ax.text(
            0.04,
            0.98,
            rf"$R^2$ = {r2:.2f}, {format_panel_p_value(float(overall_p))}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.7,
        )
        ax.text(
            0.04,
            0.88,
            rf"LOSO $R^2$ = {float(loso_r2):.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.7,
        )
        ax.text(
            0.0,
            1.0,
            spec["panel"],
            transform=ax.title.get_transform() + mpl.transforms.ScaledTranslation(
                0, 9.4 / 72, fig.dpi_scale_trans),
            fontsize=10.0 * panel_scale,
            fontweight="bold",
            ha="left",
            va=ax.title.get_verticalalignment(),
            clip_on=False,
        )
    ax_coef.text(0.04, 0.992, "e", transform=ax_coef.transAxes,
                 fontsize=10.0 * panel_scale, fontweight="bold", va="top", zorder=5)

    site_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=colors[gauge_id],
            markeredgecolor=SITE_SCATTER_EDGE_COLOR,
            markeredgewidth=SITE_SCATTER_EDGE_WIDTH,
            markersize=3.5,
            label=gauge_id,
        )
        for gauge_id in site_order
    ]
    for handle, gauge_id in zip(site_handles, site_order):
        handle.set_label(display_labels[gauge_id])
    y_base = np.arange(len(predictor_order))
    bar_height = 0.16
    series = [
        ("climate_only", "delta_doc_flux_mgC_m2_h", -0.24, "Absolute DOC export change, Climate-only", "///"),
        ("relative_flux_climate_only", RELATIVE_DOC_FLUX_COL, -0.08, "Relative DOC export change, Climate-only", "///"),
        ("selected_model", "delta_doc_flux_mgC_m2_h", 0.08, "Absolute DOC export change, Climate-plus-landscape", None),
        ("relative_flux_selected_model", RELATIVE_DOC_FLUX_COL, 0.24, "Relative DOC export change, Climate-plus-landscape", None),
    ]
    for model_name, response_col, offset, label, hatch in series:
        coeffs = (
            all_regression[(all_regression["model"] == model_name) & (all_regression["response"] == response_col)]
            .set_index("term")["standardized_coef"]
            .reindex(predictor_order)
        )
        coeffs = pd.to_numeric(coeffs, errors="coerce")
        valid = coeffs.notna()
        color = response_colors[response_col]
        ax_coef.barh(
            y_base[valid] + offset,
            coeffs[valid],
            height=bar_height,
            facecolor="white" if hatch else color,
            edgecolor=color,
            linewidth=0.95,
            hatch=hatch,
            label=label,
            alpha=1.0 if hatch else 0.62,
            zorder=3,
        )
        if model_name in {"selected_model", "relative_flux_selected_model"}:
            p_values = doc_flux_univariate_p if model_name == "selected_model" else relative_flux_univariate_p
            x_span = max(float(np.nanmax(np.abs(coeffs[valid]))), 0.15)
            for yy, predictor, coef in zip(y_base[valid] + offset, coeffs[valid].index, coeffs[valid].to_numpy(dtype=float)):
                if predictor not in p_values or not np.isfinite(coef):
                    continue
                p_label = format_panel_p_value(p_values[predictor])
                x_pad = 0.035 * x_span
                text_x = coef + x_pad if coef >= 0 else x_pad
                ax_coef.text(
                    text_x,
                    yy,
                    p_label,
                    color=color,
                    fontsize=7.7,
                    fontweight="bold",
                    ha="left",
                    va="center",
                    zorder=4,
                )

    ax_coef.axvline(0, color="#555555", lw=0.80, zorder=2)
    ax_coef.set_yticks(y_base)
    ax_coef.set_yticklabels([PREDICTOR_LABELS[p].replace("distance to channel", "distance\nto channel")
                            for p in predictor_order])
    ax_coef.invert_yaxis()
    ax_coef.set_ylim(len(predictor_order) - 0.45, -0.65)
    coef_min, coef_max = ax_coef.dataLim.intervalx
    ax_coef.set_xlim(min(-0.52, coef_min - 0.08), max(1.10, coef_max + 0.30))
    ax_coef.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax_coef.set_xlabel("Standardized regression coefficient", labelpad=3)
    ax_coef.grid(axis="x", color="#DADADA", linewidth=0.42, alpha=0.72, zorder=1)
    coef_legend = ax_coef.legend(
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(225 / width_mm, 8 / height_mm),
        bbox_transform=fig.transFigure,
        ncol=2,
        fontsize=7.6,
        handlelength=1.3,
        borderaxespad=0,
        labelspacing=0.50,
        columnspacing=0.70,
    )

    # Fixed physical baselines prevent superscripted units from shifting labels.
    for panel, ax in scatter_axes.items():
        box = ax.get_position()
        baseline_mm = 95 if panel in {"a", "b"} else 22
        ax.xaxis.label.set_verticalalignment("baseline")
        ax.xaxis.set_label_coords(box.x0 + box.width / 2, baseline_mm / height_mm,
                                 transform=fig.transFigure)
        label_x_mm = 7 if panel in {"a", "c"} else 78
        ax.yaxis.set_label_coords(label_x_mm / width_mm, box.y0 + box.height / 2,
                                 transform=fig.transFigure)
    coef_box = ax_coef.get_position()
    ax_coef.xaxis.label.set_verticalalignment("baseline")
    ax_coef.xaxis.set_label_coords(coef_box.x0 + coef_box.width / 2, 22 / height_mm,
                                  transform=fig.transFigure)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Both legend groups have two single-line rows at identical vertical positions.
    for row in range(2):
        model_box = coef_legend.get_texts()[row].get_window_extent(renderer)
        row_y = (model_box.y0 + model_box.y1) / (2 * fig.bbox.height)
        fig.legend(
            handles=site_handles[row * 6:(row + 1) * 6],
            loc="center",
            bbox_to_anchor=(50.5 / width_mm, row_y),
            ncol=6,
            frameon=False,
            handlelength=0.8,
            handletextpad=0.25,
            columnspacing=2.20,
            borderaxespad=0,
            borderpad=0,
            fontsize=7.6,
        )

    out_stem = OUTPUT_DIR / "fig4"
    fig.savefig(out_stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(out_stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(out_stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(out_stem.with_suffix(".tif"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def select_landscape_predictors(data: pd.DataFrame, response: str) -> list[str]:
    """Compare all eight landscape subsets while always retaining both climate terms."""
    candidates = []
    for size in range(len(TOPO_PREDICTORS) + 1):
        for subset in combinations(TOPO_PREDICTORS, size):
            predictors = [*PREDICTORS, *subset]
            candidates.append(ols_model_selection_row(data, response, predictors,
                                                       f"m{len(candidates):02d}"))
    ranking = pd.DataFrame(candidates).sort_values(["AICc", "BIC", "n_predictors"])
    return ranking.iloc[0]["predictors"].split(" + ")


def leave_one_site_out(data: pd.DataFrame, response: str,
                       model_specs: list[tuple[str, list[str]]]) -> pd.DataFrame:
    """Refit coefficients after withholding all four scenarios of one catchment."""
    # Predictor sets are selected on the full dataset and held fixed across folds,
    # matching the validation shown in the manuscript figure.
    rows = []
    for held_out_site in sorted(data["gauge_id"].unique()):
        train = data.loc[data["gauge_id"] != held_out_site]
        test = data.loc[data["gauge_id"] == held_out_site]
        for model, predictors in model_specs:
            predicted = predict_with_ols(train, test, response, predictors)
            for observed, fitted in zip(test[response], predicted):
                rows.append({"model": model, "observed": observed, "predicted": fitted})
    summaries = []
    for model, group in pd.DataFrame(rows).groupby("model", sort=False):
        observed = group["observed"].to_numpy(dtype=float)
        predicted = group["predicted"].to_numpy(dtype=float)
        sse = float(np.sum((observed - predicted) ** 2))
        sst = float(np.sum((observed - observed.mean()) ** 2))
        summaries.append({"response": response, "model": model, "LOSO_R2": 1.0 - sse / sst})
    return pd.DataFrame(summaries)


def analyse(data: pd.DataFrame) -> dict:
    """Fit both response scales and return the in-memory tables needed by the plot."""
    results = {}
    for response, prefix in ((ABSOLUTE_DOC_FLUX_COL, ""), (RELATIVE_DOC_FLUX_COL, "relative_flux_")):
        selected = select_landscape_predictors(data, response)
        climate = fit_ols(data, response, PREDICTORS, prefix + "climate_only")
        landscape = fit_ols(data, response, selected, prefix + "selected_model")
        loso = leave_one_site_out(data, response,
                                 [(prefix + "climate_only", PREDICTORS),
                                  (prefix + "selected_model", selected)])
        # Panel e annotates separate univariate slope tests, not partial-coefficient tests.
        univariate_p = univariate_slope_p_values(data, response, COMPARISON_PREDICTOR_ORDER)
        results[response] = {"climate": climate, "landscape": landscape, "loso": loso,
                             "univariate_p": univariate_p, "selected": selected}
    return results


def main() -> None:
    """Run the analysis, report selected models and save Figure 4 in four formats."""
    setup_style()
    data = load_analysis_data()
    results = analyse(data)
    absolute = results[ABSOLUTE_DOC_FLUX_COL]
    relative = results[RELATIVE_DOC_FLUX_COL]
    regression = pd.concat([absolute["climate"][0], absolute["landscape"][0]], ignore_index=True)
    predictions = pd.concat([absolute["climate"][1], absolute["landscape"][1]], ignore_index=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_model_comparison_figure(
        predictions, regression, absolute["loso"],
        relative["climate"][1], relative["climate"][0],
        relative["landscape"][1], relative["landscape"][0], relative["loso"],
        absolute["univariate_p"], relative["univariate_p"], data["gauge_id"].nunique(),
    )
    for response, result in results.items():
        print(f"{response}: selected predictors = {', '.join(result['selected'])}")
        for _, row in result["loso"].iterrows():
            print(f"  {row['model']}: LOSO R-squared = {row['LOSO_R2']:.4f}")
    print(f"Saved Figure 4 to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
