"""Recalculate the bootstrap and pairing tests for Extended Data Figure 2.

Read cleaned catchment-model summaries, select the landscape regressions and
run 10,000 replicates per response and test. Intermediate results stay in memory;
only PNG, PDF, SVG and TIFF figures are written to the output folder.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# All paths are relative to this script; no other figure package is required.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
N_BOOTSTRAP = 10_000
RANDOM_SEED = 20260623
CLIMATE_PREDICTORS = ["delta_t_degC", "delta_pr_mm_yr"]
LANDSCAPE_PREDICTORS = ["mean_elevation_m", "mean_distance_to_channel_m", "slope_p90_p10_deg"]
ABSOLUTE_FLUX_RESPONSE = "delta_doc_flux_mgC_m2_h"
RELATIVE_FLUX_RESPONSE = "relative_delta_doc_flux_percent"

PREDICTOR_ORDER = [
    "delta_pr_mm_yr",
    "delta_t_degC",
    "mean_elevation_m",
    "mean_distance_to_channel_m",
    "slope_p90_p10_deg",
]

PREDICTOR_LABELS = {
    "delta_pr_mm_yr": r"$\Delta$P",
    "delta_t_degC": r"$\Delta$T",
    "mean_elevation_m": "Mean elevation",
    "mean_distance_to_channel_m": "Mean distance\nto channel",
    "slope_p90_p10_deg": "Slope range",
}

MODEL_ORDER = [
    "absolute_doc_flux_topography",
    "relative_doc_flux_topography",
]

MODEL_LABELS = {
    "absolute_doc_flux_topography": "Absolute DOC export change",
    "relative_doc_flux_topography": "Relative DOC export change",
}

MODEL_TICK_LABELS = {
    "absolute_doc_flux_topography": "Change in\nDOC export",
    "relative_doc_flux_topography": "Relative change in\nDOC export",
}

MODEL_COLORS = {
    "absolute_doc_flux_topography": "#C94E3C",
    "relative_doc_flux_topography": "#7A6BAE",
}

MODEL_OFFSETS = {
    "absolute_doc_flux_topography": -0.17,
    "relative_doc_flux_topography": 0.17,
}


def load_data() -> pd.DataFrame:
    """Read all 48 records and calculate relative DOC-export change."""
    data = pd.read_csv(DATA_DIR / "regression_inputs.csv", dtype={"gauge_id": str})
    numeric = [*CLIMATE_PREDICTORS, *LANDSCAPE_PREDICTORS,
               "current_doc_flux_mgC_m2_h", ABSOLUTE_FLUX_RESPONSE]
    if not np.isfinite(data[numeric].to_numpy(float)).all():
        raise ValueError("Regression inputs must be finite.")
    if (data["current_doc_flux_mgC_m2_h"] <= 0).any():
        raise ValueError("Relative change requires positive current DOC fluxes.")
    if len(data) != 48 or data["gauge_id"].nunique() != 12:
        raise ValueError("Expected 48 records from 12 catchments.")
    if data.duplicated(["gauge_id", "future_model"]).any() or not data.groupby("gauge_id").size().eq(4).all():
        raise ValueError("Each catchment must have four distinct climate projections.")
    data[RELATIVE_FLUX_RESPONSE] = data[ABSOLUTE_FLUX_RESPONSE] / data["current_doc_flux_mgC_m2_h"] * 100.0
    return data


def select_predictors(data: pd.DataFrame, response: str) -> list[str]:
    """Select the minimum-AICc subset, keeping both climate predictors in every model."""
    candidates = []
    for size in range(len(LANDSCAPE_PREDICTORS) + 1):
        for optional in combinations(LANDSCAPE_PREDICTORS, size):
            predictors = [*CLIMATE_PREDICTORS, *optional]
            y = data[response].to_numpy(float)
            design = np.column_stack([np.ones(len(data)), data[predictors].to_numpy(float)])
            beta = np.linalg.lstsq(design, y, rcond=None)[0]
            sse = float(np.sum((y - design @ beta) ** 2))
            n = len(y)
            k = design.shape[1] + 1  # Intercept, slopes and residual variance.
            variance = max(sse / n, np.finfo(float).tiny)
            log_likelihood = float(-0.5 * n * (np.log(2.0 * np.pi) + 1.0 + np.log(variance)))
            aic = float(-2.0 * log_likelihood + 2.0 * k)
            aicc = float(aic + 2.0 * k * (k + 1) / (n - k - 1))
            bic = float(-2.0 * log_likelihood + np.log(n) * k)
            candidates.append((aicc, bic, len(predictors), predictors))
    return min(candidates, key=lambda item: item[:3])[3]


def model_specs(data: pd.DataFrame) -> list[tuple[str, str, list[str]]]:
    """Select each response model once; retain its predictors in all resamples."""
    return [(name, response, select_predictors(data, response)) for name, response in
            zip(MODEL_ORDER, [ABSOLUTE_FLUX_RESPONSE, RELATIVE_FLUX_RESPONSE])]


def fit_coefficients(data: pd.DataFrame, response: str, predictors: list[str]) -> np.ndarray:
    """Fit OLS with an intercept, rejecting rank-deficient designs."""
    y = data[response].to_numpy(dtype=float)
    x = data[predictors].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(data)), x])
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise np.linalg.LinAlgError("Rank-deficient regression design")
    return np.linalg.lstsq(design, y, rcond=None)[0]


def predict(data: pd.DataFrame, beta: np.ndarray, predictors: list[str]) -> np.ndarray:
    """Apply the fitted regression to the supplied predictor values."""
    x = data[predictors].to_numpy(dtype=float)
    return np.column_stack([np.ones(len(data)), x]) @ beta


def standardize(data: pd.DataFrame, response: str, predictors: list[str]) -> pd.DataFrame:
    """Standardize response and predictors within each resample using sample SDs."""
    columns = [response, *predictors]
    result = pd.DataFrame(index=data.index)
    for column in columns:
        values = data[column].to_numpy(dtype=float)
        sd = np.nanstd(values, ddof=1)
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"Cannot standardize {column}: zero or invalid standard deviation")
        result[column] = (values - np.nanmean(values)) / sd
    return result


def in_sample_r2(data: pd.DataFrame, response: str, predictors: list[str]) -> float:
    """Calculate ordinary fitted R-squared for a fixed predictor set."""
    beta = fit_coefficients(data, response, predictors)
    observed = data[response].to_numpy(float)
    predicted = predict(data, beta, predictors)
    sse = float(np.sum((observed - predicted) ** 2))
    sst = float(np.sum((observed - observed.mean()) ** 2))
    return float(1.0 - sse / sst)


def site_block_bootstrap(data: pd.DataFrame, specs: list) -> pd.DataFrame:
    """Resample whole catchments, retaining all four climate scenarios per draw."""
    records = []
    for spec_index, (model_name, response, predictors) in enumerate(specs):
        subset = data.dropna(subset=[response, *predictors]).copy()
        site_groups = {
            site: group.copy()
            for site, group in subset.groupby("gauge_id", sort=True)
        }
        model_sites = np.array(sorted(site_groups))
        rng = np.random.default_rng(RANDOM_SEED + spec_index)
        accepted = 0
        attempts = 0
        max_attempts = N_BOOTSTRAP * 20
        while accepted < N_BOOTSTRAP:
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(f"Too many rank-deficient site-bootstrap samples for {model_name}")
            sampled_sites = rng.choice(model_sites, size=len(model_sites), replace=True)
            sampled = pd.concat(
                [site_groups[site] for site in sampled_sites],
                ignore_index=True,
            )
            try:
                # Keep the raw-design rank check to reproduce the accepted samples.
                fit_coefficients(sampled, response, predictors)
                standardized_beta = fit_coefficients(standardize(sampled, response, predictors), response, predictors)
            except (np.linalg.LinAlgError, ValueError):
                continue

            for term, coefficient in zip(["Intercept", *predictors], standardized_beta):
                records.append(
                    {
                        "model": model_name,
                        "response": response,
                        "predictors": " + ".join(predictors),
                        "replicate": accepted + 1,
                        "term": term,
                        "coefficient_type": "standardized",
                        "coefficient": coefficient,
                    }
                )
            accepted += 1

        print(f"Bootstrap complete: {model_name} ({accepted:,} accepted; {attempts:,} draws)", flush=True)
    return pd.DataFrame.from_records(records)


def site_level_permutation(data: pd.DataFrame, specs: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Shuffle entire landscape vectors among catchments, keeping climate/response rows fixed."""
    sites = np.array(sorted(pd.unique(data["gauge_id"])))
    records, summary_rows = [], []
    for spec_index, (model_name, response, predictors) in enumerate(specs):
        subset = data.copy()
        topo_predictors = [p for p in predictors if p not in CLIMATE_PREDICTORS]
        topo = subset.groupby("gauge_id", sort=True)[topo_predictors].first().loc[sites]
        rng = np.random.default_rng(RANDOM_SEED + spec_index)
        climate_r2 = in_sample_r2(subset, response, CLIMATE_PREDICTORS)
        full_r2 = in_sample_r2(subset, response, predictors)
        for replicate in range(1, N_BOOTSTRAP + 1):
            # Move the landscape predictors together to preserve their covariance.
            shuffled = topo.iloc[rng.permutation(len(sites))].copy().reset_index(drop=True)
            shuffled.insert(0, "gauge_id", sites)
            permuted = subset.drop(columns=topo_predictors).merge(
                shuffled, on="gauge_id", how="left", validate="many_to_one")
            records.append({"model": model_name, "replicate": replicate,
                            "climate_r2": climate_r2,
                            "permuted_full_r2": in_sample_r2(permuted, response, predictors)})
        summary_rows.append({"model": model_name, "test": "incremental_R2_from_topography",
                             "observed_statistic": full_r2 - climate_r2})
        print(f"Pairing test complete: {model_name} ({N_BOOTSTRAP:,} permutations)", flush=True)
    return pd.DataFrame.from_records(records), pd.DataFrame.from_records(summary_rows)


def setup_style() -> None:
    """Set manuscript fonts and editable text in vector exports."""
    plt.rcParams.update(
        {
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
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a bold panel letter alongside the title."""
    ax.text(
        0.0,
        1.0,
        label,
        transform=ax.title.get_transform(),
        ha="left",
        va=ax.title.get_verticalalignment(),
        fontsize=10.0,
        fontweight="bold",
        clip_on=False,
    )


def plot_bootstrap(ax: plt.Axes, bootstrap: pd.DataFrame) -> tuple[pd.DataFrame, list[mpatches.Patch]]:
    """Draw standardized-coefficient distributions with percentile whiskers."""
    data = bootstrap[
        bootstrap["coefficient_type"].eq("standardized")
        & bootstrap["term"].isin(PREDICTOR_ORDER)
        & bootstrap["model"].isin(MODEL_ORDER)
    ].copy()
    if data.empty:
        raise ValueError("No standardized bootstrap coefficients were found.")

    y_base = np.arange(len(PREDICTOR_ORDER))
    source_rows: list[dict[str, float | str | int]] = []

    for model_name in MODEL_ORDER:
        model_data = data[data["model"].eq(model_name)]
        color = MODEL_COLORS[model_name]
        for y, predictor in zip(y_base, PREDICTOR_ORDER):
            values = model_data.loc[model_data["term"].eq(predictor), "coefficient"].dropna().to_numpy(dtype=float)
            if values.size == 0:
                continue
            position = y + MODEL_OFFSETS[model_name]
            box = ax.boxplot(
                values,
                positions=[position],
                vert=False,
                widths=0.25,
                patch_artist=True,
                showfliers=False,
                whis=(2.5, 97.5),
                medianprops={"color": "white", "linewidth": 1.05},
                boxprops={"facecolor": color, "edgecolor": color, "linewidth": 0.8},
                whiskerprops={"color": color, "linewidth": 0.8},
                capprops={"color": color, "linewidth": 0.8},
            )
            for artist in box["boxes"]:
                artist.set_alpha(0.90)
            median = float(np.median(values))
            ax.plot(median, position, marker="o", ms=2.4, color="white", markeredgewidth=0, zorder=4)
            source_rows.append(
                {
                    "panel": "a",
                    "model": model_name,
                    "term": predictor,
                    "n": int(values.size),
                    "median": median,
                    "ci_2_5": float(np.quantile(values, 0.025)),
                    "ci_97_5": float(np.quantile(values, 0.975)),
                }
            )

    ax.axvline(0, color="#5A5A5A", lw=0.75, ls=(0, (3, 3)), zorder=1)
    ax.set_yticks(y_base)
    ax.set_yticklabels([PREDICTOR_LABELS[predictor] for predictor in PREDICTOR_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Standardized regression coefficient")
    ax.set_title("Site-block bootstrap coefficients", pad=5)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.42, zorder=0)
    handles = [
        mpatches.Patch(facecolor=MODEL_COLORS[model], edgecolor=MODEL_COLORS[model], label=MODEL_LABELS[model])
        for model in MODEL_ORDER
    ]
    add_panel_label(ax, "a")
    return pd.DataFrame(source_rows), handles


def plot_permutation(
    ax: plt.Axes, replicates: pd.DataFrame, summary: pd.DataFrame
) -> tuple[pd.DataFrame, list[mlines.Line2D]]:
    """Compare random landscape pairings with observed and climate-only fits."""
    r2_rows = summary[summary["test"].eq("incremental_R2_from_topography")].copy()
    source_rows: list[dict[str, float | str | int]] = []

    for position, model_name in enumerate(MODEL_ORDER, start=1):
        model_reps = replicates[replicates["model"].eq(model_name)].copy()
        if model_reps.empty:
            continue
        values = model_reps["permuted_full_r2"].to_numpy(dtype=float)
        climate_r2 = float(model_reps["climate_r2"].iloc[0])
        observed_increment = float(r2_rows.loc[r2_rows["model"].eq(model_name), "observed_statistic"].iloc[0])
        observed_full_r2 = climate_r2 + observed_increment

        box = ax.boxplot(
            values,
            positions=[position],
            widths=0.48,
            patch_artist=True,
            showfliers=False,
            whis=(2.5, 97.5),
            medianprops={"color": "white", "linewidth": 1.15},
            boxprops={"facecolor": "#8FA0AE", "edgecolor": "#8FA0AE", "linewidth": 0.85},
            whiskerprops={"color": "#8FA0AE", "linewidth": 0.85},
            capprops={"color": "#8FA0AE", "linewidth": 0.85},
        )
        for artist in box["boxes"]:
            artist.set_alpha(0.82)

        ax.scatter(
            position,
            observed_full_r2,
            marker="D",
            s=30,
            color="#1F1F1F",
            edgecolor="white",
            linewidth=0.45,
            zorder=5,
        )
        ax.hlines(
            climate_r2,
            position - 0.29,
            position + 0.29,
            color="#6F6F6F",
            linewidth=1.0,
            linestyle=(0, (3, 2)),
            zorder=4,
        )
        source_rows.append(
            {
                "panel": "b",
                "model": model_name,
                "term": "full_model_R2",
                "n": int(values.size),
                "climate_only_R2": climate_r2,
                "observed_paired_topography_R2": observed_full_r2,
                "permuted_median": float(np.median(values)),
                "permuted_ci_2_5": float(np.quantile(values, 0.025)),
                "permuted_ci_97_5": float(np.quantile(values, 0.975)),
            }
        )

    ax.set_xticks(range(1, len(MODEL_ORDER) + 1))
    ax.set_xticklabels([MODEL_TICK_LABELS[model] for model in MODEL_ORDER], rotation=0)
    ax.set_ylabel(r"Model $R^2$")
    ax.set_ylim(0.35, 1.0)
    ax.set_title("Site-topography pairing test", pad=5)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.42, zorder=0)
    ax.text(
        0.04,
        0.98,
        "10,000 permutations",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.7,
        color="#3A3A3A",
    )
    handles = [
        mlines.Line2D([], [], marker="s", markersize=6, linestyle="None", color="#8FA0AE", label="Random pairing"),
        mlines.Line2D([], [], marker="D", markersize=4.8, linestyle="None", color="#1F1F1F", label="Observed pairing"),
        mlines.Line2D([], [], color="#6F6F6F", linestyle=(0, (3, 2)), linewidth=1.0, label="Climate only"),
    ]
    add_panel_label(ax, "b")
    return pd.DataFrame(source_rows), handles


def plot_robustness_figure(
    bootstrap: pd.DataFrame,
    replicates: pd.DataFrame,
    permutation_summary: pd.DataFrame,
) -> tuple[plt.Figure, pd.DataFrame]:
    # Preserve panel sizes and typography while tightening the surrounding gaps.
    """Assemble the two robustness panels with their shared legends."""
    width_mm, height_mm = 219.0, 91.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax_a = fig.add_axes([32 / width_mm, 18 / height_mm, 101 / width_mm, 64 / height_mm])
    ax_b = fig.add_axes([149 / width_mm, 18 / height_mm, 65 / width_mm, 64 / height_mm])

    bootstrap_source, bootstrap_handles = plot_bootstrap(ax_a, bootstrap)
    permutation_source, permutation_handles = plot_permutation(ax_b, replicates, permutation_summary)

    fig.legend(
        handles=bootstrap_handles,
        loc="lower center",
        bbox_to_anchor=(82.5 / width_mm, 4.0 / height_mm),
        ncol=2,
        columnspacing=1.1,
        handlelength=1.4,
        handletextpad=0.35,
        borderaxespad=0,
        borderpad=0,
        frameon=False,
    )
    fig.legend(
        handles=permutation_handles,
        loc="lower center",
        bbox_to_anchor=(181.5 / width_mm, 4.0 / height_mm),
        ncol=3,
        columnspacing=0.82,
        handlelength=1.4,
        handletextpad=0.35,
        borderaxespad=0,
        borderpad=0,
        frameon=False,
    )

    source = pd.concat([bootstrap_source, permutation_source], ignore_index=True, sort=False)
    return fig, source


def main() -> None:
    """Run both robustness tests from basic data, then export the figure."""
    setup_style()
    data = load_data()
    specs = model_specs(data)
    for name, _, predictors in specs:
        print(f"{name}: {', '.join(predictors)}", flush=True)
    bootstrap = site_block_bootstrap(data, specs)
    replicates, permutation_summary = site_level_permutation(data, specs)
    fig, _ = plot_robustness_figure(bootstrap, replicates, permutation_summary)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT_DIR / "extended_data_fig2"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".tif"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Saved Extended Data Figure 2 to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
