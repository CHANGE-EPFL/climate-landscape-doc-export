"""Reproduce Extended Data Figure 1 from two local source tables.

Show catchment climate and vegetation context alongside all 31 climate-model
configurations. Calculate the covariance ellipses and median markers directly
from the model-catchment changes; no intermediate tables are written.
"""
from __future__ import annotations

from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# File locations are independent of the terminal's working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
OUT_STEM = OUTPUT_DIR / "extended_data_fig1"

FIGURE_SIZE = (183 / 25.4, 80 / 25.4)

LABEL_SIZE = 8.5

TICK_SIZE = 7.7

LEGEND_SIZE = 8.0

PANEL_SIZE = 10.0

SCENARIO_LABELS = {
    "FGOALS-g3": "Low warming",
    "NorESM2-MM": "Moderate warming",
    "ACCESS-CM2": "High warming",
    "KACE-1-0-G": "Very high warming",
}

MODEL_ORDER = ["FGOALS-g3", "NorESM2-MM", "ACCESS-CM2", "KACE-1-0-G"]

COLORS = ["#fcb86a", "#f47646", "#ba2b0d", "#3b0e61"]

MARKERS = ["o", "s", "^", "v"]

ELLIPSE_PROBABILITY = 0.80

BIOME_COLORS = {
    "Evergreen forest": "#1F6F5B",
    "Deciduous/mixed forest": "#8AA64A",
    "Grassland": "#D5A329",
}

REGION_MARKERS = {
    "CONUS": ("o", "Continental US"),
    "Alps_CH": ("^", "Swiss Alps"),
}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read all catchments and model changes, then identify the four study models."""
    sites = pd.read_csv(DATA_DIR / "catchment_climate.csv", dtype={"gauge_id": str})
    background = pd.read_csv(DATA_DIR / "climate_changes.csv", dtype={"site_id": str})
    expected_ids = set(sites["gauge_id"])
    if len(sites) != 12 or len(expected_ids) != 12:
        raise ValueError("Expected one climate record for each of the 12 catchments.")
    if not np.isfinite(sites[["ta_mean_degC", "pr_total_mm_yr"]].to_numpy(float)).all():
        raise ValueError("Catchment climate values must be finite.")
    if not sites["biome_group"].isin(BIOME_COLORS).all() or not sites["region"].isin(REGION_MARKERS).all():
        raise ValueError("Unexpected catchment vegetation group or region.")
    if background["model_member"].nunique() != 31 or len(background) != 372:
        raise ValueError("Expected 31 climate-model configurations x 12 catchments.")
    if background.duplicated(["model_member", "site_id"]).any():
        raise ValueError("Duplicate climate-model/catchment records.")
    if not np.isfinite(background[["dta_K", "dpr_percent"]].to_numpy(float)).all():
        raise ValueError("Climate changes must be finite.")
    for _, group in background.groupby("model_member", sort=False):
        if set(group["site_id"]) != expected_ids:
            raise ValueError("Incomplete catchment coverage for a climate model.")
    # These are the four projections used in the study, not a new selection run.
    selected = background.loc[background["model_member"].isin(MODEL_ORDER)].copy()
    if selected["model_member"].nunique() != 4:
        raise ValueError("One or more of the four study projections are missing.")
    return sites, selected, background


def biome_legend_handles() -> list[Line2D]:
    """Create shared legends for vegetation colours and regional marker shapes."""
    biome_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=color,
            markeredgecolor="#222222",
            markeredgewidth=0.6,
            markersize=7.0,
            label=biome,
        )
        for biome, color in BIOME_COLORS.items()
    ]
    region_handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="None",
            markerfacecolor="white",
            markeredgecolor="#222222",
            markeredgewidth=0.9,
            markersize=7.0,
            label=label,
        )
        for marker, label in REGION_MARKERS.values()
    ]
    return [*biome_handles, *region_handles]


def draw_biome_panel(ax: plt.Axes, data: pd.DataFrame) -> None:
    """Plot contemporary temperature and precipitation for all 12 catchments."""
    for region, (marker, _) in REGION_MARKERS.items():
        for biome, color in BIOME_COLORS.items():
            subset = data[(data["region"] == region) & (data["biome_group"] == biome)]
            if subset.empty:
                continue
            ax.scatter(
                subset["ta_mean_degC"],
                subset["pr_total_mm_yr"],
                s=58,
                marker=marker,
                color=color,
                edgecolor="#222222",
                linewidth=0.6,
                alpha=0.92,
                zorder=3,
            )

    ax.set_xlabel("Mean air temperature (°C)")
    ax.set_ylabel(r"Annual precipitation (mm yr$^{-1}$)")
    ax.grid(color="#E1E1E1", linewidth=0.45, alpha=0.85)
    ax.set_xlim(-2.0, 17.4)
    ax.set_ylim(600, 3150)
    ax.set_yticks([750, 1500, 2250, 3000])


def covariance_ellipse(
    x: np.ndarray,
    y: np.ndarray,
    probability: float,
    n_points: int = 360,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate an 80% bivariate-normal contour using population covariance."""
    values = np.column_stack([x, y])
    center = np.mean(values, axis=0)
    # Divide by n, as in the original figure, rather than n - 1.
    covariance = np.cov(values, rowvar=False, bias=True)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0)
    eigenvectors = eigenvectors[:, order]

    radius = np.sqrt(-2.0 * np.log(1.0 - probability))
    theta = np.linspace(0, 2 * np.pi, n_points)
    unit_circle = np.vstack([np.cos(theta), np.sin(theta)])
    ellipse = (
        center[:, None]
        + eigenvectors
        @ np.diag(radius * np.sqrt(eigenvalues))
        @ unit_circle
    )
    return ellipse[0], ellipse[1]


def plot_ellipse(
    data: pd.DataFrame,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    background: pd.DataFrame | None = None,
    *,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Draw all model-catchment changes and the highlighted model envelopes."""
    owns_figure = ax is None
    if owns_figure:
        fig, ax = plt.subplots(figsize=(7.2, 6.1), facecolor="white")
    else:
        fig = ax.figure
    fig.patch.set_alpha(1.0)
    ax.set_facecolor("white")

    if background is not None and not background.empty:
        other = background.loc[
            ~background["model_member"].astype(str).isin(MODEL_ORDER)
        ].dropna(subset=["dta_K", "dpr_percent"])
        all_background = background.dropna(subset=["dta_K", "dpr_percent"])
        if not other.empty:
            ax.scatter(
                other["dta_K"].to_numpy(dtype=float),
                other["dpr_percent"].to_numpy(dtype=float),
                s=28,
                marker="o",
                facecolor="#8f8f8f",
                edgecolor="#5f5f5f",
                linewidth=0.22,
                alpha=0.62,
                label="Other climate models",
                zorder=0.8,
            )
            if not all_background.empty:
                other_ellipse_x, other_ellipse_y = covariance_ellipse(
                    all_background["dta_K"].to_numpy(dtype=float),
                    all_background["dpr_percent"].to_numpy(dtype=float),
                    ELLIPSE_PROBABILITY,
                )
                ax.fill(
                    other_ellipse_x,
                    other_ellipse_y,
                    facecolor="#9b9b9b",
                    edgecolor="none",
                    alpha=0.055,
                    zorder=0.5,
                )
                ax.plot(
                    other_ellipse_x,
                    other_ellipse_y,
                    color="#8c8c8c",
                    linewidth=1.45,
                    linestyle=(0, (5, 4)),
                    alpha=0.78,
                    zorder=0.7,
                )

    for model, color, marker in zip(MODEL_ORDER, COLORS, MARKERS):
        group = data.loc[data["model_member"] == model].dropna(
            subset=["dta_K", "dpr_percent"]
        )
        x = group["dta_K"].to_numpy(dtype=float)
        y = group["dpr_percent"].to_numpy(dtype=float)

        ellipse_x, ellipse_y = covariance_ellipse(x, y, ELLIPSE_PROBABILITY)
        ax.fill(
            ellipse_x,
            ellipse_y,
            facecolor=color,
            edgecolor="none",
            alpha=0.13,
            zorder=1,
        )
        ax.plot(
            ellipse_x,
            ellipse_y,
            color=color,
            linewidth=2.1,
            zorder=2,
        )

        ax.scatter(
            x,
            y,
            s=105,
            marker=marker,
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            alpha=0.90,
            label=SCENARIO_LABELS[model],
            zorder=3,
        )
        # Diamonds show coordinate-wise medians, not ellipse centres.
        ax.scatter(
            np.median(x),
            np.median(y),
            s=175,
            marker="D",
            facecolor=color,
            edgecolor="#222222",
            linewidth=1.1,
            zorder=4,
        )

    ax.axhline(
        0,
        color="#666666",
        linewidth=1.0,
        linestyle=(0, (5, 4)),
        zorder=0,
    )
    ax.set_xlabel(r"$\Delta T$ (K)")
    ax.set_ylabel(r"$\Delta P$ (%)")
    ax.grid(color="#d8d8d8", linewidth=0.65, alpha=0.72)
    ax.set_axisbelow(True)
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)

    return fig


def shared_limits(
    data: pd.DataFrame,
    background: pd.DataFrame | None = None,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Choose axis limits that include the points and four highlighted ellipses."""
    x_extent = data["dta_K"].dropna().to_list()
    y_extent = data["dpr_percent"].dropna().to_list()
    if background is not None and not background.empty:
        x_extent.extend(background["dta_K"].dropna().to_list())
        y_extent.extend(background["dpr_percent"].dropna().to_list())
    for model in MODEL_ORDER:
        group = data.loc[data["model_member"] == model].dropna(
            subset=["dta_K", "dpr_percent"]
        )
        ellipse_x, ellipse_y = covariance_ellipse(
            group["dta_K"].to_numpy(dtype=float),
            group["dpr_percent"].to_numpy(dtype=float),
            ELLIPSE_PROBABILITY,
        )
        x_extent.extend(ellipse_x.tolist())
        y_extent.extend(ellipse_y.tolist())

    x_min, x_max = min(x_extent), max(x_extent)
    y_min, y_max = min(y_extent), max(y_extent)
    x_pad = 0.06 * (x_max - x_min)
    y_pad = 0.06 * (y_max - y_min)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)


def setup_style() -> None:
    """Set the manuscript fonts and editable vector-text settings."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": TICK_SIZE,
        "axes.labelsize": LABEL_SIZE,
        "xtick.labelsize": TICK_SIZE,
        "ytick.labelsize": TICK_SIZE,
        "legend.fontsize": LEGEND_SIZE,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def draw_biome(ax: plt.Axes, data: pd.DataFrame) -> None:
    """Style panel a and place the vegetation and region legends."""
    draw_biome_panel(ax, data)
    for points in ax.collections:
        points.set_sizes([30])
        points.set_linewidth(0.45)
    ax.set_xticks([0, 5, 10, 15])
    handles = biome_legend_handles()
    for handle in handles:
        handle.set_markersize(4.7)
        handle.set_markeredgewidth(0.55)
    biome_legend = ax.legend(
        handles=handles[:3], loc="lower center", bbox_to_anchor=(0.43, 1.095),
        ncol=3, frameon=False, borderaxespad=0, handletextpad=0.25,
        columnspacing=0.40, handlelength=0.75, borderpad=0,
    )
    ax.add_artist(biome_legend)
    ax.legend(
        handles=handles[3:], loc="lower center", bbox_to_anchor=(0.43, 1.025),
        ncol=2, frameon=False, borderaxespad=0, handletextpad=0.35,
        columnspacing=1.0, handlelength=0.75, borderpad=0,
    )


def draw_climate(ax: plt.Axes, selected: pd.DataFrame, background: pd.DataFrame) -> None:
    """Style panel b and distinguish the four study projections from the ensemble."""
    limits = shared_limits(selected, background)
    plot_ellipse(selected, *limits, background, ax=ax)
    # Scale physical marks from the original single-panel canvas, not the data.
    for points in ax.collections:
        points.set_sizes(points.get_sizes() * 0.30)
        points.set_linewidths(points.get_linewidths() * 0.60)
    for line in ax.lines:
        line.set_linewidth(line.get_linewidth() * 0.60)
    ax.grid(color="#d8d8d8", linewidth=0.45, alpha=0.72)
    ax.set_xticks([2, 4, 6, 8, 10])
    ax.set_yticks([-20, 0, 20, 40, 60])
    handles = [
        Line2D([], [], marker=marker, linestyle="none", markersize=4.7,
               markerfacecolor=color, markeredgecolor="white", markeredgewidth=0.5,
               label=SCENARIO_LABELS[model])
        for model, color, marker in zip(MODEL_ORDER, COLORS, MARKERS)
    ]
    handles.extend([
        Line2D([], [], marker="o", linestyle="none", markersize=3.6,
               markerfacecolor="#8f8f8f", markeredgecolor="#5f5f5f",
               markeredgewidth=0.25, alpha=0.82, label="Other climate models"),
        Line2D([], [], color="#8c8c8c", linewidth=0.9, linestyle=(0, (5, 4)),
               alpha=0.78, label="All-model ellipse"),
    ])
    # Column-wise filling keeps selected scenarios in the first two columns.
    ax.legend(
        handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.025),
        ncol=3, frameon=False, borderaxespad=0, handletextpad=0.25,
        columnspacing=0.40, handlelength=0.75, labelspacing=0.45, borderpad=0,
    )


def build_figure() -> plt.Figure:
    """Assemble the two panels at the manuscript size of 183 x 80 mm."""
    setup_style()
    biome_data, selected, background = load_data()
    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="white")
    grid = fig.add_gridspec(1, 2, left=0.105, right=0.985, bottom=0.16,
                           top=0.85, wspace=0.28)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    draw_biome(axes[0], biome_data)
    draw_climate(axes[1], selected, background)
    for label, ax in zip("ab", axes):
        ax.text(0.025, 0.975, label, transform=ax.transAxes, ha="left", va="top",
                fontsize=PANEL_SIZE, fontweight="bold", zorder=10)
        ax.xaxis.labelpad = 4
        ax.yaxis.labelpad = 4
        ax.tick_params(pad=2)
    fig.align_xlabels(axes)
    return fig


def main() -> None:
    """Calculate the display summaries and export the two-panel figure."""
    fig = build_figure()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=600, facecolor="white", transparent=False)
    fig.savefig(OUT_STEM.with_suffix(".pdf"), facecolor="white", transparent=False)
    fig.savefig(OUT_STEM.with_suffix(".svg"), facecolor="white", transparent=False)
    fig.savefig(OUT_STEM.with_suffix(".tif"), dpi=600, facecolor="white", transparent=False,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Saved Extended Data Figure 1 to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
