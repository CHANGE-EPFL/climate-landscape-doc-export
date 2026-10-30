"""Plot current and future catchment responses for Figure 2.

Read precomputed ratios, C-Q slopes and aridity indices from data/. Only
plotting summaries (medians and boxplots) are calculated here. Save the
six-panel figure in output/ as PNG, PDF, SVG and TIFF.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Resolve paths from this file so the script can run from any directory.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
SITE_LABELS_FILE = DATA_DIR / "sites.csv"

# Figure fonts and export settings; SVG and PDF retain editable text.
plt.rcParams["font.family"] = "sans-serif"

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]

plt.rcParams["svg.fonttype"] = "none"

plt.rcParams["pdf.fonttype"] = 42

plt.rcParams.update(
    {
        "font.size": 10.5,
        "axes.labelsize": 11.5,
        "xtick.labelsize": 10.0,
        "ytick.labelsize": 9.5,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

# Fixed horizontal offsets separate the four climate models at each site.
FUTURE_JITTER = np.array([-0.065, -0.022, 0.022, 0.065])

SITE_SPACING = 0.72

# Panels a-e use these metrics; panel f shows the catchment aridity index.
MAIN_METRIC_ORDER = [
    "p_ratio",
    "Qlat_out",
    "Q_sim_mmh",
    "DOC_flux",
    "cq_slope_b",
]

PRIMARY_MODEL_LABELS = {
    "FGOALS-g3": "Low warming",
    "NorESM2-MM": "Moderate warming",
    "ACCESS-CM2": "High warming",
    "KACE-1-0-G": "Very high warming",
}

PRIMARY_METRIC_LABELS = {
    "p_ratio": "P ratio (-)",
    "Qlat_out": r"$Q_{\mathrm{lat}}$ ratio (-)",
    "Q_sim_mmh": r"$Q$ ratio (-)",
    "DOC_flux": r"$F_{\mathrm{DOC}}$ ratio (-)",
    "cq_slope_b": r"$C_{\mathrm{DOC}}$-Q slope, $b$",
}

PRIMARY_LABEL_SIZE = 8.5

PRIMARY_TICK_SIZE = 7.7

PRIMARY_LEGEND_SIZE = 7.6

PRIMARY_PANEL_SIZE = 10.0

ARIDITY_REGION_COLORS = {"US": "#73819C", "CH": "#9CA97D"}

MODEL_ORDER = ["FGOALS-g3", "NorESM2-MM", "ACCESS-CM2", "KACE-1-0-G"]

MODEL_COLORS = {
    "FGOALS-g3": "#fcb86a",
    "NorESM2-MM": "#f47646",
    "ACCESS-CM2": "#ba2b0d",
    "KACE-1-0-G": "#3b0e61",
}

FIXED_ABSOLUTE_Y_LIMITS = {"cq_slope_b": (-0.25, 1.0)}

def site_display_labels(site_ids: list[str], labels_file: Path = SITE_LABELS_FILE) -> list[str]:
    """Return US1-US8 and CH1-CH4 labels in the requested site order."""
    sites = pd.read_csv(labels_file, dtype={"gauge_id": str}).set_index("gauge_id")
    return sites.loc[site_ids, "plot_label"].tolist()


def site_aridity(site_ids: list[str], labels_file: Path = SITE_LABELS_FILE,
                 aridity_column: str | None = None) -> pd.Series:
    """Read the precomputed aridity indices in the same order as the panels."""
    sites = pd.read_csv(labels_file, dtype={"gauge_id": str}).set_index("gauge_id")
    return sites.loc[site_ids, "aridity_modis"].astype(float)


def padded_limits(values: np.ndarray) -> tuple[float, float]:
    """Add vertical space around the range of plotted values."""
    lower = float(np.nanmin(values))
    upper = float(np.nanmax(values))
    span = upper - lower
    pad = span * 0.16 if span > 0 else max(abs(upper) * 0.16, 1.0)
    return lower - pad, upper + pad

def plot_site_ensemble(
    axis: plt.Axes, x: float, values: pd.DataFrame, point_size: float = 56,
) -> None:
    """Draw four future values, their median and the current reference line."""
    current = float(values.loc[values["scenario"].eq("Current"), "value"].iloc[0])
    future = values[values["scenario_group"].eq("Future")].set_index("scenario")["value"]
    future_values = []
    for offset, model in zip(FUTURE_JITTER, MODEL_ORDER):
        value = float(future.loc[model])
        future_values.append(value)
        axis.scatter(
            x + offset,
            value,
            s=point_size,
            color=MODEL_COLORS[model],
            edgecolor="white",
            linewidth=0.45,
            zorder=3,
        )
    axis.hlines(current, x - 0.29, x + 0.29, color="#737373", linewidth=1.35, linestyle=(0, (3, 2)), zorder=1)
    axis.hlines(float(np.median(future_values)), x - 0.14, x + 0.14, color="#333333", linewidth=1.55, zorder=2)

def plot_all_site_boxplot(
    axis: plt.Axes,
    values: pd.DataFrame,
    future_position: float,
    point_size: float = 13,
) -> None:
    """Pool the 48 future site-model values for the All column."""
    current = values.loc[values["scenario_group"].eq("Current"), "value"].to_numpy(dtype=float)
    future = values.loc[values["scenario_group"].eq("Future"), "value"].to_numpy(dtype=float)
    axis.hlines(
        float(np.median(current)),
        future_position - 0.38,
        future_position + 0.38,
        color="#737373",
        linewidth=1.35,
        linestyle=(0, (3, 2)),
        zorder=1,
    )
    box = axis.boxplot(
        [future],
        positions=[future_position],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#333333", "linewidth": 1.2},
        boxprops={"edgecolor": "#4D4D4D", "linewidth": 0.75},
        whiskerprops={"color": "#4D4D4D", "linewidth": 0.75},
        capprops={"color": "#4D4D4D", "linewidth": 0.75},
    )
    for patch, color in zip(box["boxes"], ["#C94E3C"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.58)

    # Reproducible jitter keeps overlapping values visible without changing data.
    rng = np.random.default_rng(20260624 + len(future))
    axis.scatter(
        future_position + rng.uniform(-0.12, 0.12, size=len(future)),
        future,
        s=point_size,
        color="#C94E3C",
        alpha=0.82,
        linewidth=0,
        zorder=3,
    )

def plot_aridity_axis(
    axis: plt.Axes,
    site_ids: list[str],
    site_positions: np.ndarray,
    labels_file: Path = SITE_LABELS_FILE,
    aridity_column: str | None = None,
    region_colors: dict[str, str] | None = None,
) -> None:
    """Draw aridity bars aligned with the site columns and coloured by region."""
    aridity = site_aridity(site_ids, labels_file, aridity_column).to_numpy(dtype=float)
    colors = (
        [region_colors[label[:2]] for label in site_display_labels(site_ids, labels_file)]
        if region_colors is not None else "#9D9D9D"
    )
    axis.bar(
        site_positions,
        aridity,
        width=0.36,
        color=colors,
        edgecolor="#4D4D4D",
        linewidth=0.55,
        zorder=3,
    )
    axis.set_ylabel("Aridity\nindex (-)", labelpad=8)
    axis.set_ylim(0.0, max(1.35, float(np.nanmax(aridity)) * 1.08))
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.42, zorder=0)
    axis.tick_params(direction="out", length=3.0, pad=2.0)

def build_styled_grid(
    site_ids: list[str], by_site: dict[str, pd.DataFrame], all_data: pd.DataFrame,
    metric_order: list[str], metric_labels: dict[str, str],
    *, ylabel_x: float = 0.040,
) -> plt.Figure:
    """Arrange the metric rows, pooled summaries, aridity bars and legend."""
    # Consistent physical text sizes across all six panels.
    style = {
        "font.size": PRIMARY_TICK_SIZE,
        "axes.labelsize": PRIMARY_LABEL_SIZE,
        "xtick.labelsize": PRIMARY_TICK_SIZE,
        "ytick.labelsize": PRIMARY_TICK_SIZE,
        "legend.fontsize": PRIMARY_LEGEND_SIZE,
        "axes.linewidth": 0.7,
    }
    cq_row = metric_order.index("cq_slope_b") if "cq_slope_b" in metric_order else None
    n_metrics = len(metric_order)
    row_heights = [1.0] * n_metrics + [0.08, 0.30 if cq_row is not None else 1.0]
    with plt.rc_context(style):
        height_mm = 140 if cq_row is not None else 156.6
        bottom_mm = 9.1
        fig = plt.figure(figsize=(183 / 25.4, height_mm / 25.4), facecolor="white")
        site_positions = np.arange(len(site_ids), dtype=float) * SITE_SPACING
        left_limit = -0.55
        site_column_end = site_positions[-1] + 0.43
        all_column_width = 2.40
        all_position = site_column_end + all_column_width / 2
        right_limit = site_column_end + all_column_width
        # The right column holds pooled boxplots above the scenario legend.
        grid = fig.add_gridspec(
            n_metrics + 2, 2, left=0.115, right=0.985,
            bottom=bottom_mm / height_mm, top=1 - 3.5 / height_mm,
            width_ratios=[site_column_end - left_limit, all_column_width], wspace=0,
            height_ratios=row_heights, hspace=0.14,
        )
        metric_axes = [fig.add_subplot(grid[row, :] if key != "cq_slope_b" else grid[row, 0])
                       for row, key in enumerate(metric_order)]
        aridity_axis = fig.add_subplot(grid[-1, 0])
        if cq_row is not None:
            legend_axis = fig.add_subplot(grid[cq_row:, 1])
        else:
            legend_axis = fig.add_subplot(grid[-1, 1])
        legend_axis.set_label("Scenario legend")
        legend_axis.axis("off")
        for key, axis in zip(metric_order, metric_axes):
            values = all_data.loc[all_data["metric"].eq(key)]
            limits = FIXED_ABSOLUTE_Y_LIMITS.get(key)
            if limits is None:
                limits = padded_limits(values["value"].to_numpy(dtype=float))
            if key == "cq_slope_b":
                # Keep the lowest markers clear of the spine in the shorter row.
                limits = (min(limits[0], float(values["value"].min()) - 0.06), limits[1])
                # Shade the near-zero C-Q slope band used for chemostatic behaviour.
                axis.fill_between([-0.55, site_positions[-1] + 0.34], -0.1, 0.1,
                                  color="#D8D8D8", alpha=0.38, zorder=0)
                for y, label in ((0.225, "flushing"), (0.0, "chemostatic"), (-0.185, "dilution")):
                    axis.text(-0.48, y, label, ha="left", va="center", fontsize=6.5,
                              color="#555555")
            for x, site_id in zip(site_positions, site_ids):
                site_data = by_site[site_id]
                plot_site_ensemble(axis, x, site_data.loc[site_data["metric"].eq(key)], point_size=20)
            if key != "cq_slope_b":
                plot_all_site_boxplot(axis, values, all_position, point_size=7)
                axis.axvline(site_column_end, color="#C8C8C8", linewidth=0.65, zorder=0)
            axis.set_ylim(*limits)
            axis.set_ylabel(metric_labels[key], fontsize=PRIMARY_LABEL_SIZE)
            axis.grid(axis="y", color="#D9D9D9", linewidth=0.42, zorder=0)
            axis.set_xticks([])

        plot_aridity_axis(aridity_axis, site_ids, site_positions, region_colors=ARIDITY_REGION_COLORS)
        aridity_axis.set_ylabel("AI", fontsize=PRIMARY_LABEL_SIZE)
        aridity_axis.set_yticks([0, 2] if cq_row is not None else [0, 1, 2])
        aridity_axis.set_xticks(site_positions)
        aridity_axis.set_xticklabels(site_display_labels(site_ids))
        narrow_axes = [aridity_axis]
        if cq_row is not None:
            narrow_axes.append(metric_axes[cq_row])
        for letter, axis in zip("abcdef", [*metric_axes, aridity_axis]):
            # Matching data spans and grid widths preserve site-column alignment.
            end = site_column_end if axis in narrow_axes else right_limit
            axis.set_xlim(left_limit, end)
            axis.tick_params(direction="out", length=2.6, width=0.65, pad=2,
                             labelsize=PRIMARY_TICK_SIZE)
            bounds = axis.get_position()
            axis.yaxis.set_label_coords(ylabel_x, bounds.y0 + bounds.height / 2,
                                       transform=fig.transFigure)
            axis.yaxis.label.set_verticalalignment("center")
            axis.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(3, -3),
                          textcoords="offset points", ha="left", va="top",
                          fontsize=PRIMARY_PANEL_SIZE, fontweight="bold", zorder=10)

        handles = [
            mlines.Line2D([], [], color="#737373", linestyle=(0, (3, 2)), linewidth=1.35, label="Current"),
            mlines.Line2D([], [], color="#333333", linewidth=1.55, label="Future median"),
        ]
        handles.extend([
            mlines.Line2D([], [], marker="o", linestyle="none", markersize=4.4,
                          markerfacecolor=MODEL_COLORS[model], markeredgecolor="white",
                          markeredgewidth=0.45, label=PRIMARY_MODEL_LABELS[model])
            for model in MODEL_ORDER
        ])
        legend_axis.legend(
            handles=handles, loc="lower center", ncol=1, frameon=False,
            fontsize=PRIMARY_LEGEND_SIZE, handlelength=1.4,
            handletextpad=0.4, labelspacing=0.9 if cq_row is not None else 0.4,
            borderaxespad=0, borderpad=0,
        )
        last_box_row = max(row for row, key in enumerate(metric_order) if key != "cq_slope_b")
        metric_axes[last_box_row].annotate(
            "All", xy=(all_position, 0), xycoords=("data", "axes fraction"),
            xytext=(0, -5), textcoords="offset points", ha="center", va="top",
            fontsize=PRIMARY_TICK_SIZE,
        )
    return fig

def build_primary_grid(
    site_ids: list[str], by_site: dict[str, pd.DataFrame], all_data: pd.DataFrame,
) -> plt.Figure:
    """Build the manuscript layout with the five response metrics and AI."""
    return build_styled_grid(
        site_ids, by_site, all_data, MAIN_METRIC_ORDER, PRIMARY_METRIC_LABELS,
        ylabel_x=0.060,
    )

def main() -> None:
    """Load the two source tables and export the complete figure."""
    # Preserve leading zeros in gauge IDs and the site order listed in sites.csv.
    sites = pd.read_csv(SITE_LABELS_FILE, dtype={"gauge_id": str})
    site_ids = sites["gauge_id"].tolist()
    all_data = pd.read_csv(DATA_DIR / "ensemble_responses.csv", dtype={"gauge_id": str, "scenario": str})
    by_site = {site: all_data.loc[all_data["gauge_id"].eq(site)].copy() for site in site_ids}
    fig = build_primary_grid(site_ids, by_site, all_data)
    # Use the same physical dimensions and 600 dpi for all output formats.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf", ".svg"):
        fig.savefig(OUTPUT_DIR / ("fig2" + suffix), dpi=600, facecolor="white")
    fig.savefig(OUTPUT_DIR / "fig2.tif", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Saved Figure 2 to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
