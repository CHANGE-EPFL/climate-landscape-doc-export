"""Plot Figure 1: GPP distributions, validation metrics and DOC comparisons.

Run this script to save the six-panel figure in output/. Source tables are
read from data/; individual panel drawing functions are in plot_helpers.py.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Save figures without opening a graphical window.
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PathCollection
from matplotlib.patches import Rectangle
import geopandas as gpd

import plot_helpers as validation


# Resolve paths from this file so the script works from any working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
FIG1_DIR = SCRIPT_DIR / "output"
DATA_DIR = validation.DATA_DIR
OUT_STEM = FIG1_DIR / "fig1"

# Final double-column size; all panels share these physical font sizes.
FIGURE_SIZE = (183 / 25.4, 146 / 25.4)
LABEL_SIZE = 8.5
TICK_SIZE = 7.7
LEGEND_SIZE = 7.6
PANEL_SIZE = 10.0
LEFT_BLOCK_RIGHT = 0.620
MAP_LABEL_SIZE = 6.0
MAP_REGION_COLORS = {"US": "#73819C", "CH": "#9CA97D"}

# Offset crowded station labels in points; gauge IDs keep placements stable.
MAP_LABEL_OFFSETS = {
    "12447390": (2, 2, "left", "bottom"),
    "09352900": (2, -2, "left", "top"),
    "10343500": (-2, 0, "right", "center"),
    "07362587": (2, -2, "left", "top"),
    "01545600": (-2, -2, "right", "top"),
    "01434025": (2, -1, "left", "top"),
    "01054200": (2, 2, "left", "bottom"),
    "12056500": (-2, 0, "right", "center"),
    "2203": (-2, -2, "right", "top"),
    "2420": (2, -2, "left", "top"),
    "2481": (2, 2, "left", "bottom"),
    "2469": (-2, 2, "right", "bottom"),
}


def format_axis(ax: mpl.axes.Axes) -> None:
    """Apply consistent text sizes and line widths at the final figure size."""
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE, width=0.65, length=2.6, pad=2)
    ax.xaxis.label.set_fontsize(LABEL_SIZE)
    ax.yaxis.label.set_fontsize(LABEL_SIZE)
    ax.xaxis.labelpad = 3
    ax.yaxis.labelpad = 3
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
    for text in ax.texts:
        text.set_fontsize(LEGEND_SIZE)


def panel_letter(ax: mpl.axes.Axes, letter: str) -> None:
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(3, -3),
                textcoords="offset points", ha="left", va="top",
                fontsize=PANEL_SIZE, fontweight="bold", zorder=20)


def draw_gpp(ax: mpl.axes.Axes, region: str) -> None:
    """Panels a-b: plot regional GPP distributions with supplied spatial correlations."""
    data = pd.read_csv(DATA_DIR / f"gpp_distribution_{region}.csv")
    correlations = pd.read_csv(DATA_DIR / f"gpp_spatial_correlations_{region}.csv")
    validation.plot_gpp_distribution(data, None, correlations, include_slope=region == "us8", ax=ax)
    for line in ax.lines:
        line.set_linewidth(1.4)
    ymax = max(float(np.max(line.get_ydata())) for line in ax.lines)
    ax.set_ylim(0, ymax * (1.95 if region == "us8" else 1.65))
    ax.legend(loc="upper right", frameon=False, fontsize=LEGEND_SIZE, handlelength=2.2,
              handletextpad=0.45, labelspacing=0.22, borderaxespad=0.3)
    format_axis(ax)


def draw_metrics(ax: mpl.axes.Axes) -> mpl.axes.Axes:
    """Panel c: compare catchment-level relative differences and KGE values."""
    # The source table contains the precomputed validation metrics.
    data = pd.read_csv(DATA_DIR / "validation_relative_difference_kge_boxplot.csv")
    rel_names = ["MODIS_GPP", "US SLOPE_GPP", "QE", "DOC"]
    kge_names = ["Discharge SIM", "DOC"]
    rel = data[data["axis"].eq("Relative difference")]
    kge = data[data["axis"].eq("KGE")]
    rel_values = [rel.loc[rel["metric"].eq(name), "value"].dropna().to_numpy() for name in rel_names]
    kge_values = [kge.loc[kge["metric"].eq(name), "value"].dropna().to_numpy() for name in kge_names]
    rel_positions = [1, 2, 3, 4]
    kge_positions = [5.35, 6.35]
    right = ax.twinx()  # Relative differences use the left axis; KGE uses the right.
    # Fixed seed keeps the overlaid catchment points in the original positions.
    rng = np.random.default_rng(8)
    validation.draw_box_group(ax, rel_values, rel_positions, validation.RELATIVE_METRIC_COLORS, rng)
    validation.draw_box_group(right, kge_values, kge_positions, validation.KGE_METRIC_COLORS, rng)
    ax.set_xlim(0.5, 6.85)
    ax.set_ylim(-70, 125)
    right.set_ylim(-0.45, 1.2)
    ax.set_yticks([-50, 0, 50, 100])
    right.set_yticks([-0.25, 0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticks(rel_positions + kge_positions)
    ax.set_xticklabels(["MODIS\nGPP", "SLOPE\nGPP", r"$Q_E$", "DOC", "Discharge", "DOC"])
    ax.set_ylabel("Relative difference (%)")
    right.set_ylabel("KGE")
    right.spines["right"].set_visible(True)
    ax.hlines(0, 0.55, 4.45, color="#555555", lw=0.7, ls=(0, (4, 3)), zorder=0)
    ax.axvline(4.7, color="#b8b8b8", lw=0.6, zorder=0)
    for y in ax.get_yticks():
        ax.hlines(y, 0.55, 4.45, color="#dedede", lw=0.45, alpha=0.8, zorder=0)
    for y in right.get_yticks():
        right.hlines(y, 4.95, 6.8, color="#dedede", lw=0.45, alpha=0.8, zorder=0)
    for axis in (ax, right):
        format_axis(axis)
        for text in axis.texts:
            text.set_y(0.89)
            text.set_va("top")
        for collection in axis.collections:
            if isinstance(collection, PathCollection):
                collection.set_sizes([12])
                collection.set_linewidths([0.3])
        for line in axis.lines:
            line.set_linewidth(0.65)
        for patch in axis.patches:
            patch.set_linewidth(0.65)
    return right


def draw_doc_scatter(ax: mpl.axes.Axes, cax: mpl.axes.Axes | None, data: pd.DataFrame, flux: bool) -> None:
    """Panels d-e: compare matched DOC values on linear or logarithmic axes."""
    symbol = r"$F_{\mathrm{DOC}}$" if flux else r"$C_{\mathrm{DOC}}$"
    unit = r"(mg C m$^{-2}$ h$^{-1}$)" if flux else r"(mg C L$^{-1}$)"
    validation.scatter_1to1(
        data, "doc_flux_obs" if flux else "doc_obs", "doc_flux_sim" if flux else "doc_sim",
        f"Observed {symbol} {unit}", f"Simulated {symbol} {unit}", None,
        log_axes=flux, fixed_limits=None if flux else (0, 14),
        major_tick_step=None if flux else 4,
        log_major_ticks=[1e-5, 1e-3, 1e-1, 1e2] if flux else None,
        ax=ax, cax=cax, add_colorbar=cax is not None,
    )
    ax.set_aspect("auto")
    format_axis(ax)
    if cax is not None:
        format_axis(cax)
        cax.yaxis.label.set_fontsize(LEGEND_SIZE)
        cax.tick_params(labelsize=TICK_SIZE, length=2, pad=1.5)
    for text in ax.texts:
        text.set_position((0.25, 0.97))
        text.set_horizontalalignment("left")
        text.set_verticalalignment("top")
    for collection in ax.collections:
        if isinstance(collection, PathCollection):
            collection.set_sizes([10])
            collection.set_linewidths([0.18])
    for line in ax.lines:
        line.set_linewidth(0.75)


def draw_site_location_insets(ax: mpl.axes.Axes) -> None:
    """Add station maps using the supplied projected boundaries and coordinates."""
    boundaries = pd.read_csv(DATA_DIR / "map_boundaries.csv")
    stations = pd.read_csv(DATA_DIR / "map_stations.csv", dtype={"gauge_id": str}, float_precision="round_trip")
    labels = pd.read_csv(DATA_DIR / "site_labels.csv", dtype={"gauge_id": str})
    labels["gauge_id"] = labels["gauge_id"].map(validation.normalize_gauge)
    regions = (
        ("US", "EPSG:5070", [0.34, 0.140, 0.65, 0.205]),
        ("CH", "EPSG:2056", [0.38, 0.015, 0.61, 0.155]),
    )
    # One opaque background also covers the grid in the gap between the maps.
    bounds = np.asarray([region[2] for region in regions])
    lower = bounds[:, :2].min(axis=0) - 0.005
    upper = (bounds[:, :2] + bounds[:, 2:]).max(axis=0) + 0.005
    ax.add_patch(Rectangle(lower, *(upper - lower), transform=ax.transAxes,
                           facecolor="white", edgecolor="none", linewidth=0,
                           zorder=7, label="Site location background"))
    for region, crs, bounds in regions:
        polygons = boundaries.loc[boundaries["region"].eq(region)]
        country = gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_wkt(polygons["wkt"].tolist(), crs=crs))
        site_table = stations.loc[stations["region"].eq(region)].copy()
        sites = gpd.GeoDataFrame(site_table, geometry=gpd.points_from_xy(site_table["x"], site_table["y"]), crs=crs)
        expected = set(labels.loc[labels["plot_label"].str.startswith(region), "plot_label"])
        if set(sites["plot_label"]) != expected:
            raise ValueError(f"Incomplete {region} station labels for the location inset")

        inset = ax.inset_axes(bounds, zorder=8)
        inset.set_label(f"{region} site locations")
        inset.set_facecolor("white")
        country.plot(ax=inset, color="#dedede", edgecolor="#a0a0a0", linewidth=0.35,
                     aspect="equal", zorder=1)
        xmin, ymin, xmax, ymax = country.total_bounds
        width, height = xmax - xmin, ymax - ymin
        inset.set_xlim(xmin - 0.16 * width, xmax + 0.19 * width)
        inset.set_ylim(ymin - 0.12 * height, ymax + 0.26 * height)
        inset.scatter(sites.geometry.x, sites.geometry.y, s=14, color=MAP_REGION_COLORS[region],
                      edgecolor="none", linewidth=0, zorder=3)
        for site in sites.itertuples():
            dx, dy, ha, va = MAP_LABEL_OFFSETS[site.gauge_id]
            inset.annotate(site.plot_label, (site.geometry.x, site.geometry.y),
                           xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                           fontsize=MAP_LABEL_SIZE, color="#242424", zorder=4)
        inset.set_xticks([])
        inset.set_yticks([])
        for spine in inset.spines.values():
            spine.set_visible(False)


def build_figure() -> plt.Figure:
    """Assemble the six panels, shared formatting and station-location insets."""
    validation.setup_style()
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": TICK_SIZE, "axes.labelsize": LABEL_SIZE,
                         "legend.fontsize": LEGEND_SIZE, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="white")
    # Panel positions are [left, bottom, width, height] in figure coordinates.
    axes = {
        "a": fig.add_axes([0.090, 0.755, 0.218, 0.225]),
        "b": fig.add_axes([0.378, 0.755, LEFT_BLOCK_RIGHT - 0.378, 0.225]),
        "c": fig.add_axes([0.090, 0.420, 0.463, 0.250]),
        "d": fig.add_axes([0.090, 0.125, 0.180, 0.226]),
        "e": fig.add_axes([0.378, 0.125, 0.180, 0.226]),
        "f": fig.add_axes([0.680, 0.125, 0.305, 0.855]),
    }
    cax_e = fig.add_axes([0.568, 0.125, 0.009, 0.226])
    # Keep leading zeros in US station IDs; reuse these pairs in panels d-f.
    doc = pd.read_csv(DATA_DIR / "doc_1to1_12sites.csv", dtype={"gauge_id": str})
    doc["gauge_id"] = doc["gauge_id"].map(validation.normalize_gauge)
    draw_gpp(axes["a"], "us8")
    draw_gpp(axes["b"], "ch4")
    metrics_right = draw_metrics(axes["c"])
    draw_doc_scatter(axes["d"], None, doc, flux=False)
    draw_doc_scatter(axes["e"], cax_e, doc, flux=True)
    validation.plot_doc_violin(doc, None, ax=axes["f"])
    format_axis(axes["f"])
    axes["f"].set_xticks(np.arange(0, 13, 2))
    for line in axes["f"].lines:
        line.set_linewidth(1.05 if line.get_linewidth() > 1.1 else 0.65)
    for collection in axes["f"].collections:
        if isinstance(collection, PathCollection):
            collection.set_sizes([16])
    legend = axes["f"].get_legend()
    handles = legend.legend_handles
    axes["f"].legend(handles=handles, loc="upper right", ncol=2, frameon=False,
                     fontsize=LEGEND_SIZE, handlelength=1.4, handletextpad=0.4,
                     columnspacing=0.65, labelspacing=0.4, borderaxespad=0.4)
    # Single-line labels share baselines independent of tick-label extents.
    for name in ("d", "e", "f"):
        ax = axes[name]
        bounds = ax.get_position()
        ax.xaxis.set_label_coords(bounds.x0 + bounds.width / 2, 0.060,
                                 transform=fig.transFigure)
        ax.xaxis.label.set_verticalalignment("baseline")
    for names, x in (("acd", 0.038), ("be", 0.319)):
        for name in names:
            ax = axes[name]
            bounds = ax.get_position()
            ax.yaxis.set_label_coords(x, bounds.y0 + bounds.height / 2,
                                     transform=fig.transFigure)
            ax.yaxis.label.set_verticalalignment("center")
    # For rotated text, bottom alignment anchors its visible right edge.
    for ax in (metrics_right, cax_e):
        bounds = ax.get_position()
        ax.yaxis.set_label_coords(LEFT_BLOCK_RIGHT, bounds.y0 + bounds.height / 2,
                                 transform=fig.transFigure)
        ax.yaxis.label.set_rotation_mode("anchor")
        ax.yaxis.label.set_verticalalignment("bottom")
    for letter, ax in axes.items():
        panel_letter(ax, letter)
    draw_site_location_insets(axes["f"])
    return fig


def main() -> None:
    FIG1_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_figure()
    # Raster previews and editable vector formats share the same figure layout.
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(OUT_STEM.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(OUT_STEM.with_suffix(".svg"), dpi=600, facecolor="white")
    fig.savefig(OUT_STEM.with_suffix(".tif"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    print(OUT_STEM)
    plt.close(fig)


if __name__ == "__main__":
    main()
