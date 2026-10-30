"""Fit the DOC-change regression, project it across CONUS and draw Figure 5.

All inputs are local: basin predictors, a small regression training table,
display-resolution basin-ID rasters and map boundaries. Calculated predictions,
boxplot summaries, quadrant counts and covariance ellipses remain in memory.
"""
from __future__ import annotations

from pathlib import Path
import geopandas as gpd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import MultiPolygon, Polygon

# Resolve paths from this file, independently of the terminal directory.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
MAP_DIR = DATA_DIR / "maps"
OUTPUT_DIR = SCRIPT_DIR / "output"
FENNEMAN_GPKG = MAP_DIR / "boundaries.gpkg"
TIGER_CONUS_BOUNDARY_GPKG = MAP_DIR / "boundaries.gpkg"
X_COL = "climate_pressure_percent"
Y_COL = "landscape_component_percent"
FINAL_COL = "predicted_change_percent"
REGION_LABELS = {
    1: "Appalachian Highlands", 3: "Interior Highlands",
    4: "Interior Plains", 8: "Rocky Mountain System",
}
AXIS_LIMIT = 70.0
ELLIPSE_PROBABILITY = 0.50
CHI2_DF2_QUANTILE = -2.0 * np.log(1.0 - ELLIPSE_PROBABILITY)
ELLIPSE_SCALE = np.sqrt(CHI2_DF2_QUANTILE)
FOCUS_REGION_IDS = [1, 3, 4, 8]

REGION_COLORS = {
    1: "#4f63aa",  # Appalachian Highlands
    3: "#d63f4e",  # Interior Highlands
    4: "#f47646",  # Interior Plains
    8: "#9e0243",  # Rocky Mountain System
}

REGION_DRAW_ORDER = [4, 1, 3, 8]
OUTER_BOUNDARY_COLOR = "#4A4A4A"
OUTER_BOUNDARY_WIDTH = 0.45
_CONUS_OUTER_BOUNDARY = None

SCIENTIFIC_DIVERGING_COLORS = [
    (0.00, "#4b5578"),
    (0.22, "#8b8fa8"),
    (0.38, "#b7bed0"),
    (0.50, "#ffffff"),
    (0.62, "#efddc5"),
    (0.74, "#dabdad"),
    (0.84, "#d48f73"),
    (0.92, "#db7757"),
    (1.00, "#dc4f3b"),
]
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "slate_coral", SCIENTIFIC_DIVERGING_COLORS
)
SIGN_COLORS = {"positive": DIVERGING_CMAP(1.0), "negative": DIVERGING_CMAP(0.0)}



def basin_tiles() -> list[Path]:
    """Read the two spatial tiles in their original plotting order."""
    return sorted(MAP_DIR.glob("basin_ids_*.tif"))


def fit_projection_model() -> dict[str, float]:
    """Refit the selected relative-change model from 48 catchment-model records."""
    training = pd.read_csv(DATA_DIR / "regression_training.csv", dtype={"gauge_id": str},
                            float_precision="round_trip")
    predictors = ["delta_t_degC", "delta_pr_mm_yr", "mean_elevation_m",
                  "mean_distance_to_channel_m", "slope_p90_p10_deg"]
    baseline = training["current_doc_flux_mgC_m2_h"].to_numpy(float)
    if not np.isfinite(baseline).all() or (baseline <= 0).any():
        raise ValueError("Current DOC fluxes must be positive and finite.")
    response = training["delta_doc_flux_mgC_m2_h"].to_numpy(float) / baseline * 100.0
    design = np.column_stack([np.ones(len(training)), training[predictors].to_numpy(float)])
    if not np.isfinite(design).all() or not np.isfinite(response).all():
        raise ValueError("Regression training values must be finite.")
    coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
    return dict(zip(["Intercept", *predictors], coefficients))


def load_values() -> pd.DataFrame:
    """Calculate climate pressure, landscape component and their sum for each basin."""
    df = pd.read_csv(DATA_DIR / "basin_inputs.csv", float_precision="round_trip")
    if not np.isfinite(df.to_numpy(dtype=float)).all():
        raise ValueError("Basin inputs must be finite; use the supplied screened table.")
    coef = fit_projection_model()
    df[X_COL] = (coef["delta_t_degC"] * df["mean_dTa_degC_4models"]
                 + coef["delta_pr_mm_yr"] * df["mean_dPr_mm_yr_4models"])
    df[Y_COL] = (coef["Intercept"]
                 + coef["mean_distance_to_channel_m"] * df["mean_distance_to_channel_m"]
                 + coef["slope_p90_p10_deg"] * df["slope_p90_p10_deg"]
                 + coef["mean_elevation_m"] * df["mean_elevation_m"])
    df[FINAL_COL] = df[X_COL] + df[Y_COL]
    return df


def read_basin_tile(path: Path, values_by_id: np.ndarray):
    """Join calculated basin values to the supplied display-resolution ID grid."""
    with rasterio.open(path) as src:
        basin = src.read(1)
        out = np.full(basin.shape, np.nan, dtype=np.float32)
        valid = (basin > 0) & (basin < len(values_by_id))
        out[valid] = values_by_id[basin[valid]]
        extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)
    return np.ma.masked_invalid(out), extent


def conus_outer_boundary():
    """Cache the supplied CONUS outline; no download is needed."""
    global _CONUS_OUTER_BOUNDARY
    if _CONUS_OUTER_BOUNDARY is None:
        gdf = gpd.read_file(TIGER_CONUS_BOUNDARY_GPKG, layer="conus_boundary").to_crs("EPSG:5070")
        _CONUS_OUTER_BOUNDARY = gdf.geometry.union_all()
    return _CONUS_OUTER_BOUNDARY


def set_nature_style() -> None:
    """Set fonts and preserve editable text in vector exports."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.labelsize": 7.3,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.3,
            "legend.fontsize": 5.8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.65,
        }
    )


def values_array(values: pd.DataFrame, column: str) -> np.ndarray:
    """Build an ID lookup for assigning basin values to map cells."""
    max_id = int(values["coarse_basin_id"].max())
    out = np.full(max_id + 1, np.nan, dtype=np.float32)
    ids = values["coarse_basin_id"].to_numpy(dtype=np.int64)
    out[ids] = values[column].to_numpy(dtype=np.float32)
    return out


def draw_conus_outer_boundary(ax: plt.Axes) -> None:
    """Draw only the exterior coastline and national outline."""
    boundary = conus_outer_boundary()
    if boundary is None or boundary.is_empty:
        return
    polygons = boundary.geoms if isinstance(boundary, MultiPolygon) else [boundary]
    for polygon in polygons:
        if not isinstance(polygon, Polygon) or polygon.is_empty:
            continue
        x, y = polygon.exterior.xy
        ax.plot(
            x,
            y,
            color=OUTER_BOUNDARY_COLOR,
            linewidth=OUTER_BOUNDARY_WIDTH,
            alpha=0.95,
            zorder=6,
        )


def robust_abs_bound(series: pd.Series, percentile: float = 98.0) -> float:
    """Use the 98th percentile of absolute values for a symmetric colour scale."""
    vals = series.to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    return float(np.nanpercentile(np.abs(vals), percentile))


def draw_basin_map(
    ax: plt.Axes,
    values: pd.DataFrame,
    column: str,
    vmin: float,
    vmax: float,
    title: str,
    draw_regions: bool = False,
) -> mpl.image.AxesImage:
    """Draw the basin raster tiles and optional regional outlines."""
    value_lookup = values_array(values, column)
    image = None
    extents = []
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    for path in basin_tiles():
        data, extent = read_basin_tile(path, value_lookup)
        extents.append(extent)
        image = ax.imshow(
            data,
            extent=extent,
            origin="upper",
            cmap=DIVERGING_CMAP,
            norm=norm,
        )
    if image is None:
        raise RuntimeError("No map tiles were drawn.")
    ax.set_axis_off()
    ax.set_xlim(min(e[0] for e in extents), max(e[1] for e in extents))
    ax.set_ylim(min(e[2] for e in extents), max(e[3] for e in extents))
    ax.set_anchor("N")
    ax.set_title(title, pad=5)
    draw_conus_outer_boundary(ax)
    if draw_regions:
        draw_fenneman_boundaries(ax)
    return image


def add_signed_boxplot_inset(
    ax: plt.Axes,
    values: pd.DataFrame,
    column: str,
    norm: TwoSlopeNorm,
    *,
    fixed_sign_colors: bool = False,
    tick_labelsize: float = 5.6,
) -> None:
    """Summarize positive and negative values separately and label their medians."""
    vals = values[column].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    positive = vals[vals > 0]
    negative = vals[vals < 0]
    groups = [positive, negative]
    if len(positive) == 0 and len(negative) == 0:
        return

    inset = ax.inset_axes([0.095, -0.095, 0.215, 0.375])
    inset.set_facecolor("none")
    inset.patch.set_alpha(0.0)
    box = inset.boxplot(
        groups,
        positions=[1, 2],
        widths=0.50,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "white" if fixed_sign_colors else "0.05", "linewidth": 0.85},
        whiskerprops={"color": "0.18", "linewidth": 0.70},
        capprops={"color": "0.18", "linewidth": 0.70},
    )
    for patch, group, sign in zip(box["boxes"], groups, ["positive", "negative"]):
        if len(group) == 0:
            patch.set_facecolor("white")
            patch.set_edgecolor("0.40")
            continue
        color = SIGN_COLORS[sign] if fixed_sign_colors else DIVERGING_CMAP(norm(float(np.nanmean(group))))
        patch.set_facecolor(color)
        patch.set_edgecolor("0.18")
        patch.set_alpha(1.0 if fixed_sign_colors else 0.90)
        patch.set_linewidth(0.85)

    inset.axhline(0, color="0.50", linewidth=0.60, linestyle="--", zorder=0)
    inset.set_xlim(0.45, 2.55)
    lower, upper = np.nanpercentile(vals, [1, 99])
    lower = min(float(lower), 0.0)
    upper = max(float(upper), 0.0)
    span = max(upper - lower, 1.0)
    lower -= 0.08 * span
    upper += 0.08 * span
    tick_step = 5.0 if span <= 45 else 10.0
    lower_tick = tick_step * np.floor(lower / tick_step)
    upper_tick = tick_step * np.ceil(upper / tick_step)
    inset.set_ylim(lower_tick, upper_tick)
    inset.set_xticks([])
    median_ticks = []
    if len(negative) > 0:
        median_ticks.append(float(np.nanmedian(negative)))
    if len(positive) > 0:
        median_ticks.append(float(np.nanmedian(positive)))
    median_ticks = sorted(median_ticks)
    inset.set_yticks(median_ticks)
    inset.set_yticklabels([f"{tick:+.1f}" for tick in median_ticks])
    inset.tick_params(axis="y", labelsize=tick_labelsize, length=1.8, pad=1.0, width=0.45)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.spines["left"].set_color("0.30")
    inset.spines["bottom"].set_color("0.30")
    inset.spines["left"].set_linewidth(0.60)
    inset.spines["bottom"].set_linewidth(0.60)


def draw_fenneman_boundaries(ax: plt.Axes) -> None:
    """Overlay the four highlighted physiographic divisions in a fixed order."""
    divisions = gpd.read_file(FENNEMAN_GPKG, layer="fenneman_divisions")
    divisions = divisions[divisions["fenneman_division_id"].isin(FOCUS_REGION_IDS)].copy()
    divisions["draw_order"] = divisions["fenneman_division_id"].map(
        {region_id: order for order, region_id in enumerate(REGION_DRAW_ORDER)}
    )
    divisions = divisions.sort_values("draw_order")
    for _, row in divisions.iterrows():
        region_id = int(row["fenneman_division_id"])
        color = REGION_COLORS[region_id]
        geom = row.geometry
        polygons = geom.geoms if isinstance(geom, MultiPolygon) else [geom]
        for polygon in polygons:
            if not isinstance(polygon, Polygon) or polygon.is_empty:
                continue
            x, y = polygon.exterior.xy
            ax.plot(
                x,
                y,
                color=color,
                linewidth=0.90,
                alpha=0.96,
                zorder=5,
            )


def ellipse_from_points(x: np.ndarray, y: np.ndarray) -> dict[str, float] | None:
    """Calculate a 50% bivariate-normal covariance ellipse from basin values."""
    pts = np.column_stack([x, y])
    pts = pts[np.isfinite(pts).all(axis=1)]
    if len(pts) < 3:
        return None
    cov = np.cov(pts, rowvar=False)
    if not np.isfinite(cov).all():
        return None
    vals, vecs = np.linalg.eigh(cov)
    vals = np.maximum(vals, 0)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    if vals[0] <= 0:
        return None
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    center = pts.mean(axis=0)
    return {
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "width": float(2 * ELLIPSE_SCALE * np.sqrt(vals[0])),
        "height": float(2 * ELLIPSE_SCALE * np.sqrt(vals[1])),
        "angle_deg": float(angle),
    }


def make_ellipse_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize the component distribution within each highlighted division."""
    label_by_id = REGION_LABELS
    rows = []
    for region_id in FOCUS_REGION_IDS:
        group = df[df["fenneman_division_id"].eq(region_id)]
        ell = ellipse_from_points(group[X_COL].to_numpy(float), group[Y_COL].to_numpy(float))
        row = {
            "fenneman_division_id": region_id,
            "fenneman_division_label": label_by_id[region_id],
            "basin_count": len(group),
            "ellipse_probability": ELLIPSE_PROBABILITY if ell is not None else np.nan,
        }
        if ell is not None:
            row.update(ell)
        rows.append(row)
    out = pd.DataFrame(rows)
    return out


def make_quadrant_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Count basins by the signs of climate pressure and landscape component."""
    quadrants = [
        ("Q1", df[X_COL] >= 0, df[Y_COL] >= 0),
        ("Q2", df[X_COL] < 0, df[Y_COL] >= 0),
        ("Q3", df[X_COL] < 0, df[Y_COL] < 0),
        ("Q4", df[X_COL] >= 0, df[Y_COL] < 0),
    ]
    rows = []
    total_count = len(df)
    for label, x_mask, y_mask in quadrants:
        group = df[x_mask & y_mask]
        rows.append(
            {
                "quadrant": label,
                "basin_count": len(group),
                "basin_percent": len(group) / total_count * 100,
            }
        )
    out = pd.DataFrame(rows)
    return out


def add_quadrant_labels(ax: plt.Axes, quadrant_summary: pd.DataFrame) -> None:
    """Annotate basin counts and percentages for the four quadrants."""
    positions = {
        "Q1": (64, 8, "right", "bottom"),
        "Q2": (-64, 8, "left", "bottom"),
        "Q3": (-64, -8, "left", "top"),
        "Q4": (64, -8, "right", "top"),
    }
    for _, row in quadrant_summary.iterrows():
        x, y, ha, va = positions[row["quadrant"]]
        ax.text(
            x,
            y,
            f"{row['quadrant']}: {row['basin_percent']:.1f}%\n(n={int(row['basin_count']):,})",
            ha=ha,
            va=va,
            fontsize=7.7,
            color="0.12",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.70, "pad": 1.3},
            zorder=6,
        )


def add_change_fraction_labels(ax: plt.Axes, df: pd.DataFrame) -> None:
    """Show the percentages of basins with positive and negative total changes."""
    prediction_component = df[X_COL].to_numpy(dtype=float) + df[Y_COL].to_numpy(dtype=float)
    finite = np.isfinite(prediction_component)
    total = int(finite.sum())
    if total == 0:
        return
    increase_percent = float((prediction_component[finite] > 0).sum() / total * 100.0)
    decrease_percent = float((prediction_component[finite] < 0).sum() / total * 100.0)
    ax.text(
        -45,
        57,
        rf"{increase_percent:.1f}% $\uparrow$",
        ha="center",
        va="center",
        fontsize=8.5,
        fontweight="bold",
        color="#b2182b",
        rotation=-45,
        rotation_mode="anchor",
        zorder=7,
    )
    ax.text(
        -55,
        47,
        rf"{decrease_percent:.1f}% $\downarrow$",
        ha="center",
        va="center",
        fontsize=8.5,
        fontweight="bold",
        color="#2166ac",
        rotation=-45,
        rotation_mode="anchor",
        zorder=7,
    )


def draw_scatter(
    ax: plt.Axes,
    df: pd.DataFrame,
    ellipse_summary: pd.DataFrame,
    quadrant_summary: pd.DataFrame,
    vmin: float,
    vmax: float,
) -> None:
    """Plot all basin components, coloured by their predicted total change."""
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    ax.scatter(
        df[X_COL],
        df[Y_COL],
        c=df[FINAL_COL],
        s=1.35,
        alpha=0.42,
        linewidths=0,
        cmap=DIVERGING_CMAP,
        norm=norm,
        rasterized=True,
    )
    for _, row in ellipse_summary.iterrows():
        region_id = int(row["fenneman_division_id"])
        color = REGION_COLORS[region_id]
        ell = Ellipse(
            xy=(row["center_x"], row["center_y"]),
            width=row["width"],
            height=row["height"],
            angle=row["angle_deg"],
            facecolor="none",
            edgecolor=color,
            linewidth=1.05,
            alpha=0.98,
            label=row["fenneman_division_label"],
        )
        ax.add_patch(ell)
        ax.scatter(
            row["center_x"],
            row["center_y"],
            s=18,
            marker="o",
            facecolor="white",
            edgecolor=color,
            linewidth=1.05,
            zorder=5,
        )

    ax.axvline(0, color="0.25", linewidth=0.55)
    ax.axhline(0, color="0.25", linewidth=0.55)
    line_x = np.array([-AXIS_LIMIT, AXIS_LIMIT])
    ax.plot(line_x, -line_x, color="0.35", linewidth=0.6, linestyle="--", label="Reference line")
    ax.set_xlim(-AXIS_LIMIT, AXIS_LIMIT)
    ax.set_ylim(-AXIS_LIMIT, AXIS_LIMIT)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(-70, 71, 20))
    ax.set_yticks(np.arange(-70, 71, 20))
    ax.set_xlabel("Climate pressure (%)")
    ax.set_ylabel("Landscape component (%)")
    ax.grid(color="0.90", linewidth=0.30)
    add_quadrant_labels(ax, quadrant_summary)
    add_change_fraction_labels(ax, df)
    ax.legend(
        title="Fenneman divisions",
        loc="lower left",
        frameon=False,
        fontsize=7.6,
        title_fontsize=7.6,
        handlelength=1.0,
        handletextpad=0.35,
        labelspacing=0.25,
        borderpad=0,
        borderaxespad=0.2,
    )


def add_horizontal_colorbar(
    fig: plt.Figure,
    image: mpl.image.AxesImage,
    cax: plt.Axes,
    label: str,
) -> None:
    """Add a colourbar with the manuscript label sizes."""
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    cbar.set_label(label, fontsize=8.5, labelpad=2)
    cbar.ax.tick_params(labelsize=7.7, length=2, pad=1.5)


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
def build_figure(
    df: pd.DataFrame,
    ellipse_summary: pd.DataFrame,
    quadrant_summary: pd.DataFrame,
) -> plt.Figure:
    """Assemble three maps and the component scatter plot with their annotations."""
    climate_bound = robust_abs_bound(df[X_COL])
    baseline_bound = robust_abs_bound(df[Y_COL])
    final_bound = robust_abs_bound(df[FINAL_COL])
    common_map_bound = max(climate_bound, baseline_bound, final_bound)

    width_mm, height_mm = 230.0, 198.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))

    def axes_mm(left: float, bottom: float, width: float, height: float):
        return fig.add_axes([left / width_mm, bottom / height_mm,
                             width / width_mm, height / height_mm])

    ax_a = axes_mm(8, 127, 100, 60)
    ax_b = axes_mm(121, 127, 100, 60)
    ax_c = axes_mm(8, 37, 100, 60)
    ax_d = axes_mm(134, 20, 77, 77)
    cax_a = axes_mm(18, 114, 80, 4)
    cax_b = axes_mm(131, 114, 80, 4)
    cax_c = axes_mm(18, 24, 80, 4)

    image_a = draw_basin_map(
        ax_a,
        df,
        X_COL,
        -common_map_bound,
        common_map_bound,
        "",
        draw_regions=True,
    )
    image_b = draw_basin_map(
        ax_b,
        df,
        Y_COL,
        -common_map_bound,
        common_map_bound,
        "",
        draw_regions=True,
    )
    image_c = draw_basin_map(
        ax_c,
        df,
        FINAL_COL,
        -common_map_bound,
        common_map_bound,
        "",
        draw_regions=False,
    )
    for ax, column, image in [(ax_a, X_COL, image_a), (ax_b, Y_COL, image_b), (ax_c, FINAL_COL, image_c)]:
        add_signed_boxplot_inset(ax, df, column, image.norm, fixed_sign_colors=True, tick_labelsize=7.7)
    draw_scatter(ax_d, df, ellipse_summary, quadrant_summary, -common_map_bound, common_map_bound)
    ax_d.set_title("", pad=5)

    add_horizontal_colorbar(fig, image_a, cax_a, "Mean climate pressure (%)")
    add_horizontal_colorbar(fig, image_b, cax_b, "Landscape component (%)")
    add_horizontal_colorbar(fig, image_c, cax_c, r"Predicted relative change in $F_{\mathrm{DOC}}$ (%)")

    for panel, ax in zip("abcd", [ax_a, ax_b, ax_c, ax_d]):
        ax.text(
            0.0, 1.0, panel, transform=ax.title.get_transform(),
            ha="left", va=ax.title.get_verticalalignment(),
            fontsize=10.0, fontweight="bold", clip_on=False,
        )

    # Align the lower labels without changing the maps or scatter geometry.
    fig.canvas.draw()
    # Enlarge plot text while preserving all legend text and its layout.
    legend_texts = {text for legend in fig.findobj(mpl.legend.Legend)
                    for text in legend.findobj(mpl.text.Text)}
    for text in fig.findobj(mpl.text.Text):
        if text not in legend_texts:
            text.set_fontsize(text.get_fontsize() * width_mm / 183.0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_baselines = []
    for ax in [cax_c, ax_d]:
        label = ax.xaxis.label
        _, _, descent = renderer.get_text_width_height_descent(
            label.get_text(), label.get_fontproperties(), ismath="$" in label.get_text()
        )
        label_baselines.append((label.get_window_extent(renderer).y0 + descent) / fig.bbox.height)
    label_y = label_baselines[1]
    cax_c.set_position(cax_c.get_position().translated(0, label_y - label_baselines[0]))
    for ax in [cax_c, ax_d]:
        pos = ax.get_position()
        ax.xaxis.label.set_verticalalignment("baseline")
        ax.xaxis.set_label_coords((pos.x0 + pos.x1) / 2, label_y, transform=fig.transFigure)
    return fig


def main() -> None:
    """Run the projection and calculate all display summaries before exporting."""
    set_nature_style()
    values = load_values()
    ellipses = make_ellipse_summary(values)
    quadrants = make_quadrant_summary(values)
    fig = build_figure(values, ellipses, quadrants)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / "fig5.png", dpi=600, facecolor="white")
    fig.savefig(OUTPUT_DIR / "fig5.pdf", dpi=600, facecolor="white")
    fig.savefig(OUTPUT_DIR / "fig5.svg", dpi=600, facecolor="white")
    fig.savefig(OUTPUT_DIR / "fig5.tif", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"Saved Figure 5 for {len(values):,} basins to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
