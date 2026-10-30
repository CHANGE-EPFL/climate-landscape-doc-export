"""Analyse annual process anomalies and draw Figure 3 from two local CSV files.

Fit catchment-specific regressions, partition R-squared using Shapley values,
and compare multivariate and univariate coefficients. Intermediate results stay
in memory; only the finished figure is saved.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import re


# Paths are relative to this script, independent of the terminal directory.
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
ANOMALY_CSV = DATA_DIR / "annual_anomalies.csv"
SITE_CSV = DATA_DIR / "sites.csv"
PREDICTORS = ["dPlant_input", "dQlat_out", "dDOC_uptake"]
SHARE_COLUMNS = ["Plant_input_share", "Lateral_flow_share", "DOC_uptake_share"]
PREDICTOR_LABELS = ["Plant input", "Lateral flow", "DOC uptake"]
RESPONSES = [
    {
        "key": "dDOC_export",
        "name": "DOC_export_plant_input",
        "panel": "a",
        "title": "Change in DOC export",
        "unit": r"mg C m$^{-2}$ yr$^{-1}$",
    },
    {
        "key": "dDOC_concentration",
        "name": "DOC_concentration_plant_input",
        "panel": "b",
        "title": "Change in DOC concentration",
        "unit": r"mg C L$^{-1}$",
    },
]

PREDICTOR_ORDER = ["Plant input", "Lateral flow", "DOC uptake"]
UNIVARIATE_TITLES = {
    "dDOC_export": "Change in DOC export",
    "dDOC_concentration": "Change in DOC concentration",
}
FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 136
LABEL_SIZE = 8.5
TICK_SIZE = 7.7
LEGEND_SIZE = 7.6
PANEL_SIZE = 10.0
POINT_SIZE = 11.5
POINT_EDGE_WIDTH = 0.55
SCATTER_EDGE_COLOR = "#3F3F3F"
UNIVARIATE_COLOR = "#7A7A7A"
MULTIVARIATE_COLOR = "#969696"
COEFFICIENT_EDGE_COLOR = "#737373"
SIGN_FLIP_COLOR = "#B23A48"


def subset_r2(X: np.ndarray, y: np.ndarray, columns: list[int]) -> float:
    """Fit a predictor subset with an intercept and return its R-squared."""
    if columns:
        design = np.column_stack([np.ones(len(y)), X[:, columns]])
        beta = np.linalg.lstsq(design, y, rcond=None)[0]
        predicted = design @ beta
    else:
        predicted = np.full_like(y, y.mean(), dtype=float)
    return float(1.0 - np.sum((y - predicted) ** 2) / np.sum((y - y.mean()) ** 2))


def shapley_r2(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Average each predictor's incremental R-squared over all subset orders."""
    from math import factorial

    p = X.shape[1]
    contributions = np.zeros(p, dtype=float)
    for predictor in range(p):
        for mask in range(2 ** p):
            if mask & (1 << predictor):
                continue
            columns = [j for j in range(p) if mask & (1 << j)]
            size = len(columns)
            weight = factorial(size) * factorial(p - size - 1) / factorial(p)
            contributions[predictor] += weight * (
                subset_r2(X, y, sorted(columns + [predictor])) - subset_r2(X, y, columns)
            )
    return contributions


def analyse_annual_anomalies() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, float]]:
    """Calculate every panel's inputs from 24 annual records per catchment."""
    data = pd.read_csv(ANOMALY_CSV, dtype={"GaugeID": str})
    numeric = ["dDOC_export", "dDOC_concentration", *PREDICTORS]
    if not np.isfinite(data[numeric].to_numpy(dtype=float)).all():
        raise ValueError("Annual anomalies must contain finite responses and predictors.")
    if data.duplicated(["GaugeID", "Scenario", "YearIndex"]).any():
        raise ValueError("Duplicate catchment-scenario-year records in annual_anomalies.csv.")

    predictions, contributions, coefficients = [], [], []
    pooled_r2 = {}
    for response in RESPONSES:
        response_key = response["key"]
        residual_sum = 0.0
        within_site_sum = 0.0
        for gauge_id, site in data.groupby("GaugeID", sort=True):
            X = site[PREDICTORS].to_numpy(dtype=float)
            y = site[response_key].to_numpy(dtype=float)
            x_sd = X.std(axis=0, ddof=1)
            y_sd = y.std(ddof=1)
            if len(y) <= X.shape[1] + 1 or y_sd == 0 or np.any(x_sd == 0):
                raise ValueError(f"Insufficient variation or records for {gauge_id}, {response_key}.")

            # Standardize across all four models and six years within this catchment.
            Xz = (X - X.mean(axis=0)) / x_sd
            yz = (y - y.mean()) / y_sd
            design = np.column_stack([np.ones(len(y)), Xz])
            beta = np.linalg.lstsq(design, yz, rcond=None)[0]
            fitted = (design @ beta) * y_sd + y.mean()
            pred = site[["GaugeID", "Scenario", "YearIndex"]].copy()
            pred["response"] = response_key
            pred["observed_anomaly"] = y
            pred["predicted_anomaly"] = fitted
            predictions.append(pred)

            # Pool residual sums while centring each response on its own site mean.
            residual_sum += float(np.sum((y - fitted) ** 2))
            within_site_sum += float(np.sum((y - y.mean()) ** 2))
            shares = shapley_r2(Xz, yz)
            shares = shares / float(np.nansum(shares))
            contributions.append({"GaugeID": gauge_id, "response": response_key,
                                  **dict(zip(SHARE_COLUMNS, shares))})

            for index, (predictor, label) in enumerate(zip(PREDICTORS, PREDICTOR_ORDER)):
                # Fit each predictor alone using the same sample-standard-deviation convention.
                x = site[predictor].to_numpy(dtype=float)
                xz = (x - x.mean()) / x.std(ddof=1)
                single_design = np.column_stack([np.ones(len(y)), xz])
                single_beta = np.linalg.lstsq(single_design, yz, rcond=None)[0]
                coefficients.append({
                    "GaugeID": gauge_id, "response": response_key, "predictor": label,
                    "univariate_standardized_beta": float(single_beta[1]),
                    "multivariate_standardized_beta": float(beta[index + 1]),
                })
        pooled_r2[response_key] = 1.0 - residual_sum / within_site_sum

    # Preserve the drawing order of overlapping points in the manuscript figure.
    site_order = {gauge_id: index for index, gauge_id in enumerate(
        pd.read_csv(SITE_CSV, dtype={"gauge_id": str})["gauge_id"])}

    def ordered(table: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
        table["_site_order"] = table["GaugeID"].map(site_order)
        return table.sort_values(["_site_order", "GaugeID", *columns]).drop(columns="_site_order")

    return (
        ordered(pd.concat(predictions, ignore_index=True), ["response", "Scenario", "YearIndex"]),
        ordered(pd.DataFrame(contributions), ["response"]),
        pd.DataFrame(coefficients),
        pooled_r2,
    )


def setup_style() -> None:
    """Set the manuscript fonts and preserve editable vector text."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": LABEL_SIZE,
            "axes.labelsize": LABEL_SIZE,
            "axes.titlesize": LABEL_SIZE,
            "xtick.labelsize": TICK_SIZE,
            "ytick.labelsize": TICK_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 0.72,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.68,
            "ytick.major.width": 0.68,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def site_label_map(gauge_ids: list[str]) -> dict[str, str]:
    """Read display labels while preserving leading zeros in gauge IDs."""
    sites = pd.read_csv(SITE_CSV, dtype={"gauge_id": str}).set_index("gauge_id")
    return sites.loc[gauge_ids, "plot_label"].to_dict()


def load_biome_map(gauge_ids: list[str]) -> dict[str, str]:
    """Read the vegetation groups used for marker shapes."""
    sites = pd.read_csv(SITE_CSV, dtype={"gauge_id": str}).set_index("gauge_id")
    return sites.loc[gauge_ids, "biome_group"].to_dict()


def site_label_sort_key(label: str) -> tuple[int, int, str]:
    """Order US sites before Swiss sites, using the numeric label suffix."""
    match = re.fullmatch(r"(US|CH)(\d+)", str(label))
    if match:
        return (0 if match.group(1) == "US" else 1), int(match.group(2)), str(label)
    return 2, 999, str(label)


def site_colors(gauge_ids: list[str]) -> dict[str, tuple[float, float, float, float]]:
    """Assign the site palette in the same order as the shared legend."""
    cmap = plt.get_cmap("tab20")
    labels = site_label_map(gauge_ids)
    ordered = sorted(gauge_ids, key=lambda gauge_id: site_label_sort_key(labels[gauge_id]))
    return {gauge_id: cmap(ii % 20) for ii, gauge_id in enumerate(ordered)}


def biome_marker_map(biomes: list[str]) -> dict[str, str]:
    """Assign a marker shape to each vegetation group."""
    preferred = {
        "Evergreen forest": "^",
        "Deciduous/mixed forest": "s",
        "Grassland": "o",
        "Alpine sparse": "D",
        "Cropland": "P",
    }
    fallback = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]
    markers: dict[str, str] = {}
    for biome in sorted(biomes):
        markers[biome] = preferred.get(biome, fallback[len(markers) % len(fallback)])
    return markers


def barycentric_to_xy(points: np.ndarray) -> np.ndarray:
    """Map plant, flow and uptake shares to the left, right and top vertices."""
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    return points @ vertices


def add_ternary_edge_ticks(ax: plt.Axes) -> None:
    """Draw tick marks around the three edges of the triangle."""
    tick_values = [0.2, 0.4, 0.6, 0.8]
    label_kwargs = {"fontsize": 6.9, "color": "#666666", "ha": "center", "va": "center"}
    tick_kwargs = {"color": "#777777", "lw": 0.48, "zorder": 2, "solid_capstyle": "butt"}
    tick_len = 0.018
    for value in tick_values:
        label = f"{int(value * 100)}%"

        # Bottom edge: lateral-flow fraction increases from left to right.
        xy = barycentric_to_xy(np.array([[1 - value, value, 0.0]])).ravel()
        ax.plot([xy[0], xy[0]], [xy[1], xy[1] - tick_len], **tick_kwargs)
        ax.text(xy[0], xy[1] - 0.034, label, **label_kwargs)

        # Left edge: DOC-uptake fraction increases upward from the plant-input vertex.
        xy = barycentric_to_xy(np.array([[1 - value, 0.0, value]])).ravel()
        ax.plot([xy[0], xy[0] - tick_len * 0.86], [xy[1], xy[1] + tick_len * 0.50], **tick_kwargs)
        ax.text(xy[0] - 0.050, xy[1] + 0.015, label, rotation=60, **label_kwargs)

        # Right edge: DOC-uptake fraction increases upward from the lateral-flow vertex.
        xy = barycentric_to_xy(np.array([[0.0, 1 - value, value]])).ravel()
        ax.plot([xy[0], xy[0] + tick_len * 0.86], [xy[1], xy[1] + tick_len * 0.50], **tick_kwargs)
        ax.text(xy[0] + 0.050, xy[1] + 0.015, label, rotation=-60, **label_kwargs)


def draw_ternary_background(ax: plt.Axes) -> None:
    """Draw the triangle, 20% grid lines and dashed 50% reference lines."""
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    ax.plot(
        [vertices[0, 0], vertices[1, 0], vertices[2, 0], vertices[0, 0]],
        [vertices[0, 1], vertices[1, 1], vertices[2, 1], vertices[0, 1]],
        color="#2B2B2B",
        lw=0.9,
    )
    for gg in np.arange(0.2, 1.0, 0.2):
        for a, b in [
            ([gg, 1 - gg, 0], [gg, 0, 1 - gg]),
            ([1 - gg, gg, 0], [0, gg, 1 - gg]),
            ([1 - gg, 0, gg], [0, 1 - gg, gg]),
        ]:
            line = barycentric_to_xy(np.array([a, b]))
            ax.plot(line[:, 0], line[:, 1], color="#D3D3D3", lw=0.45, zorder=0)
    fifty_lines = [
        (np.array([[0.5, 0.5, 0], [0.5, 0, 0.5]]), (0.34, 0.45), "50%"),
        (np.array([[0.5, 0.5, 0], [0, 0.5, 0.5]]), (0.66, 0.45), "50%"),
        (np.array([[0.5, 0, 0.5], [0, 0.5, 0.5]]), (0.50, 0.30), "50%"),
    ]
    for bary, label_xy, label in fifty_lines:
        line = barycentric_to_xy(bary)
        ax.plot(line[:, 0], line[:, 1], color="#7A7A7A", lw=0.75, ls=(0, (3, 2)), zorder=1)
        ax.text(label_xy[0], label_xy[1], label, ha="center", va="center", fontsize=8.0, color="#555555",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.5})
    add_ternary_edge_ticks(ax)
    ax.text(0.0, -0.080, PREDICTOR_LABELS[0], ha="left", va="top", fontsize=9.4)
    ax.text(1.0, -0.080, PREDICTOR_LABELS[1], ha="right", va="top", fontsize=9.4)
    ax.text(0.5, np.sqrt(3) / 2 + 0.032, PREDICTOR_LABELS[2], ha="center", va="bottom", fontsize=9.4)
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.16, np.sqrt(3) / 2 + 0.13)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_ternary_points(
    ax: plt.Axes,
    ternary: pd.DataFrame,
    response_key: str,
    panel: str,
    title: str,
    colors: dict[str, tuple[float, float, float, float]],
    markers: dict[str, str],
    point_size: float = 58,
    point_alpha: float = 1.0,
    point_edgecolor: str = "#1F1F1F",
    point_edgewidth: float = 0.45,
) -> None:
    """Place one point per catchment using the calculated Shapley shares."""
    draw_ternary_background(ax)
    data = ternary[ternary["response"] == response_key]
    shares = data[SHARE_COLUMNS].to_numpy(dtype=float)
    xy = barycentric_to_xy(shares)
    for idx, (_, row) in enumerate(data.iterrows()):
        ax.scatter(
            xy[idx, 0],
            xy[idx, 1],
            s=point_size,
            marker=markers[row["BiomeGroup"]],
            color=colors[row["GaugeID"]],
            edgecolor=point_edgecolor,
            linewidth=point_edgewidth,
            alpha=point_alpha,
            zorder=3,
        )
    ax.set_title(title, pad=16)
    panel_label(ax, panel)


def panel_label(ax: plt.Axes, label: str) -> None:
    """Place a bold panel letter inside an axis."""
    ax.text(
        0.025,
        0.975,
        label,
        transform=ax.transAxes,
        fontsize=PANEL_SIZE,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=20,
    )


def jitter(values: np.ndarray, width: float = 0.13) -> np.ndarray:
    """Separate overlapping coefficient pairs with deterministic vertical offsets."""
    if len(values) == 0:
        return values
    return values + np.linspace(-width, width, len(values))


def sort_by_site_label(data: pd.DataFrame) -> pd.DataFrame:
    """Keep coefficient pairs in the same US1-US8 and CH1-CH4 order."""
    labels = site_label_map(data["GaugeID"].dropna().astype(str).unique().tolist())
    sorted_data = data.copy()
    sorted_data["_site_order"] = sorted_data["GaugeID"].map(
        lambda gauge_id: site_label_sort_key(labels.get(str(gauge_id), str(gauge_id)))
    )
    sorted_data = sorted_data.sort_values(["_site_order", "GaugeID"])
    return sorted_data.drop(columns="_site_order")


def pad_xlim(values: np.ndarray) -> tuple[float, float]:
    """Leave space beside the coefficients for the sign-flip counts."""
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (-1.0, 1.0)
    xmin = float(values.min())
    xmax = float(values.max())
    span = max(xmax - xmin, 0.5)
    return xmin - 0.10 * span, xmax + 0.30 * span


def deepen_light_site_colors(
    colors: dict[str, tuple[float, float, float, float]],
    luminance_threshold: float = 0.55,
    darken_factor: float = 0.78,
) -> dict[str, tuple[float, float, float, float]]:
    """Darken the palest site colours for visibility on white."""
    adjusted: dict[str, tuple[float, float, float, float]] = {}
    for gauge_id, color in colors.items():
        rgba = mpl.colors.to_rgba(color)
        rgb = np.asarray(rgba[:3], dtype=float)
        luminance = float(0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2])
        if luminance >= luminance_threshold:
            rgb = np.clip(rgb * darken_factor, 0.0, 1.0)
        adjusted[gauge_id] = (float(rgb[0]), float(rgb[1]), float(rgb[2]), float(rgba[3]))
    return adjusted


def plot_prediction_panel(
    ax: plt.Axes,
    predictions: pd.DataFrame,
    response: dict[str, str],
    colors: dict,
    markers: dict[str, str],
    pooled_r2: float,
) -> None:
    """Compare regression predictions with the simulated annual anomalies."""
    data = predictions[predictions["response"].eq(response["key"])]
    for gauge_id, site in data.groupby("GaugeID", sort=True):
        ax.scatter(
            site["observed_anomaly"], site["predicted_anomaly"],
            s=12.5, marker=markers[site["BiomeGroup"].iloc[0]],
            facecolor=mpl.colors.to_rgba(colors[gauge_id], 0.73),
            edgecolor=SCATTER_EDGE_COLOR, linewidth=0.22, zorder=3,
        )
    values = data[["observed_anomaly", "predicted_anomaly"]].to_numpy(dtype=float)
    lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
    span = max(hi - lo, 1.0)
    limits = (lo - span * 0.08, hi + span * 0.18)
    ax.plot(limits, limits, color="#555555", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    symbol = r"$F_{\mathrm{DOC}}$" if response["key"] == "dDOC_export" else r"$C_{\mathrm{DOC}}$"
    ax.set_xlabel(f"Observed {symbol} anomaly ({response['unit']})", labelpad=3)
    ax.set_ylabel(f"Predicted {symbol} anomaly ({response['unit']})", labelpad=3)
    ax.text(0.14, 0.975, rf"pooled $R^2$ = {pooled_r2:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=LEGEND_SIZE, zorder=10)
    ax.grid(color="#DADADA", linewidth=0.45, alpha=0.70, zorder=0)
    ax.tick_params(direction="out", length=2.6, pad=2)
    if response["key"] == "dDOC_export":
        ax.set_xticks([-1000, 0, 1000, 2000, 3000])
        ax.set_yticks([-1000, 0, 1000, 2000, 3000])
    else:
        ax.set_xticks([-1, 0, 1, 2, 3])
        ax.set_yticks([-1, 0, 1, 2, 3])
    panel_label(ax, response["panel"])


def plot_ternary_panel(
    ax: plt.Axes, ternary: pd.DataFrame, response_key: str,
    panel: str, colors: dict, markers: dict[str, str],
) -> mpl.text.Text:
    # Shares keep the original vertices: plant (left), flow (right), uptake (top).
    """Label the three process shares and show their direction of increase."""
    draw_ternary_points(
        ax, ternary, response_key, panel, "", colors, markers,
        point_size=27, point_alpha=0.92,
        point_edgecolor=SCATTER_EDGE_COLOR, point_edgewidth=0.45,
    )
    # Rebuild edge annotations so each slanted axis names its own fraction.
    for text in list(ax.texts):
        text.remove()
    h = np.sqrt(3) / 2
    for value in (0.2, 0.4, 0.6, 0.8):
        percent = f"{int(value * 100)}%"
        kwargs = dict(fontsize=TICK_SIZE, color="#666666", ha="center", va="center")
        ax.text(1 - value, -0.045, percent, **kwargs)
        ax.text(0.5 * value - 0.057, h * value + 0.015, percent, rotation=60, **kwargs)
        ax.text(0.5 + 0.5 * value + 0.057, h * (1 - value) + 0.015,
                percent, rotation=-60, **kwargs)
    for xy in ((0.32, 0.45), (0.68, 0.45), (0.50, 0.29)):
        ax.text(*xy, "50%", fontsize=LEGEND_SIZE, color="#555555", ha="center", va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.2})
    # Arrows follow increasing shares toward uptake (top), flow (right), plant (left).
    ax.text(0.095, 0.49, r"DOC uptake $\longrightarrow$", rotation=60,
            ha="center", va="center", fontsize=LABEL_SIZE)
    ax.text(0.905, 0.49, r"Lateral flow $\longrightarrow$", rotation=-60,
            ha="center", va="center", fontsize=LABEL_SIZE)
    plant_label = ax.text(0.5, -0.125, r"$\longleftarrow$ Plant input",
                          ha="center", va="center", fontsize=LABEL_SIZE)
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.15, h + 0.04)
    ax.set_aspect("equal", adjustable="box")
    panel_label(ax, panel)
    return plant_label


def plot_univariate_panel(
    ax: plt.Axes,
    table: pd.DataFrame,
    response_key: str,
    panel: str,
    colors: dict[str, tuple[float, float, float, float]],
    markers: dict[str, str],
    show_title: bool = True,
) -> None:
    """Compare fitted coefficients and connect pairs with opposite signs."""
    data = table[table["response"].eq(response_key)].copy()
    x_values = data[["univariate_standardized_beta", "multivariate_standardized_beta"]].to_numpy(dtype=float).ravel()
    ax.set_xlim(*pad_xlim(x_values))
    y_positions = np.arange(len(PREDICTOR_ORDER))

    for yy, predictor in zip(y_positions, PREDICTOR_ORDER):
        subset = sort_by_site_label(data[data["predictor"].eq(predictor)])
        y_site = jitter(np.full(len(subset), yy), width=0.13)
        sign_flip = np.sign(subset["univariate_standardized_beta"]) != np.sign(
            subset["multivariate_standardized_beta"]
        )
        sign_flip = (
            sign_flip
            & subset["univariate_standardized_beta"].ne(0)
            & subset["multivariate_standardized_beta"].ne(0)
        )
        for site_y, (_, row), flipped in zip(y_site, subset.iterrows(), sign_flip):
            if flipped:
                ax.plot(
                    [row["univariate_standardized_beta"], row["multivariate_standardized_beta"]],
                    [site_y, site_y],
                    color=SIGN_FLIP_COLOR,
                    lw=0.95,
                    alpha=0.82,
                    zorder=2,
                )
            marker = markers.get(row["BiomeGroup"], "o")
            ax.scatter(
                row["univariate_standardized_beta"],
                site_y,
                s=POINT_SIZE,
                marker=marker,
                facecolor="white",
                edgecolor=UNIVARIATE_COLOR,
                linewidth=POINT_EDGE_WIDTH,
                alpha=0.96,
                zorder=3,
            )
            ax.scatter(
                row["multivariate_standardized_beta"],
                site_y,
                s=POINT_SIZE,
                marker=marker,
                facecolor=MULTIVARIATE_COLOR,
                edgecolor=COEFFICIENT_EDGE_COLOR,
                linewidth=POINT_EDGE_WIDTH,
                alpha=0.92,
                zorder=4,
            )
        ax.text(
            0.985,
            yy - 0.25,
            f"flip {int(sign_flip.sum())}/{len(subset)}",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=LEGEND_SIZE,
            color=SIGN_FLIP_COLOR if sign_flip.any() else "#666666",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.45},
        )

    ax.axvline(0, color="#555555", lw=0.75, ls=(0, (3, 3)), zorder=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(PREDICTOR_ORDER)
    ax.set_ylim(2.42, -0.52)
    ax.grid(axis="x", color="#DADADA", linewidth=0.42, alpha=0.70)
    ax.tick_params(direction="out", length=3.0, pad=2.0)
    ax.set_xlabel("Standardized coefficient")
    if show_title:
        ax.set_title(UNIVARIATE_TITLES[response_key], pad=6)
    panel_label(ax, panel)


def align_panel_labels(fig: plt.Figure, axes: np.ndarray, ternary_labels: list[mpl.text.Text]) -> None:
    # Figure coordinates keep text baselines aligned despite ternary aspect correction.
    """Align axis labels and panel letters across each row."""
    fig.canvas.draw()
    for row, plant_label in zip(axes, ternary_labels):
        row_bottom = row[0].get_position().y0
        row_top = row[0].get_position().y1
        label_y = row_bottom - 8.0 / FIGURE_HEIGHT_MM
        for ax in (row[0], row[2]):
            box = ax.get_position()
            ax.xaxis.set_label_coords((box.x0 + box.x1) / 2, label_y, transform=fig.transFigure)
            ax.xaxis.label.set_verticalalignment("baseline")
        box = row[1].get_position()
        plant_label.set_transform(fig.transFigure)
        plant_label.set_position(((box.x0 + box.x1) / 2, label_y))
        plant_label.set_verticalalignment("baseline")
        row[0].yaxis.set_label_coords(7.0 / FIGURE_WIDTH_MM, (row_bottom + row_top) / 2,
                                      transform=fig.transFigure)
        row[0].yaxis.label.set_verticalalignment("center")
        for ax in row:
            for text in ax.texts:
                if text.get_text() in "abcdef":
                    text.set_transform(fig.transFigure)
                    text.set_position((ax.get_position().x0 + 1.0 / FIGURE_WIDTH_MM,
                                       row_top - 1.05 / FIGURE_HEIGHT_MM))


def make_figure() -> plt.Figure:
    """Analyse the local annual data, assemble all six panels and export them."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    predictions, ternary, univariate_table, r2_by_response = analyse_annual_anomalies()
    gauge_ids = sorted(predictions["GaugeID"].unique())
    labels = site_label_map(gauge_ids)
    biome_labels = load_biome_map(gauge_ids)
    markers = biome_marker_map(list(set(biome_labels.values())))
    colors = deepen_light_site_colors(site_colors(gauge_ids))

    predictions["BiomeGroup"] = predictions["GaugeID"].map(biome_labels)
    ternary["BiomeGroup"] = ternary["GaugeID"].map(biome_labels)
    univariate_table["BiomeGroup"] = univariate_table["GaugeID"].map(biome_labels).fillna("Unknown")

    fig = plt.figure(figsize=(FIGURE_WIDTH_MM / 25.4, FIGURE_HEIGHT_MM / 25.4))

    def axes_mm(x: float, y: float, width: float, height: float) -> plt.Axes:
        return fig.add_axes([x / FIGURE_WIDTH_MM, y / FIGURE_HEIGHT_MM,
                             width / FIGURE_WIDTH_MM, height / FIGURE_HEIGHT_MM])

    axes = np.array([
        [axes_mm(21, 73, 41, 42), axes_mm(64, 65, 60, 52), axes_mm(141, 73, 40, 42)],
        [axes_mm(21, 15, 41, 42), axes_mm(64, 7, 60, 52), axes_mm(141, 15, 40, 42)],
    ])
    ternary_labels = []
    for row, response in enumerate(RESPONSES):
        plot_prediction_panel(axes[row, 0], predictions, response, colors, markers,
                              r2_by_response[response["key"]])
        ternary_labels.append(plot_ternary_panel(axes[row, 1], ternary, response["key"],
                                                 "cd"[row], colors, markers))
        plot_univariate_panel(axes[row, 2], univariate_table, response["key"], "ef"[row],
                              colors, markers, show_title=False)
    align_panel_labels(fig, axes, ternary_labels)

    legend_gauge_ids = sorted(gauge_ids, key=lambda gauge_id: site_label_sort_key(labels[gauge_id]))
    site_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker=markers[biome_labels[gauge_id]],
            linestyle="none",
            markerfacecolor=colors[gauge_id],
            markeredgecolor=SCATTER_EDGE_COLOR,
            markeredgewidth=0.38,
            markersize=4.2,
            label=labels[gauge_id],
        )
        for gauge_id in legend_gauge_ids
    ]
    # Six columns read left-to-right in two rows, keeping all sites above a-d.
    legend_rows = (len(site_handles) + 5) // 6
    site_handles = [site_handles[row * 6 + col] for col in range(6)
                    for row in range(legend_rows) if row * 6 + col < len(site_handles)]
    fig.legend(
        handles=site_handles,
        loc="upper center",
        bbox_to_anchor=(67 / FIGURE_WIDTH_MM, 135 / FIGURE_HEIGHT_MM),
        ncol=6,
        frameon=False,
        columnspacing=0.78,
        handletextpad=0.3,
        labelspacing=0.35,
        borderaxespad=0,
    )

    method_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor=UNIVARIATE_COLOR,
            markeredgewidth=POINT_EDGE_WIDTH,
            markersize=4.0,
            label="Univariate",
        ),
        mpl.lines.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=MULTIVARIATE_COLOR,
            markeredgecolor=COEFFICIENT_EDGE_COLOR,
            markeredgewidth=POINT_EDGE_WIDTH,
            markersize=4.0,
            label="Three-variable",
        ),
        mpl.lines.Line2D([0], [0], color=SIGN_FLIP_COLOR, lw=1.0, label="Sign flip"),
    ]
    biome_handles = [
        mpl.lines.Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="none",
            markerfacecolor="#252525",
            markeredgecolor="#1F1F1F",
            markeredgewidth=0.5,
            markersize=4.2,
            label=biome,
        )
        for biome, marker in markers.items()
    ]
    fig.legend(
        handles=biome_handles,
        loc="upper center",
        bbox_to_anchor=(67 / FIGURE_WIDTH_MM, 126.6 / FIGURE_HEIGHT_MM),
        ncol=len(biome_handles),
        frameon=False,
        columnspacing=0.70,
        handletextpad=0.30,
        borderaxespad=0,
    )
    fig.legend(
        handles=method_handles,
        loc="upper center",
        bbox_to_anchor=(158 / FIGURE_WIDTH_MM, 135 / FIGURE_HEIGHT_MM),
        ncol=1,
        frameon=False,
        handletextpad=0.4,
        labelspacing=0.35,
        borderaxespad=0,
    )

    out_stem = OUTPUT_DIR / "fig3"
    fig.savefig(out_stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(out_stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(out_stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(out_stem.with_suffix(".tif"), dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"Saved Figure 3 to {OUTPUT_DIR}")
    return fig


if __name__ == "__main__":
    plt.close(make_figure())
