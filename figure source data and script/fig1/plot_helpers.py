"""Panel drawing functions for Figure 1.

Validation metrics are read from the source tables. Density curves, boxplots
and descriptive summaries are calculated from the supplied values for plotting.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import LogFormatterMathtext
from scipy.stats import gaussian_kde

DATA_DIR = Path(__file__).resolve().parent / "data"
SITE_LABELS_FILE = DATA_DIR / "site_labels.csv"

# Panel defaults; plot_fig1.py adjusts text sizes for the assembled figure.
SMALL_PANEL_LABEL_FONTSIZE = 18.0
SMALL_PANEL_TICK_FONTSIZE = 16.0
SMALL_PANEL_LEGEND_FONTSIZE = 15.0
RIGHT_PANEL_LABEL_FONTSIZE = 18.0
RIGHT_PANEL_TICK_FONTSIZE = 16.0
RIGHT_PANEL_LEGEND_FONTSIZE = 15.0
RELATIVE_METRIC_COLORS = ["#4e659b", "#8a8dc0", "#b7a8cf", "#fdcf9e"]
KGE_METRIC_COLORS = ["#efa484", "#b5756b"]
GPP_COLORS = {"MODIS": RELATIVE_METRIC_COLORS[0], "SLOPE": RELATIVE_METRIC_COLORS[1], "Simulation": "#D55E00"}
GPP_LINESTYLES = {"MODIS": "-", "SLOPE": "--", "Simulation": ":"}

def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 11,
            "axes.labelsize": 12,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "legend.fontsize": 10.5,
            "axes.linewidth": 0.9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

def normalize_gauge(x: object) -> str:
    """Keep Swiss IDs unchanged and format US IDs with eight digits."""
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    if s in {"2203", "2420", "2469", "2481"}:
        return s
    return s.zfill(8)

def q_bias_cmap() -> LinearSegmentedColormap:
    """Use green for matched discharge and blue/brown for under/overestimation."""
    colors = [
        (0.00, "#8aa3bf"),
        (0.18, "#d5e0ea"),
        (0.36, "#4f9a66"),
        (0.50, "#137a3a"),
        (0.64, "#4f9a66"),
        (0.82, "#efd1b8"),
        (1.00, "#b66a45"),
    ]
    return LinearSegmentedColormap.from_list("q_bias_center_green", colors, N=256)

def add_discharge_bias(df: pd.DataFrame) -> pd.DataFrame:
    """Derive Q_sim/Q_obs for colouring the DOC scatter plots."""
    out = df.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        # Flux / concentration gives simulated discharge in mm h^-1.
        q_sim = out["doc_flux_sim"].to_numpy(dtype=float) / out["doc_sim"].to_numpy(dtype=float)
        q_bias = q_sim / out["q_obs_mmh"].to_numpy(dtype=float)
    out["q_sim_mmh"] = q_sim
    out["q_bias"] = q_bias
    return out

def prepare_scatter_data(
    df: pd.DataFrame,
    obs_col: str,
    sim_col: str,
    log_axes: bool = False,
    min_positive: float | None = None,
) -> pd.DataFrame:
    """Keep finite matched values, requiring positive pairs on logarithmic axes."""
    cols = list(dict.fromkeys([obs_col, sim_col, "q_obs_mmh", "doc_sim", "doc_flux_sim"]))
    d = add_discharge_bias(df[cols].copy())
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    if min_positive is not None:
        d = d[(d[obs_col] >= min_positive) & (d[sim_col] >= min_positive)].copy()
    if log_axes:
        d = d[(d[obs_col] > 0) & (d[sim_col] > 0)].copy()
    return d

def scatter_1to1(
    df: pd.DataFrame,
    obs_col: str,
    sim_col: str,
    xlabel: str,
    ylabel: str,
    out: Path | None,
    log_axes: bool = False,
    min_positive: float | None = None,
    fixed_limits: tuple[float, float] | None = None,
    major_tick_step: float | None = None,
    log_major_ticks: list[float] | None = None,
    title: str | None = None,
    ax: mpl.axes.Axes | None = None,
    cax: mpl.axes.Axes | None = None,
    add_colorbar: bool = True,
) -> tuple[float, int]:
    """Draw a 1:1 comparison and annotate the precomputed R-squared and count."""
    d = prepare_scatter_data(df, obs_col, sim_col, log_axes, min_positive)
    panel = "e" if log_axes else "d"
    stats = pd.read_csv(DATA_DIR / "scatter_statistics.csv", float_precision="round_trip").set_index("panel").loc[panel]
    n, r2 = int(stats["n"]), stats["r2"]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(4.7, 4.2), constrained_layout=True)
    else:
        fig = ax.figure
    q_bias = d["q_bias"].to_numpy(dtype=float)
    color_mask = np.isfinite(q_bias) & (q_bias > 0)
    # Clip the colour scale only; point coordinates retain their supplied values.
    color_value = np.log10(np.clip(q_bias[color_mask], 0.1, 10.0))
    sc = ax.scatter(
        d.loc[color_mask, obs_col],
        d.loc[color_mask, sim_col],
        s=38,
        c=color_value,
        cmap=q_bias_cmap(),
        vmin=-1,
        vmax=1,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.25,
        zorder=3,
    )

    vals = np.r_[d[obs_col].to_numpy(), d[sim_col].to_numpy()]
    vals = vals[np.isfinite(vals)]
    if log_axes:
        vals = vals[vals > 0]
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    if fixed_limits is not None:
        lo, hi = fixed_limits
    elif log_axes:
        lo = 10 ** np.floor(np.log10(lo))
        hi = 10 ** np.ceil(np.log10(hi))
    else:
        pad = 0.055 * (hi - lo)
        lo -= pad
        hi += pad
    ax.plot([lo, hi], [lo, hi], color="black", lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    if major_tick_step is not None:
        ticks = np.arange(lo, hi + 0.5 * major_tick_step, major_tick_step)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
    if log_axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        if log_major_ticks is not None:
            ax.set_xticks(log_major_ticks)
            ax.set_yticks(log_major_ticks)
            formatter = LogFormatterMathtext(base=10)
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
    else:
        ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(xlabel, fontsize=SMALL_PANEL_LABEL_FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=SMALL_PANEL_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=SMALL_PANEL_TICK_FONTSIZE)
    if title is not None:
        ax.set_title(title, pad=6, fontsize=16.0)
    ax.grid(True, color="#d8d8d8", lw=0.5, alpha=0.8)
    ax.text(
        0.05,
        0.95,
        f"$R^2$ = {r2:.2f}\n$n$ = {n}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=SMALL_PANEL_LEGEND_FONTSIZE,
    )
    if add_colorbar:
        cb = fig.colorbar(sc, ax=ax, cax=cax, pad=0.025, fraction=0.055)
        cb.set_ticks([-1, 0, 1])
        cb.set_ticklabels(["-1", "0", "1"])
        cb.set_label(r"$\log_{10}(Q_{sim}/Q_{obs})$", fontsize=SMALL_PANEL_LABEL_FONTSIZE)
        cb.ax.tick_params(labelsize=SMALL_PANEL_TICK_FONTSIZE)
    if out is not None:
        fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
    if standalone:
        plt.close(fig)
    return r2, n

def kde_line(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Evaluate a Gaussian density estimate using all finite source values."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or np.nanstd(values) == 0:
        y = np.zeros_like(grid)
        if values.size:
            y[np.argmin(np.abs(grid - np.nanmedian(values)))] = 1.0
        return y
    return gaussian_kde(values)(grid)

def kde_for_violin(values: np.ndarray, global_xmax: float, n_grid: int = 320) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate each half-violin on its own nonnegative concentration grid."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.array([]), np.array([])
    local_max = max(np.nanmax(values) * 1.15, np.nanpercentile(values, 99) * 1.2, 0.25)
    local_max = min(max(local_max, np.nanmax(values) + 0.15), global_xmax)
    grid = np.linspace(0, local_max, n_grid)
    if values.size < 3 or np.nanstd(values) == 0:
        center = np.nanmedian(values)
        bw = max(0.08 * max(local_max, 1.0), 0.05)
        y = np.exp(-0.5 * ((grid - center) / bw) ** 2)
    else:
        kde = gaussian_kde(values, bw_method="scott")
        y = kde(grid)
    y = np.asarray(y, dtype=float)
    if y.size:
        # Close the half-violin at both ends of its plotted range.
        y[0] = 0.0
        y[-1] = 0.0
    return grid, y

def plot_gpp_distribution(
    df: pd.DataFrame,
    out: Path | None,
    correlations: pd.DataFrame,
    include_slope: bool = True,
    ax: mpl.axes.Axes | None = None,
) -> None:
    """Draw GPP density curves; legend correlations come from matched-pixel analyses."""
    order = ["MODIS", "SLOPE", "Simulation"] if include_slope else ["MODIS", "Simulation"]
    corr_map = correlations.set_index("comparison")["r"].to_dict()
    labels = {
        "MODIS": f"MODIS (r = {corr_map['SIM-MODIS']:.2f})",
        "Simulation": "Simulation",
    }
    if include_slope:
        labels["SLOPE"] = f"SLOPE (r = {corr_map['SIM-SLOPE']:.2f})"
    x = df["value"].replace([np.inf, -np.inf], np.nan).dropna()
    # Limit the displayed range, while retaining all finite values for the KDE.
    xmin, xmax = np.nanpercentile(x, [0.2, 99.8])
    grid = np.linspace(max(0, xmin), xmax, 400)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
    else:
        fig = ax.figure
    ymax = 0.0
    for src in order:
        vals = df.loc[df["source"] == src, "value"].to_numpy()
        y = kde_line(vals, grid)
        if np.any(np.isfinite(y)):
            ymax = max(ymax, float(np.nanmax(y)))
        ax.plot(grid, y, color=GPP_COLORS[src], ls=GPP_LINESTYLES[src], lw=2.2, label=labels[src])
    ax.set_xlabel("GPP (g C m$^{-2}$ d$^{-1}$)", fontsize=SMALL_PANEL_LABEL_FONTSIZE)
    ax.set_ylabel("Density", fontsize=SMALL_PANEL_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=SMALL_PANEL_TICK_FONTSIZE)
    if ymax > 0:
        ax.set_ylim(0, ymax * 1.45)
    ax.grid(True, color="#dedede", lw=0.5, alpha=0.75)
    ax.legend(
        frameon=False,
        loc="upper left",
        borderaxespad=0.25,
        ncol=1,
        fontsize=SMALL_PANEL_LEGEND_FONTSIZE,
        handlelength=1.35,
        handletextpad=0.45,
        labelspacing=0.18,
    )
    if out is not None:
        fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
    if standalone:
        plt.close(fig)

def draw_box_group(
    ax: mpl.axes.Axes,
    data: list[np.ndarray],
    positions: list[float],
    colors: list[str],
    rng: np.random.Generator,
    n_fontsize: float = 9,
) -> None:
    """Draw quartiles and 1.5-IQR whiskers, with every catchment value overlaid."""
    bp = ax.boxplot(
        data,
        positions=positions,
        patch_artist=True,
        widths=0.58,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "black", "linewidth": 1.25},
        whiskerprops={"color": "#333333", "linewidth": 1.0},
        capprops={"color": "#333333", "linewidth": 1.0},
        boxprops={"edgecolor": "#333333", "linewidth": 1.0},
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.52)
    for position, vals, color in zip(positions, data, colors):
        if vals.size == 0:
            continue
        # Horizontal jitter separates overlapping points without changing values.
        x = rng.normal(position, 0.045, size=vals.size)
        ax.scatter(
            x,
            vals,
            s=30,
            color=color,
            edgecolor="black",
            linewidth=0.35,
            alpha=0.88,
            zorder=3,
        )
        ax.text(
            position,
            1.015,
            f"n={vals.size}",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=n_fontsize,
            clip_on=False,
        )

def plot_doc_violin(df: pd.DataFrame, out: Path | None, ax: mpl.axes.Axes | None = None) -> None:
    """Panel f: paired DOC distributions with median lines and mean markers."""
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["doc_obs", "doc_sim"]).copy()
    df["gauge_label"] = df["gauge_id"].astype(str).map(lambda x: x if len(x) <= 4 else x.zfill(8))
    site_labels = pd.read_csv(
        SITE_LABELS_FILE,
        dtype={"plot_label": "string", "gauge_id": "string"},
    )
    site_labels["gauge_id"] = site_labels["gauge_id"].map(normalize_gauge)
    available = set(df["gauge_label"])
    site_labels = site_labels.loc[site_labels["gauge_id"].isin(available)].copy()
    order = site_labels["gauge_id"].tolist()
    label_map = dict(zip(site_labels["gauge_id"], site_labels["plot_label"]))

    all_vals = np.r_[df["doc_obs"].to_numpy(), df["doc_sim"].to_numpy()]
    finite_vals = all_vals[np.isfinite(all_vals)]
    xmax = np.nanmax(finite_vals) * 1.03
    xmin = max(0, np.nanmin(finite_vals) - 0.03 * (np.nanmax(finite_vals) - np.nanmin(finite_vals)))
    grid = np.linspace(xmin, xmax, 350)
    width = 0.34
    obs_color = "#4C78A8"
    sim_color = "#D55E00"

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.15, 12.45), constrained_layout=True)
    else:
        fig = ax.figure
    for yi, gid in enumerate(order):
        sub = df[df["gauge_label"] == gid]
        obs = sub["doc_obs"].to_numpy()
        sim = sub["doc_sim"].to_numpy()
        x_obs, y_obs = kde_for_violin(obs, xmax)
        x_sim, y_sim = kde_for_violin(sim, xmax)
        # Scale each half independently so shapes can be compared across sites.
        if np.nanmax(y_obs) > 0:
            y_obs = y_obs / np.nanmax(y_obs) * width
        if np.nanmax(y_sim) > 0:
            y_sim = y_sim / np.nanmax(y_sim) * width
        ax.fill_between(x_obs, yi, yi - y_obs, color=obs_color, alpha=0.34, lw=0)
        ax.fill_between(x_sim, yi, yi + y_sim, color=sim_color, alpha=0.34, lw=0)
        ax.plot(x_obs, yi - y_obs, color=obs_color, lw=1.65)
        ax.plot(x_sim, yi + y_sim, color=sim_color, lw=1.65)

        med_obs = np.nanmedian(obs)
        med_sim = np.nanmedian(sim)
        mean_obs = np.nanmean(obs)
        mean_sim = np.nanmean(sim)
        # Place summary markers on the corresponding half-violin boundary.
        med_obs_y = np.interp(med_obs, x_obs, y_obs)
        med_sim_y = np.interp(med_sim, x_sim, y_sim)
        mean_obs_y = np.interp(mean_obs, x_obs, y_obs)
        mean_sim_y = np.interp(mean_sim, x_sim, y_sim)
        ax.plot([med_obs, med_obs], [yi, yi - med_obs_y], color="#222222", lw=1.05, zorder=5)
        ax.plot([med_sim, med_sim], [yi, yi + med_sim_y], color="#222222", lw=1.05, zorder=5)
        ax.scatter([mean_obs], [yi - mean_obs_y], s=34, color=obs_color, edgecolor="white", linewidth=0.35, zorder=6)
        ax.scatter([mean_sim], [yi + mean_sim_y], s=34, color=sim_color, edgecolor="white", linewidth=0.35, zorder=6)

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([label_map[g] for g in order])
    ax.set_ylim(len(order) - 0.35, -1.35)
    ax.set_xlabel("DOC concentration (mg C L$^{-1}$)", fontsize=RIGHT_PANEL_LABEL_FONTSIZE)
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=RIGHT_PANEL_TICK_FONTSIZE)
    ax.tick_params(axis="both", width=1.15, length=4.2)
    ax.spines["left"].set_linewidth(1.15)
    ax.spines["bottom"].set_linewidth(1.15)
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", color="#dedede", lw=0.5, alpha=0.75)
    handles = [
        mpl.lines.Line2D([0], [0], color=obs_color, lw=2, label="Observed"),
        mpl.lines.Line2D([0], [0], color=sim_color, lw=2, label="Simulated"),
        mpl.lines.Line2D([0], [0], color="#222222", lw=0.9, label="Median"),
        mpl.lines.Line2D([0], [0], marker="o", color="none", markerfacecolor="#666666", markersize=5.5, label="Mean"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="upper right",
        ncol=2,
        handlelength=1.35,
        columnspacing=0.75,
        labelspacing=0.25,
        fontsize=RIGHT_PANEL_LEGEND_FONTSIZE,
    )
    if out is not None:
        fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
    if standalone:
        plt.close(fig)
