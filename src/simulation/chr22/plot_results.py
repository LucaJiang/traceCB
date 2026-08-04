"""Visualize chr22 eQTL mixture-architecture simulation summaries."""

from __future__ import annotations

import argparse
import os
import re
import textwrap
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patheffects
from matplotlib.container import BarContainer
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm


METHOD_ORDER = ["pop1_sumstat", "traceC", "traceCB", "meta", "meta_tissue"]
PLOT_METHOD_ORDER = ["pop1_sumstat", "traceC", "traceCB", "meta", "meta_tissue"]
CORE_METHODS = {"pop1_sumstat", "traceC", "traceCB"}
RE2_METHODS = {"meta", "meta_tissue"}
METHOD_LABELS = {
    "pop1_sumstat": "Original",
    "traceC": "traceC",
    "traceCB": "traceCB",
    "meta": "RE2(sc)",
    "meta_tissue": "RE2(sc+tissue)",
}
COLOR_MAP = {
    "pop1_sumstat": "#47D45A",
    "traceC": "#fb8500",
    "traceCB": "#2d00f7",
    "meta": "#fcbf49",
    "meta_tissue": "#41cef1",
}
METRIC_LABELS = {
    "power_shared": "Power: shared genes",
    "power_pop1_specific": "Power: Population 1-specific genes",
    "alpha_null": "Type I error: null genes",
    "alpha_pop2_specific": "Type I error: Population 2-specific genes",
}
METRIC_YLABELS = {
    "power_shared": "Power",
    "power_pop1_specific": "Power",
    "alpha_null": "Type I error",
    "alpha_pop2_specific": "Type I error",
}
METRIC_DEFAULT_YMAX = {
    "power_shared": 0.39,
    "power_pop1_specific": 0.39,
    "alpha_null": 0.4,
    "alpha_pop2_specific": 0.4,
}
METRIC_ALIASES = {
    "power": "power_shared",
    "alpha": "alpha_null",
}
METRIC_ORDER = [
    "alpha_null",
    "alpha_pop2_specific",
    "power_pop1_specific",
    "power_shared",
]
SETTING_KEYS = ("h1sq", "h2sq", "gc", "n1", "n2", "nt", "propt", "pcausal")
RUN_SETTING_COLS = ["run_prefix", "method"]
SCRIPT_PATH = Path(__file__).with_name("run.sh")
DEFAULT_DISPLAY_RUN_PREFIX_ORDER = [
    "baseline",
    "propt_0.01",
    "h2sq_0.2",
    "n2_1000",
    "nt_20000",
]
PARAM_LABELS = {
    "h1sq": r"$h_1^2$",
    "h2sq": r"$h_2^2$",
    "gc": r"$\rho$",
    "n1": r"$N_1$",
    "n2": r"$N_2$",
    "nt": r"$N_t$",
    "propt": r"$\pi$",
    "pcausal": r"$p_{causal}$",
}
EGENE_Z_COLUMNS = {
    "pop1_sumstat": "z1_sumstat",
    "traceC": "z1_cross",
    "traceCB": "z1_tissue",
    "meta": "z_meta",
    "meta_tissue": "z_metatissue",
}
ARCHITECTURE_ORDER = ["shared", "pop1_specific", "null", "pop2_specific"]
POWER_ARCHITECTURES = {"shared", "pop1_specific"}
ARCHITECTURE_LABELS = {
    "shared": "shared genes",
    "pop1_specific": "Population 1-specific genes",
    "null": "null genes",
    "pop2_specific": "Population 2-specific genes",
}
DEFAULT_EGENE_P_THRESHOLD = 1e-5
DEFAULT_EGENE_CHUNKSIZE = 250_000


sns.set_theme(style="whitegrid", palette="muted", color_codes=True)


def set_alpha_yticks(ax, ymin: float, ymax: float) -> None:
    yticks = [
        tick
        for tick in [0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
        if ymin <= tick <= ymax
    ]
    ax.set_yticks(yticks)
    ax.set_yticklabels(
        [f"{tick:g}" if tick == 0 else f"{tick:.2f}" for tick in yticks]
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_path", default="bench/result/chr22_eqtl_mixture")
    parser.add_argument(
        "--runname",
        default="mixture_chr22",
        help=(
            "Output name stem. If matching result directories exist, it is also "
            "used as a directory prefix filter; otherwise all runs in base_path "
            "are plotted."
        ),
    )
    parser.add_argument(
        "--metric",
        choices=[
            "all",
            "egene",
            "power",
            "alpha",
            "power_shared",
            "power_pop1_specific",
            "alpha_null",
            "alpha_pop2_specific",
        ],
        default="all",
    )
    parser.add_argument(
        "--weighting",
        choices=["gene_mean", "snp_weighted"],
        default="gene_mean",
        help="Plot weighting.",
    )
    parser.add_argument(
        "--error_unit",
        choices=["gene", "replicate"],
        default="gene",
        help=(
            "Unit used to compute 95%% CI error bars. Gene uses genes as "
            "repeat units; replicate uses simulation replicates."
        ),
    )
    parser.add_argument(
        "--errorbar",
        choices=["ci", "sd"],
        default="ci",
        help="Error bar statistic used by seaborn: 95%% CI of the mean or SD.",
    )
    parser.add_argument("--out_dir", default=os.environ.get("IMG_DIR"))
    parser.add_argument(
        "--summary_out_dir",
        default=None,
        help="Directory for summary CSV files. Defaults to --base_path.",
    )
    parser.add_argument(
        "--output_name",
        default=None,
        help=(
            "Output PDF file name when plotting one metric. With --metric all, "
            "each gene architecture is always written to a separate PDF."
        ),
    )
    parser.add_argument(
        "--egene_output_prefix",
        default="egene_count",
        help=(
            "Filename prefix for the four architecture-specific eGene count "
            "PDFs in --out_dir."
        ),
    )
    parser.add_argument(
        "--egene_p_threshold",
        type=float,
        default=DEFAULT_EGENE_P_THRESHOLD,
        help=(
            "A gene is called an eGene when at least one of its SNPs has a "
            "two-sided p-value below this threshold."
        ),
    )
    parser.add_argument(
        "--egene_chunksize",
        type=int,
        default=DEFAULT_EGENE_CHUNKSIZE,
        help="Rows read at a time from each snp_results.csv.gz file.",
    )
    parser.add_argument(
        "--egene_ymax",
        type=float,
        default=None,
        help="Optional shared y-axis maximum for all eGene count plots.",
    )
    parser.add_argument(
        "--skip_egene",
        action="store_true",
        help="Skip eGene counting and its separate count plot.",
    )
    parser.add_argument("--ymin", type=float, default=0.0)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument(
        "--power_shared_ylim",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        default=None,
        help="Y-axis limits for the power_shared panel.",
    )
    parser.add_argument(
        "--power_pop1_specific_ylim",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        default=None,
        help="Y-axis limits for the power_pop1_specific panel.",
    )
    parser.add_argument(
        "--alpha_null_ylim",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        default=None,
        help="Y-axis limits for the alpha_null panel.",
    )
    parser.add_argument(
        "--alpha_pop2_specific_ylim",
        nargs=2,
        type=float,
        metavar=("YMIN", "YMAX"),
        default=None,
        help="Y-axis limits for the alpha_pop2_specific panel.",
    )
    parser.add_argument(
        "--re2",
        action="store_true",
        help="Load RE2(sc) and RE2(sc+tissue) rows and include them in all plots.",
    )
    parser.add_argument(
        "--run_prefix_order",
        nargs="+",
        default=None,
        help=(
            "Optional explicit run_prefix order. By default, the visualizer reads "
            "the order from run.sh and then appends extra prefixes."
        ),
    )
    return parser.parse_args(argv)


def parse_setting_name(name: str) -> dict[str, float | int | str]:
    result: dict[str, float | int | str] = {}
    pattern = rf"({'|'.join(SETTING_KEYS)})_([^_]+)"
    for key, value in re.findall(pattern, name):
        try:
            number = float(value)
            result[key] = int(number) if number.is_integer() else number
        except ValueError:
            result[key] = value
    return result


def extract_run_prefix(dirname: str) -> str:
    if "_h1sq_" in dirname:
        return dirname.split("_h1sq_", maxsplit=1)[0]
    first_setting = re.search(rf"(?:^|_)({'|'.join(SETTING_KEYS)})_", dirname)
    if first_setting and first_setting.start() > 0:
        return dirname[: first_setting.start()].rstrip("_")
    return dirname


def read_run_prefix_order(script_path: Path = SCRIPT_PATH) -> list[str]:
    if not script_path.exists():
        return []
    text = script_path.read_text()
    return re.findall(r"--run_prefix\s+([^\s\\]+)", text)


def discover_summary_paths(base_path: Path, runname: str) -> list[Path]:
    matched = sorted(base_path.glob(f"{runname}*/gene_summary.csv"))
    if matched:
        return matched
    return sorted(base_path.glob("*/gene_summary.csv"))


def load_summaries(base_path: Path, runname: str, include_re2: bool) -> pd.DataFrame:
    rows = []
    for summary_path in discover_summary_paths(base_path, runname):
        run_dir = summary_path.parent.name
        df = pd.read_csv(
            summary_path,
            keep_default_na=False,
            na_values=["", "nan", "NaN"],
        )
        df["run_prefix"] = extract_run_prefix(run_dir)
        df["run_dir"] = run_dir
        params = parse_setting_name(run_dir)
        for key, value in params.items():
            df[key] = value
        rows.append(df)
    if not rows:
        raise FileNotFoundError(f"No gene_summary.csv found under {base_path}")

    df = pd.concat(rows, ignore_index=True)
    df["method"] = pd.Categorical(df["method"], METHOD_ORDER, ordered=True)
    df = df[df["method"].isin(METHOD_ORDER)]
    if not include_re2:
        df = df[~df["method"].astype(str).isin(RE2_METHODS)]
    return df


def load_egene_counts(
    base_path: Path,
    run_dirs: Sequence[str],
    include_re2: bool,
    p_threshold: float,
    chunksize: int,
) -> pd.DataFrame:
    if not 0 < p_threshold < 1:
        raise ValueError("--egene_p_threshold must be between 0 and 1.")
    if chunksize <= 0:
        raise ValueError("--egene_chunksize must be greater than 0.")

    methods = [method for method in PLOT_METHOD_ORDER if method in CORE_METHODS]
    if include_re2:
        methods.extend(method for method in PLOT_METHOD_ORDER if method in RE2_METHODS)
    z_columns = {method: EGENE_Z_COLUMNS[method] for method in methods}
    usecols = ["rep", "gene_id", "architecture", *z_columns.values()]
    z_threshold = float(norm.isf(p_threshold / 2.0))
    rows = []

    for run_dir in run_dirs:
        snp_path = base_path / run_dir / "snp_results.csv.gz"
        if not snp_path.exists():
            raise FileNotFoundError(f"Missing SNP result file: {snp_path}")
        print(f"Counting eGenes from {snp_path}")

        max_abs_z_chunks = []
        for chunk in pd.read_csv(
            snp_path,
            usecols=usecols,
            chunksize=chunksize,
            keep_default_na=False,
            na_values=["", "nan", "NaN"],
        ):
            chunk.loc[:, list(z_columns.values())] = chunk[
                list(z_columns.values())
            ].abs()
            max_abs_z_chunks.append(
                chunk.groupby(
                    ["rep", "gene_id", "architecture"],
                    sort=False,
                    observed=True,
                )[list(z_columns.values())].max()
            )
        if not max_abs_z_chunks:
            raise ValueError(f"No SNP rows found in {snp_path}")

        gene_max_abs_z = (
            pd.concat(max_abs_z_chunks)
            .groupby(
                level=["rep", "gene_id", "architecture"],
                sort=False,
                observed=True,
            )
            .max()
        )
        run_prefix = extract_run_prefix(run_dir)
        setting_params = parse_setting_name(run_dir)
        for method, z_column in z_columns.items():
            for (rep, architecture), values in gene_max_abs_z[z_column].groupby(
                level=["rep", "architecture"], sort=True
            ):
                valid = values.notna()
                egene_count = int((values[valid] > z_threshold).sum())
                valid_gene_count = int(valid.sum())
                row = {
                    "run_prefix": run_prefix,
                    "run_dir": run_dir,
                    "rep": rep,
                    "architecture": str(architecture),
                    "method": method,
                    "egene_count": egene_count,
                    "valid_gene_count": valid_gene_count,
                    "egene_fraction": (
                        egene_count / valid_gene_count
                        if valid_gene_count > 0
                        else np.nan
                    ),
                    "p_threshold": p_threshold,
                    "z_threshold": z_threshold,
                }
                row.update(setting_params)
                rows.append(row)

    if not rows:
        raise ValueError("No eGene counts could be computed.")

    counts = pd.DataFrame(rows)
    comparison_keys = ["run_prefix", "run_dir", "rep", "architecture"]
    original_counts = counts[
        counts["method"].eq("pop1_sumstat")
    ][comparison_keys + ["egene_count"]].rename(
        columns={"egene_count": "original_egene_count"}
    )
    counts = counts.merge(original_counts, on=comparison_keys, how="left")
    counts["delta_vs_original"] = (
        counts["egene_count"] - counts["original_egene_count"]
    )
    counts["percent_change_vs_original"] = np.where(
        counts["original_egene_count"] > 0,
        100.0
        * counts["delta_vs_original"]
        / counts["original_egene_count"],
        np.nan,
    )
    counts["method"] = pd.Categorical(counts["method"], METHOD_ORDER, ordered=True)
    counts["architecture"] = pd.Categorical(
        counts["architecture"], ARCHITECTURE_ORDER, ordered=True
    )
    return counts


def ordered_run_prefixes(
    df: pd.DataFrame, explicit_order: Sequence[str] | None
) -> list[str]:
    present = list(dict.fromkeys(df["run_prefix"].astype(str)))
    if explicit_order:
        requested = list(explicit_order)
    elif set(DEFAULT_DISPLAY_RUN_PREFIX_ORDER).issubset(set(present)):
        requested = DEFAULT_DISPLAY_RUN_PREFIX_ORDER
    else:
        requested = read_run_prefix_order()
    if requested:
        ordered = [prefix for prefix in requested if prefix in present]
        ordered.extend(prefix for prefix in present if prefix not in ordered)
    else:
        ordered = sorted(present)
    return ordered


def apply_run_prefix_order(
    df: pd.DataFrame, explicit_order: Sequence[str] | None
) -> pd.DataFrame:
    order = ordered_run_prefixes(df, explicit_order)
    ordered_df = df[df["run_prefix"].astype(str).isin(order)].copy()
    ordered_df["run_prefix"] = pd.Categorical(
        ordered_df["run_prefix"].astype(str),
        categories=order,
        ordered=True,
    )
    return ordered_df


def summarize_metric(
    df: pd.DataFrame, metric: str, weighting: str, error_unit: str
) -> pd.DataFrame:
    metric_df = df[df[metric].notna()].copy()
    if metric_df.empty:
        return pd.DataFrame()
    if "run_prefix" not in metric_df.columns:
        metric_df["run_prefix"] = "default"

    setting_cols = RUN_SETTING_COLS
    if error_unit == "replicate":
        unit_cols = setting_cols + ["rep"]
        if weighting == "snp_weighted":
            weighted_df = metric_df.copy()
            weighted_df["_weighted_value"] = weighted_df[metric] * weighted_df["nsnp_gene"]
            plot_df = (
                weighted_df.groupby(unit_cols, observed=True)
                .agg(
                    weighted_sum=("_weighted_value", "sum"),
                    weight_sum=("nsnp_gene", "sum"),
                )
                .reset_index()
            )
            plot_df["value"] = plot_df["weighted_sum"] / plot_df["weight_sum"]
            plot_df["weight"] = 1.0
            plot_df = plot_df[unit_cols + ["value", "weight"]]
        else:
            plot_df = (
                metric_df.groupby(unit_cols, observed=True)[metric]
                .mean()
                .rename("value")
                .reset_index()
            )
            plot_df["weight"] = 1.0
    else:
        plot_df = (
            metric_df.groupby(setting_cols + ["gene_id"], observed=True)
            .agg(value=(metric, "mean"), weight=("nsnp_gene", "first"))
            .reset_index()
        )
        if weighting == "gene_mean":
            plot_df["weight"] = 1.0

    plot_df["metric"] = metric
    plot_df["weighting"] = weighting
    plot_df["error_unit"] = error_unit
    return plot_df


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    valid = values.notna() & weights.notna()
    if not valid.any():
        return np.nan
    return float(np.average(values[valid], weights=weights[valid]))


def summarize_plot_data(plot_df: pd.DataFrame) -> pd.DataFrame:
    if plot_df.empty:
        return pd.DataFrame()

    rows = []
    for keys, group in plot_df.groupby(RUN_SETTING_COLS, observed=True):
        run_prefix, method = keys
        rows.append(
            {
                "run_prefix": run_prefix,
                "method": method,
                "mean": weighted_mean(group["value"], group["weight"]),
                "n": int(group["value"].notna().sum()),
                "metric": group["metric"].iloc[0],
                "weighting": group["weighting"].iloc[0],
                "error_unit": group["error_unit"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def filter_metric_methods(
    plot_df: pd.DataFrame, metric: str, include_re2: bool
) -> pd.DataFrame:
    if plot_df.empty or "method" not in plot_df.columns:
        return plot_df.copy()
    allowed = set(CORE_METHODS)
    if metric.startswith("alpha") and include_re2:
        allowed.update(RE2_METHODS)
    return plot_df[plot_df["method"].astype(str).isin(allowed)].copy()


def methods_for_plot(plot_df: pd.DataFrame) -> list[str]:
    if plot_df.empty or "method" not in plot_df.columns:
        return []
    present = set(plot_df["method"].astype(str))
    return [method for method in PLOT_METHOD_ORDER if method in present]


def make_legend_handles(methods: list[str]) -> list[Patch]:
    return [
        Patch(
            facecolor=COLOR_MAP[method],
            edgecolor="#333333",
            linewidth=0.5,
            label=METHOD_LABELS.get(method, method),
        )
        for method in methods
    ]


def metric_ylim(
    args: argparse.Namespace, metric: str, plot_df: pd.DataFrame
) -> tuple[float, float | None]:
    custom_ylim = getattr(args, f"{metric}_ylim")
    if custom_ylim is not None:
        ymin, ymax = custom_ylim
        if ymin >= ymax:
            raise ValueError(f"--{metric}_ylim requires YMIN < YMAX.")
        return ymin, ymax

    ymin = args.ymin
    ymax = args.ymax
    if ymax is None:
        ymax = METRIC_DEFAULT_YMAX.get(metric)
        if metric.startswith("alpha") and not plot_df.empty:
            observed_max = float(plot_df["value"].max())
            ymax = max(ymax or 0.0, min(1.0, observed_max * 1.2))
    return ymin, ymax


def format_run_prefix_value(value: str) -> str:
    try:
        number = float(value)
        return f"{number:g}"
    except ValueError:
        return value.replace("_", " ")


def format_run_prefix_label(value: object) -> str:
    prefix = str(value)
    if prefix == "baseline":
        return "Baseline"

    for key in sorted(PARAM_LABELS, key=len, reverse=True):
        marker = f"{key}_"
        if prefix.startswith(marker):
            raw_value = prefix[len(marker) :]
            if raw_value:
                return f"{PARAM_LABELS[key]} = {format_run_prefix_value(raw_value)}"
    return prefix.replace("_", " ")


def wrap_run_prefix_labels(values: Sequence[object]) -> list[str]:
    return [
        "\n".join(
            textwrap.wrap(
                format_run_prefix_label(value),
                width=18,
                break_long_words=False,
                break_on_hyphens=False,
            )
        )
        for value in values
    ]


def plot_metric_axis(
    ax,
    plot_df: pd.DataFrame,
    metric: str,
    run_prefix_order: Sequence[str],
    ymin: float,
    ymax: float | None,
    errorbar: str = "ci",
    show_xlabel: bool = True,
    show_ylabel: bool = True,
) -> None:
    ax.set_title(METRIC_LABELS.get(metric, metric), fontsize=12, pad=9)
    if plot_df.empty:
        ax.text(0.5, 0.5, f"No data for {metric}", ha="center", va="center")
        ax.set_axis_off()
        return

    methods = methods_for_plot(plot_df)
    barplot_kwargs = {
        "data": plot_df,
        "x": "run_prefix",
        "y": "value",
        "hue": "method",
        "order": list(run_prefix_order),
        "hue_order": methods,
        "estimator": "mean",
        "errorbar": ("ci", 95) if errorbar == "ci" else "sd",
        "n_boot": 1000,
        "seed": 20260430,
        "dodge": True,
        "gap": 0.12,
        "linewidth": 0.5,
        "edgecolor": "#333333",
        "palette": COLOR_MAP,
        "saturation": 1.0,
        "capsize": 0.06,
        "err_kws": {"linewidth": 1.1},
        "legend": False,
        "ax": ax,
    }
    if plot_df["weight"].nunique(dropna=False) > 1:
        barplot_kwargs["weights"] = "weight"
    sns.barplot(**barplot_kwargs)

    ax.set_xlabel("Simulation setting" if show_xlabel else "", labelpad=7)
    ax.set_ylabel(METRIC_YLABELS.get(metric, metric) if show_ylabel else "")
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(range(len(run_prefix_order)))
    ax.set_xticklabels(wrap_run_prefix_labels(run_prefix_order), rotation=0)
    ax.grid(axis="y", color="#D8DDE6", linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.01)
    if metric.startswith("alpha"):
        set_alpha_yticks(ax, ymin, ymax if ymax is not None else 1.0)
        ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=1.3)


def plot_metric(
    plot_df: pd.DataFrame,
    metric: str,
    args: argparse.Namespace,
    out_path: Path,
) -> None:
    run_prefix_order = ordered_run_prefixes(plot_df, args.run_prefix_order)
    fig_width = max(10.8, 1.25 * len(run_prefix_order))
    fig, ax = plt.subplots(figsize=(fig_width, 4.4))
    ymin, ymax = metric_ylim(args, metric, plot_df)
    plot_metric_axis(
        ax,
        plot_df,
        metric,
        run_prefix_order,
        ymin,
        ymax,
        errorbar=args.errorbar,
    )

    legend_methods = methods_for_plot(plot_df)
    handles = make_legend_handles(legend_methods)
    fig.legend(
        handles=handles,
        labels=[handle.get_label() for handle in handles],
        fontsize=11,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=max(1, len(handles)),
        frameon=False,
    )
    fig.subplots_adjust(top=0.80, bottom=0.18, left=0.09, right=0.99)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {out_path}")


def format_scientific(value: float) -> str:
    return f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e+")


def egene_barplot(
    ax,
    data: pd.DataFrame,
    x: str,
    order: Sequence,
    errorbar: str,
) -> None:
    methods = methods_for_plot(data)
    sns.barplot(
        data=data,
        x=x,
        y="egene_count",
        hue="method",
        order=list(order),
        hue_order=methods,
        estimator="mean",
        errorbar=("ci", 95) if errorbar == "ci" else "sd",
        n_boot=1000,
        seed=20260430,
        dodge=True,
        gap=0.12,
        linewidth=0.5,
        edgecolor="#333333",
        palette=COLOR_MAP,
        saturation=1.0,
        capsize=0.06,
        err_kws={"linewidth": 1.1},
        legend=False,
        ax=ax,
    )
    annotate_egene_bars(ax, methods)


def annotate_egene_bars(ax, methods: Sequence[str]) -> None:
    bar_containers = [
        container
        for container in ax.containers
        if isinstance(container, BarContainer)
    ]
    observed_max = max(
        (
            float(bar.get_height())
            for container in bar_containers
            for bar in container.patches
            if np.isfinite(bar.get_height())
        ),
        default=0.0,
    )
    rotation = 90 if observed_max >= 100 else 0

    for method, container in zip(methods, bar_containers, strict=False):
        color = COLOR_MAP[method]
        for bar in container.patches:
            height = float(bar.get_height())
            if not np.isfinite(height):
                continue
            label = ax.annotate(
                f"{int(round(height))}",
                xy=(bar.get_x() + bar.get_width() / 2, max(height, 0.0)),
                xytext=(0, 2),
                textcoords="offset points",
                ha="center",
                va="bottom",
                rotation=rotation,
                fontsize=7.5,
                color=color,
                clip_on=False,
            )
            label.set_path_effects(
                [patheffects.withStroke(linewidth=1.1, foreground="white")]
            )


def egene_ymax(counts: pd.DataFrame, requested_ymax: float | None) -> float:
    if requested_ymax is not None:
        if requested_ymax <= 0:
            raise ValueError("--egene_ymax must be greater than 0.")
        return requested_ymax
    observed_max = float(counts["egene_count"].max())
    return max(1.0, np.ceil(observed_max * 1.15))


def is_historical_grid(counts: pd.DataFrame) -> bool:
    required = {"h2sq", "n2", "propt"}
    return (
        required.issubset(counts.columns)
        and counts["run_prefix"].astype(str).nunique() == 1
        and counts["h2sq"].dropna().nunique() > 1
        and counts["n2"].dropna().nunique() > 1
        and counts["propt"].dropna().nunique() > 1
    )


def add_egene_legend(fig, counts: pd.DataFrame, y: float = 0.995) -> None:
    handles = make_legend_handles(methods_for_plot(counts))
    fig.legend(
        handles=handles,
        labels=[handle.get_label() for handle in handles],
        fontsize=11,
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        ncol=max(1, len(handles)),
        frameon=False,
    )


def plot_egene_grid(
    counts: pd.DataFrame,
    architecture: str,
    args: argparse.Namespace,
    out_path: Path,
) -> None:
    h2_values = sorted(counts["h2sq"].dropna().unique())
    n2_values = sorted(counts["n2"].dropna().unique())
    propt_values = sorted(counts["propt"].dropna().unique())
    fig, axes = plt.subplots(
        len(h2_values),
        len(n2_values),
        figsize=(4.4 * len(n2_values), 3.3 * len(h2_values)),
        squeeze=False,
        sharey=True,
    )
    x = np.arange(len(propt_values))
    ymax = egene_ymax(counts, args.egene_ymax)

    for row_i, h2sq in enumerate(h2_values):
        for col_i, n2 in enumerate(n2_values):
            ax = axes[row_i][col_i]
            facet = counts[
                counts["h2sq"].eq(h2sq) & counts["n2"].eq(n2)
            ]
            egene_barplot(ax, facet, "propt", propt_values, args.errorbar)
            ax.set_title(rf"$h_2^2$={h2sq:g}, $N_2$={int(n2)}")
            ax.set_xticks(x)
            ax.set_xticklabels([f"{value:g}" for value in propt_values])
            ax.set_xlabel(r"$\pi$")
            ax.set_ylim(0, ymax)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.grid(axis="y", color="#E5E5E5", linewidth=0.8)
            ax.grid(axis="x", visible=False)
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_ylabel("")

    label = ARCHITECTURE_LABELS.get(architecture, architecture)
    fig.supylabel(f"Number of eGenes: {label}", x=0.01)
    add_egene_legend(fig, counts, y=1.0)
    fig.tight_layout(rect=(0.035, 0, 1, 0.90))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"eGene plot saved to {out_path}")


def plot_egene_scenarios(
    counts: pd.DataFrame,
    architecture: str,
    args: argparse.Namespace,
    out_path: Path,
) -> None:
    run_prefix_order = ordered_run_prefixes(counts, args.run_prefix_order)
    fig_width = max(10.8, 1.25 * len(run_prefix_order))
    fig, ax = plt.subplots(figsize=(fig_width, 4.4))
    egene_barplot(ax, counts, "run_prefix", run_prefix_order, args.errorbar)

    p_threshold = float(counts["p_threshold"].iloc[0])
    label = ARCHITECTURE_LABELS.get(architecture, architecture)
    ax.set_title(
        f"eGenes: {label} (at least one SNP with "
        f"two-sided p < {format_scientific(p_threshold)})",
        fontsize=12,
        pad=9,
    )
    ax.set_xlabel("Simulation setting", labelpad=7)
    ax.set_ylabel("Number of eGenes")
    ax.set_ylim(0, egene_ymax(counts, args.egene_ymax))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(range(len(run_prefix_order)))
    ax.set_xticklabels(wrap_run_prefix_labels(run_prefix_order), rotation=0)
    ax.grid(axis="y", color="#D8DDE6", linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(x=0.01)

    add_egene_legend(fig, counts)
    fig.subplots_adjust(top=0.80, bottom=0.18, left=0.09, right=0.99)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"eGene plot saved to {out_path}")


def plot_egene_counts(
    counts: pd.DataFrame,
    architecture: str,
    args: argparse.Namespace,
    out_path: Path,
) -> None:
    architecture_counts = counts[
        counts["architecture"].astype(str).eq(architecture)
    ].copy()
    if architecture in POWER_ARCHITECTURES:
        architecture_counts = architecture_counts[
            ~architecture_counts["method"].astype(str).isin(RE2_METHODS)
        ].copy()
    if architecture_counts.empty:
        raise ValueError(f"No eGene counts found for architecture={architecture}.")
    if is_historical_grid(architecture_counts):
        plot_egene_grid(architecture_counts, architecture, args, out_path)
    else:
        plot_egene_scenarios(architecture_counts, architecture, args, out_path)


def metrics_to_plot(metric_arg: str) -> list[str]:
    if metric_arg == "all":
        return METRIC_ORDER
    if metric_arg == "egene":
        return []
    return [METRIC_ALIASES.get(metric_arg, metric_arg)]


def resolve_output_path(
    args: argparse.Namespace, out_dir: Path, metric: str
) -> Path:
    if args.output_name:
        if args.metric == "all":
            raise ValueError("--output_name can only be used with a single metric.")
        out_path = Path(args.output_name)
        if not out_path.suffix:
            out_path = out_path.with_suffix(".pdf")
        if not out_path.is_absolute():
            out_path = out_dir / out_path
        return out_path

    return out_dir / (
        f"{args.runname}_{metric}_{args.weighting}_{args.error_unit}.pdf"
    )


def resolve_egene_output_path(
    args: argparse.Namespace, out_dir: Path, architecture: str
) -> Path:
    prefix = Path(args.egene_output_prefix)
    if prefix.suffix:
        prefix = prefix.with_suffix("")
    out_path = prefix.parent / f"{prefix.name}_{architecture}.pdf"
    if not out_path.is_absolute():
        out_path = out_dir / out_path
    return out_path


def warn_if_replicate_ci_unavailable(summary: pd.DataFrame, warned: bool) -> bool:
    if warned or summary.empty or "rep" not in summary.columns:
        return warned
    rep_counts = summary.groupby(RUN_SETTING_COLS, observed=True)["rep"].nunique()
    if not rep_counts.empty and rep_counts.max() <= 1:
        print(
            "Warning: each setting has only one replicate; replicate CI "
            "error bars will be omitted. Re-run with --nrep > 1 to show "
            "simulation uncertainty."
        )
        return True
    return warned


def write_summary_table(
    summary: pd.DataFrame, metric: str, args: argparse.Namespace, out_dir: Path
) -> None:
    table = summarize_plot_data(summary)
    if table.empty:
        return
    csv_path = (
        out_dir
        / f"{args.runname}_{metric}_{args.weighting}_{args.error_unit}_summary.csv"
    )
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(csv_path, index=False)
    print(f"Summary saved to {csv_path}")


def write_egene_count_table(
    counts: pd.DataFrame, args: argparse.Namespace, out_dir: Path
) -> None:
    csv_path = out_dir / f"{args.runname}_egene_count_summary.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    sort_cols = [
        col
        for col in [
            "run_prefix",
            "h2sq",
            "n2",
            "propt",
            "rep",
            "architecture",
            "method",
        ]
        if col in counts.columns
    ]
    counts.sort_values(sort_cols).to_csv(csv_path, index=False)
    print(f"eGene count summary saved to {csv_path}")


def main() -> None:
    args = parse_args()
    base_path = Path(args.base_path)
    out_dir = Path(args.out_dir) if args.out_dir else base_path / "img"
    summary_out_dir = Path(args.summary_out_dir) if args.summary_out_dir else base_path
    df = load_summaries(base_path, args.runname, args.re2)
    df = apply_run_prefix_order(df, args.run_prefix_order)
    warned_no_rep_ci = False
    metrics = metrics_to_plot(args.metric)
    plotted_any_metric = False

    for metric in metrics:
        summary = summarize_metric(df, metric, args.weighting, args.error_unit)
        summary = filter_metric_methods(summary, metric, args.re2)
        if args.error_unit == "replicate":
            warned_no_rep_ci = warn_if_replicate_ci_unavailable(
                summary, warned_no_rep_ci
            )
        if summary.empty:
            print(f"Warning: no non-missing values for {metric}.")
        else:
            write_summary_table(summary, metric, args, summary_out_dir)
            plot_metric(
                summary,
                metric,
                args,
                resolve_output_path(args, out_dir, metric),
            )
            plotted_any_metric = True

    if metrics and not plotted_any_metric:
        raise ValueError("No non-missing metric values to plot.")

    if args.metric == "egene" and args.skip_egene:
        raise ValueError("--metric egene cannot be combined with --skip_egene.")

    if not args.skip_egene:
        run_dirs = list(dict.fromkeys(df["run_dir"].astype(str)))
        egene_counts = load_egene_counts(
            base_path,
            run_dirs,
            args.re2,
            args.egene_p_threshold,
            args.egene_chunksize,
        )
        egene_counts = apply_run_prefix_order(egene_counts, args.run_prefix_order)
        write_egene_count_table(egene_counts, args, summary_out_dir)
        for architecture in ARCHITECTURE_ORDER:
            plot_egene_counts(
                egene_counts,
                architecture,
                args,
                resolve_egene_output_path(args, out_dir, architecture),
            )


if __name__ == "__main__":
    main()
