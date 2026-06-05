"""Visualize chr22 eQTL mixture-architecture simulation summaries."""

from __future__ import annotations

import argparse
import re
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METHOD_ORDER = ["pop1_sumstat", "traceC", "traceCB", "meta", "meta_tissue"]
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
    "power_pop1_specific": "Power: pop1-specific genes",
    "alpha_null": "Type I error: null genes",
    "alpha_pop2_specific": "Type I error: pop2-specific genes",
}
METRIC_DEFAULT_YMAX = {
    "power_shared": 1.0,
    "power_pop1_specific": 1.0,
    "alpha_null": 0.38,
    "alpha_pop2_specific": 0.3,
}
METRIC_ALIASES = {
    "power": "power_shared",
    "alpha": "alpha_null",
}
PLOT_SETTING_COLS = ["h2sq", "n2", "nt", "propt", "method"]
FIXED_SETTING_COLS = ["h1sq", "gc", "n1", "pcausal"]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_path", default="bench/result/chr22_eqtl_mixture")
    parser.add_argument("--runname", default="mixture_chr22")
    parser.add_argument(
        "--metric",
        choices=[
            "all",
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
        help="Main plot weighting. CSV summaries are always saved for both weightings.",
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
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--ymin", type=float, default=0.0)
    parser.add_argument("--ymax", type=float, default=None)
    return parser.parse_args(argv)


def parse_setting_name(name: str) -> dict[str, float | int | str]:
    result: dict[str, float | int | str] = {}
    pattern = r"(h1sq|h2sq|gc|n1|n2|nt|propt|pcausal)_([^_]+)"
    for key, value in re.findall(pattern, name):
        try:
            number = float(value)
            result[key] = int(number) if number.is_integer() else number
        except ValueError:
            result[key] = value
    return result


def load_summaries(base_path: Path, runname: str) -> pd.DataFrame:
    rows = []
    for summary_path in sorted(base_path.glob(f"{runname}*/gene_summary.csv")):
        df = pd.read_csv(
            summary_path,
            keep_default_na=False,
            na_values=["", "nan", "NaN"],
        )
        params = parse_setting_name(summary_path.parent.name)
        for key, value in params.items():
            df[key] = value
        rows.append(df)
    if not rows:
        raise FileNotFoundError(f"No gene_summary.csv found under {base_path}/{runname}*")
    df = pd.concat(rows, ignore_index=True)
    df["method"] = pd.Categorical(df["method"], METHOD_ORDER, ordered=True)
    return df[df["method"].isin(METHOD_ORDER)]


def validate_plot_settings(df: pd.DataFrame) -> None:
    varying = []
    for col in FIXED_SETTING_COLS:
        if col in df and df[col].dropna().nunique() > 1:
            values = ", ".join(str(v) for v in sorted(df[col].dropna().unique()))
            varying.append(f"{col}=[{values}]")
    if varying:
        raise ValueError(
            "visual_chr22_eqtl_simulation.py plots only h2sq, n2, nt, and propt. "
            "Refusing to aggregate across varying settings: " + "; ".join(varying)
        )


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    valid = values.notna() & weights.notna()
    if not valid.any():
        return np.nan
    return float(np.average(values[valid], weights=weights[valid]))


def summarize_metric(df: pd.DataFrame, metric: str, weighting: str, error_unit: str) -> pd.DataFrame:
    metric_df = df[df[metric].notna()].copy()
    if metric_df.empty:
        return pd.DataFrame()

    setting_cols = PLOT_SETTING_COLS
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


def summarize_plot_data(plot_df: pd.DataFrame) -> pd.DataFrame:
    if plot_df.empty:
        return pd.DataFrame()
    rows = []
    for keys, group in plot_df.groupby(PLOT_SETTING_COLS, observed=True):
        rows.append(
            (
                *keys,
                weighted_mean(group["value"], group["weight"]),
                int(group["value"].notna().sum()),
                group["metric"].iloc[0],
                group["weighting"].iloc[0],
                group["error_unit"].iloc[0],
            )
        )
    return pd.DataFrame(
        rows,
        columns=PLOT_SETTING_COLS
        + ["mean", "n", "metric", "weighting", "error_unit"],
    )


def plot_metric(
    plot_df: pd.DataFrame,
    metric: str,
    out_path: Path,
    ymin: float,
    ymax: float | None,
    errorbar: str = "ci",
) -> None:
    if plot_df.empty:
        raise ValueError(f"No non-missing values for {metric}.")

    h2_values = sorted(plot_df["h2sq"].dropna().unique())
    n2_values = sorted(plot_df["n2"].dropna().unique())
    propt_values = sorted(plot_df["propt"].dropna().unique())
    methods = [method for method in METHOD_ORDER if method in set(plot_df["method"].astype(str))]

    fig, axes = plt.subplots(
        len(h2_values),
        len(n2_values),
        figsize=(4.4 * len(n2_values), 3.3 * len(h2_values)),
        squeeze=False,
        sharey=True,
    )
    width = min(0.14, 0.78 / max(len(methods), 1))
    x = np.arange(len(propt_values))
    legend_handles = None
    legend_labels = None

    for row_i, h2sq in enumerate(h2_values):
        for col_i, n2 in enumerate(n2_values):
            ax = axes[row_i][col_i]
            facet = plot_df[(plot_df["h2sq"] == h2sq) & (plot_df["n2"] == n2)]
            barplot_kwargs = {
                "data": facet,
                "x": "propt",
                "y": "value",
                "hue": "method",
                "order": propt_values,
                "hue_order": methods,
                "estimator": "mean",
                "errorbar": ("ci", 95) if errorbar == "ci" else "sd",
                "n_boot": 1000,
                "seed": 20260430,
                "dodge": 0.4,
                "linewidth": 0.5,
                "palette": COLOR_MAP,
                "saturation": 1.0,
                "capsize": 0.08,
                "err_kws": {"linewidth": 1.1},
                "legend": "auto" if legend_handles is None else False,
                "ax": ax,
            }
            if facet["weight"].nunique(dropna=False) > 1:
                barplot_kwargs["weights"] = "weight"
            sns.barplot(**barplot_kwargs)
            if legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()
            if ax.get_legend() is not None:
                ax.get_legend().remove()
            ax.set_title(rf"$h_2^2$={h2sq:g}, $N_2$={int(n2)}")
            ax.set_xticks(x)
            ax.set_xticklabels([str(v) for v in propt_values])
            ax.set_xlabel(r"$\pi$")
            ax.grid(axis="y", color="#E5E5E5", linewidth=0.8)
            ax.spines[["top", "right"]].set_visible(False)
            if col_i == 0:
                ax.set_ylabel(METRIC_LABELS.get(metric, metric))
            if metric.startswith("alpha"):
                ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=1.3)

    if ymax is None:
        ymax = METRIC_DEFAULT_YMAX.get(metric)
        if metric.startswith("alpha"):
            observed_max = float(plot_df["value"].max())
            ymax = max(ymax or 0.0, min(1.0, observed_max * 1.2))
    for ax in axes.flatten():
        ax.set_ylim(ymin, ymax)

    handles = legend_handles or []
    labels = [METHOD_LABELS.get(label, label) for label in (legend_labels or [])]
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=len(methods),
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {out_path}")


def metrics_to_plot(metric_arg: str) -> list[str]:
    if metric_arg == "all":
        return [
            "power_shared",
            "power_pop1_specific",
            "alpha_null",
            "alpha_pop2_specific",
        ]
    return [METRIC_ALIASES.get(metric_arg, metric_arg)]


def main() -> None:
    args = parse_args()
    base_path = Path(args.base_path)
    out_dir = Path(args.out_dir) if args.out_dir else base_path / "img"
    df = load_summaries(base_path, args.runname)
    validate_plot_settings(df)
    warned_no_rep_ci = False

    for metric in metrics_to_plot(args.metric):
        for weighting in ("gene_mean", "snp_weighted"):
            summary = summarize_metric(df, metric, weighting, args.error_unit)
            summary_table = summarize_plot_data(summary)
            if (
                args.error_unit == "replicate"
                and not warned_no_rep_ci
                and not summary.empty
                and summary_table["n"].max() <= 1
            ):
                print(
                    "Warning: each setting has only one replicate; replicate CI "
                    "error bars will be omitted. Re-run with --nrep > 1 to show "
                    "simulation uncertainty."
                )
                warned_no_rep_ci = True
            csv_path = out_dir / f"{args.runname}_{metric}_{weighting}_{args.error_unit}_summary.csv"
            csv_path.parent.mkdir(parents=True, exist_ok=True)
            summary_table.to_csv(csv_path, index=False)
            print(f"Saved summary to {csv_path}")
            if weighting == args.weighting:
                out_path = out_dir / f"{metric}.pdf"
                plot_metric(
                    summary,
                    metric,
                    out_path,
                    args.ymin,
                    args.ymax,
                    errorbar=args.errorbar,
                )


if __name__ == "__main__":
    main()
