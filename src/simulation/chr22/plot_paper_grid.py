"""Redraw the archived 16-setting chr22 paper grid as mean/CI bars with individual gene dots.

Generate exactly four PDFs for Supplementary Figures 7–10 from the recovered gene-level paper grid. The current chr22 scenario runs use a different design. No simulations are rerun and no statistics tables or preview files are written.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from seaborn.algorithms import bootstrap
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
METRIC_ORDER = [
    "alpha_null",
    "alpha_pop2_specific",
    "power_pop1_specific",
    "power_shared",
]

RESULT_DIR = Path(__file__).resolve().parents[3] / "bench" / "result"
INPUT_DIR = RESULT_DIR / "chr22_paper_grid"
OUTPUT_DIR = RESULT_DIR / "img" / "chr22_paper_grid_bars_dots"
BOOTSTRAP_SEED = 20260430
JITTER_SEED = 20260928

KEYS = ["metric", "h2sq", "n2", "propt", "method"]
PROPORTIONS = [0.01, 0.3, 0.6, 0.9]
YLIMITS = {
    "alpha_null": (0.0, 0.45),
    "alpha_pop2_specific": (0.0, 1.03),
    "power_pop1_specific": (0.0, 1.03),
    "power_shared": (0.0, 1.03),
}


def summarize(points: pd.DataFrame) -> pd.DataFrame:
    required = KEYS + ["gene_id", "value"]
    missing = set(required) - set(points.columns)
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")
    if points[required].isna().any().any():
        raise ValueError(
            "Missing gene identifiers, settings, or values in the paper grid."
        )
    if set(points.metric) != set(METRIC_ORDER):
        raise ValueError("The input must contain all four paper metrics and no others.")
    if points.duplicated(KEYS + ["gene_id"]).any():
        raise ValueError("Duplicate gene-level observations in the paper grid.")
    rows = []
    for key, group in points.groupby(KEYS, sort=True, observed=True):
        values = group.sort_values("gene_id")["value"].to_numpy()
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"Invalid proportions in {key}.")
        ci = np.percentile(
            bootstrap(values, func="mean", n_boot=1000, seed=BOOTSTRAP_SEED),
            [2.5, 97.5],
        )
        rows.append(
            dict(
                zip(KEYS, key),
                n=len(values),
                mean=values.mean(),
                ci_low=ci[0],
                ci_high=ci[1],
            )
        )
    return pd.DataFrame(rows)


def validate_grid(summary: pd.DataFrame, points: pd.DataFrame, metric: str) -> None:
    data = summary[summary.metric.eq(metric)]
    methods = METHOD_ORDER if metric.startswith("alpha") else METHOD_ORDER[:3]
    expected = pd.MultiIndex.from_product(
        [[0.1, 0.2], [100, 400], PROPORTIONS, methods], names=KEYS[1:]
    )
    actual = pd.MultiIndex.from_frame(data[KEYS[1:]])
    if actual.has_duplicates or set(actual) != set(expected):
        raise ValueError(f"Incomplete or unexpected parameter grid for {metric}.")
    ymin, ymax = YLIMITS[metric]
    metric_points = points[points.metric.eq(metric)]
    if not metric_points.value.between(ymin, ymax).all():
        raise ValueError(f"Individual values exceed the plotting range for {metric}.")
    if not (data.ci_low.ge(ymin).all() and data.ci_high.le(ymax).all()):
        raise ValueError(
            f"Confidence intervals exceed the plotting range for {metric}."
        )


def plot(
    summary: pd.DataFrame, points: pd.DataFrame, metric: str, output: Path
) -> None:
    data = summary[summary.metric.eq(metric)]
    metric_points = points[points.metric.eq(metric)]
    methods = METHOD_ORDER if metric.startswith("alpha") else METHOD_ORDER[:3]
    ymin, ymax = YLIMITS[metric]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.8), sharey=True)
    width = 0.78 / len(methods)
    offsets = (np.arange(len(methods)) - (len(methods) - 1) / 2) * width
    # Keep a gene's horizontal jitter identical across methods and settings.
    genes = sorted(metric_points.gene_id.unique())
    rng = np.random.default_rng(JITTER_SEED)
    jitter = dict(zip(genes, rng.uniform(-0.40, 0.40, len(genes))))
    dot_size, dot_alpha = (2.0, 0.16) if metric.startswith("alpha") else (2.6, 0.24)
    for row, h2sq in enumerate([0.1, 0.2]):
        for col, n2 in enumerate([100, 400]):
            ax = axes[row, col]
            facet = data[data.h2sq.eq(h2sq) & data.n2.eq(n2)]
            for method, offset in zip(methods, offsets):
                group = (
                    facet[facet.method.eq(method)].set_index("propt").loc[PROPORTIONS]
                )
                mean = group["mean"].to_numpy()
                x = np.arange(4) + offset
                color = COLOR_MAP[method]
                ax.bar(
                    x,
                    mean,
                    width=width * 0.90,
                    color=color,
                    edgecolor=color,
                    linewidth=0.75,
                    zorder=2,
                )
                for x_center, propt in zip(x, PROPORTIONS):
                    values = metric_points[
                        metric_points.h2sq.eq(h2sq)
                        & metric_points.n2.eq(n2)
                        & metric_points.propt.eq(propt)
                        & metric_points.method.eq(method)
                    ].sort_values("gene_id")
                    xs = x_center + values.gene_id.map(jitter).to_numpy() * width
                    ax.scatter(
                        xs,
                        values.value,
                        s=dot_size,
                        color="#343434",
                        alpha=dot_alpha,
                        linewidths=0,
                        zorder=3,
                        clip_on=False,
                    )
                ax.errorbar(
                    x,
                    mean,
                    yerr=np.array([mean - group.ci_low, group.ci_high - mean]),
                    fmt="none",
                    ecolor="#202020",
                    capsize=2.0,
                    elinewidth=0.9,
                    capthick=0.9,
                    zorder=5,
                )
            ax.set_title(rf"$h_2^2 = {h2sq:g},\ N_2 = {n2}$", fontsize=12)
            ax.set_xticks(range(4), [f"{p:g}" for p in PROPORTIONS])
            ax.set_xlim(-0.55, 3.55)
            ax.set_ylim(ymin, ymax)
            ax.set_xlabel(r"$\pi$")
            ax.set_ylabel(
                ("Type I error" if metric.startswith("alpha") else "Power")
                if col == 0
                else ""
            )
            ax.grid(axis="y", color="#e8e8e8", linewidth=0.7)
            ax.grid(axis="x", visible=False)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
            if metric.startswith("alpha"):
                ax.axhline(
                    0.05, color="#e63946", linestyle="--", linewidth=1.0, zorder=4
                )
            ticks = (
                [0, 0.05, 0.1, 0.2, 0.3, 0.4]
                if metric == "alpha_null"
                else [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            )
            ax.set_yticks(ticks)
    handles = [
        Patch(facecolor=COLOR_MAP[m], edgecolor=COLOR_MAP[m], label=METHOD_LABELS[m])
        for m in methods
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=len(methods),
        frameon=False,
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(output / f"{metric}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-points",
        type=Path,
        default=INPUT_DIR / "fig07_10_gene_points.tsv.gz",
        help="Recovered gene-level input table.",
    )
    parser.add_argument(
        "--reference-stats",
        type=Path,
        default=INPUT_DIR / "fig07_10_n_means_bootstrap95.tsv",
        help="Archived mean/CI table used to verify the computed statistics.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for the four PDF figures.",
    )
    args = parser.parse_args()
    sns.set_theme(
        style="whitegrid", font_scale=1.0, rc={"pdf.fonttype": 42, "ps.fonttype": 42}
    )
    points = pd.read_csv(args.gene_points, sep="\t")
    summary = summarize(points)
    reference = pd.read_csv(args.reference_stats, sep="\t")
    check = summary.merge(
        reference, on=KEYS, suffixes=("", "_reference"), validate="one_to_one"
    )
    if not len(check) == len(summary) == len(reference):
        raise ValueError(
            "Computed groups do not match the archived reference statistics."
        )
    for column in ["n", "mean", "ci_low", "ci_high"]:
        np.testing.assert_allclose(
            check[column], check[f"{column}_reference"], rtol=0, atol=1e-12
        )
    for metric in METRIC_ORDER:
        validate_grid(summary, points, metric)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for metric in METRIC_ORDER:
        plot(summary, points, metric, args.out_dir)
        print(f"Saved {args.out_dir / (metric + '.pdf')}")


if __name__ == "__main__":
    main()
