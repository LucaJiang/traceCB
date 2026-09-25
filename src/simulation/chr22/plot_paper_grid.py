"""Redraw the archived 16-setting chr22 paper grid as mean/CI line plots.

The current chr22 scenario runs use a different design. This entry point reads
the recovered, gene-level paper grid explicitly rather than mixing the runs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from seaborn.algorithms import bootstrap
import seaborn as sns

from plot_results import COLOR_MAP, METHOD_LABELS, METHOD_ORDER, METRIC_ORDER

KEYS = ["metric", "h2sq", "n2", "propt", "method"]
MARKERS = dict(zip(METHOD_ORDER, ["o", "D", "D", "s", "s"]))
PROPORTIONS = [0.01, 0.3, 0.6, 0.9]
YLIMITS = {
    "alpha_null": (0.0, 0.45),
    "alpha_pop2_specific": (0.0, 0.65),
    "power_pop1_specific": (0.2, 0.4),
    "power_shared": (0.25, 0.38),
}


def summarize(points: pd.DataFrame) -> pd.DataFrame:
    if points.duplicated(KEYS + ["gene_id"]).any():
        raise ValueError("Duplicate gene-level observations in the paper grid.")
    rows = []
    for key, group in points.groupby(KEYS, sort=True, observed=True):
        values = group.sort_values("gene_id")["value"].to_numpy()
        if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
            raise ValueError(f"Invalid proportions in {key}.")
        ci = np.percentile(bootstrap(values, func="mean", n_boot=1000, seed=20260430), [2.5, 97.5])
        rows.append(dict(zip(KEYS, key), n=len(values), mean=values.mean(), ci_low=ci[0], ci_high=ci[1]))
    return pd.DataFrame(rows)


def plot(summary: pd.DataFrame, metric: str, output: Path) -> None:
    data = summary[summary.metric.eq(metric)]
    methods = METHOD_ORDER if metric.startswith("alpha") else METHOD_ORDER[:3]
    expected = pd.MultiIndex.from_product([[0.1, 0.2], [100, 400], PROPORTIONS, methods], names=KEYS[1:])
    actual = pd.MultiIndex.from_frame(data[KEYS[1:]])
    if actual.has_duplicates or set(actual) != set(expected):
        raise ValueError(f"Incomplete or unexpected parameter grid for {metric}.")
    ymin, ymax = YLIMITS[metric]
    clipped = (data.ci_low.lt(ymin) | data.ci_high.gt(ymax)).sum()
    if clipped:
        print(f"Note: {metric}: {clipped} confidence intervals extend beyond the requested y-axis range {ymin}–{ymax}; statistics are unchanged.")
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.0), sharey=True)
    offsets = np.linspace(-0.16, 0.16, len(methods))
    for row, h2sq in enumerate([0.1, 0.2]):
        for col, n2 in enumerate([100, 400]):
            ax = axes[row, col]
            facet = data[data.h2sq.eq(h2sq) & data.n2.eq(n2)]
            for method, offset in zip(methods, offsets):
                group = facet[facet.method.eq(method)].set_index("propt").loc[PROPORTIONS]
                mean = group["mean"].to_numpy()
                ax.errorbar(np.arange(4) + offset, mean, yerr=np.array([mean - group.ci_low, group.ci_high - mean]), color=COLOR_MAP[method], marker=MARKERS[method], linestyle="--", linewidth=1.6, markersize=5, capsize=2.5, elinewidth=1.1, zorder=3)
            ax.set_title(rf"$h_2^2 = {h2sq:g},\ N_2 = {n2}$", fontsize=12)
            ax.set_xticks(range(4), [f"{p:g}" for p in PROPORTIONS])
            ax.set_xlim(-0.4, 3.4)
            ax.set_ylim(ymin, ymax)
            ax.set_xlabel(r"$\pi$")
            ax.set_ylabel(("Type I error" if metric.startswith("alpha") else "Power") if col == 0 else "")
            if metric.startswith("alpha"):
                ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=1.1, zorder=2)
                ticks = [0, 0.05, 0.1, 0.2, 0.3, 0.4] if metric == "alpha_null" else [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
                ax.set_yticks(ticks)
            else:
                ticks = [0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38] if metric == "power_shared" else [0.2, 0.25, 0.3, 0.35, 0.4]
                ax.set_yticks(ticks)
    handles = [Line2D([0], [0], color=COLOR_MAP[m], marker=MARKERS[m], linestyle="--", linewidth=1.6, label=METHOD_LABELS[m]) for m in methods]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=len(methods), frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{metric}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-points", required=True, type=Path)
    parser.add_argument("--reference-stats", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--metrics", nargs="+", choices=METRIC_ORDER, default=METRIC_ORDER)
    args = parser.parse_args()
    sns.set_theme(style="darkgrid", font_scale=1.0, rc={"pdf.fonttype": 42, "ps.fonttype": 42})
    summary = summarize(pd.read_csv(args.gene_points, sep="\t"))
    reference = pd.read_csv(args.reference_stats, sep="\t")
    check = summary.merge(reference, on=KEYS, suffixes=("", "_reference"), validate="one_to_one")
    assert len(check) == len(summary) == len(reference)
    for column in ["n", "mean", "ci_low", "ci_high"]:
        np.testing.assert_allclose(check[column], check[f"{column}_reference"], rtol=0, atol=1e-12)
    summary["source"] = "archived_gene_values_5f72d47"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary[summary.metric.isin(args.metrics)].to_csv(args.out_dir / "plotted_statistics.tsv", sep="\t", index=False)
    for metric in args.metrics:
        plot(summary, metric, args.out_dir)
        print(f"Saved {args.out_dir / (metric + '.pdf')}")


if __name__ == "__main__":
    main()
