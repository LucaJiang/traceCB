"""Plot masked-omega simulation comparisons.

This visualizer pairs with ``simulate_masked_omega.py`` and the shell wrapper
``run_masked_omega.sh``. It can combine multiple runnames so one figure can
facet over genetic correlation values while comparing original traceC/traceCB
against masked-input variants.

Outputs are saved to ``<img_dir>`` as PDF.
"""

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

sns.set_theme(style="darkgrid", palette="muted", color_codes=True)

METHODS = [
    "original_pop1",
    "original_pop2",
    "traceC_original",
    "traceCB_original",
    "pop1sc_pop2bulk",
    "pop2sc_pop2bulk",
]

color_map = {
    "original_pop1": "#47D45A",
    "original_pop2": "#0091FF",
    "traceC_original": "#fb8500",
    "traceCB_original": "#2d00f7",
    "pop1sc_pop2bulk": "#fcbf49",
    "pop2sc_pop2bulk": "#41cef1",
}

marker_map = {
    "original_pop1": "o",
    "original_pop2": "o",
    "traceC_original": "D",
    "traceCB_original": "D",
    "pop1sc_pop2bulk": "s",
    "pop2sc_pop2bulk": "s",
}

legend_mapping = {
    "original_pop1": "Original pop1",
    "original_pop2": "Original pop2",
    "traceC_original": "traceC",
    "traceCB_original": "traceCB",
    "pop1sc_pop2bulk": "pop1 sc + bulk",
    "pop2sc_pop2bulk": "pop2 sc + bulk",
}

all_param_mapping = {
    "h1sq": r"$h_1^2$",
    "h2sq": r"$h_2^2$",
    "propt": r"$\pi$",
    "gc": r"$\rho$",
    "n1": r"$N_1$",
    "n2": r"$N_2$",
    "nt": r"$N_t$",
    "nsnp": r"$N_{snp}$",
    "pcausal": r"$p_{causal}$",
}


def set_alpha_yticks(ax, ymin, ymax):
    yticks = [
        tick
        for tick in [0, 0.05, 0.10, 0.20, 0.30, 0.40]
        if ymin <= tick <= ymax
    ]
    ax.set_yticks(yticks)
    ax.set_yticklabels(
        [f"{tick:g}" if tick == 0 else f"{tick:.2f}" for tick in yticks]
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize masked-omega traceC/traceCB simulation results."
    )
    parser.add_argument(
        "--base_path",
        "-b",
        default="bench/result/masked_omega_compare",
        help="Base output path used by compare_tracec_tracecb_masked_omega.py.",
    )
    parser.add_argument("--runname", "-r", nargs="+", default=["masked_omega_compare"])
    parser.add_argument(
        "--save_prefix",
        default=None,
        help="Output filename prefix under --img_dir. Defaults to runname.",
    )
    parser.add_argument(
        "--img_dir",
        default=os.environ.get("IMG_DIR"),
        help="Directory for figure PDFs. Defaults to $IMG_DIR or <base_path>/img.",
    )
    parser.add_argument("--x", default="propt")
    parser.add_argument("--row", default="h1sq")
    parser.add_argument("--col", default="h2sq")
    parser.add_argument("--power_ymax", type=float, default=None)
    parser.add_argument("--alpha_ymax", type=float, default=None)
    parser.add_argument("--ymin", type=float, default=0.0)
    parser.add_argument(
        "--metric",
        choices=["power", "alpha", "both"],
        default="both",
        help="Which metric to plot.",
    )
    return parser.parse_args()


def plot_metric(
    result_df,
    metric,
    row="h1sq",
    col="h2sq",
    x="propt",
    ymin=0.0,
    ymax=None,
    save_name="./",
):
    plot_df = result_df.copy()
    method_order = [
        method for method in METHODS if method in plot_df["method"].unique()
    ]
    plot_df["source"] = pd.Categorical(
        plot_df["method"], categories=method_order, ordered=True
    )
    plot_df = plot_df.sort_values(["source", x])
    param_mapping = {
        row: all_param_mapping.get(row, row),
        col: all_param_mapping.get(col, col),
        x: all_param_mapping.get(x, x),
    }
    plot_df = plot_df.rename(columns=param_mapping)
    row_name = param_mapping.get(row, row)
    col_name = param_mapping.get(col, col)
    x_name = param_mapping.get(x, x)
    y_name = "Power" if metric == "power" else "Type I error"

    g = sns.FacetGrid(
        plot_df,
        row=row_name,
        col=col_name,
        margin_titles=True,
        height=2.8,
        aspect=1.35,
        despine=False,
    )
    g.map_dataframe(
        sns.pointplot,
        x=x_name,
        y=metric,
        hue="source",
        hue_order=method_order,
        estimator="mean",
        errorbar=("ci", 95),
        dodge=0.45,
        linestyles="--",
        linewidth=1.6,
        markers=[marker_map[method] for method in method_order],
        palette=color_map,
    )

    handles, labels = [], []
    for method in method_order:
        handles.append(
            Line2D(
                [0],
                [0],
                marker=marker_map[method],
                color=color_map[method],
                linestyle="--",
                markersize=7,
                label=legend_mapping[method],
            )
        )
        labels.append(legend_mapping[method])

    g.figure.legend(
        handles=handles,
        labels=labels,
        title="",
        title_fontsize=0.1,
        fontsize=11,
        bbox_to_anchor=(0.53, 0.94),
        loc="center",
        ncol=min(len(method_order), 6),
        frameon=False,
    )

    for ax in g.axes.flatten():
        ax.set_ylim(ymin, ymax)
        ax.grid(True, alpha=0.8)
        ax.set_ylabel(y_name)
        if metric == "alpha":
            set_alpha_yticks(ax, ymin, ymax)
            ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=1.6)

    g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.88,
        right=0.95,
        left=0.1,
        bottom=0.12,
        wspace=0.02,
        hspace=0.08,
    )
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()
    print(f"{y_name} plot saved to {save_name}.pdf")


def main():
    args = parse_args()
    result_dfs = []
    for runname in args.runname:
        result_path = os.path.join(args.base_path, runname)
        result_csv = os.path.join(result_path, "replicate_metrics.csv")
        if not os.path.exists(result_csv):
            print(f"Warning: skip missing result file {result_csv}")
            continue
        run_df = pd.read_csv(result_csv)
        run_df["runname"] = runname
        result_dfs.append(run_df)
    if not result_dfs:
        raise FileNotFoundError(
            "No replicate_metrics.csv found for requested runname(s)"
        )
    result_df = pd.concat(result_dfs, ignore_index=True)
    img_dir = args.img_dir or os.path.join(args.base_path, "img")
    os.makedirs(img_dir, exist_ok=True)
    save_prefix = args.save_prefix or "_".join(args.runname)

    if args.metric in ("power", "both"):
        plot_metric(
            result_df,
            "power",
            row=args.row,
            col=args.col,
            x=args.x,
            ymin=args.ymin,
            ymax=args.power_ymax,
            save_name=os.path.join(img_dir, f"{save_prefix}_power"),
        )
    if args.metric in ("alpha", "both"):
        plot_metric(
            result_df,
            "alpha",
            row=args.row,
            col=args.col,
            x=args.x,
            ymin=args.ymin,
            ymax=args.alpha_ymax,
            save_name=os.path.join(img_dir, f"{save_prefix}_alpha"),
        )


if __name__ == "__main__":
    main()
