import numpy as np
import pandas as pd
import os, argparse, glob
from scipy.stats import norm
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

## Usage:
# python src/simulation/visual_simulation.py --metric power --runname nt_n2_propt --ymax 0.55
# python src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.48
# python src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.38
# python src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.28

## true omega
# python src/simulation/visual_simulation.py --metric power --runname nt_n2_propt --ymax 0.87
# python src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.62
# python src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.44
# python src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.28

## convert p-value to z-score and vice versa
p2z = lambda p: np.abs(norm.ppf(p / 2))
z2p = lambda z: 2 * norm.sf(abs(z))
P_THRESHOLD = 0.05
sns.set_theme(style="darkgrid", palette="muted", color_codes=True)

base_path = "/Users/lucajiang/learn/CityU/traceCB/bench/result/"


def calculate_metrics(gt, pred, eval_method):
    tp = np.sum(gt & pred)
    tn = np.sum(~gt & ~pred)
    fp = np.sum(~gt & pred)
    fn = np.sum(gt & ~pred)
    if eval_method == "alpha":
        alpha = fp / (tn + fp + 1e-12)
        return alpha
    elif eval_method == "power" or eval_method == "sen":
        power = tp / (tp + fn + 1e-12)
        return power
    elif eval_method == "spe":
        specificity = tn / (tn + fp + 1e-12)
        return specificity
    raise ValueError("eval_method should be alpha, power, or specificity")


result_param_names = [
    "h1sq",
    "h2sq",
    "gc",
    "n1",
    "n2",
    "nt",
    "nsnp",
    "propt",
    "pcausal",
    "omega",
]

methods_power = ["sumstat", "cross", "tissue", "meta"]
# methods_power = ["sumstat", "cross", "tissue", "meta", "metatissue"]
methods_alpha = ["sumstat", "cross", "tissue", "meta", "metatissue"]


def get_result_table(result_path, eval_method, target_id=1):
    if eval_method == "power":
        methods = methods_power
    elif eval_method == "alpha":
        methods = methods_alpha
    result_df = pd.DataFrame(columns=result_param_names + methods)
    all_folder = glob.glob(result_path + "/*")

    for result_name in all_folder:
        if not os.path.isdir(result_name):
            continue
        result_full_path = result_name
        result_row = pd.DataFrame(columns=result_df.columns)

        # get parameter values from folder name
        folder_name = os.path.basename(result_name)
        param_pairs = folder_name.split("_")
        i = 0

        while i < len(param_pairs):
            if param_pairs[i] in result_param_names:
                param_name = param_pairs[i]
                if i + 1 < len(param_pairs):
                    if param_pairs[i + 1].lower() == "true":
                        param_value = True
                    elif param_pairs[i + 1].lower() == "false":
                        param_value = False
                    else:
                        try:
                            param_value = float(param_pairs[i + 1])
                        except ValueError:
                            param_value = param_pairs[i + 1]
                    result_row.loc[0, param_name] = param_value
                    i += 2
                else:
                    i += 1
            else:
                i += 1

        csv_files = glob.glob(f"{result_full_path}/simulation_*.csv")
        if not csv_files:
            print(f"Warning: in {result_full_path} not found csv files")
            continue

        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            # causal,sign1,sign2,sign_t,z1_sumstat,z1_cross,z1_tissue,z2_sumstat,z2_cross,z2_tissue,zt_sumstat,z_meta,z_metatissue
            if eval_method == "alpha":
                gt = np.zeros_like(df["causal"].values, dtype=bool)
            elif eval_method == "power":
                gt = df["causal"].values == 1
            sumstat_pred = z2p(df[f"z{target_id}_sumstat"].values) < P_THRESHOLD
            cross_pred = z2p(df[f"z{target_id}_cross"].values) < P_THRESHOLD
            tissue_pred = z2p(df[f"z{target_id}_tissue"].values) < P_THRESHOLD
            meta_pred = z2p(df["z_meta"].values) < P_THRESHOLD
            metatissue_pred = z2p(df["z_metatissue"].values) < P_THRESHOLD
            result_row.loc[0, "sumstat"] = calculate_metrics(
                gt, sumstat_pred, eval_method
            )
            result_row.loc[0, "cross"] = calculate_metrics(gt, cross_pred, eval_method)
            result_row.loc[0, "tissue"] = calculate_metrics(
                gt, tissue_pred, eval_method
            )
            result_row.loc[0, "meta"] = calculate_metrics(gt, meta_pred, eval_method)
            result_row.loc[0, "metatissue"] = calculate_metrics(
                gt, metatissue_pred, eval_method
            )
            result_df = pd.concat([result_df, result_row], ignore_index=True)
    # if N1, N2, Nt are in columns, convert them to int
    for param in ["n1", "n2", "nt"]:
        if param in result_df.columns:
            result_df[param] = result_df[param].astype(int)
    return result_df


color_map = {
    "sumstat": "#47D45A",
    "cross": "#fb8500",
    "tissue": "#2d00f7",
    "meta": "#fcbf49",
    "metatissue": "#41cef1",
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
legend_mapping = {
    "sumstat": "Original",
    "cross": "traceC",
    "tissue": "traceCB",
    "meta": "RE2(sc)",
    "metatissue": "RE2(sc+tissue)",
}


def plot_power_analysis(
    result_df, row="h1sq", col="h2sq", x="propt", ymax=0.48, save_name="./"
):
    # convert to pivot table
    melted_df = pd.melt(
        result_df.loc[:, [row, col, x] + methods_power],
        id_vars=[row, col, x],
        value_vars=methods_power,
        var_name="source",
        value_name="power",
    )

    param_mapping = {
        row: all_param_mapping.get(row, row),
        col: all_param_mapping.get(col, col),
        x: all_param_mapping.get(x, x),
    }
    renamed_df = melted_df.rename(columns=param_mapping)

    # Define plot order from bottom to top to ensure desired layering
    plot_order = [
        m
        for m in ["meta", "metatissue", "sumstat", "cross", "tissue"]
        if m in methods_power
    ]

    g = sns.FacetGrid(
        renamed_df,
        row=param_mapping.get(row, row),
        col=param_mapping.get(col, col),
        margin_titles=True,
        height=2.1,
        aspect=1.2,
        despine=False,
    )
    g.map_dataframe(
        sns.pointplot,
        x=param_mapping.get(x, x),
        y="power",
        hue="source",
        hue_order=plot_order,
        estimator="mean",
        # errorbar="se",
        errorbar=("ci", 95),
        dodge=0.4,
        linestyles="--",
        linewidth=1.6,
        markers="D",
        palette=color_map,
    )

    handles, labels = [], []
    for method in methods_power:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="D",
                color=color_map[method],
                linestyle="--",
                markersize=8,
                label=legend_mapping[method],
            )
        )
        labels.append(legend_mapping[method])

    g.figure.legend(
        handles=handles,
        labels=labels,
        title="",
        title_fontsize=0.1,
        fontsize=12,
        bbox_to_anchor=(0.53, 0.94),
        loc="center",
        ncol=5,
        frameon=False,
    )

    for ax in g.axes.flatten():
        ax.set_ylim(0, ymax)
        # ax.set_yticks([0, 0.2, 0.4])
        # ax.set_yticklabels(["0", "0.2", "0.4"])
        ax.grid(True, alpha=0.8)
        ax.set_ylabel("Power")

    g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.85,
        right=0.95,
        left=0.1,
        bottom=0.1,
        wspace=0.02,
        hspace=0.08,
    )

    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()
    print(f"Power analysis plot saved to {save_name}.pdf")


def plot_alpha_analysis(
    result_df, row="pcau", col="h22sq", x="pi", ymax=0.3, save_name="./"
):
    # convert to pivot table
    melted_df = pd.melt(
        result_df.loc[:, [row, col, x] + methods_alpha],
        id_vars=[row, col, x],
        value_vars=methods_alpha,
        var_name="source",
        value_name="Type 1 error",
    )

    param_mapping = {
        row: all_param_mapping.get(row, row),
        col: all_param_mapping.get(col, col),
        x: all_param_mapping.get(x, x),
    }

    # Define plot order from bottom to top to ensure desired layering
    plot_order = [
        m
        for m in ["meta", "metatissue", "sumstat", "cross", "tissue"]
        if m in methods_alpha
    ]

    g = sns.FacetGrid(
        melted_df.rename(columns=param_mapping),
        row=param_mapping[row],
        col=param_mapping[col],
        margin_titles=True,
        height=2.1,
        aspect=1.2,
    )
    g.map_dataframe(
        sns.pointplot,
        x=param_mapping[x],
        y="Type 1 error",
        hue="source",
        hue_order=plot_order,
        # order=sorted(melted_df[x].unique()),
        # hue_order=methods,
        palette=color_map,
        estimator="mean",
        # errorbar="se",
        errorbar=("ci", 95),
        dodge=0.4,
        markers="D",
        linestyles="--",
        linewidth=1.6,
    )

    for ax in g.axes.flatten():
        ax.set_ylim(0, ymax)
        ax.set_yticks([0, 0.05, 0.10, 0.2])
        ax.set_yticklabels(["0", "0.05", "0.1", "0.2"])
        ax.set_ylabel("Type I error")
        ax.grid(True, alpha=0.8)
        ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=2)

    handles, labels = [], []
    for method in methods_alpha:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="D",
                color=color_map[method],
                linestyle="--",
                markersize=8,
                label=legend_mapping[method],
            )
        )
        labels.append(legend_mapping[method])

    g.figure.legend(
        handles=handles,
        labels=labels,
        title="",
        title_fontsize=0.1,
        fontsize=12,
        bbox_to_anchor=(0.53, 0.94),
        loc="center",
        ncol=5,
        frameon=False,
    )

    g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.85,
        right=0.95,
        left=0.1,
        bottom=0.1,
        wspace=0.02,
        hspace=0.08,
    )
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()
    print(f"Power analysis plot saved to {save_name}.pdf")


def add_parser(parser):
    parser.add_argument(
        "--runname",
        "-r",
        type=str,
        default="h1sq_h2sq_gc",
        help="run name for the simulation",
    )
    parser.add_argument(
        "--metric",
        "-m",
        type=str,
        default="power",
        help="type of metric to evaluate, alpha or power",
    )
    parser.add_argument(
        "--target",
        "-t",
        type=int,
        default=1,
        help="target is population 1 or 2, type 1 or 2",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=0.48,
        help="ymax for the plot",
    )
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = add_parser(parser)
    args = parser.parse_args()
    metric = args.metric
    target = args.target
    runname = args.runname
    ymax = args.ymax

    result_path = os.path.join(base_path, runname)
    result_df = get_result_table(result_path, metric, target)
    result_df.to_csv(os.path.join(result_path, "result_df.csv"), index=False)

    if metric == "power":
        plot_func = plot_power_analysis
    elif metric == "alpha":
        plot_func = plot_alpha_analysis
    else:
        raise ValueError("metric should be power or alpha")
    os.makedirs(os.path.join(base_path, "img"), exist_ok=True)
    tmp = runname.split("_")  # h1sq_h2sq_gc -> ['h1sq', 'h2sq', 'gc']
    plot_func(
        result_df=result_df,
        row=tmp[-3],
        col=tmp[-2],
        x=tmp[-1],
        ymax=ymax,
        save_name=os.path.join(base_path, "img", runname),
    )
    print("Visual done")
