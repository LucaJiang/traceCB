"""Aggregate and plot small-window traceCB simulation replicates.

This is the paired visualizer for ``simulation.py``. It reads
``simulation_<rep>.csv`` files from ``<base_path>/<runname>/<setting>/``,
computes power or type I error for each method, writes
``<base_path>/<runname>/result_df.csv``, and saves figures in
``<base_path>/img``.

The production plotting commands live in ``src/simulation/run_simulation.sh``
immediately after the simulation commands that generate each runname.
"""

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
    "gmmproptmode",
    "gmmproptmodescale",
    "gmmproptnormalvar",
    "gmmpropt",
    "pcausal",
    "causaloverlap",
    "causalmaxabscor",
    "causalpartition",
    "nullregionprop",
    "omega",
]

methods_power = ["sumstat", "cross", "tissue"]
# methods_power = ["sumstat", "cross", "tissue", "meta", "metatissue"]
methods_alpha = ["sumstat", "cross", "tissue", "meta", "metatissue"]
GMM_PROPT_MODE_ORDER = [
    "exact",
    "underestimate",
    "normal",
    "overestimate",
]
GMM_PROPT_MODE_LABELS = {
    "exact": "exact",
    "underestimate": "underestimate",
    "overestimate": "overestimate",
    "normal": "normal",
}

GMM_PROPT_PANEL_ORDER = [
    "exact",
    "underestimate (0.1)",
    "underestimate (0.2)",
    "normal (var=0.1)",
    "overestimate (0.1)",
    "overestimate (0.2)",
]


def build_gmm_propt_panel(mode, scale, normal_var=np.nan):
    if pd.isna(mode) or mode == "exact":
        return GMM_PROPT_MODE_LABELS["exact"]
    if mode in ("underestimate", "overestimate") and pd.notna(scale):
        return f"{GMM_PROPT_MODE_LABELS[mode]} ({float(scale):.1f})"
    if mode == "normal" and pd.notna(normal_var):
        return f"{GMM_PROPT_MODE_LABELS['normal']} (var={float(normal_var):.3g})"
    return GMM_PROPT_MODE_LABELS.get(mode, str(mode))


def get_gmm_panel_order(df):
    if "gmmproptpanel" in df.columns:
        present = set(df["gmmproptpanel"].dropna().astype(str))
        ordered = [panel for panel in GMM_PROPT_PANEL_ORDER if panel in present]
        extras = sorted(present - set(ordered))
        if ordered:
            return ordered + extras

    order = []
    for mode in GMM_PROPT_MODE_ORDER:
        if mode == "exact":
            label = build_gmm_propt_panel(mode, np.nan)
            if label in df["gmmproptpanel"].values:
                order.append(label)
            continue
        if mode in ("underestimate", "overestimate"):
            scales = (
                df.loc[df["gmmproptmode"] == mode, "gmmproptmodescale"]
                .dropna()
                .astype(float)
                .sort_values()
                .unique()
            )
            for scale in scales:
                label = build_gmm_propt_panel(mode, scale)
                if label in df["gmmproptpanel"].values:
                    order.append(label)
            continue
        if mode == "normal" and "gmmproptnormalvar" in df.columns:
            nvars = (
                df.loc[df["gmmproptmode"] == mode, "gmmproptnormalvar"]
                .dropna()
                .astype(float)
                .sort_values()
                .unique()
            )
            for nv in nvars:
                label = build_gmm_propt_panel(mode, np.nan, nv)
                if label in df["gmmproptpanel"].values:
                    order.append(label)
            label = build_gmm_propt_panel(mode, np.nan, np.nan)
            if label in df["gmmproptpanel"].values:
                order.append(label)
            continue
        label = build_gmm_propt_panel(mode, np.nan)
        if label in df["gmmproptpanel"].values:
            order.append(label)
    return order


def parse_param_value(raw_value):
    if raw_value.lower() == "true":
        return True
    if raw_value.lower() == "false":
        return False
    try:
        return float(raw_value)
    except ValueError:
        return raw_value


def parse_result_folder(folder_name):
    """Parse parameter values from a simulation setting folder name."""
    param_set = set(result_param_names)
    tokens = folder_name.split("_")
    params = {}
    i = 0
    while i < len(tokens):
        if tokens[i] not in param_set:
            i += 1
            continue
        param_name = tokens[i]
        j = i + 1
        value_tokens = []
        while j < len(tokens) and tokens[j] not in param_set:
            value_tokens.append(tokens[j])
            j += 1
        if value_tokens:
            params[param_name] = parse_param_value("_".join(value_tokens))
        i = j
    return params


def get_result_table(result_path, eval_method, target_id=1):
    if eval_method == "power":
        methods = methods_power
    elif eval_method in ("alpha", "alpha_a"):
        methods = methods_alpha
    else:
        raise ValueError("eval_method should be alpha, alpha_a, or power")
    result_df = pd.DataFrame(columns=result_param_names + methods)
    all_folder = glob.glob(result_path + "/*")

    for result_name in all_folder:
        if not os.path.isdir(result_name):
            continue
        result_full_path = result_name
        result_row = pd.DataFrame(columns=result_df.columns)

        # get parameter values from folder name
        for param_name, param_value in parse_result_folder(
            os.path.basename(result_name)
        ).items():
            result_row.loc[0, param_name] = param_value

        csv_files = glob.glob(f"{result_full_path}/simulation_*.csv")
        if not csv_files:
            print(f"Warning: in {result_full_path} not found csv files")
            continue

        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            if "gmm_propt" in df.columns:
                result_row.loc[0, "gmmpropt"] = df["gmm_propt"].iloc[0]
            # causal,sign1,sign2,sign_t,z1_sumstat,z1_cross,z1_tissue,z2_sumstat,z2_cross,z2_tissue,zt_sumstat,z_meta,z_metatissue
            if eval_method == "alpha":
                eval_mask = np.ones(df.shape[0], dtype=bool)
                gt = np.zeros(np.sum(eval_mask), dtype=bool)
                metric_method = "alpha"
            elif eval_method == "alpha_a":
                if "region_a" not in df.columns:
                    raise ValueError(
                        f"{csv_file} does not contain region_a for alpha_a."
                    )
                eval_mask = df["region_a"].values == 1
                gt = np.zeros(np.sum(eval_mask), dtype=bool)
                metric_method = "alpha"
            elif eval_method == "power":
                eval_mask = np.ones(df.shape[0], dtype=bool)
                gt = df["causal"].values == 1
                metric_method = "power"
            sumstat_pred = z2p(df[f"z{target_id}_sumstat"].values) < P_THRESHOLD
            cross_pred = z2p(df[f"z{target_id}_cross"].values) < P_THRESHOLD
            tissue_pred = z2p(df[f"z{target_id}_tissue"].values) < P_THRESHOLD
            meta_pred = z2p(df["z_meta"].values) < P_THRESHOLD
            metatissue_pred = z2p(df["z_metatissue"].values) < P_THRESHOLD
            result_row.loc[0, "sumstat"] = calculate_metrics(
                gt, sumstat_pred[eval_mask], metric_method
            )
            result_row.loc[0, "cross"] = calculate_metrics(
                gt, cross_pred[eval_mask], metric_method
            )
            result_row.loc[0, "tissue"] = calculate_metrics(
                gt, tissue_pred[eval_mask], metric_method
            )
            result_row.loc[0, "meta"] = calculate_metrics(
                gt, meta_pred[eval_mask], metric_method
            )
            result_row.loc[0, "metatissue"] = calculate_metrics(
                gt, metatissue_pred[eval_mask], metric_method
            )
            result_df = pd.concat([result_df, result_row], ignore_index=True)
    # if N1, N2, Nt are in columns, convert them to int
    for param in ["n1", "n2", "nt"]:
        if param in result_df.columns:
            result_df[param] = result_df[param].astype(int)
    if "gmmproptmode" in result_df.columns:
        result_df["gmmproptmode"] = result_df["gmmproptmode"].fillna("exact")
        result_df["gmmproptmode"] = pd.Categorical(
            result_df["gmmproptmode"],
            categories=GMM_PROPT_MODE_ORDER,
            ordered=True,
        )
    if "gmmproptmodescale" in result_df.columns:
        result_df["gmmproptmodescale"] = pd.to_numeric(
            result_df["gmmproptmodescale"], errors="coerce"
        )
    if "gmmpropt" in result_df.columns and "propt" in result_df.columns:
        result_df["gmmpropt"] = result_df["gmmpropt"].combine_first(result_df["propt"])
    if "propt" in result_df.columns:
        result_df["propt"] = result_df["propt"].astype(float)
    if "gmmproptnormalvar" in result_df.columns:
        result_df["gmmproptnormalvar"] = pd.to_numeric(
            result_df["gmmproptnormalvar"], errors="coerce"
        )
    if "gmmproptmode" in result_df.columns:
        result_df["gmmproptpanel"] = result_df.apply(
            lambda row: build_gmm_propt_panel(
                row["gmmproptmode"],
                row.get("gmmproptmodescale", np.nan),
                row.get("gmmproptnormalvar", np.nan),
            ),
            axis=1,
        )
        panel_order = get_gmm_panel_order(result_df)
        result_df["gmmproptpanel"] = pd.Categorical(
            result_df["gmmproptpanel"],
            categories=panel_order,
            ordered=True,
        )
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
    "gmmproptmode": r"$\hat{\pi}$ mode",
    "gmmproptpanel": r"$\hat{\pi}$ mode",
    "gmmpropt": r"$\hat{\pi}$",
    "gc": r"$\rho$",
    "n1": r"$N_1$",
    "n2": r"$N_2$",
    "nt": r"$N_t$",
    "nsnp": r"$N_{snp}$",
    "pcausal": r"$p_{causal}$",
    "causaloverlap": r"$p_{shared}$",
    "causalmaxabscor": r"$r_{max}$",
    "causalpartition": "partition",
    "nullregionprop": r"$p_A$",
}
legend_mapping = {
    "sumstat": "Original",
    "cross": "traceC",
    "tissue": "traceCB",
    "meta": "RE2(sc)",
    "metatissue": "RE2(sc+tissue)",
}

marker_map = {
    "sumstat": "o",
    "cross": "D",
    "tissue": "D",
}


def get_marker(method):
    return marker_map.get(method, "s")


def infer_plot_axes(runname):
    tmp = runname.split("_")
    while tmp[-1] in ("nocormax",):
        tmp = tmp[:-1]
    row, col, x = tmp[-3], tmp[-2], tmp[-1]
    if x == "gmmproptmode":
        x = "propt"
        col = "gmmproptpanel"
    return row, col, x


def prepare_plot_df(df, row, col, x, value_name, methods):
    melted_df = pd.melt(
        df.loc[:, [row, col, x] + methods],
        id_vars=[row, col, x],
        value_vars=methods,
        var_name="source",
        value_name=value_name,
    )
    param_mapping = {
        row: all_param_mapping.get(row, row),
        col: all_param_mapping.get(col, col),
        x: all_param_mapping.get(x, x),
    }
    renamed_df = melted_df.rename(columns=param_mapping)
    facet_col = param_mapping.get(col, col)
    return renamed_df, param_mapping


def plot_power_analysis(
    result_df, row="h1sq", col="h2sq", x="propt", ymin=0.0, ymax=0.48, save_name="./"
):
    renamed_df, param_mapping = prepare_plot_df(
        result_df, row, col, x, "power", methods_power
    )

    # Define plot order from bottom to top to ensure desired layering
    plot_order = [
        m
        for m in ["meta", "metatissue", "sumstat", "cross", "tissue"]
        if m in methods_power
    ]

    facet_kwargs = {
        "margin_titles": True,
        "height": 2.1,
        "aspect": 1.2,
        "despine": False,
    }
    if col == "gmmproptpanel":
        g = sns.FacetGrid(
            renamed_df,
            col=param_mapping.get(col, col),
            col_order=get_gmm_panel_order(result_df),
            col_wrap=3,
            **facet_kwargs,
        )
    else:
        g = sns.FacetGrid(
            renamed_df,
            row=param_mapping.get(row, row),
            col=param_mapping.get(col, col),
            **facet_kwargs,
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
        markers=[get_marker(method) for method in plot_order],
        palette=color_map,
    )

    handles, labels = [], []
    for method in methods_power:
        handles.append(
            Line2D(
                [0],
                [0],
                marker=get_marker(method),
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
        ax.set_ylim(ymin, ymax)
        # ax.set_yticks([0, 0.2, 0.4])
        # ax.set_yticklabels(["0", "0.2", "0.4"])
        ax.grid(True, alpha=0.8)
        ax.set_ylabel("Power")
        ax.set_xlabel(param_mapping.get(x, x), labelpad=10)
        ax.set_title(ax.get_title(), pad=12)

    if col == "gmmproptpanel":
        g.set_titles(template="{col_name}")
    else:
        g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.85,
        right=0.95,
        left=0.1,
        bottom=0.1,
        wspace=0.02,
        hspace=0.28,
    )

    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()
    print(f"Plot saved to {save_name}.pdf")


def plot_alpha_analysis(
    result_df, row="pcau", col="h22sq", x="pi", ymin=0.0, ymax=0.3, save_name="./"
):
    renamed_df, param_mapping = prepare_plot_df(
        result_df, row, col, x, "Type 1 error", methods_alpha
    )

    # Define plot order from bottom to top to ensure desired layering
    plot_order = [
        m
        for m in ["meta", "metatissue", "sumstat", "cross", "tissue"]
        if m in methods_alpha
    ]

    facet_kwargs = {
        "margin_titles": True,
        "height": 2.1,
        "aspect": 1.2,
    }
    if col == "gmmproptpanel":
        g = sns.FacetGrid(
            renamed_df,
            col=param_mapping[col],
            col_order=get_gmm_panel_order(result_df),
            col_wrap=3,
            **facet_kwargs,
        )
    else:
        g = sns.FacetGrid(
            renamed_df,
            row=param_mapping[row],
            col=param_mapping[col],
            **facet_kwargs,
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
        markers=[get_marker(method) for method in plot_order],
        linestyles="--",
        linewidth=1.6,
    )

    for ax in g.axes.flatten():
        ax.set_ylim(ymin, ymax)
        yticks = [
            tick
            for tick in [0, 0.05, 0.10, 0.20, 0.30, 0.40]
            if ymin <= tick <= ymax
        ]
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            [f"{tick:g}" if tick == 0 else f"{tick:.2f}" for tick in yticks]
        )
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel("Type I error")
        ax.set_xlabel(param_mapping.get(x, x), labelpad=10)
        ax.grid(True, alpha=0.8)
        ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=2)
        ax.set_title(ax.get_title(), pad=12)

    handles, labels = [], []
    for method in methods_alpha:
        handles.append(
            Line2D(
                [0],
                [0],
                marker=get_marker(method),
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

    if col == "gmmproptpanel":
        g.set_titles(template="{col_name}")
    else:
        g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.85,
        right=0.95,
        left=0.1,
        bottom=0.1,
        wspace=0.02,
        hspace=0.28,
    )
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()
    print(f"Plot saved to {save_name}.pdf")


def add_parser(parser):
    parser.add_argument(
        "--base_path",
        "-b",
        type=str,
        default="./bench/result/",
        help="base path for the simulation results",
    )
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
        help="type of metric to evaluate, alpha, alpha_a, or power",
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
    parser.add_argument(
        "--ymin",
        type=float,
        default=0.0,
        help="ymin for the plot",
    )
    parser.add_argument(
        "--row",
        type=str,
        default=None,
        help="Facet row variable. Defaults to values inferred from --runname.",
    )
    parser.add_argument(
        "--col",
        type=str,
        default=None,
        help="Facet column variable. Defaults to values inferred from --runname.",
    )
    parser.add_argument(
        "--x",
        type=str,
        default=None,
        help="X-axis variable. Defaults to values inferred from --runname.",
    )
    parser.add_argument(
        "--omega",
        type=str,
        choices=["all", "true", "false"],
        default="all",
        help="Filter by omega mode before plotting.",
    )
    parser.add_argument(
        "--save_suffix",
        type=str,
        default=None,
        help="Optional suffix appended to output figure and result table names.",
    )
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = add_parser(parser)
    args = parser.parse_args()
    base_path = args.base_path
    metric = args.metric
    target = args.target
    runname = args.runname
    ymax = args.ymax
    ymin = args.ymin

    result_path = os.path.join(base_path, runname)
    if not os.path.isdir(result_path):
        raise FileNotFoundError(
            f"Simulation result directory does not exist: {result_path}"
        )
    result_df = get_result_table(result_path, metric, target)
    omega_filter = args.omega.lower()
    if omega_filter != "all":
        if "omega" not in result_df.columns:
            raise ValueError("--omega was provided, but result table has no omega column.")
        result_df = result_df[
            result_df["omega"].astype(str).str.lower() == omega_filter
        ].copy()
        if result_df.empty:
            raise ValueError(f"No rows left after --omega {omega_filter} filter.")
    save_suffix = args.save_suffix
    if save_suffix is None:
        if omega_filter == "true":
            save_suffix = "_trueomega"
        elif omega_filter == "false":
            save_suffix = "_estomega"
        else:
            save_suffix = ""
    result_table_name = f"result_df{save_suffix}.csv"
    result_df.to_csv(os.path.join(result_path, result_table_name), index=False)

    if metric == "power":
        plot_func = plot_power_analysis
    elif metric in ("alpha", "alpha_a"):
        plot_func = plot_alpha_analysis
    else:
        raise ValueError("metric should be power, alpha, or alpha_a")
    os.makedirs(os.path.join(base_path, "img"), exist_ok=True)
    inferred_row, inferred_col, inferred_x = infer_plot_axes(runname)
    row = args.row or inferred_row
    col = args.col or inferred_col
    x = args.x or inferred_x
    plot_func(
        result_df=result_df,
        row=row,
        col=col,
        x=x,
        ymin=ymin,
        ymax=ymax,
        save_name=os.path.join(base_path, "img", f"{runname}{save_suffix}"),
    )
