"""Aggregate and plot small-window traceCB simulation replicates.

This is the paired visualizer for ``simulation.py``. It reads
``simulation_<rep>.csv`` files from ``<base_path>/<runname>/<setting>/``,
computes power or type I error for each method, writes
``<base_path>/<runname>/result_df.csv``, and saves figures in ``--img_dir``
or ``$IMG_DIR``.

The production plotting commands live in ``src/simulation/run_main.sh``
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
# python src/simulation/plot_results.py --metric power --runname nt_n2_propt --ymax 0.55
# python src/simulation/plot_results.py --metric power --runname h2sq_gc_propt --ymax 0.48
# python src/simulation/plot_results.py --metric power --runname n1_pcausal_propt --ymax 0.38
# python src/simulation/plot_results.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.28

## true omega
# python src/simulation/plot_results.py --metric power --runname nt_n2_propt --ymax 0.87
# python src/simulation/plot_results.py --metric power --runname h2sq_gc_propt --ymax 0.62
# python src/simulation/plot_results.py --metric power --runname n1_pcausal_propt --ymax 0.44
# python src/simulation/plot_results.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.28

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
    "nt1",
    "nt2",
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
            tracecb2_columns = {
                "z1_tracecb2",
                "z2_tracecb2",
                "zt1_sumstat",
                "zt2_sumstat",
            }
            if tracecb2_columns.intersection(df.columns):
                raise ValueError(
                    f"{csv_file} looks like traceCB^2 output. "
                    "Rename that run with a tracecb2_ prefix and use "
                    "src/simulation/experiments/plot_tracecb2.py for plotting."
                )
            required_columns = [
                "causal",
                f"z{target_id}_sumstat",
                f"z{target_id}_cross",
                f"z{target_id}_tissue",
                "z_meta",
                "z_metatissue",
            ]
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                raise ValueError(
                    f"{csv_file} is missing required simulation.py columns: "
                    + ", ".join(missing_columns)
                )
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
    # Keep integer-valued parameters nullable because some result grids omit
    # one or more sample-size dimensions.
    for param in ["n1", "n2", "nt", "nt1", "nt2"]:
        if param in result_df.columns:
            result_df[param] = pd.to_numeric(
                result_df[param], errors="raise"
            ).astype("Int64")
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
    "h2sq_nullregionprop": r"$h_2^2$, $|A|$",
    "propt": r"$\pi$",
    "gmmproptmode": r"$\hat{\pi}$ mode",
    "gmmproptpanel": r"$\hat{\pi}$ mode",
    "gmmpropt": r"$\hat{\pi}$",
    "gc": r"$\rho$",
    "n1": r"$N_1$",
    "n2": r"$N_2$",
    "nt": r"$N_t$",
    "nt1": r"$N_{t1}$",
    "nt2": r"$N_{t2}$",
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


def format_plot_value(value):
    if pd.isna(value):
        return "NA"
    if isinstance(value, (float, np.floating, int, np.integer)):
        return f"{value:g}"
    return str(value)


def add_plot_derived_columns(df, variables):
    df = df.copy()
    if "h2sq_nullregionprop" in variables:
        required_cols = ["h2sq", "nullregionprop"]
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            raise ValueError(
                "h2sq_nullregionprop requires columns: "
                + ", ".join(required_cols)
                + f". Missing: {', '.join(missing)}"
            )

        def build_row_label(row):
            h2sq = format_plot_value(row["h2sq"])
            nullregionprop = format_plot_value(row["nullregionprop"])
            return (
                rf"$h_2^2 = {h2sq}$"
                + "\n"
                + rf"$p_{{\mathcal{{A}}}} = {nullregionprop}$"
            )

        order_df = (
            df.loc[:, required_cols]
            .drop_duplicates()
            .assign(
                _nullregionprop_order=lambda x: pd.to_numeric(
                    x["nullregionprop"], errors="coerce"
                ),
                _h2sq_order=lambda x: pd.to_numeric(x["h2sq"], errors="coerce"),
            )
            .sort_values(
                ["_nullregionprop_order", "_h2sq_order", "nullregionprop", "h2sq"]
            )
        )
        categories = [build_row_label(row) for _, row in order_df.iterrows()]
        df["h2sq_nullregionprop"] = df.apply(build_row_label, axis=1)
        df["h2sq_nullregionprop"] = pd.Categorical(
            df["h2sq_nullregionprop"],
            categories=categories,
            ordered=True,
        )
    return df


def get_facet_order(df, column):
    series = df[column]
    if isinstance(series.dtype, pd.CategoricalDtype):
        present = set(series.dropna())
        return [value for value in series.cat.categories if value in present]

    values = list(series.dropna().unique())
    numeric_values = pd.to_numeric(pd.Series(values), errors="coerce")
    if numeric_values.notna().all():
        return [
            value
            for _, value in sorted(
                zip(numeric_values.astype(float).tolist(), values),
                key=lambda item: item[0],
            )
        ]
    return sorted(values, key=str)


def prepare_plot_df(df, row, col, x, value_name, methods):
    df = add_plot_derived_columns(df, [row, col, x])
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
    facet_orders = {
        "row": get_facet_order(renamed_df, param_mapping.get(row, row)),
        "col": get_facet_order(renamed_df, param_mapping.get(col, col)),
    }
    return renamed_df, param_mapping, facet_orders


def apply_plot_filters(df, filter_specs):
    for filter_spec in filter_specs:
        if "=" not in filter_spec:
            raise ValueError(
                f"Invalid --plot_filter {filter_spec!r}; expected column=value1,value2."
            )
        column, raw_values = filter_spec.split("=", 1)
        column = column.strip()
        values = [
            parse_param_value(value.strip())
            for value in raw_values.split(",")
            if value.strip()
        ]
        if not column or not values:
            raise ValueError(
                f"Invalid --plot_filter {filter_spec!r}; expected column=value1,value2."
            )
        if column not in df.columns:
            raise ValueError(f"--plot_filter column not found in result table: {column}")

        numeric_values = []
        for value in values:
            try:
                numeric_values.append(float(value))
            except (TypeError, ValueError):
                numeric_values = []
                break

        if numeric_values:
            series = pd.to_numeric(df[column], errors="coerce")
            mask = np.logical_or.reduce(
                [
                    np.isclose(series.astype(float), value, rtol=0, atol=1e-12)
                    for value in numeric_values
                ]
            )
        else:
            mask = df[column].astype(str).isin([str(value) for value in values])

        df = df.loc[mask].copy()
        if df.empty:
            raise ValueError(f"No rows left after --plot_filter {filter_spec!r}.")
    return df


def plot_power_analysis(
    result_df, row="h1sq", col="h2sq", x="propt", ymin=0.0, ymax=0.48, save_name="./"
):
    renamed_df, param_mapping, facet_orders = prepare_plot_df(
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
            row_order=facet_orders["row"],
            col_order=facet_orders["col"],
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
    elif row == "h2sq_nullregionprop":
        g.set_titles(row_template="{row_name}", col_template="{col_var} = {col_name}")
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
    renamed_df, param_mapping, facet_orders = prepare_plot_df(
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
            row_order=facet_orders["row"],
            col_order=facet_orders["col"],
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
    elif row == "h2sq_nullregionprop":
        g.set_titles(row_template="{row_name}", col_template="{col_var} = {col_name}")
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
    parser.add_argument(
        "--plot_filter",
        action="append",
        default=[],
        help=(
            "Optional result-table filter as column=value1,value2. "
            "Can be repeated."
        ),
    )
    parser.add_argument(
        "--img_dir",
        type=str,
        default=os.environ.get("IMG_DIR"),
        help="Directory for figure PDFs. Defaults to $IMG_DIR or <base_path>/img.",
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
    result_df = apply_plot_filters(result_df, args.plot_filter)
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
    img_dir = args.img_dir or os.path.join(base_path, "img")
    os.makedirs(img_dir, exist_ok=True)
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
        save_name=os.path.join(img_dir, f"{runname}{save_suffix}"),
    )
