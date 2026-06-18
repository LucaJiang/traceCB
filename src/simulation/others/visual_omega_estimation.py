import argparse
import glob
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_theme(style="darkgrid", palette="colorblind", color_codes=True)

PANEL_COLOR_CYCLE = [
    "#1b1b1b",
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
]
GMM_PROPT_MODE_ORDER = ["exact", "underestimate", "overestimate", "normal"]

RESULT_PARAM_NAMES = [
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
    "pcausal",
    "causaloverlap",
    "causalmaxabscor",
    "causalpartition",
    "nullregionprop",
    "omega",
]

COMPONENT_LABELS = {
    "omega11": r"$\omega_{11}$",
    "omega22": r"$\omega_{22}$",
    "omega12": r"$\omega_{12}$",
}


def build_gmm_panel(row):
    mode = row.get("gmmproptmode", "exact")
    if pd.isna(mode) or mode == "exact":
        return "exact"
    if mode in ("underestimate", "overestimate"):
        scale = row.get("gmmproptmodescale")
        if pd.notna(scale):
            return f"{mode} ({float(scale):g})"
    if mode == "normal":
        normal_var = row.get("gmmproptnormalvar")
        if pd.notna(normal_var):
            return f"normal (var={float(normal_var):g})"
    return str(mode)


def get_gmm_panel_order(df):
    order = []
    for mode in GMM_PROPT_MODE_ORDER:
        mode_df = df[df["gmmproptmode"] == mode]
        if mode == "exact":
            if "exact" in df["gmmproptpanel"].values:
                order.append("exact")
            continue
        if mode in ("underestimate", "overestimate"):
            if "gmmproptmodescale" not in mode_df:
                continue
            for scale in sorted(mode_df["gmmproptmodescale"].dropna().astype(float).unique()):
                label = f"{mode} ({scale:g})"
                if label in df["gmmproptpanel"].values:
                    order.append(label)
            continue
        if mode == "normal":
            if "gmmproptnormalvar" not in mode_df:
                continue
            for normal_var in sorted(
                mode_df["gmmproptnormalvar"].dropna().astype(float).unique()
            ):
                label = f"normal (var={normal_var:g})"
                if label in df["gmmproptpanel"].values:
                    order.append(label)
    extras = sorted(set(df["gmmproptpanel"].dropna()) - set(order))
    return order + extras


def get_panel_palette(hue_order):
    colors = PANEL_COLOR_CYCLE
    if len(hue_order) > len(colors):
        colors = sns.color_palette("husl", len(hue_order)).as_hex()
    return dict(zip(hue_order, colors[: len(hue_order)]))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize estimated omega versus realized true omega."
    )
    parser.add_argument("--base_path", "-b", default="bench/result_estOmega_omegaCompare")
    parser.add_argument("--runname", "-r", default="alpha_pcausal_propt_gmmproptmode")
    parser.add_argument("--x", default="propt")
    parser.add_argument(
        "--img_dir",
        default=os.environ.get("IMG_DIR"),
        help="Directory for figure PDFs. Defaults to $IMG_DIR or <base_path>/img.",
    )
    return parser.parse_args()


def parse_result_folder(folder_name):
    values = {}
    parts = os.path.basename(folder_name).split("_")
    i = 0
    while i < len(parts):
        if parts[i] in RESULT_PARAM_NAMES and i + 1 < len(parts):
            key = parts[i]
            raw_value = parts[i + 1]
            if raw_value.lower() == "true":
                value = True
            elif raw_value.lower() == "false":
                value = False
            else:
                try:
                    value = float(raw_value)
                except ValueError:
                    value = raw_value
            values[key] = value
            i += 2
        else:
            i += 1
    return values


def read_omega_summaries(result_path):
    rows = []
    for folder in glob.glob(os.path.join(result_path, "*")):
        if not os.path.isdir(folder):
            continue
        params = parse_result_folder(folder)
        for summary_file in glob.glob(os.path.join(folder, "omega_summary_*.csv")):
            df = pd.read_csv(summary_file)
            for key, value in params.items():
                df[key] = value
            rows.append(df)
    if not rows:
        raise FileNotFoundError(
            f"No omega_summary_*.csv files found under {result_path}"
        )
    result_df = pd.concat(rows, ignore_index=True)
    result_df["gmmproptmode"] = result_df.get(
        "gmmproptmode", pd.Series(index=result_df.index)
    ).fillna("exact")
    result_df["gmmproptpanel"] = result_df.apply(build_gmm_panel, axis=1)
    return result_df


def build_long_df(result_df):
    components = ["omega11", "omega22", "omega12"]
    rows = []
    for component in components:
        true_col = f"{component}_true"
        est_col = f"{component}_est"
        tmp = result_df.copy()
        tmp["component"] = COMPONENT_LABELS[component]
        tmp["true"] = tmp[true_col]
        tmp["estimated"] = tmp[est_col]
        tmp["diff"] = tmp["estimated"] - tmp["true"]
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def plot_diff(long_df, x, save_name):
    hue_order = get_gmm_panel_order(long_df)
    palette = get_panel_palette(hue_order)
    g = sns.FacetGrid(
        long_df,
        col="component",
        col_wrap=2,
        sharey=False,
        height=2.8,
        aspect=1.25,
        despine=False,
    )
    g.map_dataframe(
        sns.pointplot,
        x=x,
        y="diff",
        hue="gmmproptpanel",
        hue_order=hue_order,
        palette=palette,
        estimator="mean",
        errorbar=("ci", 95),
        markers="D",
        linestyles="--",
        linewidth=1.5,
        alpha=0.88,
        err_kws={"alpha": 0.42, "linewidth": 1.1},
    )
    for ax in g.axes.flatten():
        ax.axhline(0, color="#e63946", linestyle="--", linewidth=1.4)
        ax.grid(True, alpha=0.8)
        ax.set_ylabel("Estimated - true")
    g.set_titles(template="{col_name}")
    g.add_legend(title="")
    g.figure.subplots_adjust(top=0.9, hspace=0.3)
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()


def plot_scatter(long_df, save_name):
    hue_order = get_gmm_panel_order(long_df)
    palette = get_panel_palette(hue_order)
    g = sns.FacetGrid(
        long_df,
        col="component",
        col_wrap=2,
        sharex=False,
        sharey=False,
        height=2.8,
        aspect=1.25,
        despine=False,
    )
    g.map_dataframe(
        sns.scatterplot,
        x="true",
        y="estimated",
        hue="gmmproptpanel",
        hue_order=hue_order,
        palette=palette,
        alpha=0.42,
        s=20,
        edgecolor=None,
    )
    for ax in g.axes.flatten():
        low = min(ax.get_xlim()[0], ax.get_ylim()[0])
        high = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([low, high], [low, high], color="#e63946", linestyle="--", linewidth=1.2)
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.grid(True, alpha=0.8)
    g.set_titles(template="{col_name}")
    g.add_legend(title="")
    g.figure.subplots_adjust(top=0.9, hspace=0.3)
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()
    result_path = os.path.join(args.base_path, args.runname)
    result_df = read_omega_summaries(result_path)
    img_dir = args.img_dir or os.path.join(args.base_path, "img")
    os.makedirs(img_dir, exist_ok=True)
    summary_path = os.path.join(result_path, "omega_summary_df.csv")
    result_df.to_csv(summary_path, index=False)
    long_df = build_long_df(result_df)
    plot_diff(
        long_df,
        args.x,
        os.path.join(img_dir, f"{args.runname}_omega_diff"),
    )
    plot_scatter(
        long_df,
        os.path.join(img_dir, f"{args.runname}_omega_est_vs_true"),
    )
    print(f"Omega summary saved to {summary_path}")


if __name__ == "__main__":
    main()
