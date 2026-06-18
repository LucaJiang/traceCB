"""Aggregate and plot traceCB^2 simulation replicates.

This is the paired visualizer for ``simulation_tracecb2.py``. It reads
replicate CSV files from ``<base_path>/<runname>/<setting>/``, writes
``result_df_tracecb2.csv`` or ``result_df_tracecb2_trueomega.csv`` in the
runname directory, and saves paper figures under ``--img_dir`` or ``$IMG_DIR``.

The curated plotting commands are at the end of ``run_tracecb2.sh``.
"""

import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import norm


P_THRESHOLD = 0.05
z2p = lambda z: 2 * norm.sf(abs(z))

PARAM_NAMES = [
    "h1sq",
    "h2sq",
    "gc",
    "n1",
    "n2",
    "nt1",
    "nt2",
    "nsnp",
    "propt",
    "pcausal",
    "omega",
]
METHODS = ["sumstat", "cross", "tissue", "tracecb2"]
METHOD_LABELS = {
    "sumstat": "Original",
    "cross": "traceC",
    "tissue": "traceCB",
    "tracecb2": r"traceCB$^2$",
}
COLORS = {
    "sumstat": "#47D45A",
    "cross": "#fb8500",
    "tissue": "#2d00f7",
    "tracecb2": "#00a6a6",
}
MARKERS = {
    "sumstat": "o",
    "cross": "D",
    "tissue": "D",
    "tracecb2": "^",
}
PARAM_LABELS = {
    "h1sq": r"$h_1^2$",
    "h2sq": r"$h_2^2$",
    "gc": r"$\rho$",
    "n1": r"$N_1$",
    "n2": r"$N_2$",
    "nt1": r"$N_{t1}$",
    "nt2": r"$N_{t2}$",
    "nsnp": r"$N_{snp}$",
    "propt": r"$\pi$",
    "pcausal": r"$p_{causal}$",
}

sns.set_theme(style="darkgrid", palette="muted", color_codes=True)


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


def calculate_metric(gt, pred, metric):
    tp = np.sum(gt & pred)
    tn = np.sum(~gt & ~pred)
    fp = np.sum(~gt & pred)
    fn = np.sum(gt & ~pred)
    if metric == "alpha":
        return fp / (tn + fp + 1e-12)
    if metric == "power":
        return tp / (tp + fn + 1e-12)
    raise ValueError("metric must be alpha or power")


def parse_folder(folder_name):
    parts = folder_name.split("_")
    row = {}
    i = 0
    while i < len(parts):
        if parts[i] in PARAM_NAMES and i + 1 < len(parts):
            key = parts[i]
            val = parts[i + 1]
            if val.lower() == "true":
                row[key] = True
            elif val.lower() == "false":
                row[key] = False
            else:
                try:
                    row[key] = float(val)
                except ValueError:
                    row[key] = val
            i += 2
        else:
            i += 1
    return row


def get_result_table(result_path, metric, target_id=1):
    rows = []
    for result_dir in glob.glob(os.path.join(result_path, "*")):
        if not os.path.isdir(result_dir):
            continue
        base = parse_folder(os.path.basename(result_dir))
        for csv_file in glob.glob(os.path.join(result_dir, "simulation_*.csv")):
            df = pd.read_csv(csv_file)
            required_cols = [f"z{target_id}_{method}" for method in METHODS]
            if not set(required_cols).issubset(df.columns):
                continue
            if metric == "alpha":
                gt = np.zeros(df.shape[0], dtype=bool)
            elif metric == "power":
                gt = df["causal"].values == 1
            else:
                raise ValueError("metric must be alpha or power")
            row = dict(base)
            for method in METHODS:
                pred = z2p(df[f"z{target_id}_{method}"].values) < P_THRESHOLD
                row[method] = calculate_metric(gt, pred, metric)
            rows.append(row)
    result_df = pd.DataFrame(rows)
    if result_df.empty:
        raise ValueError(f"No traceCB^2 simulation CSV files found in {result_path}")
    for col in ("n1", "n2", "nt1", "nt2"):
        if col in result_df.columns:
            result_df[col] = result_df[col].astype(int)
    if "propt" in result_df.columns:
        result_df["propt"] = result_df["propt"].astype(float)
    return result_df


def infer_plot_axes(runname):
    if runname.endswith("nt1_nt2_propt"):
        return "nt2", "nt1", "propt"
    parts = runname.split("_")
    if len(parts) < 3:
        raise ValueError("runname must end with row_col_x, e.g. nt1_nt2_propt")
    return parts[-3], parts[-2], parts[-1]


def auto_ylim(df, metric, row, col, x, ymin=None, ymax=None):
    # The plot shows point estimates, so base the automatic axis on panel means
    # instead of raw replicate extrema. Raw power replicates often include 0/1.
    values = (
        df.groupby([row, col, x], observed=True)[METHODS]
        .mean()
        .astype(float)
        .to_numpy()
        .ravel()
    )
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    low = float(np.min(values))
    high = float(np.max(values))
    if metric == "alpha":
        low = min(low, 0.05)
        high = max(high, 0.05)
        pad = max(0.006, (high - low) * 0.35)
    else:
        pad = max(0.02, (high - low) * 0.18)
    y0 = max(0.0, low - pad) if ymin is None else ymin
    y1 = min(1.0, high + pad) if ymax is None else ymax
    if y1 <= y0:
        y1 = y0 + 0.1
    return y0, y1


def plot(
    result_df,
    runname,
    metric,
    ymin=None,
    ymax=None,
    base_path="bench/result",
    suffix="",
    img_dir=None,
):
    row, col, x = infer_plot_axes(runname)
    y_label = "Type I error" if metric == "alpha" else "Power"
    y0, y1 = auto_ylim(result_df, metric, row, col, x, ymin, ymax)
    melted = pd.melt(
        result_df[[row, col, x] + METHODS],
        id_vars=[row, col, x],
        value_vars=METHODS,
        var_name="source",
        value_name=y_label,
    )
    label_map = {
        row: PARAM_LABELS.get(row, row),
        col: PARAM_LABELS.get(col, col),
        x: PARAM_LABELS.get(x, x),
    }
    melted = melted.rename(columns=label_map)

    g = sns.FacetGrid(
        melted,
        row=label_map[row],
        col=label_map[col],
        margin_titles=True,
        height=2.15,
        aspect=1.2,
        despine=False,
    )
    g.map_dataframe(
        sns.pointplot,
        x=label_map[x],
        y=y_label,
        hue="source",
        hue_order=METHODS,
        estimator="mean",
        errorbar=("ci", 95),
        dodge=0.45,
        linestyles="--",
        linewidth=1.6,
        markers=[MARKERS[m] for m in METHODS],
        palette=COLORS,
    )
    for ax in g.axes.flatten():
        ax.set_ylim(y0, y1)
        ax.set_ylabel(y_label)
        ax.set_xlabel(label_map[x], labelpad=10)
        ax.grid(True, alpha=0.8)
        if metric == "alpha":
            set_alpha_yticks(ax, y0, y1)
            ax.axhline(0.05, color="#e63946", linestyle="--", linewidth=1.8)
        ax.set_title(ax.get_title(), pad=12)

    handles = [
        Line2D(
            [0],
            [0],
            marker=MARKERS[m],
            color=COLORS[m],
            linestyle="--",
            markersize=8,
            label=METHOD_LABELS[m],
        )
        for m in METHODS
    ]
    g.figure.legend(
        handles=handles,
        title="",
        fontsize=12,
        bbox_to_anchor=(0.53, 0.94),
        loc="center",
        ncol=4,
        frameon=False,
    )
    g.set_titles(template="{row_name} | {col_name}")
    g.figure.subplots_adjust(
        top=0.85,
        right=0.95,
        left=0.1,
        bottom=0.1,
        wspace=0.02,
        hspace=0.28,
    )
    img_dir = img_dir or os.path.join(base_path, "img")
    os.makedirs(img_dir, exist_ok=True)
    save_name = os.path.join(img_dir, f"{runname}_{metric}_tracecb2{suffix}.pdf")
    g.savefig(save_name, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to {save_name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_path", default="bench/result")
    parser.add_argument(
        "--img_dir",
        default=os.environ.get("IMG_DIR"),
        help="Directory for figure PDFs. Defaults to $IMG_DIR or <base_path>/img.",
    )
    parser.add_argument("--runname", required=True)
    parser.add_argument("--metric", choices=("power", "alpha"), required=True)
    parser.add_argument("--target", default=1, type=int)
    parser.add_argument("--ymin", default=None, type=float)
    parser.add_argument("--ymax", default=None, type=float)
    parser.add_argument("--drop_pcausal", default=None, type=float, nargs="+")
    parser.add_argument(
        "--omega",
        choices=("true", "false"),
        default=None,
        help="Filter to true omega or estimated omega result folders.",
    )
    args = parser.parse_args()

    result_path = os.path.join(args.base_path, args.runname)
    if not os.path.isdir(result_path):
        raise FileNotFoundError(result_path)
    result_df = get_result_table(result_path, args.metric, args.target)
    suffix = ""
    if args.omega is not None and "omega" in result_df.columns:
        omega_value = args.omega == "true"
        result_df = result_df[result_df["omega"] == omega_value]
        suffix = "_trueomega" if omega_value else ""
    if args.drop_pcausal is not None and "pcausal" in result_df.columns:
        result_df = result_df[~result_df["pcausal"].isin(args.drop_pcausal)]
    if result_df.empty:
        raise ValueError(f"No rows left after filtering {result_path}")
    result_csv = f"result_df_tracecb2{suffix}.csv"
    result_df.to_csv(os.path.join(result_path, result_csv), index=False)
    plot(
        result_df,
        args.runname,
        args.metric,
        args.ymin,
        args.ymax,
        args.base_path,
        suffix,
        args.img_dir,
    )


if __name__ == "__main__":
    main()
