"""Plot mashr power and realized sign errors using saved true betas."""

import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PARAM_NAMES = [
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

METHOD_SPECS = {
    "mashr_sc": ("mashr_sc", "mashr(sc)"),
    "mashr_sc_bulk": ("mashr_sc_bulk", "mashr(sc+bulk)"),
}

COLOR_MAP = {
    "mashr(sc)": "#d62828",
    "mashr(sc+bulk)": "#6a4c93",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot mashr true FSR curves.")
    parser.add_argument("--base_path", default="bench/result_mashr")
    parser.add_argument("--runname", default="nt_n2_propt_mashr")
    parser.add_argument("--n2", type=int, default=400)
    parser.add_argument("--nt", type=int, default=5000)
    parser.add_argument("--propt", type=float, default=0.4)
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[i / 100 for i in range(1, 16)],
    )
    parser.add_argument("--save_suffix", default="_mashr_true_fsr_curve")
    return parser.parse_args()


def parse_setting(folder):
    tokens = os.path.basename(folder).split("_")
    row = {}
    i = 0
    while i < len(tokens):
        if tokens[i] in PARAM_NAMES and i + 1 < len(tokens):
            key = tokens[i]
            raw_value = tokens[i + 1]
            if raw_value.lower() == "true":
                value = True
            elif raw_value.lower() == "false":
                value = False
            else:
                value = float(raw_value)
            row[key] = value
            i += 2
        else:
            i += 1
    return row


def setting_matches(params, n2, nt, propt):
    return (
        int(params.get("n2", -1)) == n2
        and int(params.get("nt", -1)) == nt
        and abs(float(params.get("propt", -1.0)) - propt) < 1e-12
    )


def collect_setting(result_path, n2, nt, propt):
    matched = []
    for setting_dir in sorted(glob.glob(os.path.join(result_path, "*"))):
        if not os.path.isdir(setting_dir):
            continue
        params = parse_setting(setting_dir)
        if setting_matches(params, n2, nt, propt):
            matched.append((setting_dir, params))
    if not matched:
        raise FileNotFoundError(
            f"No setting found for n2={n2}, nt={nt}, propt={propt} under {result_path}"
        )
    if len(matched) > 1:
        raise ValueError(f"Expected one matching setting, found {len(matched)}")
    return matched[0]


def compute_curve(setting_dir, params, thresholds):
    rows = []
    sim_files = sorted(glob.glob(os.path.join(setting_dir, "simulation_*.csv")))
    for _, (prefix, label) in METHOD_SPECS.items():
        paired = []
        for sim_file in sim_files:
            rep_id = os.path.splitext(os.path.basename(sim_file))[0].split("_")[-1]
            mash_file = os.path.join(setting_dir, f"{prefix}_{rep_id}.csv")
            if os.path.exists(mash_file):
                paired.append((sim_file, mash_file))
        if not paired:
            raise FileNotFoundError(f"No paired files found for {prefix}")

        total_nonzero = 0
        threshold_counts = {
            threshold: {
                "tp": 0,
                "called_nonzero": 0,
                "sign_errors": 0,
                "calls": 0,
                "sum_lfsr": 0.0,
            }
            for threshold in thresholds
        }
        for sim_file, mash_file in paired:
            sim_df = pd.read_csv(sim_file, usecols=["beta1_true"])
            mash_df = pd.read_csv(
                mash_file, usecols=["mash_lfsr_pop1sc", "z1_mashr_pm"]
            )
            beta = sim_df["beta1_true"].to_numpy()
            nonzero = beta != 0
            lfsr = mash_df["mash_lfsr_pop1sc"].to_numpy()
            estimate = mash_df["z1_mashr_pm"].to_numpy()
            total_nonzero += int(nonzero.sum())
            for threshold in thresholds:
                pred = lfsr < threshold
                called_nonzero = pred & nonzero
                sign_error = called_nonzero & (np.sign(estimate) != np.sign(beta))
                threshold_counts[threshold]["tp"] += int(called_nonzero.sum())
                threshold_counts[threshold]["called_nonzero"] += int(
                    called_nonzero.sum()
                )
                threshold_counts[threshold]["sign_errors"] += int(sign_error.sum())
                threshold_counts[threshold]["calls"] += int(pred.sum())
                threshold_counts[threshold]["sum_lfsr"] += float(lfsr[pred].sum())

        for threshold in thresholds:
            counts = threshold_counts[threshold]
            called_nonzero = counts["called_nonzero"]
            calls = counts["calls"]
            rows.append(
                {
                    **params,
                    "method": label,
                    "threshold": threshold,
                    "power": counts["tp"] / total_nonzero if total_nonzero else 0.0,
                    "actual_fsr_nonzero": (
                        counts["sign_errors"] / called_nonzero
                        if called_nonzero
                        else 0.0
                    ),
                    "mean_lfsr": counts["sum_lfsr"] / calls if calls else 0.0,
                    "discoveries": calls,
                    "discoveries_per_rep": calls / len(paired),
                    "called_nonzero": called_nonzero,
                    "sign_errors": counts["sign_errors"],
                    "n_reps": len(paired),
                    "total_nonzero": total_nonzero,
                }
            )
    return pd.DataFrame(rows)


def plot_curve(result_df, save_name):
    plot_df = result_df.melt(
        id_vars=["method", "threshold"],
        value_vars=["power", "actual_fsr_nonzero"],
        var_name="metric",
        value_name="value",
    )
    plot_df["metric"] = plot_df["metric"].map(
        {
            "power": "Power",
            "actual_fsr_nonzero": "Actual sign error rate",
        }
    )

    sns.set_theme(style="darkgrid", palette="muted", color_codes=True)
    g = sns.FacetGrid(
        plot_df,
        col="metric",
        hue="method",
        hue_order=["mashr(sc)", "mashr(sc+bulk)"],
        palette=COLOR_MAP,
        height=2.8,
        aspect=1.25,
        sharey=False,
    )
    g.map_dataframe(
        sns.lineplot,
        x="threshold",
        y="value",
        marker="o",
        linewidth=2.0,
    )
    thresholds = sorted(result_df["threshold"].unique())
    for ax, metric in zip(
        g.axes.flat,
        ["Power", "Actual sign error rate"],
    ):
        ax.set_xlabel("lfsr threshold")
        ax.set_ylabel(metric)
        ax.set_xticks(thresholds)
        ax.set_xlim(0, max(thresholds))
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.set_ylim(0, max(0.55, ax.get_ylim()[1]))
        if metric != "Power":
            ax.plot(
                [0, max(thresholds)],
                [0, max(thresholds)],
                color="#444444",
                linestyle="--",
                linewidth=1.5,
            )
    g.add_legend(title="")
    g.set_titles("{col_name}")
    g.figure.subplots_adjust(top=0.78, wspace=0.28)
    setting = result_df.iloc[0]
    title = (
        f"n2={int(setting['n2'])}, nt={int(setting['nt'])}, "
        f"propt={setting['propt']}, reps={int(setting['n_reps'])}"
    )
    g.figure.suptitle(title, fontsize=13)
    g.savefig(f"{save_name}.pdf", bbox_inches="tight")
    g.savefig(f"{save_name}.png", dpi=220, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to {save_name}.pdf")


def main():
    args = parse_args()
    result_path = os.path.join(args.base_path, args.runname)
    setting_dir, params = collect_setting(result_path, args.n2, args.nt, args.propt)
    result_df = compute_curve(setting_dir, params, args.thresholds)
    os.makedirs(os.path.join(args.base_path, "img"), exist_ok=True)
    tag = f"n2_{args.n2}_nt_{args.nt}_propt_{args.propt:g}"
    result_csv = os.path.join(result_path, f"result_df{args.save_suffix}_{tag}.csv")
    result_df.to_csv(result_csv, index=False)
    save_name = os.path.join(args.base_path, "img", f"{args.runname}{args.save_suffix}_{tag}")
    plot_curve(result_df, save_name)
    print(f"Result table saved to {result_csv}")


if __name__ == "__main__":
    main()
