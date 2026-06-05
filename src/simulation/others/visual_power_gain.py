"""Plot two rebuttal power-gain panels from new true-omega simulations."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SIMULATION_DIR = Path(__file__).resolve().parents[1]
if str(SIMULATION_DIR) not in sys.path:
    sys.path.insert(0, str(SIMULATION_DIR))

from visual_simulation import color_map

METHOD_COLUMNS = ["sumstat", "cross", "tissue"]
CATEGORY_COLORS = [
    color_map["cross"],
    color_map["tissue"],
    color_map["sumstat"],
    color_map["meta"],
    color_map["metatissue"],
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_path",
        type=Path,
        default=Path("bench/result_power_gain"),
        help="Base directory containing simulation result folders.",
    )
    parser.add_argument(
        "--tracec_runname",
        default="power_gain_tracec_n2_rho",
        help="Run name for the traceC-vs-Original N2/rho grid.",
    )
    parser.add_argument(
        "--tracecb_runname",
        default="power_gain_tracecb_propt_nt",
        help="Run name for the traceCB-vs-traceC propt/Nt grid.",
    )
    parser.add_argument(
        "--output_prefix",
        type=Path,
        default=Path("bench/result_power_gain/img/power_gain"),
        help="Output path without extension.",
    )
    return parser.parse_args()


def read_result_df(base_path, runname):
    run_dir = base_path / runname
    candidates = [
        run_dir / "result_df_trueomega.csv",
        run_dir / "result_df.csv",
    ]
    for result_file in candidates:
        if result_file.exists():
            df = pd.read_csv(result_file)
            missing = [col for col in METHOD_COLUMNS if col not in df.columns]
            if missing:
                raise ValueError(
                    f"{result_file} is missing required columns: {', '.join(missing)}"
                )
            for col in METHOD_COLUMNS:
                df[col] = pd.to_numeric(df[col], errors="coerce")
            return df.dropna(subset=METHOD_COLUMNS).copy()
    raise FileNotFoundError(
        f"Missing result_df_trueomega.csv/result_df.csv under {run_dir}. "
        "Run visual_simulation.py aggregation first."
    )


def summarize_gain(df, x, hue, numerator, denominator, gain_label, panel):
    plot_df = df.loc[:, [x, hue, numerator, denominator]].copy()
    plot_df[gain_label] = plot_df[numerator] - plot_df[denominator]
    summary = (
        plot_df.groupby([x, hue], observed=True)[gain_label]
        .agg(["mean", "sem", "count"])
        .reset_index()
        .rename(columns={"count": "n"})
    )
    summary["ci95"] = 1.96 * summary["sem"].fillna(0)
    summary["panel"] = panel
    summary["x_name"] = x
    summary["hue_name"] = hue
    summary["gain_name"] = gain_label
    return summary


def format_value(value):
    value = float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def numeric_order(values):
    return sorted(values, key=lambda value: float(value))


def set_ylim(ax, summary):
    low = (summary["mean"] - summary["ci95"]).min()
    high = (summary["mean"] + summary["ci95"]).max()
    low = min(0, low)
    high = max(0, high)
    pad = max((high - low) * 0.18, 0.01)
    ax.set_ylim(low - pad, high + pad)


def plot_panel(ax, summary, x, hue, y_label, title, x_label, legend_title):
    x_order = numeric_order(summary[x].dropna().unique())
    hue_order = numeric_order(summary[hue].dropna().unique())
    x_positions = {value: idx for idx, value in enumerate(x_order)}

    for idx, hue_value in enumerate(hue_order):
        line_df = summary.loc[summary[hue] == hue_value].sort_values(x)
        positions = [x_positions[value] for value in line_df[x]]
        ax.errorbar(
            positions,
            line_df["mean"],
            yerr=line_df["ci95"],
            color=CATEGORY_COLORS[idx % len(CATEGORY_COLORS)],
            marker="o",
            markersize=5.5,
            linewidth=1.8,
            linestyle="-",
            capsize=3,
            label=format_value(hue_value),
        )

    ax.axhline(0, color="#6c757d", linestyle=":", linewidth=1.2)
    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels([format_value(value) for value in x_order])
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title, pad=8)
    ax.grid(True, alpha=0.75)
    set_ylim(ax, summary)
    ax.legend(title=legend_title, frameon=False, fontsize=8.5, title_fontsize=9)


def plot_figure(tracec_summary, tracecb_summary, output_prefix):
    sns.set_theme(
        style="darkgrid",
        palette="muted",
        color_codes=True,
        context="paper",
        font_scale=1.25,
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), constrained_layout=True)
    plot_panel(
        axes[0],
        tracec_summary,
        x="n2",
        hue="gc",
        y_label="Power gain (traceC - Original)",
        title=r"traceC gain by $N_2$",
        x_label=r"$N_2$",
        legend_title=r"$\rho$",
    )
    plot_panel(
        axes[1],
        tracecb_summary,
        x="propt",
        hue="nt",
        y_label="Power gain (traceCB - traceC)",
        title=r"traceCB gain by $\pi$",
        x_label=r"$\pi$",
        legend_title=r"$N_t$",
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    tracec_df = read_result_df(args.base_path, args.tracec_runname)
    tracecb_df = read_result_df(args.base_path, args.tracecb_runname)
    tracec_summary = summarize_gain(
        tracec_df,
        x="n2",
        hue="gc",
        numerator="cross",
        denominator="sumstat",
        gain_label="traceC - Original",
        panel="traceC",
    )
    tracecb_summary = summarize_gain(
        tracecb_df,
        x="propt",
        hue="nt",
        numerator="tissue",
        denominator="cross",
        gain_label="traceCB - traceC",
        panel="traceCB",
    )

    combined = pd.concat([tracec_summary, tracecb_summary], ignore_index=True)
    output_prefix = args.output_prefix
    plot_figure(tracec_summary, tracecb_summary, output_prefix)
    combined.to_csv(f"{output_prefix}_data.csv", index=False)
    print(f"Saved {output_prefix}.pdf")
    print(f"Saved {output_prefix}.png")
    print(f"Saved {output_prefix}_data.csv")


if __name__ == "__main__":
    main()
