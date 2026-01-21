import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import json
import os
import sys


def main():
    # Define paths
    base_path = "/home/wjiang49/traceCB"
    data_path = os.path.join(base_path, "tmp/timing/summary_timing.csv")
    meta_path = os.path.join(base_path, "src/visual/metadata.json")
    output_plot = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/timing_heatmap.pdf"

    # Check if files exist
    if not os.path.exists(data_path):
        print(f"Error: Data file not found at {data_path}")
        return
    if not os.path.exists(meta_path):
        print(f"Error: Metadata file not found at {meta_path}")
        return

    # Load data
    df = pd.read_csv(data_path)
    with open(meta_path, "r") as f:
        meta = json.load(f)

    id2name = meta.get("id2name", {})

    # Map study_id to study_name
    df["study_name"] = df["study_id"].map(lambda x: id2name.get(x, x))

    # Pivot the data: index=study, columns=chromosome, values=elapsed_seconds
    # Ensure chromosome is int for proper sorting
    df["chromosome"] = df["chromosome"].astype(int)
    pivot_df = df.pivot(
        index="study_name", columns="chromosome", values="elapsed_seconds"
    )

    # Sort columns (Chromosomes 1-22)
    sorted_cols = sorted(pivot_df.columns)
    pivot_df = pivot_df[sorted_cols]

    # Sort rows (Studies) alphabetically by study name
    ordered_studies = sorted(pivot_df.index)
    pivot_df = pivot_df.reindex(ordered_studies)

    # Calculate stats
    chr_means = pivot_df.mean(axis=0)  # Top bar (Average)
    study_totals = pivot_df.sum(axis=1)  # Right bar (Total)
    print("Study average time (s):", study_totals.mean())
    # Study average time (s): 3201.7

    # Set up the plot grid and style
    sns.set_theme(style="white", font_scale=1.2)
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(
        2, 2, width_ratios=[6, 1], height_ratios=[1, 6], wspace=0.03, hspace=0.03
    )

    # Define colors
    bar_color = "#a8dcb1"
    bar_edge_color = "#457b9d"
    heatmap_cmap = "RdYlBu_r"  # Red (slow) to Blue (fast), or rev.
    # Actually YlOrRd is good for time (darker = more time). Let's stick to a nice one.
    heatmap_cmap = sns.cubehelix_palette(
        start=0.5, rot=-0.5, as_cmap=True
    )  # nice custom map
    # Or just "Viridis" or "Rocket"
    heatmap_cmap = "Blues"  # lighter is less time, darker is more time

    # 1. Top Bar Plot (Average time per Chromosome)
    ax_top = plt.subplot(gs[0, 0])
    ax_top.bar(
        range(len(chr_means)),
        chr_means.values,
        color=bar_color,
        edgecolor=bar_edge_color,
        width=0.8,
    )
    ax_top.set_xlim(-0.5, len(chr_means) - 0.5)
    ax_top.set_xticks([])
    ax_top.set_ylabel("Avg Time (s)", fontsize=14, labelpad=10)
    ax_top.set_title("Running Time Distribution", fontsize=20, pad=20)
    sns.despine(ax=ax_top, bottom=True)  # Remove bottom line since it touches heatmap

    # Add value labels
    max_val_top = chr_means.max()
    for i, v in enumerate(chr_means.values):
        ax_top.text(
            i,
            v + (max_val_top * 0.02),
            f"{int(v)}",
            ha="center",
            va="bottom",
            fontsize=10,
            rotation=45,
            color="#1d3557",
        )

    # 2. Right Bar Plot (Total time per Study)
    ax_right = plt.subplot(gs[1, 1])
    y_pos = range(len(study_totals))
    ax_right.barh(
        y_pos,
        study_totals.values,
        color=bar_color,
        edgecolor=bar_edge_color,
        height=0.8,
    )
    ax_right.set_ylim(len(study_totals) - 0.5, -0.5)
    ax_right.set_yticks([])
    ax_right.set_xlabel("Total CPU Time (s)", fontsize=14, labelpad=10)
    sns.despine(ax=ax_right, left=True)

    # Add value labels
    max_val_right = study_totals.max()
    for i, v in enumerate(study_totals.values):
        ax_right.text(
            v + (max_val_right * 0.02),
            i,
            f" {int(v)}",
            va="center",
            fontsize=11,
            color="#1d3557",
        )

    # 3. Heatmap
    ax_main = plt.subplot(gs[1, 0])
    sns.heatmap(
        pivot_df,
        annot=True,
        fmt=".0f",
        cmap=heatmap_cmap,
        ax=ax_main,
        cbar=False,
        annot_kws={"size": 10},
        linewidths=0.5,
        linecolor="white",
    )

    ax_main.set_xlabel("Chromosome", fontsize=16, labelpad=10)
    ax_main.set_ylabel("Study", fontsize=16, labelpad=10)

    # Rotate y-axis labels to be horizontal
    plt.setp(ax_main.get_yticklabels(), rotation=0)

    plt.tight_layout()

    # Create directory if it doesn't exist (though script assumes timing dir exists)
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)

    plt.savefig(output_plot, dpi=300, bbox_inches="tight", pad_inches=0.1)
    print(f"Visualization saved to: {output_plot}")


if __name__ == "__main__":
    main()
