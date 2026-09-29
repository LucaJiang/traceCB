import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
from figures.paths import FIGURE_DIR, METADATA_FILE, TIMING_FILE, require_file


def main():
    # Define paths
    data_path = require_file(TIMING_FILE, "TRACECB_TIMING_FILE")
    meta_path = require_file(METADATA_FILE, "TRACECB_FIGURE_METADATA")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    output_plot = FIGURE_DIR / "runtime_comparison.pdf"

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

    # Layout follows the statistics version, with room for individual records.
    sns.set_theme(style="white", font="DejaVu Sans")
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[6, 1.6], height_ratios=[2.75, 6],
        wspace=0.03, hspace=0.025,
    )
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.18, top=0.96)

    # Define colors
    bar_color = "#a8dcb1"
    bar_edge_color = "#457b9d"
    heatmap_cmap = "Blues"  # lighter is less time, darker is more time

    # 1. Top Bar Plot (Average time per Chromosome)
    ax_top = fig.add_subplot(gs[0, 0])
    chromosome_positions = np.arange(len(chr_means))
    ax_top.bar(
        chromosome_positions,
        chr_means.values,
        color=bar_color,
        edgecolor=bar_edge_color,
        width=0.8,
        zorder=1,
    )
    # Each dot is one study-chromosome elapsed record, not a benchmark replicate.
    # Fixed offsets separate records without changing their observed durations.
    palette = sns.color_palette("tab10", n_colors=len(ordered_studies))
    offsets = (
        np.linspace(-0.27, 0.27, len(ordered_studies))
        if len(ordered_studies) > 1 else np.zeros(1)
    )
    for index, (study, row) in enumerate(pivot_df.iterrows()):
        ax_top.scatter(
            chromosome_positions + offsets[index],
            row.to_numpy(),
            s=17,
            color=palette[index],
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
            label=study,
        )
    ax_top.set_xlim(-0.5, len(chr_means) - 0.5)
    ax_top.set_ylim(bottom=0)
    ax_top.set_xticks([])
    ax_top.set_ylabel("Elapsed time (s)", fontsize=12, labelpad=10)
    ax_top.legend(
        loc="upper left", bbox_to_anchor=(1.01, 1.02),
        fontsize=10, frameon=False, borderaxespad=0,
    )

    # 2. Right Bar Plot (Total time per Study)
    ax_right = fig.add_subplot(gs[1, 1])
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
    ax_right.set_xlabel("Sum of chromosome\nelapsed times (s)", fontsize=12)
    ax_right.set_xlim(0, study_totals.max() * 1.2)

    # Add value labels
    max_val_right = study_totals.max()
    for i, v in enumerate(study_totals.values):
        ax_right.text(
            v + (max_val_right * 0.01),
            i,
            f"{int(v)}",
            va="center",
            fontsize=10,
        )

    # 3. Heatmap
    ax_main = fig.add_subplot(gs[1, 0])
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

    ax_main.set_xlabel("Chromosome", fontsize=12)
    ax_main.set_ylabel("Study", fontsize=12)

    # Rotate y-axis labels to be horizontal
    plt.setp(ax_main.get_yticklabels(), rotation=0)

    # Separate colorbar keeps the top bars and heatmap columns aligned.
    heatmap_position = ax_main.get_position()
    colorbar_ax = fig.add_axes([
        heatmap_position.x0 + 0.49 * heatmap_position.width,
        0.065,
        0.42 * heatmap_position.width,
        0.012,
    ])
    colorbar = fig.colorbar(
        ax_main.collections[0], cax=colorbar_ax, orientation="horizontal"
    )
    colorbar.set_label("Elapsed time (s)")
    colorbar.outline.set_visible(False)

    fig.savefig(output_plot, dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print(f"Visualization saved to: {output_plot}")


if __name__ == "__main__":
    main()
