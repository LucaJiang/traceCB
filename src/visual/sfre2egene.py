# plot number of eGenes/eSNPs for each method from RE2 results
from utils import *
from matplotlib.patches import Rectangle


def load_all_re2_summary(study_path_main=study_path_main):
    """
    Load all summary data from RE2 results in the main study directory.
    Returns a DataFrame with all results.
    columns: QTDid,GENE,NSNP,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SeSNP,TISSUE_SeSNP,CHR
    """
    all_qtdids = meta_data["QTDids"]
    all_summary_df = pd.DataFrame()
    for qtdid in all_qtdids:
        study_dir = os.path.join(study_path_main, qtdid)
        summary_dirs = glob.glob(os.path.join(study_dir, "RE2", "chr*", "summary.csv"))
        if not summary_dirs:
            continue

        summary_df = pd.DataFrame()
        for summary_file in summary_dirs:
            try:
                df = pd.read_csv(summary_file)
                df["CHR"] = os.path.basename(os.path.dirname(summary_file)).replace(
                    "chr", ""
                )
                summary_df = (
                    df
                    if summary_df.empty
                    else pd.concat([summary_df, df], ignore_index=True)
                )
            except (pd.errors.EmptyDataError, FileNotFoundError):
                print(f"Warning: Could not read or found empty file: {summary_file}")
                continue

        if not summary_df.empty:
            summary_df.loc[:, "QTDid"] = qtdid
            all_summary_df = (
                summary_df
                if all_summary_df.empty
                else pd.concat([all_summary_df, summary_df], ignore_index=True)
            )

    # Fill NA values for SNP counts with 0
    for col in ["TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"]:
        if col in all_summary_df.columns:
            all_summary_df[col] = all_summary_df[col].fillna(0)

    return all_summary_df


def count_re2_egene(df, replicate_egenes):
    """
    Count the number of eGenes for each qtdid and method from RE2 results.
    Count the number of replicate eGenes from each qtdid and method.
    Returns a DataFrame with the count of eGenes.
    """
    df = df.copy()
    df.loc[:, "is_replicate"] = df.GENE.isin(replicate_egenes)

    # S: Summary, C: Cross-study RE2, T: Cross-study+Tissue RE2
    methods = {
        "S": meta_data["method_name"][0],
        "C": meta_data["method_name"][1],
        "T": meta_data["method_name"][2],
    }

    for m_short, m_long in methods.items():
        df.loc[:, f"{m_long}"] = (df.loc[:, f"TAR_{m_short}eSNP"] > 0).astype(int)
        df.loc[:, f"{m_long}_replicate"] = (
            df.loc[:, "is_replicate"] & (df.loc[:, f"TAR_{m_short}eSNP"] > 0)
        ).astype(int)

    agg_dict = {m_long: (m_long, "sum") for m_long in methods.values()}
    agg_dict.update(
        {
            f"{m_long}_replicate": (f"{m_long}_replicate", "sum")
            for m_long in methods.values()
        }
    )

    count_df = df.groupby("QTDid").agg(**agg_dict)

    count_df.loc[:, "name"] = count_df.index.map(meta_data["id2name"])
    count_df.name = pd.Categorical(
        count_df.name, categories=meta_data["Names"], ordered=True
    )
    count_df = count_df.reset_index().sort_values("name")

    return count_df


def f3egene_re2(plot_df):
    """
    Plot eGene counts from RE2 results.
    """
    method_names = meta_data["method_name"]
    legend_order = meta_data["method_name"] + [m + " (replicate)" for m in method_names]
    fig, ax = plt.subplots(figsize=(8, 5))

    # Set x-axis positions
    studies = plot_df.name.unique()
    x = np.arange(len(studies))
    width = 0.35  # bar width

    # Plot bars for each method
    for i, m in enumerate(reversed(method_names)):
        # Original data
        ax.bar(x - width / 2, plot_df[m], width, label=m, color=meta_data["Colors"][m])

        # Replicate data
        ax.bar(
            x + width / 2,
            plot_df[m + "_replicate"],
            width,
            label=m + " (replicate)",
            color=meta_data["Colors"][m],
            hatch="\\\\",
        )

    ax.set_ylabel("Number of eGenes")
    ax.set_xlabel("Study")
    ax.set_xticks(x)
    ax.set_xticklabels(studies, rotation=20, ha="right")

    # Add celltype background rectangles
    celltype_colors = meta_data["celltype_colors"]
    celltype_ranges = {
        "Monocytes": (-0.5, 3),
        "CD4+T_cells": (2.5, 3),
        "CD8+T_cells": (5.5, 2),
        "B_cells": (7.5, 1),
        "NK_cells": (8.5, 1),
    }

    y_min, y_max = ax.get_ylim()
    y_max_plot_area = y_max - (y_max - y_min) * 0.05
    margin = 0.05
    for celltype, (x_start, rect_width) in celltype_ranges.items():
        ax.add_patch(
            Rectangle(
                (x_start + margin, y_min),
                rect_width - margin * 2,
                y_max_plot_area - y_min,
                facecolor=celltype_colors[celltype],
                edgecolor="white",
                linewidth=0.2,
                alpha=0.3,
                zorder=0,
            )
        )

        # Add celltype annotation lines
        x_end = x_start + rect_width
        ax.hlines(
            y=y_max_plot_area + (y_max - y_min) * 0.01,
            xmin=x_start + margin,
            xmax=x_end - margin,
            colors=celltype_colors[celltype],
            linewidth=5,
            linestyle="-",
            clip_on=False,
            zorder=5,
        )

    # Add celltype labels
    for celltype, (x_start, rect_width) in celltype_ranges.items():
        x_center = x_start + rect_width / 2
        ax.text(
            x_center,
            y_max_plot_area - (y_max - y_min) * 0.05,
            label_name_shorten[celltype],
            color="black",
            fontsize=10,
            ha="center",
            va="bottom",
            clip_on=False,
        )

    # Legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))
    ordered_handles = [
        label_to_handle[label] for label in legend_order if label in label_to_handle
    ]
    ordered_labels = [label for label in legend_order if label in label_to_handle]
    ax.legend(
        handles=ordered_handles,
        labels=ordered_labels,
        loc="upper right",
        bbox_to_anchor=(1.0, 0.9),
        title="Method",
    )

    plt.tight_layout()
    save_name = "f3egene_re2.png"
    plt.savefig(os.path.join(save_path, save_name), dpi=300)
    print(f"Plot saved to {os.path.join(save_path, save_name)}")


if __name__ == "__main__":
    # Load summary data from RE2 results
    all_summary_df = load_all_re2_summary()

    # Load replicate eGenes list
    replicate_df = pd.read_csv("/gpfs1/home/wjiang49/xpmm/data/hum0304_eGene.csv")
    replicate_egenes = replicate_df.gene

    # Count eGenes and replicate eGenes
    plot_df = count_re2_egene(all_summary_df, replicate_egenes)

    # Generate and save the plot
    f3egene_re2(plot_df)
