# plot f5 AFR+GTEx violin plot for effective sample size
from utils import *
from f3violin import remove_outliers
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import os
import warnings

warnings.filterwarnings("ignore")

study_path_main = "/home/group1/wjiang49/data/traceCB/AFR_eQTLGen"
save_path = "/home/group1/wjiang49/data/traceCB/AFR_eQTLGen/results"

TAR_SAMPLE_SIZE = 80


def f5violin(summary_sign_df):
    """
    Plot Fig. 5 violin plot to visualize the effective sample size for AFR+GTEx
    """
    # prepare data for plotting
    melt_df = (
        summary_sign_df.loc[:, ["NAME", "TAR_CNEFF", "TAR_TNEFF"]]
        .rename(
            columns={
                "TAR_CNEFF": meta_data["method_name"][1],
                "TAR_TNEFF": meta_data["method_name"][2],
            }
        )
        .melt(
            id_vars=["NAME"],
            value_vars=[meta_data["method_name"][1], meta_data["method_name"][2]],
            var_name="Method",
            value_name="Effective Sample Size",
        )
        .dropna(subset=["Effective Sample Size"])
    )
    # sort the x-axis by meta_data["Names"]
    melt_df["NAME"] = pd.Categorical(
        melt_df["NAME"], categories=meta_data["Names"], ordered=True
    )
    melt_df.sort_values("NAME", inplace=True)

    # plot violin plot
    plt.figure(figsize=(8, 5))
    ax = sns.boxenplot(
        x="NAME",
        y="Effective Sample Size",
        hue="Method",
        data=melt_df,
        width=0.7,
        dodge=True,
        showfliers=False,
        palette=[
            meta_data["Colors"][meta_data["method_name"][1]],
            meta_data["Colors"][meta_data["method_name"][2]],
        ],
    )
    # Add horizontal line
    plt.axhline(TAR_SAMPLE_SIZE, color="red", linestyle="--")

    # Create a single legend with both elements
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color="red", linestyle="--", label="AFR (80)"))
    labels.append("AFR (80)")
    plt.legend(
        handles=handles,
        labels=labels,
        title="Method",
        bbox_to_anchor=(0.65, 1),
        loc="upper left",
    )

    y_limits = plt.ylim()
    plt.ylim(0, y_limits[1])
    plt.xticks(rotation=30, ha="right")
    plt.xlabel("Study")
    plt.ylabel("Effective Sample Size")

    # Add celltype background rectangles
    celltype_colors = meta_data["celltype_colors"]
    celltype_ranges = {  # (起始索引, 覆盖宽度)
        "Monocytes": (-0.5, 3),
        "CD4+T_cells": (2.5, 3),
        "CD8+T_cells": (5.5, 2),
        "B_cells": (7.5, 1),
        "NK_cells": (8.5, 1),
    }
    y_min, y_max = 0, y_limits[1]
    margin = 0.05
    for celltype, (x_start, width) in celltype_ranges.items():
        ax.add_patch(
            Rectangle(
                (x_start + margin, y_min),  # 左下角坐标
                width - margin * 2,  # 宽度
                y_max - y_min,  # 高度（覆盖整个y轴范围）
                facecolor=celltype_colors[celltype],
                edgecolor="white",
                linewidth=0.2,
                alpha=0.3,
                zorder=0,  # 确保在最底层
            )
        )

        # 添加celltype标注线
        x_end = x_start + width
        ax.hlines(
            y=y_min + (y_max - y_min) * 0.02,  # 位于x轴上方
            xmin=x_start + margin,
            xmax=x_end - margin,
            colors=celltype_colors[celltype],
            linewidth=6,
            linestyle="-",
            clip_on=False,
            zorder=5,
        )

    for celltype, (x_start, width) in celltype_ranges.items():
        x_center = x_start + width / 2
        ax.text(
            x_center,
            y_min + (y_max - y_min) * 0.07,  # 位于标注线上方
            label_name_shorten[celltype],
            color="black",
            # fontsize=10,
            ha="center",
            va="top",
            clip_on=False,
        )
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "f5violin.pdf"), dpi=300, bbox_inches="tight")
    plt.close()  # 添加close以释放内存
    print(f"Violin plot saved to {os.path.join(save_path, 'f5violin.pdf')}")


if __name__ == "__main__":
    summary_sign_df, _ = load_all_summary(study_path_main)
    summary_sign_df = remove_outliers(summary_sign_df)
    summary_sign_df.loc[:, "NAME"] = summary_sign_df.QTDid.map(meta_data["id2name"])
    f5violin(summary_sign_df)
