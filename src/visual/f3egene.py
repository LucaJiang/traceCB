# plot number of eGenes/eSNPs for each method
from visual.utils import *
from matplotlib.patches import Rectangle


def count_egene(df, replicate_egenes):
    """
    Count the number of eGenes for each qtdid and method.
    Count the number of replicate eGenes from each qtdid and method.
    Returns a DataFrame with the count of eGenes.
    """
    df = df.copy()
    df.loc[:, "is_replicate"] = df.GENE.isin(replicate_egenes)
    for m in ["S", "C", "T"]:
        df.loc[:, f"TAR_{m}eGene"] = (df.loc[:, f"TAR_{m}eSNP"] > 0).astype(int)
        df.loc[:, f"TAR_{m}eGene_replicate"] = (
            df.loc[:, "is_replicate"] & df.loc[:, f"TAR_{m}eGene"]
        ).astype(int)

    # traceC 相对于 Original 的新 eGene
    df["new_in_C"] = (df["TAR_CeGene"] > 0) & (df["TAR_SeGene"] == 0)
    df["new_in_C_replicate"] = df["new_in_C"] & df["is_replicate"]

    # traceCB 相对于 traceC 的新 eGene
    df["new_in_T"] = (df["TAR_TeGene"] > 0) & (df["TAR_CeGene"] == 0)
    df["new_in_T_replicate"] = df["new_in_T"] & df["is_replicate"]

    count_df = df.groupby("QTDid").agg(
        S=("TAR_SeGene", "sum"),
        C=("TAR_CeGene", "sum"),
        T=("TAR_TeGene", "sum"),
        S_replicate=("TAR_SeGene_replicate", "sum"),
        C_replicate=("TAR_CeGene_replicate", "sum"),
        T_replicate=("TAR_TeGene_replicate", "sum"),
        new_in_C=("new_in_C", "sum"),
        new_in_C_replicate=("new_in_C_replicate", "sum"),
        new_in_T=("new_in_T", "sum"),
        new_in_T_replicate=("new_in_T_replicate", "sum"),
    )
    count_df = count_df.rename(
        columns={
            "S": meta_data["method_name"][0],
            "C": meta_data["method_name"][1],
            "T": meta_data["method_name"][2],
        }
    )
    count_df = count_df.rename(
        columns={
            "S_replicate": meta_data["method_name"][0] + "_replicate",
            "C_replicate": meta_data["method_name"][1] + "_replicate",
            "T_replicate": meta_data["method_name"][2] + "_replicate",
        }
    )
    count_df.loc[:, "name"] = count_df.index.map(meta_data["id2name"])
    # count_df.name = pd.Categorical(
    #     count_df.name, categories=meta_data["Names"], ordered=False
    # )
    count_df = count_df.reset_index().sort_values("name")

    # find average replicate rate
    print(count_df)
    replicate_rate = (
        count_df.loc[:, [m + "_replicate" for m in meta_data["method_name"]]]
        .sum()
        .values
        / count_df.loc[:, [m for m in meta_data["method_name"]]].sum().values
    )
    print("Average replicate rate:")
    print(replicate_rate)

    # 计算并打印新发现 eGene 的平均复现率
    print("\nAverage replicate rate for newly discovered eGenes:")
    new_C_total = count_df["new_in_C"].sum()
    new_C_replicate_total = count_df["new_in_C_replicate"].sum()
    new_C_rate = new_C_replicate_total / new_C_total if new_C_total > 0 else 0
    print(f"New in traceC (vs Original): {new_C_rate:.2%}")

    new_T_total = count_df["new_in_T"].sum()
    new_T_replicate_total = count_df["new_in_T_replicate"].sum()
    new_T_rate = new_T_replicate_total / new_T_total if new_T_total > 0 else 0
    print(f"New in traceCB (vs traceC): {new_T_rate:.2%}")

    return count_df


def f3egene(plot_df):
    method_names = meta_data["method_name"]
    legend_order = meta_data["method_name"] + [m + " (replicate)" for m in method_names]
    fig, ax = plt.subplots(figsize=(8, 5))  # 增加高度以容纳上方标注

    # 设置x轴位置
    plot_df.QTDid = pd.Categorical(
        plot_df.QTDid, categories=meta_data["QTDids"], ordered=True
    )
    plot_df = plot_df.sort_values("QTDid")
    studies = plot_df.name.unique()
    x = np.arange(len(studies))
    width = 0.35  # 条形宽度

    # 绘制每个方法的条形图
    for i, m in enumerate(reversed(method_names)):
        # 原始数据
        ax.bar(x - width / 2, plot_df[m], width, label=m, color=meta_data["Colors"][m])

        # replicate数据
        ax.bar(
            x + width / 2,
            plot_df[m + "_replicate"],
            width,
            label=m + " (replicate)",
            color=meta_data["Colors"][m],
            # alpha=0.7,
            hatch="\\\\",
        )

    ax.set_ylabel("Number of eGenes")
    # ax.set_xlabel("Study")
    ax.set_xticks(x)
    ax.set_xticklabels(studies, rotation=20, ha="right")

    # 添加celltype背景矩形
    celltype_colors = meta_data["celltype_colors"]
    celltype_ranges = {  # (起始索引, 覆盖宽度)
        "Monocytes": (-0.5, 3),
        "CD4+T_cells": (2.5, 3),
        "CD8+T_cells": (5.5, 2),
        "B_cells": (7.5, 1),
        "NK_cells": (8.5, 1),
    }
    ax.set_xlim(-0.5, len(meta_data["Names"]) - 0.5)
    # 设置y轴范围以便添加标注线
    y_min, y_max = ax.get_ylim()
    y_max += (y_max - y_min) * 0.1
    ax.set_ylim(y_min, y_max)  # 留出空间用于标注
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

        # 添加celltype标注线 - 在图形上边缘上方
        x_end = x_start + width
        ax.hlines(
            y=y_max - 110,  # 位于图形上边缘上方
            xmin=x_start + margin,
            xmax=x_end - margin,
            colors=celltype_colors[celltype],
            linewidth=13,
            linestyle="-",
            clip_on=False,
            zorder=5,
        )

    # 添加celltype标签 - 在标注线下方
    for celltype, (x_start, width) in celltype_ranges.items():
        x_center = x_start + width / 2
        ax.text(
            x_center,
            y_max - (y_max - y_min) * 0.05,  # 位于标注线上方
            cell_label_name[celltype],
            color="black",
            fontsize=11,
            ha="center",
            va="bottom",
            clip_on=False,
            zorder=8,
        )

    # 图例handles和labels
    handles, labels = ax.get_legend_handles_labels()
    # 创建label到handle的映射
    label_to_handle = dict(zip(labels, handles))
    # 按照legend_order重新排序
    ordered_handles = [
        label_to_handle[label] for label in legend_order if label in label_to_handle
    ]
    ordered_labels = [label for label in legend_order if label in label_to_handle]
    ax.legend(
        handles=ordered_handles,
        labels=ordered_labels,
        loc="upper right",
        bbox_to_anchor=(1.0, 0.93),
        title="Method",
    )

    plt.tight_layout()
    save_name = "f3egene.pdf"
    plt.savefig(os.path.join(save_path, save_name))
    print(f"Figure saved to {os.path.join(save_path, save_name)}")


if __name__ == "__main__":
    _, summary_df = load_all_summary()
    replicate_df = pd.read_csv("/home/group1/wjiang49/data/hum0343/hum0343_eGene.csv")
    replicate_egenes = replicate_df.gene
    plot_df = count_egene(summary_df, replicate_egenes)
    f3egene(plot_df)
