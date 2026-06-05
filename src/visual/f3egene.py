# plot number of eGenes/eSNPs for each method
from visual.utils import *
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from scipy import stats


CELLTYPE_RANGES = {
    "Monocytes": (-0.5, 3),
    "CD4+T_cells": (2.5, 3),
    "CD8+T_cells": (5.5, 2),
    "B_cells": (7.5, 1),
    "NK_cells": (8.5, 1),
}


def prediction_interval(x, y, new_x, confidence=0.95):
    """Calculate the regression confidence interval."""
    n = len(x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    sxx = np.sum((x - x_mean) ** 2)
    sxy = np.sum((x - x_mean) * (y - y_mean))
    syy = np.sum((y - y_mean) ** 2)
    s = np.sqrt((syy - sxy**2 / sxx) / (n - 2))

    t_val = stats.t.ppf((1 + confidence) / 2, n - 2)
    se = s * np.sqrt(1 / n + (new_x - x_mean) ** 2 / sxx)
    return t_val * se


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


def prepare_growth_df(plot_df):
    growth_df = plot_df.copy()
    growth_df.loc[:, "CELL_TYPE"] = growth_df.QTDid.map(meta_data["id2celltype"])
    growth_df.loc[:, "SAMPLE_SIZE"] = (
        growth_df.name.str.extract(r"\((\d+)\)").astype(int)
    )
    growth_df.loc[:, "CELL_TYPE_PROP"] = growth_df.CELL_TYPE.map(
        meta_data["celltype_proportion"]
    )
    growth_df.loc[:, "traceC_growth_ratio"] = (
        growth_df[meta_data["method_name"][1]] / growth_df[meta_data["method_name"][0]]
    )
    growth_df.loc[:, "traceCB_growth_ratio"] = (
        growth_df[meta_data["method_name"][2]] / growth_df[meta_data["method_name"][1]]
    )
    growth_df.loc[:, "traceC_growth_rate"] = (
        growth_df["new_in_C"] / growth_df[meta_data["method_name"][0]]
    )
    growth_df.loc[:, "traceCB_growth_rate"] = (
        growth_df["new_in_T"] / growth_df[meta_data["method_name"][1]]
    )
    print("eGene growth summary:")
    print(
        growth_df.loc[
            :,
            [
                "QTDid",
                "name",
                "SAMPLE_SIZE",
                "CELL_TYPE",
                "CELL_TYPE_PROP",
                "traceC_growth_ratio",
                "traceCB_growth_ratio",
                "traceC_growth_rate",
                "traceCB_growth_rate",
                "new_in_C",
                "new_in_T",
            ],
        ]
    )
    return growth_df


def add_celltype_background(ax):
    celltype_colors = meta_data["celltype_colors"]
    ax.set_xlim(-0.5, len(meta_data["Names"]) - 0.5)
    y_min, y_max = ax.get_ylim()
    y_max += (y_max - y_min) * 0.1
    ax.set_ylim(y_min, y_max)
    margin = 0.05

    for celltype, (x_start, width) in CELLTYPE_RANGES.items():
        ax.add_patch(
            Rectangle(
                (x_start + margin, y_min),
                width - margin * 2,
                y_max - y_min,
                facecolor=celltype_colors[celltype],
                edgecolor="white",
                linewidth=0.2,
                alpha=0.3,
                zorder=0,
            )
        )

        x_end = x_start + width
        ax.hlines(
            y=y_max - 110,
            xmin=x_start + margin,
            xmax=x_end - margin,
            colors=celltype_colors[celltype],
            linewidth=13,
            linestyle="-",
            clip_on=False,
            zorder=5,
        )

    for celltype, (x_start, width) in CELLTYPE_RANGES.items():
        x_center = x_start + width / 2
        ax.text(
            x_center,
            y_max - (y_max - y_min) * 0.05,
            cell_label_name[celltype],
            color="black",
            fontsize=11,
            ha="center",
            va="bottom",
            clip_on=False,
            zorder=8,
        )


def plot_growth_scatter(
    plot_df,
    x_col,
    y_col,
    x_label,
    y_label,
    save_name,
    exclude_qtdids=None,
):
    fig, ax = plt.subplots(figsize=(6, 4))
    exclude_qtdids = exclude_qtdids or []
    regression_mask = ~plot_df.QTDid.isin(exclude_qtdids)
    x = plot_df.loc[regression_mask, x_col]
    y = plot_df.loc[regression_mask, y_col]

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.linspace(plot_df[x_col].min(), plot_df[x_col].max(), 100)
    line_y = slope * line_x + intercept
    margin = prediction_interval(x, y, line_x, confidence=0.95)
    print(
        f"{save_name} regression: y = {slope:.4f} * x + {intercept:.4f}, "
        f"R^2 = {r_value**2:.4f}, p = {p_value:.4g}"
    )

    ax.plot(
        line_x,
        line_y,
        color="black",
        linewidth=2,
        linestyle="--",
        alpha=0.5,
    )
    ax.fill_between(
        line_x,
        line_y - margin,
        line_y + margin,
        color="gray",
        alpha=0.1,
    )
    sns.scatterplot(
        x=x_col,
        y=y_col,
        hue="CELL_TYPE",
        data=plot_df,
        ax=ax,
        palette=meta_data["celltype_colors"],
        s=50,
        alpha=1,
        edgecolor="gray",
        linewidth=0.5,
    )

    texts = []
    for _, row in plot_df.iterrows():
        texts.append(
            ax.text(
                row[x_col],
                row[y_col],
                row["name"],
                fontsize=10,
                color="black",
                ha="right",
                va="bottom",
            )
        )

    adjust_text(
        texts,
        ax=ax,
        ha="center",
        va="center",
        force_text=(4.0, 4.0),
        force_points=(0.5, 0.5),
        force_objects=(0.5, 0.5),
        expand_text=(1.5, 1.5),
        expand_points=(1.5, 1.5),
        expand_objects=(1.5, 1.5),
        min_arrow_len=16,
        lim=8000,
        precision=0.0001,
        arrowprops=dict(
            arrowstyle="->",
            color="grey",
            lw=0.5,
            shrinkA=4,
            shrinkB=2,
        ),
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    handles, labels = ax.get_legend_handles_labels()
    shortened_labels = [label_name_shorten.get(label, label) for label in labels]
    ax.legend(
        handles=handles,
        labels=shortened_labels,
        title="Cell Type",
        loc="lower right",
    )

    plt.tight_layout()
    plt.savefig(os.path.join(save_path, save_name), bbox_inches="tight")
    print(f"Figure saved to {os.path.join(save_path, save_name)}")


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

    add_celltype_background(ax)

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
        bbox_to_anchor=(1.0, 0.94),
        title="Method",
    )

    plt.tight_layout()
    save_name = "f3egene.pdf"
    plt.savefig(os.path.join(save_path, save_name))
    print(f"Figure saved to {os.path.join(save_path, save_name)}")


def f3egene_increment(plot_df):
    fig, ax = plt.subplots(figsize=(8.5, 5))

    plot_df.QTDid = pd.Categorical(
        plot_df.QTDid, categories=meta_data["QTDids"], ordered=True
    )
    plot_df = plot_df.sort_values("QTDid")
    studies = plot_df.name.unique()
    x = np.arange(len(studies))
    width = 0.18

    ax.bar(
        x - 1.5 * width,
        plot_df["new_in_C"],
        width,
        label=f"New in {meta_data['method_name'][1]}",
        color=meta_data["Colors"][meta_data["method_name"][1]],
    )
    ax.bar(
        x - 0.5 * width,
        plot_df["new_in_C_replicate"],
        width,
        label=f"New in {meta_data['method_name'][1]} (replicate)",
        color=meta_data["Colors"][meta_data["method_name"][1]],
        hatch="\\\\",
    )
    ax.bar(
        x + 0.5 * width,
        plot_df["new_in_T"],
        width,
        label=f"New in {meta_data['method_name'][2]}",
        color=meta_data["Colors"][meta_data["method_name"][2]],
    )
    ax.bar(
        x + 1.5 * width,
        plot_df["new_in_T_replicate"],
        width,
        label=f"New in {meta_data['method_name'][2]} (replicate)",
        color=meta_data["Colors"][meta_data["method_name"][2]],
        hatch="\\\\",
    )

    ax.set_ylabel("Number of Newly Discovered eGenes")
    ax.set_xticks(x)
    ax.set_xticklabels(studies, rotation=20, ha="right")
    add_celltype_background(ax)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.94), title="Discovery")

    plt.tight_layout()
    save_name = "f3egene_increment.pdf"
    plt.savefig(os.path.join(save_path, save_name))
    print(f"Figure saved to {os.path.join(save_path, save_name)}")


def f3egene_growth_samplesize(plot_df):
    plot_growth_scatter(
        plot_df=plot_df,
        x_col="SAMPLE_SIZE",
        y_col="traceC_growth_ratio",
        x_label="Sample Size of Study",
        y_label=f"eGene Growth Ratio ({meta_data['method_name'][1]} / {meta_data['method_name'][0]})",
        save_name="f3egene_growth_samplesize.pdf",
        exclude_qtdids=["QTD000021", "QTD000031"],
    )


def f3egene_growth_celltype_proportion(plot_df):
    plot_growth_scatter(
        plot_df=plot_df,
        x_col="CELL_TYPE_PROP",
        y_col="traceCB_growth_ratio",
        x_label="Cell Type Proportion (Mean)",
        y_label=f"eGene Growth Ratio ({meta_data['method_name'][2]} / {meta_data['method_name'][1]})",
        save_name="f3egene_growth_celltype_proportion.pdf",
    )


if __name__ == "__main__":
    _, summary_df = load_all_summary()
    replicate_df = pd.read_csv("/home/group1/wjiang49/data/hum0343/hum0343_eGene.csv")
    replicate_egenes = replicate_df.gene
    plot_df = count_egene(summary_df, replicate_egenes)
    f3egene(plot_df)
    f3egene_increment(plot_df)
    growth_df = prepare_growth_df(plot_df)
    f3egene_growth_samplesize(growth_df)
    f3egene_growth_celltype_proportion(growth_df)
