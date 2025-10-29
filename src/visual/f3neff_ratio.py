from utils import *
from f3violin import remove_outliers
from adjustText import adjust_text
from scipy import stats

# 计算regression置信区间
def prediction_interval(x, y, new_x, confidence=0.95):
    """计算预测区间"""
    n = len(x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    # 计算标准误差
    sxx = np.sum((x - x_mean) ** 2)
    sxy = np.sum((x - x_mean) * (y - y_mean))
    syy = np.sum((y - y_mean) ** 2)

    s = np.sqrt((syy - sxy**2 / sxx) / (n - 2))

    # t分布临界值
    t_val = stats.t.ppf((1 + confidence) / 2, n - 2)

    # 置信区间
    se = s * np.sqrt(1 / n + (new_x - x_mean) ** 2 / sxx)
    margin = t_val * se

    return margin


def f3neff_ratio_samplesize(summary_sign_df):
    """
    x=true sample size of the AUX, y=TAR_CNEFF / TAR_SNEFF, color=cell type
    """
    print("Beginning to plot sample size vs TAR_CNEFF_SNEFF_RATIO...")
    ## group by QTDid and get the mean of TAR_CNEFF_SNEFF_RATIO
    samplesize_df = (
        summary_sign_df.groupby("QTDid")
        .agg(
            SAMPLE_SIZE=("SAMPLE_SIZE", "mean"),
            TAR_CNEFF_SNEFF_RATIO=("TAR_CNEFF_SNEFF_RATIO", "median"),
        )
        .reset_index()
        .copy()
    )
    samplesize_df.loc[:, "CELL_TYPE"] = samplesize_df.QTDid.map(
        meta_data["id2celltype"]
    )
    samplesize_df.loc[:, "NAME"] = samplesize_df.QTDid.map(meta_data["id2name"])

    print("Sample size plot data:")
    print(samplesize_df)

    # plot
    fig, ax = plt.subplots(figsize=(6, 4))
    # 添加回归线和置信区间, exclude qtdid == "BLUEPRINT(191)" or "BLUEPRINT(167)"
    blueprint_rows = samplesize_df.QTDid.isin(["QTD000021", "QTD000031"])
    x = samplesize_df.loc[~blueprint_rows, "SAMPLE_SIZE"]
    y = samplesize_df.loc[~blueprint_rows, "TAR_CNEFF_SNEFF_RATIO"]

    # 计算回归线
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.linspace(
        samplesize_df.SAMPLE_SIZE.min(), samplesize_df.SAMPLE_SIZE.max(), 100
    )
    line_y = slope * line_x + intercept
    margin = prediction_interval(x, y, line_x, confidence=0.95)
    print(f"Regression line: y = {slope:.4f} * x + {intercept:.4f}, R^2 = {r_value**2:.4f}")

    # 绘制回归线
    ax.plot(
        line_x,
        line_y,
        color="black",
        linewidth=2,
        linestyle="--",
        alpha=0.5,
    )

    # 绘制置信区间
    ax.fill_between(
        line_x,
        line_y - margin,
        line_y + margin,
        color="gray",
        alpha=0.1,
    )
    sns.scatterplot(
        x="SAMPLE_SIZE",
        y="TAR_CNEFF_SNEFF_RATIO",
        hue="CELL_TYPE",
        data=samplesize_df,
        ax=ax,
        palette=meta_data["celltype_colors"],
        s=50,  # 稍微增大点的大小
        alpha=1,
        edgecolor="gray",
        linewidth=0.5,
    )

    # 使用adjustText自动避让
    texts_upper = []
    for i, row in samplesize_df.iterrows():
        texts_upper.append(
            ax.text(
                row.SAMPLE_SIZE,
                row.TAR_CNEFF_SNEFF_RATIO,
                row.NAME,
                fontsize=10,
                color="black",
                ha="right",
                va="bottom",
            )
        )

    # 设置避让参数
    adjust_text(
        texts_upper,
        ax=ax,
        ha="center",
        va="center",
        force_text=(4.0, 4.0),  # 增加文本之间的排斥力
        force_points=(0.5, 0.5),  # 增加文本与点之间的排斥力
        force_objects=(0.5, 0.5),
        expand_text=(1.5, 1.5),  # 增加文本周围的扩展区域
        expand_points=(1.5, 1.5),  # 增加点周围的扩展区域
        expand_objects=(1.5, 1.5),  # 增加对象周围的扩展区域
        min_arrow_len=16,
        lim=8000,  # 增加迭代次数
        precision=0.0001,  # 提高精度
        arrowprops=dict(
            arrowstyle="->",
            color="grey",
            lw=0.5,
            shrinkA=4,  # 增加箭头与文本的距离
            shrinkB=2,  # 增加箭头与点的距离
        ),
    )
    ax.set_xlabel("Sample Size of Study")
    ax.set_ylabel(
        f"{meta_data['method_name'][1]} / {meta_data['method_name'][0]} (median)"
    )

    handles, labels = ax.get_legend_handles_labels()
    shortened_labels = [label_name_shorten.get(label, label) for label in labels]
    ax.legend(
        handles=handles,
        labels=shortened_labels,
        title="Cell Type",
        loc="lower right",
    )

    plt.tight_layout()
    plt.savefig(
        f"{save_path}/f3neff_ratio_samplesize.pdf", bbox_inches="tight"
    )
    print(f"Sample size plot saved to: {save_path}/f3neff_ratio_samplesize.pdf")


def f3neff_ratio_celltype_proportion(summary_sign_df):
    """
    x=cell type proportion, y=TAR_TNEFF / TAR_CNEFF, color=cell type
    """
    print("Beginning to plot cell type proportion vs TAR_TNEFF_CNEFF_RATIO...")
    ## group by QTDid and get the mean of TAR_TNEFF_CNEFF_RATIO
    celltype_proportion_df = (
        summary_sign_df.groupby("QTDid")
        .agg(
            TAR_TNEFF_CNEFF_RATIO=("TAR_TNEFF_CNEFF_RATIO", "median"),
        )
        .reset_index()
        .copy()
    )
    celltype_proportion_df.loc[:, "CELL_TYPE"] = celltype_proportion_df.QTDid.map(
        meta_data["id2celltype"]
    )
    celltype_proportion_df.loc[:, "CELL_TYPE_PROP"] = (
        celltype_proportion_df.CELL_TYPE.map(meta_data["celltype_proportion"])
    )
    celltype_proportion_df.loc[:, "NAME"] = celltype_proportion_df.QTDid.map(
        meta_data["id2name"]
    )

    print("Celltype proportion plot data:")
    print(celltype_proportion_df)
    
    # plot
    fig, ax = plt.subplots(figsize=(6, 4))
    # 添加回归线和置信区间
    x_lower = celltype_proportion_df["CELL_TYPE_PROP"]
    y_lower = celltype_proportion_df["TAR_TNEFF_CNEFF_RATIO"]

    # 计算回归线
    slope_lower, intercept_lower, r_value_lower, p_value_lower, std_err_lower = (
        stats.linregress(x_lower, y_lower)
    )
    line_x_lower = np.linspace(x_lower.min(), x_lower.max(), 100)
    line_y_lower = slope_lower * line_x_lower + intercept_lower
    margin_lower = prediction_interval(x_lower, y_lower, line_x_lower, confidence=0.95)
    print(f"Regression line: y = {slope_lower:.4f} * x + {intercept_lower:.4f}, R^2 = {r_value_lower**2:.4f}")

    # 绘制回归线
    ax.plot(
        line_x_lower,
        line_y_lower,
        color="black",
        linewidth=2,
        linestyle="--",
        alpha=0.5,
    )

    # 绘制置信区间
    ax.fill_between(
        line_x_lower,
        line_y_lower - margin_lower,
        line_y_lower + margin_lower,
        color="gray",
        alpha=0.1,
    )
    sns.scatterplot(
        x="CELL_TYPE_PROP",
        y="TAR_TNEFF_CNEFF_RATIO",
        hue="CELL_TYPE",
        data=celltype_proportion_df,
        ax=ax,
        palette=meta_data["celltype_colors"],
        s=50,  # 稍微增大点的大小
        alpha=1,
        edgecolor="gray",
        linewidth=0.5,
    )
    
    # 采用斜向偏移策略
    texts_lower = []
    for i, row in celltype_proportion_df.iterrows():
        texts_lower.append(
            ax.text(
                row.CELL_TYPE_PROP,
                row.TAR_TNEFF_CNEFF_RATIO,
                row.NAME,
                fontsize=10,
                color="black",
                ha="right",
                va="bottom",
            )
        )
    adjust_text(
        texts_lower,
        ax=ax,
        ha="center",
        va="center",
        force_text=(4.0, 4.0),  # 增加文本之间的排斥力
        force_points=(0.5, 0.5),  # 增加文本与点之间的排斥力
        force_objects=(0.5, 0.5),
        expand_text=(1.5, 1.5),  # 增加文本周围的扩展区域
        expand_points=(1.5, 1.5),  # 增加点周围的扩展区域
        expand_objects=(1.5, 1.5),  # 增加对象周围的扩展区域
        lim=4000,  # 增加迭代次数
        precision=0.0001,  # 提高精度
        min_arrow_len=12,
        arrowprops=dict(
            arrowstyle="->",
            color="grey",
            lw=0.5,
            shrinkA=4,  # 增加箭头与文本的距离
            shrinkB=2,  # 增加箭头与点的距离
        ),
    )
    ax.set_xlabel("Cell Type Proportion (Mean)")
    ax.set_ylabel(
        f"{meta_data['method_name'][2]} / {meta_data['method_name'][1]} (median)"
    )

    handles, labels = ax.get_legend_handles_labels()
    shortened_labels = [label_name_shorten.get(label, label) for label in labels]
    ax.legend(
        handles=handles,
        labels=shortened_labels,
        title="Cell Type",
        loc="lower right",
    )

    plt.tight_layout()
    plt.savefig(
        f"{save_path}/f3neff_ratio_celltype_proportion.pdf",
        bbox_inches="tight",
    )
    print(f"Cell type proportion plot saved to: {save_path}/f3neff_ratio_celltype_proportion.pdf")


if __name__ == "__main__":
    summary_sign_df, _ = load_all_summary()    
    summary_sign_df = remove_outliers(summary_sign_df)
    summary_sign_df.loc[:, "NAME"] = summary_sign_df.QTDid.map(meta_data["id2name"])

    # prepare data
    summary_sign_df.loc[:, "SAMPLE_SIZE"] = summary_sign_df.NAME.str.extract(
        r"\((\d+)\)"  # NAME: BLUEPRINT(191)
    ).astype(int)
    summary_sign_df.loc[:, "TAR_CNEFF_SNEFF_RATIO"] = (
        summary_sign_df.TAR_CNEFF / summary_sign_df.TAR_SNEFF
    )
    summary_sign_df.loc[:, "TAR_TNEFF_CNEFF_RATIO"] = (
        summary_sign_df.TAR_TNEFF / summary_sign_df.TAR_CNEFF
    )
    f3neff_ratio_samplesize(summary_sign_df)
    f3neff_ratio_celltype_proportion(summary_sign_df)
