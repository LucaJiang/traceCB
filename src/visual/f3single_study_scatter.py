# figure 4: case study of single qtd study
# 1. plot scatter plot of effective sample size
# 2. plot by cor group
# 2.1 num of egene by method, proportion of egene by method
# 2.2 ratio of effective sample size: GMM/Original, GMM+/Original
from utils import *
from adjustText import adjust_text

target_qtdids = ["QTD000067", "QTD000371", "QTD000069", "QTD000081"]


def f3_neff_scatter(summary_df, target_qtdid, text_annot=True):
    # plot scatter plot of effective sample size
    ## p1: x=TAR_SNEFF, y=TAR_CNEFF, hue=COV_PVAL<0.05
    ## p2: x=TAR_SNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    ## p3: x=TAR_CNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    fig, ax = plt.subplots(1, 3, figsize=(12, 3))
    ## find significant and not significant xpop cov
    palette = {
        "Significant": "#db1717f8",
        "Not Significant": "#120eff",
    }
    summary_df.loc[:, "COV_SIGN"] = summary_df.COV_PVAL.apply(
        lambda x: "Significant" if x < 0.05 else "Not Significant"
    )
    ## init annotation class
    annot_class = onek1k_geneid2name()
    ## define method name
    column2method = {
        "TAR_SNEFF": meta_data["method_name"][0],
        "TAR_CNEFF": meta_data["method_name"][1],
        "TAR_TNEFF": meta_data["method_name"][2],
    }
    ## plot each subplot
    for i, (x_col, y_col) in enumerate(
        [
            ("TAR_SNEFF", "TAR_CNEFF"),
            ("TAR_SNEFF", "TAR_TNEFF"),
            ("TAR_CNEFF", "TAR_TNEFF"),
        ]
    ):
        sns.scatterplot(
            data=summary_df,
            x=x_col,
            y=y_col,
            hue="COV_SIGN",
            palette=palette,
            ax=ax[i],
            s=4,
        )
        ## add text annotation
        ### define gene to annotate
        summary_df.loc[:, "ANNOTATE_ORDER"] = summary_df.loc[
            :, y_col
        ]  # - 1.5*summary_df.loc[:, x_col]
        if text_annot:
            num_annotate_genes = 10
            row_to_annotate = summary_df.nlargest(num_annotate_genes, "ANNOTATE_ORDER")
            texts = []
            for _, row in row_to_annotate.iterrows():
                gene_name = annot_class.get_gene_name(row.GENE)
                if gene_name is not None:
                    texts.append(
                        ax[i].text(
                            row[x_col],
                            row[y_col],
                            gene_name,
                            fontsize=8,
                            fontstyle="italic",
                        )
                    )
            # Looks like you are using a tranform that doesn't support FancyArrowPatch, using ax.annotate instead. The arrows might strike through texts. Increasing shrinkA in arrowprops might help.
            adjust_text(
                texts,
                ax=ax[i],
                ha="center",
                va="center",
                force_text=(4.0, 4.0),  # 增加文本之间的排斥力
                force_points=(2, 2),  # 增加文本与点之间的排斥力
                expand_text=(4, 4),  # 增加文本周围的扩展区域
                expand_points=(1.5, 1.5),  # 增加点周围的扩展区域
                lim=10000,  # 增加迭代次数
                precision=0.0001,  # 提高精度
                arrowprops=dict(
                    arrowstyle="->",
                    color="grey",
                    lw=0.5,
                    shrinkA=4,  # 增加箭头与文本的距离
                    shrinkB=2,  # 增加箭头与点的距离
                ),
            )
        ## add baseline
        xmax = summary_df[x_col].max()
        ymax = summary_df[y_col].max()
        xymax = max(xmax, ymax)
        ax[i].plot(
            [0, xymax],
            [0, xymax],
            color="grey",
            linestyle="--",
            linewidth=0.5,
        )
        ## set title and labels
        ax[i].set_xlabel(f"{column2method[x_col]}")
        ax[i].set_ylabel(f"{column2method[y_col]}")
        ## set x and y ticks
        ax[i].ticklabel_format(
            style="sci",  # 使用科学记数法
            axis="both",  # 应用于x和y轴
            scilimits=(3, 3),  # 当数字达到10^3时触发
            useMathText=True,  # 使用数学文本格式（例如 10³ 而不是 1e3）
        )
        if i != 2:
            ax[i].get_legend().remove()
        # ax[i].get_legend().set_title("Correlation Significance")
    ax[2].legend(
        loc="lower right",
        title="Correlation Significance",
        # fontsize=6,
        # title_fontsize=6,
    )
    plt.suptitle(
        f"Effective Sample Size Comparison for {label_name_shorten[meta_data['id2celltype'][target_qtdid]]}: {meta_data['id2name'][target_qtdid]}",
    )
    plt.tight_layout()
    if text_annot:
        plt.savefig(
            f"{save_path}/f3_neff_scatter_{target_qtdid}_annot.png",
            dpi=300,
            bbox_inches="tight",
        )
    else:
        plt.savefig(
            f"{save_path}/f3_neff_scatter_{target_qtdid}_no_annot.png",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close()


def f3_cor(summary_sign_df, target_qtdid):
    # plot by cor group
    ## 1. num of egene by method, proportion of egene by method
    ## 2. ratio of effective sample size: GMM/Original, GMM+/Original
    ## define cor group
    cor_group_cut = [0.4, 0.6, 0.8]
    cor_group_labels = ["0.0-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]
    summary_sign_df.loc[:, "COR_GROUP"] = pd.cut(
        summary_sign_df.COR,
        bins=[-1] + cor_group_cut + [1.1],
        labels=cor_group_labels,
        include_lowest=True,
    )
    summary_sign_df.COR_GROUP = pd.Categorical(
        summary_sign_df.COR_GROUP, categories=cor_group_labels, ordered=True
    )

    ## define ratio columns
    c2o_ration_name = f"{meta_data['method_name'][1]}/{meta_data['method_name'][0]}"
    t2o_ration_name = f"{meta_data['method_name'][2]}/{meta_data['method_name'][0]}"
    summary_sign_df.loc[:, c2o_ration_name] = (
        summary_sign_df.TAR_CNEFF / summary_sign_df.TAR_SNEFF
    )
    summary_sign_df.loc[:, t2o_ration_name] = (
        summary_sign_df.TAR_TNEFF / summary_sign_df.TAR_SNEFF
    )

    # plot
    fig, ax = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
    ## 1. num of egene by method, proportion of egene by method
    egene_data = pd.DataFrame(
        {
            "GENE": summary_sign_df["GENE"],
            "COR_GROUP": summary_sign_df["COR_GROUP"],
            meta_data["method_name"][0]: summary_sign_df["TAR_SeSNP"] > 0,
            meta_data["method_name"][1]: summary_sign_df["TAR_CeSNP"] > 0,
            meta_data["method_name"][2]: summary_sign_df["TAR_TeSNP"] > 0,
        }
    )

    # 转换为长格式便于绘图
    melted_df = pd.melt(
        egene_data,
        id_vars=["GENE", "COR_GROUP"],
        value_vars=meta_data["method_name"],
        var_name="Method",
        value_name="is_eGene",
    )

    # 只保留是eGene的记录
    melted_df = melted_df[melted_df["is_eGene"]]

    # 统计每个相关性组中每种方法检测到的eGene数量
    plot_df = (
        melted_df.groupby(["COR_GROUP", "Method"], observed=False)
        .size()
        .reset_index(name="Count")
    )

    # 计算每个相关性组的总基因数
    cor_group_counts = (
        summary_sign_df["COR_GROUP"].value_counts().sort_index().reset_index()
    )
    cor_group_counts.columns = ["COR_GROUP", "Total"]

    # 合并总基因数到plot_df
    plot_df = plot_df.merge(cor_group_counts, on="COR_GROUP", how="left")

    # 计算比例
    plot_df["Proportion"] = plot_df["Count"] / plot_df["Total"]

    # 绘制条形图
    sns.barplot(
        x="COR_GROUP",
        y="Count",
        hue="Method",
        data=plot_df,
        ax=ax[0],
        palette=meta_data["Colors"],
        order=cor_group_labels,  # 确保x轴顺序与cor group一致
        hue_order=meta_data["method_name"],
    )

    ax[0].set_xlabel("")  # 不显示x轴标签，因为有共享的x轴
    ax[0].set_ylabel("Number of eGenes")

    # 添加第二个y轴显示比例
    ax0_twin = ax[0].twinx()
    sns.lineplot(
        x="COR_GROUP",
        y="Proportion",
        hue="Method",
        data=plot_df,
        marker="o",
        ax=ax0_twin,
        palette=meta_data["Colors"],
        hue_order=meta_data["method_name"],
        legend=False,  # 避免重复图例
    )

    ax0_twin.set_ylabel("Proportion of eGenes")
    ax0_twin.set_ylim(0, 1)

    # 添加图例
    handles1, labels1 = ax[0].get_legend_handles_labels()
    ax[0].legend(handles1, labels1, title="Method", loc="upper left")

    ## 2. 比例箱线图：展示GMM和GMM+相对于Original的有效样本量比率
    c2o_ration_name = f"{meta_data['method_name'][1]}/{meta_data['method_name'][0]}"
    t2o_ration_name = f"{meta_data['method_name'][2]}/{meta_data['method_name'][0]}"
    this_palette = {
        c2o_ration_name: meta_data["Colors"][meta_data["method_name"][1]],
        t2o_ration_name: meta_data["Colors"][meta_data["method_name"][2]],
    }
    # 转换比率数据为长格式用于绘图
    ratio_melted = pd.melt(
        summary_sign_df,
        id_vars=["COR_GROUP"],
        value_vars=[c2o_ration_name, t2o_ration_name],
        var_name="Ratio",
        value_name="Value",
    )

    # 绘制箱线图，不显示异常值点
    sns.boxplot(
        x="COR_GROUP",
        y="Value",
        hue="Ratio",
        data=ratio_melted,
        fliersize=0,  # 不显示异常值
        ax=ax[1],
        palette=this_palette,
    )

    # 设置y轴范围和标签
    ax[1].set_ylim([-0.5, 5.5])  # 根据实际数据调整上限
    ax[1].set_ylabel("Effective Sample Size Fold Change Ratio")
    ax[1].set_xlabel("Correlation Group")
    ax[1].legend(title="Ratio", loc="upper left")

    # 设置整体标题
    plt.suptitle(
        f"eGene Detection and Effective Sample Size\nin {label_name_shorten[meta_data['id2celltype'][target_qtdid]]}: {meta_data['id2name'][target_qtdid]}",
        fontsize=12,
        y=0.98,
    )

    # 调整布局并保存
    plt.tight_layout()
    plt.savefig(
        os.path.join(save_path, f"f3_cor_{target_qtdid}.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def f3_cor_density(summary_df, target_qtdid):
    # plot cor density and significant cor density
    target_summary_df = summary_df[summary_df.QTDid == target_qtdid].copy()

    # 筛选出显著相关的数据
    significant_df = target_summary_df[target_summary_df.COV_PVAL < 0.05].copy()

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # 绘制所有COR的分布曲线
    sns.kdeplot(
        data=target_summary_df,
        x="COR",
        color="#120eff",  # 蓝色
        fill=True,
        ax=ax,
        alpha=0.5,
        linewidth=1.5,
        label="All genes",
    )

    # 绘制仅显著COR的分布曲线
    sns.kdeplot(
        data=significant_df,
        x="COR",
        color="#db1717f8",  # 红色
        fill=True,
        ax=ax,
        alpha=0.5,
        linewidth=1.5,
        label="Significant genes",
    )

    ax.set_xlabel("Genetic Correlation")
    ax.set_ylabel("Density")
    ax.set_title(
        f"Genetic Correlation Density for {label_name_shorten[meta_data['id2celltype'][target_qtdid]]}: {meta_data['id2name'][target_qtdid]}"
    )

    # 设置图例
    ax.legend(
        title="Type",
        fontsize=10,
        title_fontsize=10,
        loc="upper left",
    )

    plt.tight_layout()
    plt.savefig(
        os.path.join(save_path, f"f3_cor_density_{target_qtdid}.png"),
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    summary_sign_df_all, summary_df_all = load_all_summary()
    for target_qtdid in target_qtdids:
        summary_df = summary_df_all[summary_df_all.QTDid == target_qtdid].copy()
        # QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
        ## filter h2 not significant
        f3_cor_density(summary_df, target_qtdid)
        z_threshold = p2z(0.05)
        summary_df.loc[:, "H1_SIGN"] = (
            summary_df.H1SQ / summary_df.H1SQSE
        ) > z_threshold
        summary_df.loc[:, "H2_SIGN"] = (
            summary_df.H2SQ / summary_df.H2SQSE
        ) > z_threshold
        both_sign_index = summary_df.H1_SIGN & summary_df.H2_SIGN
        summary_df = summary_df.loc[both_sign_index, :]
        ## scatter plot
        f3_neff_scatter(summary_df, target_qtdid)
        f3_neff_scatter(summary_df, target_qtdid, text_annot=False)
        # ## cor group plot
        summary_sign_df = summary_sign_df_all[
            summary_sign_df_all.QTDid == target_qtdid
        ].copy()
        f3_cor(summary_sign_df, target_qtdid)
        print(
            f"Processed {target_qtdid}: {meta_data['id2celltype'][target_qtdid]} - {meta_data['id2name'][target_qtdid]}"
        )
