# figure 3: case study of single qtd study
# 1. plot scatter plot of effective sample size
# 2. plot by cor group
# 3. plot cor density
# 4. plot eGene upset
from figures.utils import *
from adjustText import adjust_text
from upsetplot import UpSet, plot
from figures.neff_labels import label_neff_scatter
import warnings

warnings.filterwarnings("ignore")

# target_qtdids = ["QTD000081"]
target_qtdids = [
    "QTD000021",
    "QTD000031",
    "QTD000066",
    "QTD000067",
    "QTD000069",
    "QTD000073",
    "QTD000081",
    "QTD000115",
    "QTD000371",
    "QTD000372",
]
summary_output_dir = f"{save_path}/single_study"


def plot_neff_scatter(summary_df, target_qtdid, text_annot=True):
    """Use first_two label ranking only for configured EAS + GTEx analyses."""
    # plot scatter plot of effective sample size
    ## p1: x=TAR_SNEFF, y=TAR_CNEFF, hue=COV_PVAL<0.05
    ## p2: x=TAR_SNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    ## p3: x=TAR_CNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    palette = {
        "Significant": "#db1717f8",
        "Not Significant": "#120eff",
    }
    summary_df.loc[:, "COV_SIGN"] = summary_df.COV_PVAL.apply(
        lambda x: "Significant" if x < 0.10 else "Not Significant"
    )
    ## init annotation class
    annot_class = geneid2name()
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
        # ## add text annotation
        # ### define gene to annotate
        # summary_df.loc[:, "ANNOTATE_ORDER"] = summary_df.loc[
        #     :, y_col
        # ]  # - 1.5*summary_df.loc[:, x_col]
        # if text_annot:
        #     num_annotate_genes = 10
        #     row_to_annotate = summary_df.nlargest(num_annotate_genes, "ANNOTATE_ORDER")
        #     texts = []
        #     for _, row in row_to_annotate.iterrows():
        #         gene_name = annot_class.get_gene_name(row.GENE)
        #         if gene_name is not None and row[x_col] != row[y_col]:
        #             texts.append(
        #                 ax[i].text(
        #                     row[x_col],
        #                     row[y_col],
        #                     gene_name,
        #                     fontsize=7,
        #                     fontstyle="italic",
        #                 )
        #             )
        #     # Looks like you are using a tranform that doesn't support FancyArrowPatch, using ax.annotate instead. The arrows might strike through texts. Increasing shrinkA in arrowprops might help.
        #     adjust_text(
        #         texts,
        #         ax=ax[i],
        #         ha="center",
        #         va="center",
        #         force_text=(4.0, 4.0),  # 增加文本之间的排斥力
        #         force_points=(4, 4),  # 增加文本与点之间的排斥力
        #         expand_text=(4, 4),  # 增加文本周围的扩展区域
        #         expand_points=(1.5, 1.5),  # 增加点周围的扩展区域
        #         lim=10000,  # 增加迭代次数
        #         precision=0.0001,  # 提高精度
        #         # arrowprops=dict(
        #         #     arrowstyle="->",
        #         #     color="grey",
        #         #     lw=0.5,
        #         #     shrinkA=4,  # 增加箭头与文本的距离
        #         #     shrinkB=2,  # 增加箭头与点的距离
        #         # ),
        #     )
        ## add baseline
        xmax = summary_df[x_col].max()
        ymax = summary_df[y_col].max()
        xymax = min(xmax, ymax)
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
        # fontsize=8,
        # title_fontsize=8,
    )
    plt.suptitle("Effective Sample Size Comparison", y=0.92)
    plt.tight_layout()
    if text_annot:
        target_population = os.environ.get("TARGET_POPULATION", "EAS")
        tissue_source = os.environ.get("TISSUE_SOURCE", "eQTLGen")
        ranking_mode = (
            "first_two"
            if (target_population, tissue_source) == ("EAS", "GTEx")
            else "all_three"
        )
        labeled_genes = label_neff_scatter(
            ax,                                   # 原来的三个子图 axes
            summary_df,                           # 原来的完整散点数据
            gene_name=annot_class.get_gene_name,   # 保留原来的基因名映射
            max_genes=5,                          # 最多标注 5 个基因
            label_fontsize=9.0,                    # 先试 9，再试 10、11
            gap_pt=1.5,                            # 标签周围的避让余量
            ranking_mode=ranking_mode,
        )

        # 保存所选基因及三个比值，方便你检查。
        label_dir = os.path.join(save_path, "single_study")
        os.makedirs(label_dir, exist_ok=True)

        labeled_genes.to_csv(
            os.path.join(
                label_dir,
                f"f3neff_scatter_{target_qtdid}_labels.tsv",
            ),
            sep="\t",
            index=False,
        )

        print(
            f"{target_qtdid} ({target_population} + {tissue_source}): "
            f"ranking_mode={ranking_mode}; labeled {len(labeled_genes)} genes"
        )
        print(
            labeled_genes[
                [
                    "GENE_NAME",
                    "ranking_mode",
                    "selection_score",
                    "traceC_Original",
                    "traceCB_Original",
                    "traceCB_traceC",
                    "min_ratio",
                ]
            ].to_string(index=False)
        )
    if text_annot:
        plt.savefig(
            f"{save_path}/single_study/f3neff_scatter_{target_qtdid}_annot.pdf",
            bbox_inches="tight",
        )
    else:
        plt.savefig(
            f"{save_path}/single_study/f3neff_scatter_{target_qtdid}_no_annot.pdf",
            bbox_inches="tight",
        )
    plt.close()


def original_correlations(summary_df):
    """Read finite correlations recorded before clipping the LDSC covariance."""
    if "COR_X_ORI" not in summary_df.columns:
        raise ValueError(
            "Unclipped correlations require COR_X_ORI in summary.csv. "
            "Use results that retain this column or regenerate the summaries; "
            "the original values cannot be recovered from clipped COR_X."
        )
    cor = pd.to_numeric(summary_df["COR_X_ORI"], errors="coerce")
    return cor.where(np.isfinite(cor))


def correlation_groups(summary_df, target_population):
    """Group original correlations, retaining finite values outside [-1, 1]."""
    cuts_by_population = {
        "EAS": [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8],
        "AFR": [-0.8, -0.5, 0.0, 0.5, 0.8],
    }
    if target_population not in cuts_by_population:
        raise ValueError("TARGET_POPULATION must be EAS or AFR")
    cuts = cuts_by_population[target_population]
    cor_group_labels = (
        [f"≤{cuts[0]:.1f}"]
        + [f"{left:.1f}~{right:.1f}" for left, right in zip(cuts, cuts[1:])]
        + [f">{cuts[-1]:.1f}"]
    )
    return pd.cut(
        original_correlations(summary_df),
        bins=[-np.inf, *cuts, np.inf],
        labels=cor_group_labels,
        right=True,
    )


def plot_cor_box(summary_sign_df, target_qtdid, target_population=None):
    if target_population is None:
        target_population = os.environ.get("TARGET_POPULATION", "EAS")
    summary_sign_df = summary_sign_df.copy()
    summary_sign_df["COR_GROUP"] = correlation_groups(
        summary_sign_df, target_population
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
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
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
    ## set negtive value to nan
    ratio_melted.loc[ratio_melted["Value"] < 0, "Value"] = np.nan

    # 绘制箱线图，不显示异常值点
    sns.boxplot(
        x="COR_GROUP",
        y="Value",
        hue="Ratio",
        data=ratio_melted,
        order=summary_sign_df["COR_GROUP"].cat.categories,
        fliersize=0,  # 不显示异常值
        ax=ax,
        palette=this_palette,
    )
    # x轴标签旋转15度
    ax.tick_params(axis="x", rotation=15)

    # 设置y轴范围和标签
    # ax.set_ylim([-0.5, 5.5])  # 根据实际数据调整上限
    # ax.set_ylim([-0.5, 60.5])  # 根据实际数据调整上限
    ax.set_ylim([-0.2, 5.5])  # 根据实际数据调整上限
    ax.set_ylabel("Effective Sample Size Fold Change Ratio")
    ax.set_xlabel("Correlation Group")
    ax.legend(title="Ratio", loc="upper left")

    # 设置整体标题
    # plt.suptitle(
    #     f"Effective Sample Size Fold Change Ratio\nin {label_name_shorten[meta_data['id2celltype'][target_qtdid]]}: {meta_data['id2name'][target_qtdid]}",
    #     fontsize=12,
    #     y=0.98,
    # )

    # 调整布局并保存
    plt.tight_layout()
    plt.savefig(
        os.path.join(save_path, "single_study", f"f3cor_{target_qtdid}.pdf"),
        bbox_inches="tight",
    )
    plt.close()


def plot_cor_density(summary_df, summary_sign_df, target_qtdid):
    # plot cor density and significant cor density
    # QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
    target_summary_df = summary_df.copy()
    target_summary_sign_df = summary_sign_df.copy()

    # 添加分组标签
    target_summary_df.loc[:, "Type"] = "Heritability Significant"
    target_summary_sign_df.loc[:, "Type"] = "Correlation Significant"
    # remove negtive cor
    # target_summary_sign_df = target_summary_sign_df[target_summary_sign_df.COR > 0]

    # 合并数据
    combined_df = pd.concat(
        [target_summary_df, target_summary_sign_df], ignore_index=True
    )
    combined_df["COR_X_ORI"] = original_correlations(combined_df)

    # print(f"Study {target_qtdid}:")
    # print(target_summary_df.COR_X_ORI.describe())
    # print(target_summary_sign_df.COR_X_ORI.describe())

    # 使用 displot 并添加 rug
    g = sns.displot(
        data=combined_df,
        x="COR_X_ORI",
        hue="Type",
        kind="hist",
        bins=200,
        stat="count",
        kde=True,
        # rug=True,
        alpha=0.5,
        palette={
            "Heritability Significant": "#120eff",
            "Correlation Significant": "#db1717f8",
        },
        hue_order=["Heritability Significant", "Correlation Significant"],
        height=4,
        # aspect=0.7,
        aspect=1.0,
    )
    # 设置图例
    sns.move_legend(
        g,
        "upper left",
        bbox_to_anchor=(0.16, 0.85),
        title="Significance Type",
        frameon=True,
        fontsize=8,
        title_fontsize=8,
    )

    g.set_axis_labels("Genetic Correlation", "Count")

    study_name = meta_data["id2name"][target_qtdid]
    celltype = meta_data["id2celltype"][target_qtdid]
    celltype_label = cell_label_name[celltype]
    title = f"{study_name} {celltype_label}"
    g.despine(top=False, right=False)
    g.fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(
        os.path.join(save_path, "single_study", f"f3cor_density_{target_qtdid}.pdf"),
        bbox_inches="tight",
    )
    plt.close()


def summarize_cor_density_by_study(
    summary_df_all, summary_df_h2_sign, summary_sign_df, target_qtdid
):
    cor_series = original_correlations(summary_df_h2_sign)

    return {
        "QTDid": target_qtdid,
        "Study": meta_data["id2name"][target_qtdid],
        "CellType": meta_data["id2celltype"][target_qtdid],
        "Total_Genes": summary_df_all["GENE"].nunique(),
        "Heritability_Significant_Genes": summary_df_h2_sign["GENE"].nunique(),
        "Correlation_Significant_Genes": summary_sign_df["GENE"].nunique(),
        "Median_cor": cor_series.median(),
        "Median_abs_cor": cor_series.abs().median(),
    }


def format_eGene_df(df):
    gene_indicators = pd.DataFrame(
        {
            meta_data["method_name"][0]: df["TAR_SeSNP"] > 0,
            meta_data["method_name"][1]: df["TAR_CeSNP"] > 0,
            meta_data["method_name"][2]: df["TAR_TeSNP"] > 0,
        }
    )
    result = gene_indicators.value_counts().sort_index()
    # remove rows with all False
    result = result.drop((False, False, False))
    return result.to_frame(name="value")


def plot_upset(eGene_df, target_qtdid):
    plt.figure(figsize=(4, 4))
    upset = UpSet(
        eGene_df["value"],
        sort_categories_by="-input",
        min_subset_size=10,
        show_counts=True,
    )
    for i in range(3):
        upset.style_categories(
            meta_data["method_name"][i],
            bar_facecolor=meta_data["Colors"][meta_data["method_name"][i]],
        )
    upset.plot()
    plt.suptitle(f"Number of eGenes")
    plt.savefig(
        f"{save_path}/single_study/f3egene_upset_{target_qtdid}.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


# def plot_cor_density(summary_df, summary_sign_df, target_qtdid):
#     # plot cor density and significant cor density
#     # QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP

#     target_summary_df = summary_df.copy()

#     # 剪切COR值以获得更好的可视化效果
#     target_summary_df.loc[:, "COR"] = target_summary_df["COR"].clip(lower=-10, upper=10)
#     # print(target_summary_df.COR.describe())

#     # 使用 displot
#     g = sns.displot(
#         data=target_summary_df,
#         x="COR",
#         kind="hist",
#         bins=100,
#         stat="count",
#         kde=True,
#         alpha=0.5,
#         color="#120eff",  # 直接使用蓝色
#         height=4,
#         aspect=1.2,
#     )

#     g.set_axis_labels("Genetic Correlation", "Count")
#     g.fig.suptitle(
#         f"Genetic Correlation Density \n{label_name_shorten[meta_data['id2celltype'][target_qtdid]]}: {meta_data['id2name'][target_qtdid]}"
#     )
#     g.despine(top=False, right=False)

#     plt.tight_layout()
#     plt.savefig(
#         os.path.join(save_path, f"single_study/f3_cor_density_h2_sign_only_{target_qtdid}.pdf"),
#         bbox_inches="tight",
#     )
#     plt.close()


if __name__ == "__main__":
    if not os.path.exists(f"{save_path}/single_study/"):
        os.makedirs(f"{save_path}/single_study/")
    if not os.path.exists(summary_output_dir):
        os.makedirs(summary_output_dir)
    summary_sign_df_all, summary_df_all = load_all_summary()
    summary_rows = []
    for target_qtdid in target_qtdids:
        summary_df = summary_df_all[summary_df_all.QTDid == target_qtdid].copy()
        summary_df_all_study = summary_df.copy()
        summary_sign_df = summary_sign_df_all[
            summary_sign_df_all.QTDid == target_qtdid
        ].copy()
        # QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
        ## filter h2 not significant
        z_threshold = p2z(0.1)
        summary_df.loc[:, "H1_SIGN"] = (
            summary_df.H1SQ / summary_df.H1SQSE
        ) > z_threshold
        summary_df.loc[:, "H2_SIGN"] = (
            summary_df.H2SQ / summary_df.H2SQSE
        ) > z_threshold
        both_sign_index = (
            summary_df.H1_SIGN
            & summary_df.H2_SIGN
            & (summary_df.H1SQ > 1e-12)
            & (summary_df.H2SQ > 1e-12)
        )
        summary_df = summary_df.loc[both_sign_index, :]
        summary_rows.append(
            summarize_cor_density_by_study(
                summary_df_all_study, summary_df, summary_sign_df, target_qtdid
            )
        )
        plot_cor_density(summary_df, summary_sign_df, target_qtdid)
        ## scatter plot
        plot_neff_scatter(summary_df, target_qtdid)
        plot_neff_scatter(summary_df, target_qtdid, text_annot=False)
        # ## cor group plot
        plot_cor_box(summary_sign_df, target_qtdid)
        eGene_df = format_eGene_df(summary_df)
        plot_upset(eGene_df, target_qtdid)
        print(
            f"Processed {target_qtdid}: {meta_data['id2celltype'][target_qtdid]} - {meta_data['id2name'][target_qtdid]}"
        )
    summary_table = pd.DataFrame(summary_rows)
    summary_table.to_csv(
        os.path.join(summary_output_dir, "f3cor_density_summary_by_study.tsv"),
        sep="\t",
        index=False,
    )
