# fig 5 plot scatter plot for effective sample size of AFR+GTEx
from utils import *
from adjustText import adjust_text

target_qtdids = ["QTD000067", "QTD000371", "QTD000069", "QTD000081"]
study_path_main = "/gpfs1/scratch/wjiang49/xpmm/AFR_GTEx"
save_path = "/gpfs1/scratch/wjiang49/xpmm/AFR_GTEx/results"


def f6_neff_scatter(summary_df, target_qtdid):
    # plot scatter plot of effective sample size
    ## p1: x=TAR_SNEFF, y=TAR_CNEFF, hue=COV_PVAL<0.05
    ## p2: x=TAR_SNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    ## p3: x=TAR_CNEFF, y=TAR_TNEFF, hue=COV_PVAL<0.05
    fig, ax = plt.subplots(1, 3, figsize=(8, 4))
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
        num_annotate_genes = 5
        row_to_annotate = summary_df.loc[
            summary_df.COV_SIGN == "Significant", :
        ].nlargest(num_annotate_genes, "ANNOTATE_ORDER")
        texts = []
        for _, row in row_to_annotate.iterrows():
            gene_name = annot_class.get_gene_name(row.GENE)
            if gene_name is not None:
                # Only annotate if gene name is found and p-value is significant
                texts.append(
                    ax[i].text(
                        row[x_col],
                        row[y_col],
                        gene_name,
                        fontsize=5,
                        fontstyle="italic",  # Add italic style
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
            min_arrow_len=8,  # 增加箭头最小长度
            lim=10000,  # 增加迭代次数
            precision=0.0001,  # 提高精度
            arrowprops=dict(
                arrowstyle="-",
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
        if i != 2:
            ax[i].get_legend().remove()
        # ax[i].get_legend().set_title("Correlation Significance")
    ax[2].legend(
        loc="lower right",
        title="Covariance Significance",
        fontsize=6,
        title_fontsize=6,
    )
    plt.suptitle(
        f"Effective Sample Size Comparison for {meta_data['id2celltype'][target_qtdid]}: {meta_data['id2name'][target_qtdid]}",
    )
    plt.tight_layout()
    plt.savefig(
        f"{save_path}/f6_neff_scatter_{target_qtdid}.png",
        dpi=600,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    _, summary_df_all = load_all_summary(study_path_main)
    for target_qtdid in target_qtdids:
        summary_df = summary_df_all[summary_df_all.QTDid == target_qtdid].copy()
        z_threshold = p2z(0.05)
        summary_df.loc[:, "H1_SIGN"] = (
            summary_df.H1SQ / summary_df.H1SQSE
        ) > z_threshold
        summary_df.loc[:, "H2_SIGN"] = (
            summary_df.H2SQ / summary_df.H2SQSE
        ) > z_threshold
        both_sign_index = summary_df.H1_SIGN & summary_df.H2_SIGN
        summary_df = summary_df.loc[both_sign_index, :]
        f6_neff_scatter(summary_df, target_qtdid)
