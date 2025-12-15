# plot manhattan plot for given gene
from utils import *

# plot_gene_info structure: (gene_id, qtdid, chromosome, gene_name)
plot_gene_info = (
    # ## coloc for bcx_mon to mon
    # ### h4
    # ("ENSG00000159840", "QTD000081", 7, "ZYX"),
    # ("ENSG00000162104", "QTD000081", 16, "ADCY9"),
    # ("ENSG00000076641", "QTD000021", 8, "PAG1"),
    # ("ENSG00000076641", "QTD000081", 8, "PAG1"),
    # ("ENSG00000104783", "QTD000021", 19, "KCNN4"),
    # ("ENSG00000104783", "QTD000081", 19, "KCNN4"),
    # ("ENSG00000178878", "QTD000021", 12, "APOLD1"),
    # ("ENSG00000115919", "QTD000081", 2, "KYNU"),
    # ("ENSG00000076641", "QTD000069", 8, "PAG1"),
    # ("ENSG00000135272", "QTD000021", 7, "MDFIC"),
    # ("ENSG00000135272", "QTD000081", 7, "MDFIC"),
    # ### h3
    # ("ENSG00000166928", "QTD000081", 11, "MS4A14"),
    # ("ENSG00000234745", "QTD000021", 6, "HLA-B"),
    # ("ENSG00000166928", "QTD000069", 11, "MS4A14"),
    # ("ENSG00000171840", "QTD000069", 12, "NINJ2"),
    # ("ENSG00000171840", "QTD000021", 12, "NINJ2"),
    # ("ENSG00000123636", "QTD000021", 2, "BAZ2B"),
    # ("ENSG00000113273", "QTD000021", 5, "ARSB"),
    # ("ENSG00000104518", "QTD000021", 8, "GSDMD"),
    # ("ENSG00000104518", "QTD000069", 8, "GSDMD"),
    # ("ENSG00000104518", "QTD000081", 8, "GSDMD"),
    # ("ENSG00000166928", "QTD000021", 11, "MS4A14"),
    ## heatmap bcx mon
    # ("ENSG00000172543", "QTD000021", 11, "CTSW"),
    # ("ENSG00000172543", "QTD000081", 11, "CTSW"),
    # ("ENSG00000172543", "QTD000069", 11, "CTSW"),
    # ("ENSG00000011454", "QTD000021", 9, "RABGAP1"),
    # ("ENSG00000011454", "QTD000069", 9, "RABGAP1"),
    # ("ENSG00000011454", "QTD000081", 9, "RABGAP1"),
    # ("ENSG00000148700", "QTD000021", 10, "ADD3"),
    # ("ENSG00000148700", "QTD000069", 10, "ADD3"),
    # ("ENSG00000148700", "QTD000081", 10, "ADD3"),
    # ("ENSG00000118369", "QTD000021", 11, "USP35"),
    # ("ENSG00000118369", "QTD000081", 11, "USP35"),
    # ("ENSG00000118369", "QTD000069", 11, "USP35"),
    ### test WDR48   "ENSG00000114742, 3: 39,052,013-39,096,671"
    # ("ENSG00000114742", "QTD000021", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000031", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000066", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000067", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000069", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000073", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000081", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000115", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000371", 3, "WDR48"),
    # ("ENSG00000114742", "QTD000372", 3, "WDR48"),
    # ("ENSG00000079277", "QTD000081", 1, "GeneName"),
    ("ENSG00000138964", "QTD000081", 22, "test"),
)

gwas_path = "/home/wjiang49/group/wjiang49/data/xpmm/coloc/bcx/bcx_mon_GWAS.csv"

MIN_PVAL = 1e-20
MAX_RANGE = 200_000
lookup_table = get_gtex_lookup_table().set_index("RSID")  # RSID, POS


def plot_manhattan(gene_id, qtdid, chromosome, gene_name):
    """
    Plot Manhattan plot for a given gene.
    GWAS+3 method
    """
    # Load eqtl
    eqtl_path = f"{study_path_main}/{qtdid}/GMM/chr{chromosome}/{gene_id}.csv"
    if not os.path.exists(eqtl_path):
        print(f"eQTL file not found: {eqtl_path}")
        return None
    eqtl_df = pd.read_csv(eqtl_path)

    # Load GWAS data
    gwas_df = pd.read_csv(gwas_path)  # SNP,Z,CHR,POS,MAF,BETA,SE

    # Merge data on RSID/SNP
    merged_df = pd.merge(eqtl_df, gwas_df, left_on="RSID", right_on="SNP", how="inner")
    merged_df.loc[:, "PVAL"] = merged_df.Z.apply(z2p)
    # # clip pos if too large
    min_pos = merged_df.POS.min()
    max_pos = merged_df.POS.max()
    if max_pos - min_pos > MAX_RANGE:
        mid_pos = (min_pos + max_pos) / 2
        new_min_pos = mid_pos - MAX_RANGE / 2
        new_max_pos = mid_pos + MAX_RANGE / 2
        merged_df = merged_df[
            (merged_df["POS"] >= new_min_pos) & (merged_df["POS"] <= new_max_pos)
        ]
        print(
            f"Clipped position range for {gene_name} in {qtdid}: {new_min_pos} - {new_max_pos} from {min_pos} - {max_pos}"
        )

    # # convert pos to b38
    # merged_df.loc[:, "POS"] = merged_df.RSID.map(lookup_table.POS)
    # sort by POS
    merged_df.sort_values(by="POS", inplace=True)

    save_name = f"{save_path}/{qtdid}_{gene_id}_{gene_name}_manhattan.png"

    # plot 4 manhattan plots (GWAS + 3 eQTL methods)
    target_methods = ["GWAS", "TAR_SPVAL", "TAR_CPVAL", "TAR_TPVAL"]
    method_labels = [
        "GWAS",
        meta_data["method_name"][0],
        meta_data["method_name"][1],
        meta_data["method_name"][2],
    ]
    # Clip PVAL to avoid log10(0)
    for col in target_methods[1:] + ["PVAL"]:
        merged_df[col] = merged_df[col].clip(lower=MIN_PVAL)

        ## verbose output
        ### print SNP info which has PVAL < 1e-8
        if (merged_df[col] < 1e-8).any():
            print(f"Gene: {gene_name}, SNPs with PVAL < 1e-8 in {col}:")
            for idx, row in merged_df[merged_df[col] < 1e-8].iterrows():
                print(f"  RSID: {row['RSID']}, PVAL: {row[col]}, POS: {row['POS']}")

    fig, axes = plt.subplots(4, 1, figsize=(8, 4), sharex=True, sharey=True)

    # 添加总标题
    fig.suptitle(
        f"Manhattan Plot of {gene_name} in GWAS and {meta_data["id2name"][qtdid]}",
        fontsize=14,
        fontweight="bold",
    )

    for i, method in enumerate(target_methods):
        if method == "GWAS":
            # Plot GWAS data
            log10p = -np.log10(merged_df["PVAL"])
            pos = merged_df["POS"] / 1e6  # 转换为Mb
        else:
            # Plot eQTL data
            log10p = -np.log10(merged_df[method])
            pos = merged_df["POS"] / 1e6  # 转换为Mb

        # 创建不同显著性水平的掩码
        mask_8 = log10p > 8  # p < 1e-8
        mask_5 = (log10p > 5) & (log10p <= 8)  # 1e-8 <= p < 1e-5
        mask_ns = log10p <= 5  # p >= 1e-5

        # 使用不同颜色绘制不同显著性水平的点
        axes[i].scatter(
            pos[mask_ns], log10p[mask_ns], c="gray", s=2, label="p ≥ 1e-5", alpha=0.6
        )
        axes[i].scatter(
            pos[mask_5],
            log10p[mask_5],
            c="blue",
            s=2,
            label="1e-8 ≤ p < 1e-5",
        )
        axes[i].scatter(
            pos[mask_8],
            log10p[mask_8],
            c="red",
            s=2,
            label="p < 1e-8",
        )

        # 添加显著性阈值线
        max_log10p = max(log10p) if len(log10p) > 0 else 0
        if max_log10p > 8:
            axes[i].axhline(y=8, color="red", linestyle="--", alpha=0.5, linewidth=1)
        if max_log10p > 5:
            axes[i].axhline(y=5, color="blue", linestyle="--", alpha=0.5, linewidth=1)

        # 在左上角添加方法标签
        axes[i].text(
            0.02,
            0.90,
            method_labels[i],
            transform=axes[i].transAxes,
            fontsize=8,
            fontweight="bold",
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        axes[i].set_ylabel("-log10(P)")
        axes[i].grid(True, alpha=0.3)
        # 去掉右边和上边的边框
        axes[i].spines["right"].set_visible(False)
        axes[i].spines["top"].set_visible(False)

        # # 只在第一个子图显示图例
        # if i == 0:
        #     axes[i].legend(loc="upper right")

    plt.xlabel(f"Chr {chromosome} Position (Mb)")
    plt.subplots_adjust(hspace=0.01)  # 减少垂直间距
    plt.tight_layout()
    plt.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close()

    return save_name


def plot_aux_manhattan(gene_id, qtdid, chromosome, gene_name):
    """
    Plot auxiliary population Manhattan plot for a given gene.
    """
    # Load eqtl
    celltype = meta_data["id2celltype"][qtdid]
    celltype_eqtl_path = f"{study_path_main}/{qtdid}/AUX_{celltype}/chr{chromosome}.csv"
    tissue_eqtl_path = f"{study_path_main}/{qtdid}/Tissue/chr{chromosome}.csv"
    if not os.path.exists(celltype_eqtl_path):
        print(f"eQTL file not found: {celltype_eqtl_path}")
        return None
    if not os.path.exists(tissue_eqtl_path):
        print(f"Tissue eQTL file not found: {tissue_eqtl_path}")
        return None
    celltype_eqtl_df = pd.read_csv(celltype_eqtl_path)  # RSID,GENE,BETA,SE,Z,PVAL,N
    celltype_eqtl_df = celltype_eqtl_df[celltype_eqtl_df["GENE"] == gene_id]
    tissue_eqtl_df = pd.read_csv(tissue_eqtl_path)
    tissue_eqtl_df = tissue_eqtl_df[tissue_eqtl_df["GENE"] == gene_id]
    if celltype_eqtl_df.empty or tissue_eqtl_df.empty:
        print(
            f"No eQTL data found for gene {gene_id} in {qtdid} at chromosome {chromosome}."
        )
        return None

    # Load GWAS data
    gwas_df = pd.read_csv(gwas_path)  # SNP,Z,CHR,POS,MAF,BETA,SE
    gwas_df = gwas_df[gwas_df["CHR"] == chromosome]
    gwas_df.loc[:, "PVAL"] = gwas_df.Z.apply(z2p)

    # Merge data on RSID/SNP
    merged_df = pd.merge(
        celltype_eqtl_df, tissue_eqtl_df, on="RSID", suffixes=("_celltype", "_tissue")
    )
    merged_df = pd.merge(
        merged_df, gwas_df, left_on="RSID", right_on="SNP", how="inner"
    )
    merged_df.rename(columns={"PVAL": "PVAL_gwas"}, inplace=True)

    save_name = f"{save_path}/{qtdid}_{gene_id}_{gene_name}_aux_manhattan.png"
    # Clip PVAL to avoid log10(0)
    merged_df["PVAL_celltype"] = merged_df["PVAL_celltype"].clip(lower=MIN_PVAL)
    merged_df["PVAL_tissue"] = merged_df["PVAL_tissue"].clip(lower=MIN_PVAL)
    merged_df["PVAL_gwas"] = merged_df["PVAL_gwas"].clip(lower=MIN_PVAL)
    # clip pos if too large
    min_pos = merged_df.POS.min()
    max_pos = merged_df.POS.max()
    if max_pos - min_pos > MAX_RANGE:
        mid_pos = (min_pos + max_pos) / 2
        new_min_pos = mid_pos - MAX_RANGE / 2
        new_max_pos = mid_pos + MAX_RANGE / 2
        merged_df = merged_df[
            (merged_df["POS"] >= new_min_pos) & (merged_df["POS"] <= new_max_pos)
        ]
        print(
            f"Clipped position range for {gene_name} in {qtdid}: {new_min_pos} - {new_max_pos} from {min_pos} - {max_pos}"
        )

    # # convert pos to b38
    # merged_df.loc[:, "POS"] = merged_df.RSID.map(lookup_table.POS)

    # plot 3 manhattan plots (GWAS + celltype + tissue)
    target_methods = ["PVAL_gwas", "PVAL_celltype", "PVAL_tissue"]
    method_labels = [
        "GWAS",
        f"{celltype} eQTL",
        "Tissue eQTL",
    ]
    # Create subplots
    fig, axes = plt.subplots(3, 1, figsize=(8, 4), sharex=True, sharey=True)
    # 添加总标题
    fig.suptitle(
        f"Auxiliary Manhattan Plot of {gene_name} in {meta_data['id2name'][qtdid]}",
        fontsize=14,
        fontweight="bold",
    )
    for i, method in enumerate(target_methods):
        log10p = -np.log10(merged_df[method])
        pos = merged_df["POS"] / 1e6  # 转换为Mb
        # 创建不同显著性水平的掩码
        mask_8 = log10p > 8  # p < 1e-8
        mask_5 = (log10p > 5) & (log10p <= 8)  # 1e-8 <= p < 1e-5
        mask_ns = log10p <= 5  # p >= 1e-5
        # 使用不同颜色绘制不同显著性水平的点
        axes[i].scatter(
            pos[mask_ns], log10p[mask_ns], c="gray", s=2, label="p ≥ 1e-5", alpha=0.6
        )
        axes[i].scatter(
            pos[mask_5],
            log10p[mask_5],
            c="blue",
            s=2,
            label="1e-8 ≤ p < 1e-5",
        )
        axes[i].scatter(
            pos[mask_8],
            log10p[mask_8],
            c="red",
            s=2,
            label="p < 1e-8",
        )
        # 添加显著性阈值线
        max_log10p = max(log10p) if len(log10p) > 0 else 0
        if max_log10p > 8:
            axes[i].axhline(y=8, color="red", linestyle="--", alpha=0.5, linewidth=1)
        if max_log10p > 5:
            axes[i].axhline(y=5, color="blue", linestyle="--", alpha=0.5, linewidth=1)
        # 在左上角添加方法标签
        axes[i].text(
            0.02,
            0.90,
            method_labels[i],
            transform=axes[i].transAxes,
            fontsize=8,
            fontweight="bold",
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
        axes[i].set_ylabel("-log10(P)")
        axes[i].grid(True, alpha=0.3)
        # 去掉右边和上边的边框
        axes[i].spines["right"].set_visible(False)
        axes[i].spines["top"].set_visible(False)
        # # 只在第一个子图显示图例
        # if i == 0:
        #     axes[i].legend(loc="upper right")
    plt.xlabel(f"Chr {chromosome} Position (Mb)")
    plt.subplots_adjust(hspace=0.01)  # 减少垂直
    plt.tight_layout()
    plt.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    for gene_id, qtdid, chromosome, gene_name in plot_gene_info:
        save_name = plot_manhattan(gene_id, qtdid, chromosome, gene_name)
        # plot_aux_manhattan(gene_id, qtdid, chromosome, gene_name)
        print(f"Manhattan plot saved to: {save_name}")
