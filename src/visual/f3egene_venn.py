# plot egene replicate venn between QTDids
## ! Attention: eGene harmonized across QTDids
from utils import *
from matplotlib_venn import venn2, venn3
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
from matplotlib_venn.layout.venn2 import DefaultLayoutAlgorithm as Venn2LayoutAlgorithm


harmonized = True  # whether to harmonize eGenes across QTDids


def plot_venn3(Abc, aBc, abC, ABc, AbC, aBC, ABC, ax, labels):
    v = venn3(
        subsets=(Abc, aBc, abC, ABc, AbC, aBC, ABC),
        set_labels=tuple(labels),
        ax=ax,
        layout_algorithm=DefaultLayoutAlgorithm(
            fixed_subset_sizes=(1, 1, 1, 1, 1, 1, 1)
        ),
    )
    # print(v.set_labels)
    if v.set_labels[0]:  # 左上角标签 - 右移
        pos = v.set_labels[0].get_position()
        v.set_labels[0].set_position((pos[0] + 0.2, pos[1]))  # x坐标右移0.05
    if v.set_labels[1]:  # 右上角标签 - 左移
        pos = v.set_labels[1].get_position()
        v.set_labels[1].set_position((pos[0] - 0.2, pos[1]))  # x坐标左移0.05


def plot_venn2(Ab, aB, AB, ax, labels):
    venn2(
        subsets=(Ab, aB, AB),
        set_labels=tuple(labels),
        ax=ax,
        layout_algorithm=Venn2LayoutAlgorithm(fixed_subset_sizes=(1, 1, 1)),
    )


def get_venn3_count_plot(ct_df, ct_qid, harmonized=True):
    ## get common genes between three QTDids
    common_genes = (
        set(ct_df.loc[ct_df.QTDid == ct_qid[0], "GENE"])
        & set(ct_df.loc[ct_df.QTDid == ct_qid[1], "GENE"])
        & set(ct_df.loc[ct_df.QTDid == ct_qid[2], "GENE"])
    )
    all_genes = set(ct_df.GENE)
    if harmonized:
        unique_genes = all_genes - common_genes
    else:
        unique_genes = set()
    fig, ax = plt.subplots(1, 3, figsize=(10, 3))
    for i, m in enumerate(["S", "C", "T"]):
        # count eGenes for each method
        A_m_eGenes = ct_df.loc[
            (ct_df[f"is_{m}eGene"]) & (ct_df.QTDid == ct_qid[0]), "GENE"
        ]
        B_m_eGenes = ct_df.loc[
            (ct_df[f"is_{m}eGene"]) & (ct_df.QTDid == ct_qid[1]), "GENE"
        ]
        C_m_eGenes = ct_df.loc[
            (ct_df[f"is_{m}eGene"]) & (ct_df.QTDid == ct_qid[2]), "GENE"
        ]
        # get the intersection of eGenes
        Abc = len(
            set(A_m_eGenes) - set(B_m_eGenes) - set(C_m_eGenes) - set(unique_genes)
        )
        aBc = len(
            set(B_m_eGenes) - set(A_m_eGenes) - set(C_m_eGenes) - set(unique_genes)
        )
        abC = len(
            set(C_m_eGenes) - set(A_m_eGenes) - set(B_m_eGenes) - set(unique_genes)
        )
        ABc = len(
            set(A_m_eGenes) & set(B_m_eGenes) - set(C_m_eGenes) - set(unique_genes)
        )
        AbC = len(
            set(A_m_eGenes) & set(C_m_eGenes) - set(B_m_eGenes) - set(unique_genes)
        )
        aBC = len(
            set(B_m_eGenes) & set(C_m_eGenes) - set(A_m_eGenes) - set(unique_genes)
        )
        ABC = len(set(A_m_eGenes) & set(B_m_eGenes) & set(C_m_eGenes))
        # plot the venn diagram
        plot_venn3(
            Abc,
            aBc,
            abC,
            ABc,
            AbC,
            aBC,
            ABC,
            ax[i],
            labels=[meta_data["id2name"][qid] for qid in ct_qid],
        )
        ax[i].set_title(f"{meta_data['method_name'][i]}")
    plt.suptitle(
        f"eGene Replicate between {cell_label_name[meta_data['id2celltype'][ct_qid[0]]]}",
        fontsize=15,
    )
    plt.tight_layout()
    harmonized_suffix = "_harmonized" if harmonized else ""
    plt.savefig(
        os.path.join(
            save_path,
            f"f3egene_replicate_{meta_data['id2celltype'][ct_qid[0]]}{harmonized_suffix}.pdf",
        ),
    )
    print(
        f"Saved: {save_path}/f3egene_replicate_{meta_data['id2celltype'][ct_qid[0]]}{harmonized_suffix}.pdf"
    )


def get_venn2_count_plot(ct_df, ct_qid, harmonized=True):
    # common genes between two QTDids
    common_genes = set(ct_df.loc[ct_df.QTDid == ct_qid[0], "GENE"]) & set(
        ct_df.loc[ct_df.QTDid == ct_qid[1], "GENE"]
    )
    all_genes = set(ct_df.GENE)
    if harmonized:
        unique_genes = all_genes - common_genes
    else:
        unique_genes = set()
    fig, ax = plt.subplots(1, 3, figsize=(10, 3))
    for i, m in enumerate(["S", "C", "T"]):
        # count eGenes for each method
        A_m_eGenes = ct_df.loc[
            (ct_df[f"is_{m}eGene"]) & (ct_df.QTDid == ct_qid[0]), "GENE"
        ]
        B_m_eGenes = ct_df.loc[
            (ct_df[f"is_{m}eGene"]) & (ct_df.QTDid == ct_qid[1]), "GENE"
        ]
        # get the intersection of eGenes
        Ab = len(set(A_m_eGenes) - set(B_m_eGenes) - set(unique_genes))
        aB = len(set(B_m_eGenes) - set(A_m_eGenes) - set(unique_genes))
        AB = len(set(A_m_eGenes) & set(B_m_eGenes))
        # plot the venn diagram
        plot_venn2(
            Ab,
            aB,
            AB,
            ax[i],
            labels=[meta_data["id2name"][qid] for qid in ct_qid],
        )
        ax[i].set_title(f"{meta_data['method_name'][i]}", pad=-20)
    plt.suptitle(
        f"eGene Replicate between {cell_label_name[meta_data['id2celltype'][ct_qid[0]]]}",
        fontsize=15,
    )
    plt.tight_layout()
    harmonized_suffix = "_harmonized" if harmonized else ""
    plt.savefig(
        os.path.join(
            save_path,
            f"f3egene_replicate_{meta_data['id2celltype'][ct_qid[0]]}{harmonized_suffix}.pdf",
        ),
    )
    print(
        f"Saved: {save_path}/f3egene_replicate_{meta_data['id2celltype'][ct_qid[0]]}{harmonized_suffix}.pdf"
    )


def f3egene_replicate():
    """
    Plot Fig. 3 eGene Replicate venn between QTDids
    """
    _, summary_df = load_all_summary()  # ! use significant xpop cor genes
    # QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
    for m in ["S", "C", "T"]:
        summary_df.loc[:, f"is_{m}eGene"] = summary_df.loc[:, f"TAR_{m}eSNP"] > 0
    summary_df.loc[:, "CellType"] = summary_df.QTDid.map(meta_data["id2celltype"])
    # perpare data for venn diagrams: get Abc, aBc, abC, ABc, AbC, aBC, ABC
    ## Monocytes: QTD000021, QTD000069, QTD000081
    ct_df = summary_df.loc[summary_df.CellType == "Monocytes", :]
    ct_qid = sorted(ct_df.QTDid.unique().tolist())
    get_venn3_count_plot(ct_df, ct_qid, harmonized=harmonized)
    ## CD4+T_cells: QTD000031, QTD000067, QTD000371
    cd4_df = summary_df.loc[summary_df.CellType == "CD4+T_cells", :]
    ct_qid = sorted(cd4_df.QTDid.unique().tolist())
    get_venn3_count_plot(cd4_df, ct_qid, harmonized=harmonized)
    ## CD8+T_cells: QTD000066, QTD000372
    cd8_df = summary_df.loc[summary_df.CellType == "CD8+T_cells", :]
    ct_qid = sorted(cd8_df.QTDid.unique().tolist())
    get_venn2_count_plot(cd8_df, ct_qid, harmonized=harmonized)
    return


if __name__ == "__main__":
    f3egene_replicate()
