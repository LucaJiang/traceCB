# plot egene of target pathway: Surface or membrane & Nuclear, cytoplasm or ER
# IBSEP Fig3F
## this code is used to collect data for plotting
from utils import *


onek1k_esnp_path = "/gpfs1/scratch/wjiang49/XeQTL/other/onek1k_esnp.csv"
# Cell_type,Gene_ID,Gene_Ensembl_ID,SNP,Chromosome,Position,SNP_assessed_allele,eSNP_rank,rho_correlation_coefficient,S-statistics,pvalue,qvalue,FDR
Nuclear_gene_list = [
    "ACTA2",
    "AHI1",
    "BACH2",
    "BATF3",
    "CCDC85B",
    "CENPU",
    "CENPW",
    "CTSW",
    "DDX6",
    "DGKQ",
    "DRAP1",
    "ETS1",
    "ETV7",
    "FIBP",
    "GATA3",
    "GPX1",
    "HHEX",
    "IRF7",
    "JAZF1",
    "LBH",
    "LYST",
    "MPHOSPH9",
    "NCKIPSD",
    "NUTF2",
    "PHF5A",
    "PLCL1",
    "PPP5C",
    "PRKCB",
    "RGS14",
    "RPS26",
    "SESN3",
    "SHMT1",
    "SKAP2",
    "SLC2A4RG",
    "SNRPC",
    "SP140",
    "UBASH3A",
    "UBE2L3",
    "ULK3",
    "XBP1",
    "ZFP36L1",
    "ZFP90",
    "ZNF652",
]  # Nuclear, cytoplasm or ER
Membrane_genename_list = [
    "BLK",
    "BTN3A1",
    "CCR6",
    "CD247",
    "CD27",
    "CD37",
    "CD6",
    "CD63",
    "CD83",
    "CLEC2D",
    "CLECL1",
    "CRHR1",
    "CTLA4",
    "DSE",
    "FCRL3",
    "GNG8",
    "GPR18",
    "IFNGR2",
    "IL12RB2",
    "IL18R1",
    "IL2RA",
    "ITGA4",
    "LRRC37A2",
    "LY9",
    "MMEL1",
    "PTGIR",
    "RGS1",
    "SCAMP3",
    "SLC15A2",
    "SLC44A2",
    "TMEM258",
    "TNFRSF14",
    "UBE2D3",
]  # Membrane or surface

qtd_methods = []
for qtdid in meta_data["QTDids"]:
    for method in ["S", "C", "T"]:
        qtd_methods.append(f"{qtdid}_{method}eSNP")
methods = qtd_methods + ["onek1k"]


def genename2geneid(onek1k_esnp_df, genename):
    gene_id = onek1k_esnp_df.loc[
        onek1k_esnp_df["Gene_ID"] == genename, "Gene_Ensembl_ID"
    ]
    if gene_id.empty:
        print(f"Gene {genename} not found in onek1k_esnp_df")
        return None
    return gene_id.values[0]


def get_onek1k_esnp_celltype(onek1k_esnp_df, gene):
    """
    Get cell type of eSNPs for a given gene from the onek1k_esnp_df DataFrame.
    Returns a list of cell types where eSNPs are found for the gene.
    """
    esnp_celltype_subgroup = (
        onek1k_esnp_df.loc[onek1k_esnp_df["Gene_ID"] == gene, "Cell_type"]
        .unique()
        .tolist()
    )
    esnp_cell_type = []
    for cell_type in esnp_celltype_subgroup:
        cell_type_name = meta_data["onek1k2celltype"].get(cell_type, None)
        if cell_type_name and cell_type_name not in esnp_cell_type:
            esnp_cell_type.append(cell_type_name)
    return esnp_cell_type


def get_egene_count_df(summary_df, onek1k_esnp_df, gene_list):
    """
    Generate a DataFrame with eGene counts for each gene in gene_list across all studies and methods.
    """
    all_count_df = pd.DataFrame(0, index=gene_list, columns=methods)
    # set last column as list
    all_count_df["onek1k"] = all_count_df["onek1k"].astype(object)
    for gene in gene_list:
        geneid = genename2geneid(onek1k_esnp_df, gene)
        summary_results = summary_df.loc[
            summary_df["GENE"] == geneid,
            ["QTDid", "TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"],
        ]
        if (
            summary_results.empty or summary_results.iloc[:, [1, 2, 3]].sum().sum() == 0
        ):  # all qtdid eSNP count is 0
            all_count_df = all_count_df.drop(gene)
            continue
        for qtdid in summary_results.QTDid:
            for method in ["S", "C", "T"]:
                count = summary_results.loc[
                    summary_results["QTDid"] == qtdid, f"TAR_{method}eSNP"
                ].item()
                all_count_df.loc[gene, f"{qtdid}_{method}eSNP"] = count
        # Add onek1k eSNP
        onek1k_esnp_celltype = get_onek1k_esnp_celltype(onek1k_esnp_df, gene)
        all_count_df.loc[gene, "onek1k"] = (
            ",".join(onek1k_esnp_celltype) if onek1k_esnp_celltype else 0
        )
    return all_count_df


def plot_egene_count(all_count_df, title, save_prefix):
    # create plot dataframe
    plot_df = all_count_df.copy()
    plot_df.onek1k = plot_df.onek1k.apply(lambda x: str(x).count(",") + 1 if x else 0)
    # set as int, fill NA values with 0 first
    plot_df = plot_df.fillna(0).astype(int)
    ## notes: 0 could be unexist gene or no eSNP
    fig, ax = plt.subplots(figsize=(16, 12))
    sns.heatmap(
        plot_df,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar_kws={"label": "eSNP count"},
        linewidths=0.5,
        linecolor="black",
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("Study and method")
    ax.set_ylabel("Gene")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(save_prefix + ".png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    onek1k_esnp_df = pd.read_csv(onek1k_esnp_path)
    gene_categories = [
        {
            "name": "Nuclear",
            "genes": Nuclear_gene_list,
            "title": "Nuclear, Cytoplasm or ER: eSNP count",
        },
        {
            "name": "Membrane",
            "genes": Membrane_genename_list,
            "title": "Surface or Membrane: eSNP count",
        },
    ]
    for category in gene_categories:
        _, summary_df = load_all_summary()
        all_count_df = get_egene_count_df(summary_df, onek1k_esnp_df, category["genes"])
        csv_path = os.path.join(save_path, f"{category['name']}_egene.csv")
        all_count_df.to_csv(csv_path, index=True)
        # Plot eSNP count
        plot_path = os.path.join(save_path, f"{category['name']}_egene")
        plot_egene_count(all_count_df, category["title"], plot_path)

        print(f"完成处理 {category['name']} 基因组")
