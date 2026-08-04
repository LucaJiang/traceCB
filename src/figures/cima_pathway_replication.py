# %% plot eGene in Nuclear and Membrane Pathway with CIMA replication
from figures.utils import *
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

save_path = os.path.join(save_path, "pathway_CIMA")
os.makedirs(save_path, exist_ok=True)

CIMA_root = os.environ.get("TRACECB_CIMA_DIR", str(REPO_ROOT / "data/replication/CIMA"))
CIMA_lead_cis_xqtl_path = os.path.join(CIMA_root, "xQTL", "CIMA_Lead_cis-xQTL.csv")
CIMA_celltype_level_path = os.path.join(
    CIMA_root, "Cell_Atlas", "CIMA_Cell_Type_Level_and_Marker.xlsx"
)

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
]
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
]
pathways_gene_dict = {
    "Nuclear, Cytoplasm or ER": Nuclear_gene_list,
    "Membrane or Surface": Membrane_genename_list,
}

manual_annot_gene_dict = {
    "CD247": "ENSG00000198821",
    "CRHR1": "ENSG00000120088",
    "CTLA4": "ENSG00000163599",
    "GNG8": "ENST00000300873",
    "LRRC37A2": "ENSG00000277221",
    "BACH2": "ENSG00000112182",
    "CCDC85B": "ENSG00000175602",
    "CENPW": "ENSG00000203760",
    "DRAP1": "ENSG00000175550",
    "ETS1": "ENSG00000134954",
    "GATA3": "ENSG00000107485",
    "GPX1": "ENSG00000233276",
    "NUTF2": "ENSG00000102898",
    "PHF5A": "ENSG00000100410",
    "ZNF652": "ENSG00000198740",
}


def map_cima_l4_to_major_celltype(celltype):
    """
    Map CIMA L4 labels to the 5 broad immune groups used in the OASIS plot.

    We keep only subtypes comparable to the purified populations in the original
    replication panel and exclude mixed / adjacent populations such as NKT,
    MAIT, gamma-delta T, dendritic cells, plasma cells, megakaryocytes, HSPCs,
    and ILCs.
    """
    if not isinstance(celltype, str):
        return None
    if celltype in {"CD4", "CD4T"} or celltype.startswith("CD4_"):
        return "CD4+T_cells"
    if celltype in {"CD8", "CD8T"} or celltype.startswith("CD8_"):
        return "CD8+T_cells"
    if celltype in {"B"} or celltype.startswith(
        (
            "Bn_",
            "Transitional_B_",
            "Switched_Bm_",
            "Unswitched_Bm_",
            "pre-Switched_Bm_",
            "Atypical_Bm_",
        )
    ):
        return "B_cells"
    if celltype in {"Mono", "Monocyte"} or celltype.startswith(
        ("cMono_", "ncMono_", "intMono_")
    ):
        return "Monocytes"
    if celltype in {"NK"} or celltype.startswith(
        (
            "Mature_NK_",
            "Terminal_NK_",
            "Transitional_NK_",
            "NK_bright_",
            "Inflamed_NK_",
            "Proliferative_NK_",
        )
    ):
        return "NK_cells"
    return None


def validate_cima_celltype_hierarchy():
    hierarchy_df = pd.read_excel(CIMA_celltype_level_path)
    print("CIMA L1 hierarchy:", sorted(hierarchy_df["L1"].dropna().unique()))


def load_cima_lead_eqtl():
    if not os.path.exists(CIMA_lead_cis_xqtl_path):
        raise FileNotFoundError(f"CIMA lead cis-xQTL file not found: {CIMA_lead_cis_xqtl_path}")
    cols = ["phenotype_id", "celltype", "pval_nominal", "analysis"]
    cima_df = pd.read_csv(CIMA_lead_cis_xqtl_path, usecols=cols)
    cima_df = cima_df[cima_df["analysis"] == "cis-eQTL"].copy()
    cima_df.loc[:, "broad_celltype"] = cima_df["celltype"].map(map_cima_l4_to_major_celltype)
    print("CIMA broad cell type counts:")
    print(cima_df["broad_celltype"].value_counts(dropna=False).to_string())
    return cima_df[cima_df["broad_celltype"].notna()].copy()


def get_egene_type(ns, nc, nt):
    """
    1: SeSNP>0
    2: SeSNP=0, CeSNP>0
    3: SeSNP=0, CeSNP=0, TeSNP>0
    0: SeSNP=0, CeSNP=0, TeSNP=0
    """
    if ns > 0:
        return 1
    if nc > 0:
        return 2
    if nt > 0:
        return 3
    return 0


validate_cima_celltype_hierarchy()
cima_eqtl_df = load_cima_lead_eqtl()
gene_converter = geneid2name()
gene_name_to_id = {
    gene: manual_annot_gene_dict.get(gene, gene_converter.get_gene_id(gene))
    for genes in pathways_gene_dict.values()
    for gene in genes
}

_, all_summary_df = load_all_summary()
name2id = {v: k for k, v in meta_data["id2name"].items()}

# %%
for target_pathway, pathway_genes in pathways_gene_dict.items():

    def get_egene_CIMA(egene_pval_threshold=5e-3):
        result_df = pd.DataFrame(
            0,
            index=pathway_genes,
            columns=["geneid"] + list(OASIS_celltype_dict.keys()),
        )
        result_df["geneid"] = result_df.index.map(gene_name_to_id)
        pathway_df = cima_eqtl_df[cima_eqtl_df["phenotype_id"].isin(pathway_genes)].copy()
        for celltype in OASIS_celltype_dict.keys():
            ct_df = pathway_df[pathway_df["broad_celltype"] == celltype]
            if ct_df.empty:
                continue
            min_pvals = ct_df.groupby("phenotype_id")["pval_nominal"].min()
            replicated_genes = min_pvals[min_pvals < egene_pval_threshold].index
            result_df.loc[result_df.index.isin(replicated_genes), celltype] = 1
        return result_df

    def get_cima_gene_availability():
        """
        Mark genes that appear in the CIMA replicate summary for each broad cell type.

        Important: with the currently downloaded CIMA resource, "missing" means the gene
        is absent from the available lead cis-eQTL summary for that cell type. This is a
        practical proxy for "not available in replicate data" rather than a strict raw
        sequencing-status call.
        """
        result_df = pd.DataFrame(
            0,
            index=pathway_genes,
            columns=list(OASIS_celltype_dict.keys()),
        )
        pathway_df = cima_eqtl_df[cima_eqtl_df["phenotype_id"].isin(pathway_genes)].copy()
        for celltype in OASIS_celltype_dict.keys():
            ct_df = pathway_df[pathway_df["broad_celltype"] == celltype]
            if ct_df.empty:
                continue
            available_genes = set(ct_df["phenotype_id"].unique())
            result_df.loc[result_df.index.isin(available_genes), celltype] = 1
        return result_df

    cima_egene_df_1e5 = get_egene_CIMA(egene_pval_threshold=1e-5)
    cima_egene_df_5e3 = get_egene_CIMA(egene_pval_threshold=5e-3)
    cima_gene_available_df = get_cima_gene_availability()

    def get_egene_ours():
        result_df = pd.DataFrame(
            0, index=cima_egene_df_1e5.index, columns=["geneid"] + meta_data["QTDids"]
        )
        result_df["geneid"] = cima_egene_df_1e5["geneid"]
        for qtdid in meta_data["QTDids"]:
            qtd_df = all_summary_df[all_summary_df["QTDid"] == qtdid]
            for row in cima_egene_df_1e5.itertuples():
                gene = row.geneid
                if not isinstance(gene, str):
                    continue
                if gene in qtd_df["GENE"].values:
                    ns = qtd_df[qtd_df["GENE"] == gene]["TAR_SeSNP"].values[0]
                    nc = qtd_df[qtd_df["GENE"] == gene]["TAR_CeSNP"].values[0]
                    nt = qtd_df[qtd_df["GENE"] == gene]["TAR_TeSNP"].values[0]
                    result_df.loc[row.Index, qtdid] = get_egene_type(ns, nc, nt)
        return result_df

    gene_qtd_df = get_egene_ours()
    for column in gene_qtd_df.columns[1:]:
        gene_qtd_df.rename(columns={column: meta_data["id2name"][column]}, inplace=True)
    gene_qtd_df.drop(columns=["geneid"], inplace=True)
    cima_egene_df_1e5.drop(columns=["geneid"], inplace=True)
    cima_egene_df_5e3.drop(columns=["geneid"], inplace=True)

    remove_row = gene_qtd_df.max(axis=1) <= 1
    gene_qtd_df = gene_qtd_df[~remove_row]
    gene_qtd_df = gene_qtd_df.loc[
        gene_qtd_df.sum(axis=1).sort_values(ascending=False).index
    ]

    replicate_df = pd.DataFrame(0, index=gene_qtd_df.index, columns=gene_qtd_df.columns)
    replicate_available_df = pd.DataFrame(
        0, index=gene_qtd_df.index, columns=gene_qtd_df.columns
    )
    for study_name in replicate_df.columns:
        study_id = name2id[study_name]
        cell_type = meta_data["id2celltype"][study_id]
        replicate_df[study_name] = cima_egene_df_5e3[cell_type]
        replicate_df.loc[cima_egene_df_1e5[cell_type] > 0, study_name] = 2
        replicate_available_df[study_name] = cima_gene_available_df[cell_type]

    colors_list = ["#edede9"] + list(meta_data["Colors"].values())
    replicate_symbol_color = "#982536"
    missing_symbol_color = "#000000"
    title = f"eGenes of {target_pathway} (CIMA)"
    cmap = mcolors.ListedColormap(colors_list)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    num_genes, num_qtds = gene_qtd_df.shape
    cell_size = 0.5
    plt.figure(figsize=(num_qtds * cell_size + 3, num_genes * cell_size + 1))

    ax = sns.heatmap(
        gene_qtd_df,
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=6,
        square=True,
        linecolor="white",
    )

    for y, gene in enumerate(gene_qtd_df.index):
        for x, study_name in enumerate(gene_qtd_df.columns):
            if gene_qtd_df.iat[y, x] <= 0:
                continue
            if replicate_available_df.iat[y, x] == 0:
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    "x",
                    ha="center",
                    va="center",
                    color=missing_symbol_color,
                    fontsize=10,
                    fontweight=800,
                )
            elif replicate_df.iat[y, x] == 2:
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    "++",
                    ha="center",
                    va="center",
                    color=replicate_symbol_color,
                    fontsize=8,
                    fontweight=800,
                )
            elif replicate_df.iat[y, x] == 1:
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    "+",
                    ha="center",
                    va="center",
                    color=replicate_symbol_color,
                    fontsize=8,
                    fontweight=800,
                )

    celltype_ranges = {
        "Monocytes": (0, 3),
        "CD4+T_cells": (3, 3),
        "CD8+T_cells": (6, 2),
        "B_cells": (8, 1),
        "NK_cells": (9, 1),
    }

    margin = 0.10
    for celltype, (x_start, width) in celltype_ranges.items():
        ax.add_patch(
            Rectangle(
                (x_start, 0),
                width,
                num_genes,
                facecolor=meta_data["celltype_colors"][celltype],
                edgecolor="none",
                alpha=0.15,
            )
        )
        ax.hlines(
            y=-0.05,
            xmin=x_start + margin,
            xmax=x_start + width - margin,
            colors=meta_data["celltype_colors"][celltype],
            linewidth=6,
            clip_on=False,
            zorder=5,
        )
        ax.text(
            x_start + width / 2,
            -0.15,
            label_name_shorten[celltype],
            color="#000000",
            fontsize=14,
            ha="center",
            va="bottom",
        )

    class SymbolPatch(Patch):
        def __init__(self, symbol="", symbol_color="#000000", **kwargs):
            super().__init__(**kwargs)
            self.symbol = symbol
            self.symbol_color = symbol_color

    class SquareSymbolHandler(HandlerPatch):
        def create_artists(
            self,
            legend,
            orig_handle,
            xdescent,
            ydescent,
            width,
            height,
            fontsize,
            trans,
        ):
            center = 0.5 * width, 0.5 * height
            size = min(width, height)
            rect = plt.Rectangle(
                (center[0] - size / 2, center[1] - size / 2),
                size,
                size,
                facecolor=orig_handle.get_facecolor(),
                edgecolor=orig_handle.get_edgecolor(),
                linewidth=orig_handle.get_linewidth(),
                transform=trans,
            )
            artists = [rect]
            if orig_handle.symbol:
                text_artist = plt.Text(
                    center[0],
                    center[1],
                    orig_handle.symbol,
                    ha="center",
                    va="center",
                    color=orig_handle.symbol_color,
                    fontsize=12,
                    fontweight=800,
                    transform=trans,
                )
                artists.append(text_artist)
            return artists

    legend_elements = [
        SymbolPatch(facecolor=colors_list[1], label="Identified by BBJ"),
        SymbolPatch(
            facecolor=colors_list[2],
            label=f"Newly Identified by {meta_data['method_name'][1]}",
        ),
        SymbolPatch(
            facecolor=colors_list[3],
            label=f"Newly Identified by {meta_data['method_name'][2]}",
        ),
        SymbolPatch(
            symbol="+",
            facecolor="white",
            edgecolor="white",
            label="Replicated in CIMA (p<5e-3)",
            symbol_color=replicate_symbol_color,
        ),
        SymbolPatch(
            symbol="++",
            facecolor="white",
            edgecolor="white",
            label="Replicated in CIMA (p<1e-5)",
            symbol_color=replicate_symbol_color,
        ),
        SymbolPatch(
            symbol="x",
            facecolor="white",
            edgecolor="white",
            label="Not present in CIMA replicate summary",
            symbol_color=missing_symbol_color,
        ),
    ]

    plt.legend(
        handles=legend_elements,
        handler_map={Patch: SquareSymbolHandler(), SymbolPatch: SquareSymbolHandler()},
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        title="eGene Type",
        handlelength=1.5,
        handleheight=1.5,
        frameon=False,
    )

    plt.ylabel("Gene", fontsize=12)
    ax.set_yticklabels(ax.get_yticklabels(), style="italic")
    plt.title(title, fontsize=14, y=1.1)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=10,
        rotation_mode="anchor",
    )
    ax.tick_params(axis="both", which="both", length=0)
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    out_path = f"{save_path}/{title.replace(' ', '_').replace(',', '')}.pdf"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"img save to {out_path}")
    plt.show()

# %%
