# %% plot eGene in Nuclear and Membrane Pathway
from utils import *
from matplotlib.patches import Rectangle

# config, gene list obtain from onek1k paper
save_path = "/home/wjiang49/group/wjiang49/data/xpmm/EAS_GTEx/results"
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
pathways_gene_dict = {
    "Nuclear, Cytoplasm or ER": Nuclear_gene_list,
    "Membrane or Surface": Membrane_genename_list
}
OASIS_path = "/home/wjiang49/group/wjiang49/data/hum0197/eQTL_summary_statistics"
OASIS_celltype_dict = {
  "Monocytes": ["CD14_Mono", "CD16_Mono", "Mono", "cMono_IL1B", "cMono_S100A", "intMono", "ncMono"],
  "CD4+T_cells": ["CD4T", "CD4_CTL", "CD4_Naive", "CD4_TCM", "CD4_TEM", "Treg"],
  "CD8+T_cells": ["CD8T", "CD8_CTL", "CD8_Naive", "CD8_TCM", "CD8_TEM"],
  "B_cells": ["B_Activated", "B_Intermediate", "B_Memory", "B_Naive1", "B_Naive2", "B_Naive_OneK1K", "B"],
  "NK_cells": ["NK_CD56bright", "NK_CD56dim", "NK_HLA", "NK_cytokine"]
}
# %%
# target_pathway = "Nuclear, Cytoplasm or ER"
target_pathway = "Membrane or Surface"
pathway_genes = pathways_gene_dict[target_pathway]
# B_Activated_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz
# phenotype_id	gene	variant_id	tss_distance	qval	pval_nominal	slope	slope_se	af	pval_beta	erm	   pval_nominal_threshold	ar	shape1	   beta_shape2	   true_df	true_df	   ma_samples	   ma_count
# ENSG00000237491	LINC01409	chr1_1597761_T_C	819014	0.6232720207072241	0.0001447026366666335	0.5944198	4522	   0.10204082	   0.07681332607755824	   0.07729227077292271	   8.072481827377521e-06	   1.0474911	.99605	      140.8349	      0.000609261769720873		
def get_egene_OASIS(egene_pval_threshold=1e-5):
    result_df = pd.DataFrame(0, index=pathway_genes, columns=["geneid"] + list(OASIS_celltype_dict.keys()))
    result_df["geneid"] = result_df["geneid"].astype(str)
    for celltype, aliases in OASIS_celltype_dict.items():
        for alias in aliases:
            file_path = f"{OASIS_path}/{alias}_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz"
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                continue
            df = pd.read_csv(file_path, sep="\t")
            df = df[df["gene"].isin(pathway_genes)]
            sig_egene = df[df["pval_nominal"] < egene_pval_threshold]["gene"].unique()
            for gene in sig_egene:
                result_df.loc[gene, celltype] = 1
                result_df.loc[gene, "geneid"] = df[df["gene"] == gene]["phenotype_id"].values[0]
    return result_df

oasis_egene_df = get_egene_OASIS()
oasis_egene_df = oasis_egene_df[oasis_egene_df.geneid != "0"]
_, all_summary_df = load_all_summary()
# columns: QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP

def get_egene_type(ns, nc, nt):
    """
    Get eGene type based on SeSNP, CeSNP, TeSNP counts.
    Returns a string indicating the eGene type.
    1: SeSNP>0
    2: SeSNP=0, CeSNP>0
    3: SeSNP=0, CeSNP=0, TeSNP>0
    0: SeSNP=0, CeSNP=0, TeSNP=0
    """
    if ns > 0:
        return 1
    elif nc > 0:
        return 2
    elif nt > 0:
        return 3
    return 0

def get_egene_ours():
    result_df = pd.DataFrame(0, index=oasis_egene_df.index, columns=["geneid"]+meta_data["QTDids"])
    result_df.geneid = oasis_egene_df.geneid
    for qtdid in meta_data["QTDids"]:
        qtd_df = all_summary_df[all_summary_df["QTDid"] == qtdid]
        for row in oasis_egene_df.itertuples():
            gene = row.geneid
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
oasis_egene_df.drop(columns=["geneid"], inplace=True)
# %% remove no info row
remove_row = gene_qtd_df.max(axis=1) <= 1
gene_qtd_df = gene_qtd_df[~remove_row]
# sort row by its sum
gene_qtd_df = gene_qtd_df.loc[gene_qtd_df.sum(axis=1).sort_values(ascending=False).index]

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


name2id = {v: k for k, v in meta_data["id2name"].items()}

# 2. 创建一个与 gene_qtd_df 形状相同的复制标记 DataFrame
replicate_df = pd.DataFrame(0, index=gene_qtd_df.index, columns=gene_qtd_df.columns)
for study_name in replicate_df.columns:
    study_id = name2id[study_name]
    cell_type = meta_data["id2celltype"][study_id]
    if cell_type in oasis_egene_df.columns:
        # 将 oasis_egene_df 中对应细胞类型列的值广播到 replicate_df 的列
        replicate_df[study_name] = oasis_egene_df[cell_type]

# 3. 准备绘图参数
colors_list = ["#edede9"] + list(meta_data["Colors"].values())
label_name_shorten = {
    "Monocytes": "Mono", "CD4+T_cells": "CD4T", "CD8+T_cells": "CD8T",
    "B_cells": "B", "NK_cells": "NK"
}
title = f"eGenes in {target_pathway}"


cmap = mcolors.ListedColormap(colors_list)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# %%
# 4. 绘图
num_genes, num_qtds = gene_qtd_df.shape
cell_size = 0.5
plt.figure(figsize=(num_qtds * cell_size + 3, num_genes * cell_size + 1))

# 设置标记：在OASIS中为1且在我们的研究中>0
markers = np.where((replicate_df == 1) & (gene_qtd_df > 0), "+", "")

ax = sns.heatmap(
    gene_qtd_df, cmap=cmap, norm=norm, cbar=False, linewidths=8, square=True,
    linecolor="white", annot=markers, fmt="",
    annot_kws={"size": 14, "color": "#474646", "ha": "center", "va": "center", "weight": "bold"}
)

# 定义细胞类型分组的范围
celltype_ranges = {
    "Monocytes": (0, 3), "CD4+T_cells": (3, 3), "CD8+T_cells": (6, 2),
    "B_cells": (8, 1), "NK_cells": (9, 1)
}

# 绘制细胞类型分组的背景和顶线
margin = 0.1
for celltype, (x_start, width) in celltype_ranges.items():
    ax.add_patch(
        Rectangle(
            (x_start, 0), width, num_genes,
            facecolor=meta_data["celltype_colors"][celltype], edgecolor="none", alpha=0.15
        )
    )
    ax.hlines(
        y=-0.05, xmin=x_start + margin, xmax=x_start + width - margin,
        colors=meta_data["celltype_colors"][celltype], linewidth=6, clip_on=False, zorder=5
    )
    # 添加细胞类型标签
    ax.text(
        x_start + width / 2, -0.15, label_name_shorten[celltype],
        color="black", fontsize=14, ha="center", va="bottom"
    )

# 5. 创建自定义图例
class SymbolPatch(Patch):
    def __init__(self, symbol="", **kwargs):
        super().__init__(**kwargs)
        self.symbol = symbol

class SquareSymbolHandler(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width, 0.5 * height
        size = min(width, height)
        rect = plt.Rectangle(
            (center[0] - size / 2, center[1] - size / 2), size, size,
            facecolor=orig_handle.get_facecolor(), edgecolor=orig_handle.get_edgecolor(),
            linewidth=orig_handle.get_linewidth(), transform=trans
        )
        symbol = Line2D(
            [center[0]], [center[1]], marker=orig_handle.symbol, color="#474646",
            markersize=8, linestyle="None", transform=trans, markeredgewidth=2
        )
        return [rect, symbol]

legend_elements = [
    SymbolPatch(facecolor=colors_list[1], label=f"Newly Identified by {meta_data['method_name'][0]}"),
    SymbolPatch(facecolor=colors_list[2], label=f"Newly Identified by {meta_data['method_name'][1]}"),
    SymbolPatch(facecolor=colors_list[3], label=f"Newly Identified by {meta_data['method_name'][2]}"),
    SymbolPatch(symbol="+", facecolor="white", edgecolor="white", label="Replicated in OASIS"),
]

plt.legend(
    handles=legend_elements,
    handler_map={Patch: SquareSymbolHandler(), SymbolPatch: SquareSymbolHandler()},
    loc="upper left", bbox_to_anchor=(1.02, 1), title="eGene Type",
    handlelength=1.5, handleheight=1.5, frameon=False
)

# 6. 调整并保存图像
plt.xlabel("Study", fontsize=12)
plt.ylabel("Gene", fontsize=12)
ax.set_yticklabels(ax.get_yticklabels(), style="italic")
plt.title(title, fontsize=16, y=1.15)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10, rotation_mode="anchor")
ax.tick_params(axis='both', which='both', length=0) # 移除ticks
plt.tight_layout(rect=[0, 0, 0.85, 1]) # 为图例留出空间
plt.savefig(f"{save_path}/{title.replace(' ', '_').replace(',', '')}.png", dpi=300, bbox_inches="tight")
print(f"img save to {save_path}/{title.replace(' ', '_').replace(',', '')}.png")
plt.show()

# %%
