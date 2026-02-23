# %% plot fig 4 coloc heatmap for bcx mon
from visual.utils import *

COLOC_THRESHOLD = 0.7
# get annotation of gene id 2 name
cloest_protein_path = "/home/wjiang49/group/wjiang49/data/traceCB/coloc/bcx_mon.closest.protein_coding.bed"
# chr1	2980277	2980277	rs2072732	0	chr1	2985732	3355185	PRDM16_ENSG00000142611	5455
# chr1	9166344	9166344	rs6693258	0	chr1	9160364	9189161	GPR157_ENSG00000180758	0
coloc_results_path_base = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc/"  # bcx_mon_QTD000021_Monocytes_coloc.csv
# chr, gene, start, end, nsnp_eqtl, nsnp_gwas, n_snp_coloc, p_original, p_traceC, p_traceCB
# 1, ENSG00000000460, 169648341, 170163703, 2548, 4992, 550, 0.0950362338748985, 0.0950362338748985, 0.0950362338748985
replicate_path = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc/replicate/bcx_mon_replicate_tissue_coloc.csv"
save_path = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/coloc"
# get gene annotation
annot_df = pd.read_csv(cloest_protein_path, sep="\t", header=None)
annot_df.columns = [
    "chrom",
    "start",
    "end",
    "snp",
    "snp_pos",
    "closest_chrom",
    "closest_start",
    "closest_end",
    "gene_name_id",
    "_",
]
annot_df[["gene_name", "gene_id"]] = annot_df["gene_name_id"].str.split(
    "_", expand=True
)
annot_df = annot_df[["gene_name", "gene_id"]].drop_duplicates()
# %% combine celltype coloc results
combined_df = pd.DataFrame()
for i, qtdid in enumerate(meta_data["QTDids"]):
    celltype = meta_data["Celltypes"][i]
    coloc_results_path = (
        coloc_results_path_base + f"bcx_mon_eQTLGen_{qtdid}_{celltype}_coloc.csv"
    )
    coloc_df = pd.read_csv(coloc_results_path, sep=",")

    for target in ["original", "traceC", "traceCB"]:
        coloc_df[f"{target}"] = coloc_df[f"p_{target}"].apply(
            lambda x: 1 if x > COLOC_THRESHOLD else 0
        )
    coloced_df = coloc_df[coloc_df[["original", "traceC", "traceCB"]].sum(axis=1) > 0]
    annoted_coloced_df = coloced_df[
        ["chr", "gene", "start", "end", "original", "traceC", "traceCB"]
    ].merge(
        annot_df,
        left_on="gene",
        right_on="gene_id",
        how="left",
    )
    annoted_coloced_df["QTDid"] = qtdid
    annoted_coloced_df["Celltype"] = celltype
    annoted_coloced_df.drop(columns=["gene_id"], inplace=True)
    combined_df = (
        annoted_coloced_df
        if combined_df.empty
        else pd.concat([combined_df, annoted_coloced_df], axis=0)
    )
    # print(annoted_coloced_df.shape)
    # print(annoted_coloced_df.head())
    # check number of unique genes
    # print(
    #     "No duplicated genes:",
    #     len(annoted_coloced_df["gene_name"].unique()) == annoted_coloced_df.shape[0],
    # )
print(combined_df.head())
# find unique genes in the combined_df
unique_genes = combined_df["gene_name"].unique()
print("number of unique gene: ", len(unique_genes))  # 59

# summarize newly detected coloced genes for each method across all studies
num_traceC_newly_coloced = combined_df[
    (combined_df["traceC"] == 1) & (combined_df["original"] == 0)
].shape[0]
print(f"traceC newly coloced genes: {num_traceC_newly_coloced}")
num_traceCB_newly_coloced = combined_df[
    (combined_df["traceCB"] == 1) & (combined_df["original"] == 0)
].shape[0]
print(f"traceCB newly coloced genes: {num_traceCB_newly_coloced}")

# get replicate tissue coloced results
replicate_df = pd.read_csv(replicate_path, sep=",")
replicate_coloced = replicate_df[replicate_df["p_h4"] > COLOC_THRESHOLD]
replicate_coloced_annot = replicate_coloced[["gene"]].merge(
    annot_df,
    left_on="gene",
    right_on="gene_id",
    how="left",
)


# %% Define plot setting
# define 3 type of coloc
# 1: original & traceC & traceCB
# 2: traceC & traceCB & !original
# 3: traceCB & !traceC & !original
# Encode coloc types
def get_coloc_type(row):
    """Determine the coloc type for each gene
    1: original & traceC & traceCB
    2: traceC & traceCB & !original
    3: traceCB & !traceC & !original
    0: Other cases or not present
    """
    if row["original"] == 1:
        return 1
    elif row["traceC"] == 1 and row["traceCB"] == 1:
        return 2
    elif row["traceC"] == 0 and row["traceCB"] == 1:
        return 3
    else:
        return 0


# Apply function to get coloc type
combined_df["coloc_type"] = combined_df.apply(get_coloc_type, axis=1)

# Create gene-QTD matrix
gene_qtd_matrix = pd.DataFrame(0, index=unique_genes, columns=meta_data["QTDids"])

# Fill the matrix
for _, row in combined_df.iterrows():
    gene = row["gene_name"]
    qtd = row["QTDid"]
    coloc_type = row["coloc_type"]
    gene_qtd_matrix.loc[gene, qtd] = coloc_type

# Check results
print(gene_qtd_matrix.shape)
print(gene_qtd_matrix.head())
# # %% remove rows that only 1 or 2 no zero
# remove_rows = (gene_qtd_matrix != 0).sum(axis=1) <= 2
# gene_qtd_matrix = gene_qtd_matrix[~remove_rows]

# %% Create custom color map and prepare final plotting data
colors_list = [
    "#edede9",  # grey for background
    meta_data["Colors"][meta_data["method_name"][0]],  # Original
    meta_data["Colors"][meta_data["method_name"][1]],  # Newly Identified
    meta_data["Colors"][meta_data["method_name"][2]],  # Newly Identified
]

celltype_colors = meta_data["celltype_colors"]

# Create custom color map
import matplotlib.colors as mcolors

cmap = mcolors.ListedColormap(colors_list)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Remove rows with mixed 0/1 entirely
heatmap_df = gene_qtd_matrix.loc[
    gene_qtd_matrix.apply(lambda x: sum(x > 1) > 0, axis=1)
]

# Calculate the sum of colocalization types for each gene, for sorting
gene_sums = (heatmap_df > 0).sum(axis=1).sort_values(ascending=False)
sorted_genes = gene_sums.index
heatmap_df = heatmap_df.loc[sorted_genes]

# %% plot
# Get heatmap dimensions
num_genes = len(sorted_genes)
num_qtds = len(meta_data["QTDids"])
# Set figure size
cell_size = 0.5  # Size of each cell (inches)
plt.figure(figsize=(num_qtds * cell_size + 2, num_genes * cell_size + 1))
# Set markers
markers = np.where(heatmap_df == 2, "", np.where(heatmap_df == 3, "", ""))
# Plot heatmap
ax = sns.heatmap(
    heatmap_df,
    cmap=cmap,
    norm=norm,
    cbar=False,  # Do not show color bar
    linewidths=8,
    square=True,
    linecolor="white",
    annot=markers,
    fmt="",
    annot_kws={
        "size": 12,
        "color": "black",
        "ha": "center",
        "va": "center",
    },
)

# Add markers for replicated genes
for i, gene in enumerate(sorted_genes):
    if gene in replicate_coloced_annot["gene_name"].to_list():
        ax.text(
            -0.4, i + 0.5, "\u2605", fontsize=12, va="center", ha="center", color="red"
        )


# Draw background celltype rectangle annotations -------------------------
from matplotlib.patches import Rectangle

# Define start position and width for each celltype (full width coverage without gaps)
celltype_ranges = {  # (start index, coverage width)
    "Monocytes": (0, 3),
    "CD4+T_cells": (3, 3),
    "CD8+T_cells": (6, 2),
    "B_cells": (8, 1),
    "NK_cells": (9, 1),
}

margin = 0.1
for celltype, (x_start, width) in celltype_ranges.items():
    ax.add_patch(
        Rectangle(
            (x_start + margin, 0),  # Bottom-left coordinate
            width - margin * 2,
            num_genes,  # Cover all gene rows
            facecolor=celltype_colors[celltype],
            edgecolor="white",  # Remove border
            linewidth=0.2,
            alpha=0.15,
        )
    )
    x_end = x_start + width  # Calculate group right boundary
    ax.hlines(
        y=-0.05,  # 0.3 units below annotation text (adjusted based on actual layout)
        xmin=x_start + 0.1,  # 0.1 padding on the left
        xmax=x_end - 0.1,  # 0.1 padding on the right
        colors=celltype_colors[celltype],
        linewidth=6,
        linestyle="-",
        clip_on=False,  # Allow exceeding plot area
        zorder=5,  # Ensure it is at the top layer
    )

# Calculate annotation position (center x-coordinate for each group)
label_positions = {
    "Monocytes": (0 + 2.8 / 2, num_genes + 0.5),  # Centered above the rectangle
    "CD4+T_cells": (3.2 + 2.8 / 2, num_genes + 0.5),
    "CD8+T_cells": (6.2 + 1.6 / 2, num_genes + 0.5),
    "B_cells": (8 + 0.8 / 2, num_genes + 0.5),
    "NK_cells": (9 + 0.8 / 2, num_genes + 0.5),
}

for celltype, (x_center, y) in label_positions.items():
    ax.text(
        x_center,
        -0.15,
        label_name_shorten[celltype],
        color="black",
        fontsize=14,
        ha="center",
        va="bottom",
        # rotation=45,  # Italic to prevent overlap (optional)
    )

# plt.subplots_adjust(top=0.9)  # Leave space for top annotation

## # Create custom legend -----------------
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


class SymbolPatch(Patch):
    def __init__(self, symbol="", **kwargs):
        super().__init__(**kwargs)
        self.symbol = symbol  # Add symbol attribute to store marker type


class SquareSymbolHandler(HandlerPatch):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        # Correct color block positioning logic (refer to SquareHandler)
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        size = min(width, height)  # Directly take the minimum of total width/height
        x_left = center[0] - size / 2
        y_bottom = center[1] - size / 2

        # Generate color block
        rect = plt.Rectangle(
            (x_left, y_bottom),
            size,
            size,
            facecolor=orig_handle.get_facecolor(),
            edgecolor=orig_handle.get_edgecolor(),
            linewidth=orig_handle.get_linewidth(),
            transform=trans,
        )

        # Correct symbol positioning logic (based on new color block center)
        symbol_center_x = x_left + size / 2
        symbol_center_y = y_bottom + size / 2

        symbol = Line2D(
            [symbol_center_x],
            [symbol_center_y],
            marker=orig_handle.symbol,
            color="black",
            markersize=4,  # Associated with font size
            linestyle="None",
            transform=trans,
            # markeredgewidth=1.5,  # Enhance symbol visibility
        )
        return [rect, symbol]


# Create custom legend
legend_elements = [
    SymbolPatch(
        symbol="",
        facecolor=colors_list[1],
        label=meta_data["method_name"][0],  # Original
    ),
    SymbolPatch(
        symbol="",
        facecolor=colors_list[2],
        label=f"Newly Identified by {meta_data['method_name'][1]}",
    ),
    SymbolPatch(
        symbol="",
        facecolor=colors_list[3],
        label=f"Newly Identified by {meta_data['method_name'][2]}",
    ),
    Line2D(
        [0],
        [0],
        marker="*",
        color="w",
        markerfacecolor="red",
        markeredgecolor="red",
        markersize=6,
        linestyle="None",
        label="Replicated by Tissue eQTL",
    ),
]

# Adjust labels
plt.xlabel("Study", fontsize=12)
plt.ylabel("Gene", fontsize=12)
plt.title("Colocalization with BCX Monocyte Counts", fontsize=16, y=1.03, x=0.4)
plt.legend(
    handles=legend_elements,
    handler_map={Patch: SquareSymbolHandler()},  # Use custom handler for Patch class
    loc="upper left",
    bbox_to_anchor=(1, 1),
    title="Coloc Types",
    handlelength=1.5,  # Control symbol display length (square side length)
    handleheight=1.5,  # Keep equal to handlelength to ensure square
    frameon=False,  # Remove legend border
)

# x tick labels
locs = plt.xticks()[0]
plt.xticks(
    locs,
    meta_data["Names"],
    rotation=45,
    ha="right",
)
# y label move left
ax.tick_params(axis="y", pad=13)
# Make y tick labels italic
ax.set_yticklabels(ax.get_yticklabels(), fontstyle="italic")
plt.tight_layout()
plt.savefig(f"{save_path}/f4coloc_bcx_mon_heatmap.pdf", dpi=300, bbox_inches="tight")
print(f"Saved: {save_path}/f4coloc_bcx_mon_heatmap.pdf")
# %%
