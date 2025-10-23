# %% plot eGene in Nuclear and Membrane Pathway
from utils import *
from matplotlib.patches import Rectangle

# change the data path and title correspondingly
save_path = "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/other"
egene_path = "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/other/Membrane_egene.csv"
title = "Newly discovered eGenes in\nSurface or Membrane"
# egene_path = "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/other/Nuclear_egene.csv"
# title = "Newly discovered eGenes in\nNuclear, Cytoplasm or ER"


# find newly discovered eGenes
def find_new_egene(egene_df):
    """
    Find newly discovered eGenes in the given DataFrame.
    Returns gene list
    """
    new_egene = []
    for gene in egene_df.index:
        for qtdid in meta_data["QTDids"]:
            qtd_SeSNP = egene_df.loc[gene, f"{qtdid}_SeSNP"]
            qtd_CeSNP = egene_df.loc[gene, f"{qtdid}_CeSNP"]
            qtd_TeSNP = egene_df.loc[gene, f"{qtdid}_TeSNP"]
            if qtd_SeSNP == 0 and (qtd_CeSNP > 0 or qtd_TeSNP > 0):
                new_egene.append(gene)
                break
    return new_egene


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


egene_df = pd.read_csv(egene_path, index_col=0)  # index is gene name
# ,QTD000021_SeSNP,QTD000021_CeSNP,QTD000021_TeSNP,QTD000069_SeSNP,QTD000069_CeSNP,QTD000069_TeSNP,QTD000081_SeSNP,QTD000081_CeSNP,QTD000081_TeSNP,QTD000031_SeSNP,QTD000031_CeSNP,QTD000031_TeSNP,QTD000067_SeSNP,QTD000067_CeSNP,QTD000067_TeSNP,QTD000371_SeSNP,QTD000371_CeSNP,QTD000371_TeSNP,QTD000066_SeSNP,QTD000066_CeSNP,QTD000066_TeSNP,QTD000372_SeSNP,QTD000372_CeSNP,QTD000372_TeSNP,QTD000073_SeSNP,QTD000073_CeSNP,QTD000073_TeSNP,QTD000115_SeSNP,QTD000115_CeSNP,QTD000115_TeSNP,onek1k
# ACTA2, 93, 97, 108, 93, 97, 106, 93, 98, 107, 93, 99, 109, 93, 108, 108, 93, 108, 109, 93, 106, 108, 93.0, 107.0, 108.0, 59, 84, 83, 93, 109, 109, "CD4+T_cells,CD8+T_cells,Monocytes,NK_cells"
new_egene = find_new_egene(egene_df)

print(f"{title}: {len(new_egene)}")
gene_qtd_df = pd.DataFrame(0, index=new_egene, columns=meta_data["QTDids"])
onek1k_replicate_df = gene_qtd_df.copy()
for gene in new_egene:
    onek1k_replicate_celltype = egene_df.loc[gene, "onek1k"].split(",")

    for qtdid in meta_data["QTDids"]:
        gene_qtd_df.loc[gene, qtdid] = get_egene_type(
            egene_df.loc[gene, f"{qtdid}_SeSNP"],
            egene_df.loc[gene, f"{qtdid}_CeSNP"],
            egene_df.loc[gene, f"{qtdid}_TeSNP"],
        )
        qtd_celltype = meta_data["id2celltype"][qtdid]
        if qtd_celltype in onek1k_replicate_celltype:
            onek1k_replicate_df.loc[gene, qtdid] = 1

# rename columns
gene_qtd_df.columns = meta_data["Names"]
onek1k_replicate_df.columns = meta_data["Names"]

colors_list = ["#edede9", "#78F38C", "#ffbf75", "#6c89f9"]
# colors_list = ["#edede9"].append(meta_data["Colors"])
celltype_colors = meta_data["celltype_colors"]
import matplotlib.colors as mcolors

cmap = mcolors.ListedColormap(colors_list)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# 获取热图尺寸
num_genes = len(new_egene)
num_qtds = len(meta_data["QTDids"])
# 设置图形大小
cell_size = 0.5  # 每个格子的大小(英寸)
plt.figure(figsize=(num_qtds * cell_size + 2, num_genes * cell_size + 1))
# 设置标记
markers = np.where((onek1k_replicate_df == 1) & (gene_qtd_df > 0), "+", "")
# 绘制热图
ax = sns.heatmap(
    gene_qtd_df,
    cmap=cmap,
    norm=norm,
    cbar=False,  # 不显示颜色条
    linewidths=8,
    square=True,
    linecolor="white",
    annot=markers,
    fmt="",
    annot_kws={
        "size": 14,
        "color": "#474646",
        "ha": "center",
        "va": "center",
        "weight": "bold",
    },
)

# 定义每个celltype的起始位置和宽度（无间隙全宽度覆盖）
celltype_ranges = {  # (起始索引, 覆盖宽度)
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
            (x_start + margin, 0),  # 左下角坐标
            width - margin * 2,
            num_genes,  # 覆盖全部基因行
            facecolor=celltype_colors[celltype],
            edgecolor="white",  # 去除边框
            linewidth=0.2,
            alpha=0.15,
        )
    )
    x_end = x_start + width  # 计算分组右边界
    ax.hlines(
        y=-0.05,  # 位于标注文字下方0.3单位（根据实际布局调整）
        xmin=x_start + 0.1,  # 左侧留白0.1
        xmax=x_end - 0.1,  # 右侧留白0.1
        colors=celltype_colors[celltype],
        linewidth=6,
        linestyle="-",
        clip_on=False,  # 允许超出绘图区域
        zorder=5,  # 确保在最上层
    )

# 计算标注位置（每个分组的中间x坐标）
label_positions = {
    "Monocytes": (0 + 2.8 / 2, num_genes + 0.5),  # 居中在矩形上方
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
        # rotation=45,  # 斜体防止重叠（可选）
    )

## # 创建自定义图例 -----------------
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


class SymbolPatch(Patch):
    def __init__(self, symbol="", **kwargs):
        super().__init__(**kwargs)
        self.symbol = symbol  # 新增symbol属性存储标记类型


class SquareSymbolHandler(HandlerPatch):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        # 修正颜色块定位逻辑（参照SquareHandler）
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        size = min(width, height)  # 直接取总宽高最小值
        x_left = center[0] - size / 2
        y_bottom = center[1] - size / 2

        # 生成颜色块
        rect = plt.Rectangle(
            (x_left, y_bottom),
            size,
            size,
            facecolor=orig_handle.get_facecolor(),
            edgecolor=orig_handle.get_edgecolor(),
            linewidth=orig_handle.get_linewidth(),
            transform=trans,
        )

        # 修正符号定位逻辑（基于新的颜色块中心）
        symbol_center_x = x_left + size / 2
        symbol_center_y = y_bottom + size / 2

        symbol = Line2D(
            [symbol_center_x],
            [symbol_center_y],
            marker=orig_handle.symbol,
            color="#474646",
            markersize=8,  # 与字体大小关联
            linestyle="None",
            transform=trans,
            markeredgewidth=2,  # 增强符号可见性
        )
        return [rect, symbol]


# 创建自定义图例
legend_elements = [
    SymbolPatch(
        symbol="",
        facecolor=colors_list[1],
        label="Original Study",
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
    SymbolPatch(
        symbol="+",
        facecolor="white",
        edgecolor="white",
        label="Replicated in OneK1K",
    ),
]

# 调整标签
plt.xlabel("Study", fontsize=12)
plt.ylabel("Gene", fontsize=12)
ax.set_yticklabels(ax.get_yticklabels(), style="italic")
plt.title(title, fontsize=16, y=1.1)
plt.legend(
    handles=legend_elements,
    handler_map={Patch: SquareSymbolHandler()},  # 对Patch类使用自定义处理器
    loc="upper left",
    bbox_to_anchor=(1, 1),
    # ncol=2,
    title="Method",
    handlelength=1.5,  # 控制符号的显示长度（正方形边长）
    handleheight=1.5,  # 保持与handlelength相等以确保正方形
    frameon=False,  # 去除图例边框
)

# x tick rotation
ax.set_xticklabels(
    ax.get_xticklabels(), rotation=30, ha="right", fontsize=10, rotation_mode="anchor"
)
ax.tick_params(axis="y")
plt.tight_layout()
plt.savefig(f"{save_path}/{title}.png", dpi=300, bbox_inches="tight")
# plt.show()

# %%
