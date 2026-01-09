# %% visual improved coloc results
from utils import *
import matplotlib.ticker as ticker

oridata = pd.read_csv(
    "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc/coloc_statistics_summary.csv"
)
trait = "coloc_improved"
trait_name = "Other Significant Colocalization Results"
save_path = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc/"
# filename,total_genes,p_original_gt0.7,p_traceC_gt0.7,p_traceCB_gt0.7,pp_h3_single_gt0.7,pp_h3_gmm_cross_gt0.7,pp_h3_gmm_cross_tissue_gt0.7
# filter p_gmm_tissue_gt0 or p_traceC_gt0.7 - p_original_gt0.7 >= 2
oridata.loc[:, "improve"] = oridata.apply(
    lambda x: np.max([x["p_traceCB_gt0.7"], x["p_traceC_gt0.7"]])
    - x["p_original_gt0.7"],
    axis=1,
)
improved_data = oridata[oridata["improve"] >= 3]
print(improved_data)
# improved_data.to_csv(
#     "/Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/coloc/coloc_improved.csv",
#     index=False,
# )

exclude = []
# exclude = ["bcx_mon", "bcx_mcv"]
bcx_index = improved_data["filename"].str.startswith("bcx")
for bcx in exclude:
    bcx_index = bcx_index & ~improved_data["filename"].str.startswith(bcx)
bcx_data = improved_data[bcx_index]

# %%
df_trait = bcx_data.loc[
    :,
    [
        "filename",
        "p_original_gt0.7",
        "p_traceC_gt0.7",
        "p_traceCB_gt0.7",
    ],
].rename(
    columns={
        "p_original_gt0.7": meta_data["method_name"][0],
        "p_traceC_gt0.7": meta_data["method_name"][1],
        "p_traceCB_gt0.7": meta_data["method_name"][2],
    }
)
bcx_study_name_dict = {
    "rbc": "Red Blood Cell Count",
    "hgb": "Hemoglobin Concentration",
    "hct": "Hematocrit",
    "mch": "Mean Corpuscular Hemoglobin",
    "mcv": "Mean Corpuscular Volume",
    "mchc": "Mean Corpuscular Hemoglobin Concentration",
    "rdw": "RBC Distribution Width",
    "wbc": "Total White Blood Cell Count",
    "neu": "Neutrophil Count",
    "lym": "Lymphocyte Count",
    "mon": "Monocyte Count",
    "bas": "Basophil Count",
    "eos": "Eosinophil Count",
    "plt": "Platelet Count",
    "mpv": "Mean Platelet Volume",
}
# red blood cell count (RBC count) hemoglobin concentration (HGB) hematocrit (HCT) mean corpuscular hemoglobin (MCH) mean corpuscular volume (MCV) mean corpuscular hemoglobin concentration (MCHC) RBC distribution width (RDW) total white blood cell count (WBC count) neutrophil count (Neutro) lymphocyte count (Lympho) monocyte count (Mono) basophil count (Baso) eosinophil count (Eosin) platelet count (PLT count) mean platelet volume (MPV)

# filename: bcx_hct_QTD000073_B_cells
df_trait.loc[:, "study"] = (
    df_trait["filename"].str.split("_").str[1].map(bcx_study_name_dict)
)
df_trait.loc[:, "qtdid"] = df_trait["filename"].str.split("_").str[2]
df_trait.loc[:, "qtd_name"] = df_trait.qtdid.map(meta_data["id2name"])
df_trait.loc[:, "celltype"] = df_trait.qtdid.map(meta_data["id2celltype"])
df_trait.loc[:, "coloc_name"] = df_trait.apply(
    lambda x: f"{x['study']} - {label_name_shorten[x['celltype']]}: {x['qtd_name']}",
    axis=1,
)
df_trait.set_index("coloc_name", inplace=True)

# 按最大值排序
df_trait = (
    df_trait.assign(max_value=df_trait[meta_data["method_name"]].max(axis=1))
    .sort_values(by="max_value", ascending=False)
    .drop(columns=["max_value"])
)

# %% 计算堆叠数据
df_stacked = df_trait[meta_data["method_name"]].copy()
df_stacked[meta_data["method_name"][1]] = (
    df_stacked[meta_data["method_name"][1]] - df_stacked[meta_data["method_name"][0]]
)
df_stacked[meta_data["method_name"][2]] = (
    df_stacked[meta_data["method_name"][2]]
    - df_stacked[meta_data["method_name"][1]]
    - df_stacked[meta_data["method_name"][0]]
)
# 设置负值为 0
df_stacked[df_stacked < 0] = 0

# 绘图 - 使用 meta_data 的配色方案
fig, ax = plt.subplots(figsize=(8, 3))
df_stacked.plot(
    kind="barh",
    stacked=True,
    ax=ax,
    color=[meta_data["Colors"][method] for method in meta_data["method_name"]],
    width=0.8,
)

# 设置标题和标签
# ax.set_title(f"Colocalization Results", fontsize=14)
ax.set_xlabel("Number of Colocalized Loci")
ax.set_ylabel("BCX trait - Coloc cell type: Study")
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.legend(title="Method", loc="upper right")

plt.tight_layout()
plt.savefig(
    os.path.join(save_path, f"f5coloc_improved_bar.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()

# %%
