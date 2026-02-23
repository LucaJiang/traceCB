#  plot boxplot of cell type proportion
# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from visual.utils import label_name_shorten, meta_data

map_dict = {
    "B_cells": ["B cells naive", "B cells memory"],
    "CD4+T_cells": [
        "T cells CD4 naive",
        "T cells CD4 memory resting",
        "T cells CD4 memory activated",
        "T cells follicular helper",
        "T cells regulatory (Tregs)",
    ],
    "CD8+T_cells": ["T cells CD8"],
    "NK_cells": ["NK cells resting", "NK cells activated"],
    "Monocytes": ["Monocytes"],
}


data_path = "/Users/lucajiang/learn/CityU/xpmm/data/decon/ind_celltype_proportion.csv"
raw_df = pd.read_csv(data_path, header=0)
raw_df.columns = [x.replace(".", " ") for x in raw_df.columns]
raw_df.head()

# %% get celltype level proportion
ct_proportion_df = pd.DataFrame(columns=["Mixture"] + list(map_dict.keys()))
ct_proportion_df["Mixture"] = raw_df["Mixture"]
for ct, ct_list in map_dict.items():
    ct_proportion_df[ct] = raw_df[ct_list].sum(axis=1)
ct_proportion_df = ct_proportion_df.set_index("Mixture")
ct_proportion_df


# %%
# Convert data to long format for plotting
df_long = ct_proportion_df.melt(var_name="Cell Type", value_name="Proportion")
# Use new labels
df_long["Cell Type"] = df_long["Cell Type"].map(label_name_shorten)


# Create boxplot
plt.figure(figsize=(6, 6), dpi=300)
# Update palette to match new short labels
short_name_colors = {
    label_name_shorten[k]: v for k, v in meta_data["celltype_colors"].items()
}
# Create violin plot
ax = sns.violinplot(
    x="Cell Type",
    y="Proportion",
    data=df_long,
    palette=short_name_colors,
    hue="Cell Type",
    cut=0,  # Set the extension degree of the violin plot
    linewidth=1,
)

# Add title and labels
ax.set_title("Cell Type Proportions of GTEx whole blood", fontsize=16)
ax.set_xlabel("Cell Type", fontsize=12)
ax.set_ylabel("Proportion", fontsize=12)

# Rotate x-axis labels to prevent overlap
# plt.xticks(rotation=30, ha="center")
plt.tight_layout()
plt.savefig(
    "/Users/lucajiang/learn/CityU/xpmm/docs/supplymentary/ct_proportion.pdf",
    bbox_inches="tight",
)
plt.show()
# %%
