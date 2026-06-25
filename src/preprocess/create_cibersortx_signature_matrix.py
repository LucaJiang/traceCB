# format gtex blood gct file to the input format of CIBERSORTx
import pandas as pd
import numpy as np
import time
import gzip
from sklearn.feature_selection import f_classif

base_path = "~/traceCB/data/deconv/"
rsem_file = "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
annot_file = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
prob_map_file = "gencode.v23.annotation.gene.probemap"
lm22full_file = "LM22full.txt"
save_file = "GTEx_blood_signature_matrix.txt"
# load
annot = pd.read_csv(base_path + annot_file, sep="\t")  # SAMPID  SMNABTCHT
lm22_full = pd.read_csv(lm22full_file, sep="\t")
map_df = pd.read_csv(base_path + prob_map_file, sep="\t", header=0)

target_samples = annot.loc[
    annot.SMNABTCHT == "RNA isolation_PAXgene Blood RNA (Manual)", "SAMPID"
].values
target_columns = ["transcript_id", "gene_id"] + target_samples.tolist()
# 先确定文件中的列名
with gzip.open(base_path + rsem_file, "rt") as f:
    # 跳过前两行
    next(f)
    next(f)
    # 读取标题行
    header = next(f).strip().split("\t")

# 确定要读取的列索引
usecols = []
for i, col in enumerate(header):
    if col in target_samples:
        usecols.append(i)

# 只读取需要的列
rsem = pd.read_csv(
    base_path + rsem_file, sep="\t", compression="gzip", skiprows=2, usecols=usecols
)
# #1.2
# 199324	17382
# transcript_id	gene_id	GTEX-1117F-0226-SM-5GZZ7	GTEX-1117F-0426-SM-5EGHI	GTEX-1117F-0526-SM-5EGHJ	GTEX-1117F-0626-SM-5N9CS	GTEX-1117F-0726-SM-5GIEN
# ENST00000373020.8	ENSG00000000003.14	26.32
print("shape of target RSEM file:", rsem.shape)
print("head of target RSEM file:\n", rsem.head())
#        transcript_id             gene_id  GTEX-111YS-0006-SM-5NQBE  GTEX-1122O-0005-SM-5O99J
# 0  ENST00000373020.8  ENSG00000000003.14                      0.03                      0.03
# 1  ENST00000494424.1  ENSG00000000003.14                      0.00                      0.00
# 2  ENST00000496771.5  ENSG00000000003.14                      0.00

# id	gene	chrom	chromStart	chromEnd	strand
# ENSG00000223972.5	DDX11L1	chr1	11869	14409	+

# annot and merge
rsem_annot = pd.merge(
    rsem, map_df.iloc[:, :2], left_on="gene_id", right_on="id", how="left"
)

lm22full_genes = lm22_full.genesinput.values
rsem_annot_lm22 = rsem_annot[rsem_annot["gene"].isin(lm22full_genes)].reset_index(
    drop=True
)
matrix = lm22_full.values.astype("float64")  # (基因数, 细胞类型数)

# find DE genes for each (sub) cell type and select top genes
num_top_gene_each_ct = 30  # 每个细胞类型选择的基因数量
selected_genes = []
for ct_idx, ct in enumerate(lm22_full.columns):
    mask = np.zeros(matrix.shape[1], dtype=bool)
    mask[ct_idx] = True
    f_vals, _ = f_classif(matrix.T, mask)
    top_indices = np.argsort(f_vals)[-num_top_gene_each_ct:]
    selected_genes.extend(lm22_full.index[top_indices])

selected_genes = list(set(selected_genes))
print(f"Selected genes: {len(selected_genes)}")
filtered_lm22 = lm22_full[lm22_full.index.isin(selected_genes)]
filtered_lm22.to_csv(
    save_file,
    sep="\t",
)
