# Preprocess eQTLCatalogue data

run following python script

```python
import pandas as pd
import os
import numpy as np
from scipy.stats import norm

data_samplesize_dict = {
    "QTD000021": 191,
    "QTD000031": 167,
    "QTD000066": 277,
    "QTD000067": 290,
    "QTD000069": 286,
    "QTD000070": 280,
    "QTD000073": 262,
    "QTD000081": 420,
    "QTD000115": 247,
    "QTD000371": 280,
    "QTD000372": 269,
}

p2z = lambda p: np.abs(norm.ppf(p / 2))
base_path = "/gpfs1/scratch/ResearchGroups/bios_mingxcai/data/eQTLCatalogue/"
save_path = base_path + "by_celltype_chr/"
os.makedirs(save_path, exist_ok=True)

for file_id, sample_size in data_samplesize_dict.items():
    file_path = f"{base_path}{file_id}.all.tsv.gz"
    dtype_spec = {"chromosome": str}
    print(f"Processing {file_id} with sample size {sample_size}")
    os.makedirs(os.path.join(save_path, file_id), exist_ok=True)
    df = pd.read_csv(
        file_path, sep="\t", compression="gzip", header=0, dtype=dtype_spec
    )
    print(f"Loaded {file_id} with shape {df.shape}")
    df = df.loc[df.chromosome != "X", :]
    df.chromosome = df.chromosome.astype(int)
    df = df.loc[df.type == "SNP", :]
    df.loc[:, "n"] = sample_size
    for chr in range(1, 23):
        print(f"Processing chr{chr}")
        df_chr = df.loc[df.chromosome == chr, :]
        # drop duplicates molecular_trait_id
        df_chr = (
            df_chr.sort_values("molecular_trait_id")
            .reset_index(drop=True)
            .drop_duplicates(subset=["rsid", "gene_id"])
        )
        df_chr.loc[:, "z"] = df_chr.beta / df_chr.se
        df_new = df_chr.loc[
            :,
            [
                "chromosome",
                "rsid",
                "gene_id",
                "position",
                "alt",
                "ref",
                "beta",
                "se",
                "pvalue",
                "z",
                "n",
            ],
        ]
        df_new.columns = [
            "CHR",
            "RSID",
            "GENE",
            "POS",
            "A1",
            "A2",
            "BETA",
            "SE",
            "PVAL",
            "Z",
            "N",
        ]
        save_file_name = f"{save_path}{file_id}/chr{chr}.csv"
        df_new.to_csv(save_file_name, index=False)

```
