import os
import glob
import json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import norm

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 1000)
plt.rcParams["font.family"] = "DejaVu Sans"

# study_path_main = "/home/group1/wjiang49/data/traceCB/EAS_GTEx"
# save_path = "/home/group1/wjiang49/data/traceCB/EAS_GTEx/results"
# study_path_main = "/home/group1/wjiang49/data/traceCB/AFR_GTEx"
# save_path = "/home/group1/wjiang49/data/traceCB/AFR_GTEx/results"
study_path_main = "/home/group1/wjiang49/data/traceCB/EAS_eQTLGen"
save_path = "/home/group1/wjiang49/data/traceCB/EAS_eQTLGen/results"

onek1k_path = "/home/wjiang49/group/wjiang49/data/traceCB/onek1k_supp/onek1k_esnp.csv"
gtex_lookup_table_path = "/home/wjiang49/group/wjiang49/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz"

## convert p-value to z-score and vice versa
p2z = lambda p: np.abs(norm.ppf(p / 2))
z2p = lambda z: 2 * norm.sf(abs(z))

label_name_shorten = {
    "Monocytes": "Monocyte",
    "CD4+T_cells": "CD4T",
    "CD8+T_cells": "CD8T",
    "B_cells": "B",
    "NK_cells": "NK",
}

OASIS_path = "/home/wjiang49/group/wjiang49/data/hum0197/eQTL_summary_statistics"
OASIS_celltype_dict = {
    "Monocytes": ["Mono"],
    "CD4+T_cells": ["CD4T"],
    "CD8+T_cells": ["CD8T"],
    "B_cells": ["B"],
    "NK_cells": ["NK"]
}

def load_json(json_file="metadata.json"):
    """
    Load a JSON file and return its content.
    """
    with open(json_file, "r") as f:
        data = json.load(f)
    return data


# json path is same as utils.py path
json_file_path = os.path.join(os.path.dirname(__file__), "metadata.json")
meta_data = load_json(json_file_path)


def load_summary(study_dir):
    """
    Load summary data from GMM results.
    returns summary_sign_df: DataFrame with significant correlations
    returns summary_df: all results DataFrame
    """
    #     <save_path_main>/QTD@/GMM/chr@/summary.csv
    # GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
    summary_dirs = glob.glob(os.path.join(study_dir, "GMM", "chr*", "summary.csv"))
    summary_df = pd.DataFrame()
    for summary_file in summary_dirs:
        df = pd.read_csv(summary_file, header=0, index_col=None)
        df.loc[:, "CHR"] = os.path.basename(os.path.dirname(summary_file)).replace("chr", "")
        summary_df = (
            df if summary_df.empty else pd.concat([summary_df, df], ignore_index=True)
        )
    summary_df.loc[:, "COR_ORI"] = (
        summary_df.COV / summary_df.H1SQ.apply(np.sqrt) / summary_df.H2SQ.apply(np.sqrt)
    )  # convert correlation to genetic correlation, original cor
    summary_df.loc[:, "COR"] = summary_df.COR_ORI.clip(-1, 1)  # clip cor to [-1,1]
    summary_df.loc[summary_df.COV_PVAL > 0.05, "COR"] = 1e-12  # if COR_PVAL>0.05, set COR=0
    summary_sign_df = summary_df.dropna()
    na_raw = summary_df.TAR_CNEFF.isna()
    summary_df.loc[na_raw, "TAR_CNEFF"] = summary_df.loc[na_raw, "TAR_SNEFF"]
    summary_df.loc[na_raw, "TAR_TNEFF"] = summary_df.loc[na_raw, "TAR_SNEFF"]
    summary_df.loc[na_raw, "AUX_CNEFF"] = summary_df.loc[na_raw, "AUX_SNEFF"]
    summary_df.loc[na_raw, "AUX_TNEFF"] = summary_df.loc[na_raw, "AUX_SNEFF"]
    summary_df.loc[na_raw, "TAR_CeSNP"] = summary_df.loc[na_raw, "TAR_SeSNP"]
    summary_df.loc[na_raw, "TAR_TeSNP"] = summary_df.loc[na_raw, "TAR_SeSNP"]
    summary_df.loc[na_raw, "AUX_CeSNP"] = summary_df.loc[na_raw, "AUX_SeSNP"]
    summary_df.loc[na_raw, "AUX_TeSNP"] = summary_df.loc[na_raw, "AUX_SeSNP"]
    return summary_sign_df.copy(), summary_df.copy()



def load_all_summary(study_path_main=study_path_main):
    """
    Load all summary data from GMM results in the main study directory.
    Returns two DataFrames with significant correlations and all results.
    columns: QTDid,GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
    """
    all_qtdids = meta_data["QTDids"]
    all_summary_sign_df = pd.DataFrame()
    all_summary_df = pd.DataFrame()
    for qtdid in all_qtdids:
        study_dir = os.path.join(study_path_main, qtdid)
        summary_sign_df, summary_df = load_summary(study_dir)
        summary_df.loc[:, "QTDid"] = qtdid
        summary_sign_df.loc[:, "QTDid"] = qtdid
        all_summary_sign_df = (
            summary_sign_df
            if all_summary_sign_df.empty
            else pd.concat([all_summary_sign_df, summary_sign_df], ignore_index=True)
        )
        all_summary_df = (
            summary_df
            if all_summary_df.empty
            else pd.concat([all_summary_df, summary_df], ignore_index=True)
        )
    return all_summary_sign_df, all_summary_df

def get_oasis_egene_by_celltype(egene_pval_threshold=5e-3):
    """
    Loads eGenes from OASIS dataset, grouped by cell type.
    Returns a dictionary mapping cell types to a set of eGene IDs.
    """
    oasis_egenes_by_celltype = {cell: set() for cell in OASIS_celltype_dict.keys()}
    for celltype, aliases in OASIS_celltype_dict.items():
        for alias in aliases:
            file_path = f"{OASIS_path}/{alias}_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz"
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                continue
            df = pd.read_csv(file_path, sep="\t")
            # OASIS文件使用'phenotype_id'作为gene id列
            replicated_genes = set(df.loc[df["pval_nominal"] < egene_pval_threshold, "phenotype_id"])
            oasis_egenes_by_celltype[celltype].update(replicated_genes)
    return oasis_egenes_by_celltype

class onek1k_geneid2name(object):

    def __init__(self):
        # Cell_type	Gene_ID	Gene_Ensembl_ID	SNP	Chromosome
        self.onek1k_df = pd.read_csv(
            onek1k_path, usecols=["Gene_ID", "Gene_Ensembl_ID"]
        )
        self.onek1k_df = self.onek1k_df.drop_duplicates(subset=["Gene_ID"])
        self.onek1k_df.rename(
            columns={"Gene_ID": "GENE_NAME", "Gene_Ensembl_ID": "GENE_ID"}, inplace=True
        )

    def get_gene_name(self, gene_id):
        """
        Get gene name from gene ID.
        :param gene_id: Gene ID to look up.
        :return: Gene name if found, otherwise None.
        """
        if isinstance(gene_id, str):
            gene_id = gene_id.strip()
        if not isinstance(gene_id, str):
            return None
        gene_name = self.onek1k_df.loc[self.onek1k_df.GENE_ID == gene_id, "GENE_NAME"]
        return gene_name.iloc[0] if not gene_name.empty else None

    def get_gene_id(self, gene_name):
        """
        Get gene ID from gene name.
        :param gene_name: Gene name to look up.
        :return: Gene ID if found, otherwise None.
        """
        if isinstance(gene_name, str):
            gene_name = gene_name.strip()
        if not isinstance(gene_name, str):
            return None
        gene_id = self.onek1k_df.loc[self.onek1k_df.GENE_NAME == gene_name, "GENE_ID"]
        return gene_id.iloc[0] if not gene_id.empty else None


def get_gtex_lookup_table():
    """
    Load the GTEx lookup table for variant IDs and rsIDs.
    Returns a DataFrame with columns: POS, RSID
    POS based on b38
    """
    lookup_df = pd.read_csv(
        gtex_lookup_table_path,
        sep="\t",
        compression="gzip",
        usecols=["variant_id", "rs_id_dbSNP151_GRCh38p7"],
    )
    # variant_id chr variant_pos ref alt num_alt_per_site rs_id_dbSNP151_GRCh38p7 variant_id_b37
    # chr1_13526_C_T_b38 chr1 13526 C T 1 rs1209314672 1_13526_C_T_b37
    # chr1_13550_G_A_b38 chr1 13550 G A 1 rs554008981 1_13550_G_A_b37

    # variant_id: chr{chromosome}_{position}_{ref}_{alt}_b38
    lookup_df = lookup_df.drop_duplicates(subset=["rs_id_dbSNP151_GRCh38p7"])
    variant_parts = lookup_df["variant_id"].str.split("_", expand=True)
    lookup_df.loc[:, "POS"] = variant_parts[1].astype("int32")  # 位置转为整数
    lookup_df.rename(columns={"rs_id_dbSNP151_GRCh38p7": "RSID"}, inplace=True)
    lookup_df = lookup_df[["POS", "RSID"]]
    return lookup_df


if __name__ == "__main__":
    json_file = load_json(json_file_path)
    for key, value in json_file.items():
        print(f"{key}: {value}")
