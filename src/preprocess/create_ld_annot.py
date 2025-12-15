# Create annotation file for calculating LD
## output: chr.print_snps.txt and chr.annot.gz

# chr.annot.gz
# CHR	BP	SNP	CM	base    ENSG00000273443	ENSG00000273487
# 1	731718	rs142557973	0.44130709	1   0   1

# chr.print_snps.txt
# rs142557973

import pandas as pd
import os, argparse, glob


def get_args():
    parser = argparse.ArgumentParser(description="Create LD annotation files")
    parser.add_argument(
        "-s",
        "--studyid",
        type=str,
        required=True,
        help="Study ID for the eQTL catalogue",
    )
    parser.add_argument(
        "-sd",
        "--save_dir_main",
        help="save directory",
        type=str,
    )
    parser.add_argument(
        "-apld",
        "--auxiliary_ld_dir",
        help="auxiliary population LD file path",
        type=str,
    )
    parser.add_argument(
        "-tpld",
        "--target_ld_dir",
        help="target population LD file path",
        type=str,
    )
    return parser.parse_args()


def create_ld_annotation(args, chr):
    studyid = args.studyid
    save_dir_main = args.save_dir_main
    tpld_dir = args.target_ld_dir
    apld_dir = args.auxiliary_ld_dir
    save_path = os.path.join(save_dir_main, studyid, "LDSC", "LD_annotation")
    os.makedirs(save_path, exist_ok=True)

    gene_snp_path = os.path.join(save_dir_main, studyid, "INFO", f"chr{chr}.csv")
    gene_snp = pd.read_csv(gene_snp_path, header=0)

    tpref1000G = pd.read_csv(
        glob.glob(os.path.join(tpld_dir, f"1000G.*.QC.maf.{chr}.bim"))[0],
        header=None,
        sep="\t",
    )
    tpref1000G.columns = ["CHR", "RSID", "CM", "BP", "A1", "A2"]
    apref1000G = pd.read_csv(
        glob.glob(os.path.join(apld_dir, f"1000G.*.QC.maf.{chr}.bim"))[0],
        header=None,
        sep="\t",
    )
    apref1000G.columns = ["CHR", "RSID", "CM", "BP", "A1", "A2"]
    is_common_snp = tpref1000G.RSID.isin(apref1000G.RSID) & tpref1000G.RSID.isin(
        gene_snp.RSID
    )
    is_common_snp = tpref1000G.RSID.isin(apref1000G.RSID) & tpref1000G.RSID.isin(
        gene_snp.RSID
    )

    tp_duplicated_snps = tpref1000G.RSID.duplicated(keep=False)
    tp_dup_ids = set(tpref1000G.loc[tp_duplicated_snps, "RSID"])
    ap_duplicated_snps = apref1000G.RSID.duplicated(keep=False)
    ap_dup_ids = set(apref1000G.loc[ap_duplicated_snps, "RSID"])
    is_not_duplicated = ~tpref1000G.RSID.isin(tp_dup_ids | ap_dup_ids)
    if len(tp_dup_ids) > 0 or len(ap_dup_ids) > 0:
        print(
            f"Warning: Remove {len(tp_dup_ids)} duplicated SNPs in target LD and {len(ap_dup_ids)} in auxiliary LD for chr{chr}."
        )
    is_common_snp_filtered = is_common_snp & is_not_duplicated

    annot_snp = tpref1000G.loc[is_common_snp_filtered, :].sort_values(by="BP")
    unique_gene = gene_snp.GENE.unique()
    snp_gene_matrix = pd.crosstab(gene_snp.RSID, gene_snp.GENE).reindex(
        index=annot_snp.RSID, columns=unique_gene, fill_value=0
    )
    annot = annot_snp.merge(
        snp_gene_matrix, left_on="RSID", right_index=True, how="inner"
    )
    annot.rename(columns={"RSID": "SNP"}, inplace=True)
    annot.insert(loc=4, column="base", value=1)
    annot.loc[:, ["CHR", "BP", "SNP", "CM", "base"] + list(unique_gene)].to_csv(
        os.path.join(save_path, f"{chr}.annot.gz"),
        sep="\t",
        index=False,
        compression="gzip",
    )

    snp_list = pd.DataFrame(annot_snp.RSID)
    snp_list.to_csv(
        os.path.join(save_path, f"{chr}.print_snps.txt"), header=False, index=False
    )

    return snp_list


if __name__ == "__main__":
    args = get_args()
    all_snp_list = pd.DataFrame()
    print(f"Processing {args.studyid}...")
    for chr in range(1, 23):
        chr_snp_list = create_ld_annotation(args, chr)
        all_snp_list = pd.concat([all_snp_list, chr_snp_list], ignore_index=True)

    all_save_path = os.path.join(
        args.save_dir_main,
        args.studyid,
        "LDSC",
        "LD_annotation",
        "all_chroms.print_snps.txt",
    )
    all_snp_list.to_csv(all_save_path, header=False, index=False)
