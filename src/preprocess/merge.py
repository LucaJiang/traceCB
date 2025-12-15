# Merge data in chr level
# remove palindromic, mismatched snps
# find common snps, genes
# align other to target pop(bbj/afb)

#  cmd: python -m pdb ~/xpmm/src/mergy.py -tpd /gpfs1/scratch/ResearchGroups/bios_mingxcai/data/BBJ_eQTL/by_celltype_chr -tpn 98 -apd /gpfs1/scratch/ResearchGroups/bios_mingxcai/data/eQTLCatalogue/by_celltype_chr -tid /gpfs1/scratch/ResearchGroups/bios_mingxcai/data/GTEx/GTEx_Whole_Blood_by_chr -tin 982 -tpld /gpfs1/scratch/ResearchGroups/bios_mingxcai/data/1000G/1000G_EAS_EUR/EAS -apld /gpfs1/scratch/ResearchGroups/bios_mingxcai/data/1000G/1000G_EAS_EUR/EUR -s QTD000066 -t CD8+T_cells -c 22 -sd /gpfs1/scratch/wjiang49/xpmm/EAS_GTEx

import numpy as np
import pandas as pd
import os, glob, time
import argparse, warnings
from scipy.stats import norm

warnings.filterwarnings("ignore")

z2p = lambda z: 2 * norm.sf(abs(z))
is_palindromic = lambda A11, A12: A11 + A12 in ["AT", "TA", "CG", "GC"]
is_palindromic_vectorized = np.vectorize(is_palindromic)


def arg_parser():
    parser = argparse.ArgumentParser(description="Prepare and align data")
    parser.add_argument(
        "-tpd",
        "--target_dir",
        help="target population file path",
        type=str,
    )
    parser.add_argument(
        "-tpn",
        "--target_samplesize",
        help="target population sample size",
        type=int,
    )
    parser.add_argument(
        "-apd",
        "--auxiliary_dir",
        help="auxiliary population file path",
        type=str,
    )
    parser.add_argument(
        "-tid",
        "--tissue_dir",
        help="tissue file path",
        type=str,
    )
    parser.add_argument(
        "-tin",
        "--tissue_samplesize",
        help="tissue population sample size",
        type=int,
    )
    parser.add_argument(
        "-tpld",
        "--target_ld_dir",
        help="target population LD file path",
        type=str,
    )
    parser.add_argument(
        "-apld",
        "--auxiliary_ld_dir",
        help="auxiliary population LD file path",
        type=str,
    )
    parser.add_argument(
        "-s",
        "--studyid",
        help="eqtlCatalogue study id",
        type=str,
        default="QTD000066",
    )
    parser.add_argument(
        "-t",
        "--celltype",
        help="cell type name",
        type=str,
        default="CD8+T_cells",
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        nargs="+",
        help="chromosome number",
        type=int,
        default=22,
    )
    parser.add_argument(
        "-sd",
        "--save_dir_main",
        help="save directory",
        type=str,
    )
    return parser.parse_args()


# NAME FORMAT: CHR, SNP, GENE, POS, A1, A2, BETA, SE, Z, PVAL, N
def load_eQTLGen(path, _nouse, chr):
    # INPUT: chr : int
    # OUTPUT: df, unique_snps, unique_genes
    raw_df = pd.read_csv(
        os.path.join(
            path, "chr" + str(chr) + ".tsv.gz"
        ),  # e.g., /path/to/eQTLGen/chr1.tsv.gz
        sep="\t",
        header=0,
        compression="gzip",
    )
    # Pvalue	SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR	BonferroniP
    # 3.2717E-310	rs649539	1	109706393	A	G	-150.1966	ENSG00000116299	KIAA1324	1	109702851	35	30606	.0	4.1662E-302
    # 3.2717E-310	rs647294	1	109706880	A	G	-149.8928	ENSG00000116299	KIAA1324	1	109702851	34	30484	.0	4.1662E-302
    raw_df.rename(
        columns={
            "SNPChr": "CHR",
            "SNPPos": "POS",
            "AssessedAllele": "A1",
            "OtherAllele": "A2",
            "Zscore": "Z",
            "Pvalue": "PVAL",
            "NrSamples": "N",
            "Gene": "GENE",
            "SNP": "RSID",
        },
        inplace=True,
    )
    # filter out palindromic snps
    df = raw_df.loc[
        ~is_palindromic_vectorized(raw_df["A1"], raw_df["A2"]),
        ["CHR", "POS", "RSID", "GENE", "A1", "A2", "Z", "PVAL", "N"],
    ]
    df.loc[:, "SE"] = 1 / np.sqrt(df.N)
    df.loc[:, "BETA"] = df.Z * df.SE
    unique_snps = df.RSID.unique()
    unique_genes = df.GENE.unique()
    return df, unique_snps, unique_genes


# NAME FORMAT: CELL, CHR, SNP, GENE, POS, A1, A2, BETA, SE, Z, PVAL, N
def load_gtex(path, samplesize, chr):
    # INPUT: chr : int
    # OUTPUT: df, unique_snps, unique_genes
    raw_df = pd.read_csv(
        os.path.join(path, "chr" + str(chr) + ".csv"), sep=",", header=0
    )
    raw_df.loc[:, "N"] = samplesize  # 670 samples in GTEx peripheral blood
    raw_df.loc[:, "Z"] = raw_df.BETA / raw_df.SE
    df = raw_df.loc[
        :, ["CHR", "RSID", "GENE", "POS", "A1", "A2", "BETA", "SE", "Z", "PVAL", "N"]
    ]
    # filter out palindromic snps
    df = df.loc[~is_palindromic_vectorized(df.A1, df.A2)]
    return df, df.RSID.unique(), df.GENE.unique()


def load_eQTLCatalogue(path, study, chr):
    raw_df = pd.read_csv(
        os.path.join(path, study, "chr" + str(chr) + ".csv"), sep=",", header=0
    )
    # CHR, SNP, GENE, POS, A1, A2, BETA, SE, PVAL, Z, N
    # 1, rs12238997, ENSG00000215790, 758351, A, G, -0.323775, 0.227148, 0.155233, -1.421288174311356, 277
    df = raw_df.loc[~is_palindromic_vectorized(raw_df["A1"], raw_df["A2"]), :]
    return df, df["RSID"].unique(), df["GENE"].unique()


def load_bbj(path, celltype, samplesize, chr):
    # CHR,RSID,POS,REF,ALT,Gene,BETA,Z,PVALUE
    # 20,rs6057421,29442702,C,T,ENSG00000205611,0.343792288853712,2.02104953504628,0.0458946843901897
    df = pd.read_csv(
        os.path.join(path, celltype, "chr" + str(chr) + ".csv"), sep=",", header=0
    )
    df.loc[:, "N"] = samplesize
    df.loc[:, "SE"] = df.BETA / df.Z
    return df, df.RSID.unique(), df.GENE.unique()


def load_afb(path, celltype, samplesize, chr):
    # GENE	POS	RSID	REF	ALT	Z	P	BETA	R2	SE	CisDist
    # ENSG00000048740	10812856	rs1000135	T	C	-1.62686963508815	0.108459224056508	-0.149765695781306	0.0380018669484576	0.0920575887281783	0
    # ENSG00000148429	10812856	rs1000135	T	C	-1.61020393435347	0.112055141388626	-0.199888920684129	0.0372561288382921	0.12413888478318	647654
    CELLTYPE_AFB = {
        "B_cells": "B",
        "NK_cells": "NK",
        "CD4+T_cells": "T.CD4",
        "CD8+T_cells": "T.CD8",
        "Monocytes": "MONO",
    }
    df = pd.read_csv(
        os.path.join(
            path,
            CELLTYPE_AFB[celltype] + "__NS",
            "eQTL_AFB_chr" + str(chr) + "_assoc.txt.gz",
        ),
        sep="\t",
        header=0,
        compression="gzip",
    )
    df.rename(
        columns={
            "ALT": "A1",
            "REF": "A2",
            "P": "PVAL",
        },
        inplace=True,
    )
    df.loc[:, "N"] = samplesize
    df.loc[:, "CHR"] = chr
    df.POS = df.POS.astype(int)
    # filter out palindromic snps
    df = df.loc[~is_palindromic_vectorized(df["A1"], df["A2"]), :]
    return df, df.RSID.unique(), df.GENE.unique()


def load_ld(tpld_path, apld_path, chr):
    # OUTPUT: ld_df, unique_snps
    chr = str(chr)
    match_file_tar = glob.glob(os.path.join(tpld_path, f"1000G.*.QC.maf.{chr}.bim"))[0]
    raw_target_population_df = pd.read_csv(
        match_file_tar,
        header=None,
        sep="\t",
    )
    raw_target_population_df.columns = ["CHR", "RSID", "CM", "POS", "A1", "A2"]
    match_file_aux = glob.glob(os.path.join(apld_path, f"1000G.*.QC.maf.{chr}.bim"))[0]
    raw_aux_population_df = pd.read_csv(
        match_file_aux,
        header=None,
        sep="\t",
    )
    raw_aux_population_df.columns = ["CHR", "RSID", "CM", "POS", "A1", "A2"]
    # merge two populations
    common_snp = set(raw_target_population_df["RSID"]).intersection(
        set(raw_aux_population_df["RSID"])
    )
    return list(common_snp)


def remove_from_df(df, snp_list=None, gene_list=None):
    # keep rows with snp in snp_list and gene in gene_list
    if snp_list:
        df = df[df["RSID"].isin(snp_list)]
    if gene_list:
        df = df[df["GENE"].isin(gene_list)]
    return df


def record_gene_snp_count(save_dir, chr, num_snps, num_genes):
    # write num_snps and num_genes into gene_snp_count/count.csv
    new_row = pd.DataFrame(
        {"CHR": [chr], "NUM_SNPS": [num_snps], "NUM_GENES": [num_genes]}
    )
    if os.path.exists(os.path.join(save_dir, "gene_snp_count", "count.csv")):
        old_df = pd.read_csv(
            os.path.join(save_dir, "gene_snp_count", "count.csv"),
            sep=",",
            header=0,
        )
        new_df = pd.concat([old_df, new_row], ignore_index=True)
        new_df["CHR"] = new_df["CHR"].astype(str)
        new_df.sort_values(by="CHR", inplace=True)
        new_df.to_csv(
            os.path.join(save_dir, "gene_snp_count", "count.csv"),
            index=False,
        )
    else:
        new_row.to_csv(
            os.path.join(save_dir, "gene_snp_count", "count.csv"),
            index=False,
        )


def main(args, chr, save_dir):
    study = args.studyid
    celltype = args.celltype
    print(f"\nStart merging chr {chr} of {study} and {celltype} at {time.ctime()}")
    ## load ld
    tpld_path = args.target_ld_dir
    apld_path = args.auxiliary_ld_dir
    ld_snps = load_ld(tpld_path, apld_path, chr)
    ## load target population eQTL data
    target_population_dir = args.target_dir
    target_samplesize = args.target_samplesize
    if "bbj" in target_population_dir.lower():
        load_target_population_function = load_bbj
    elif "af" in target_population_dir.lower():
        load_target_population_function = load_afb
    else:
        raise ValueError(
            f"Unsupported target population directory: {target_population_dir}"
        )
    tar_cell_all_df, tar_cell_snps, tar_cell_genes = load_target_population_function(
        target_population_dir, celltype, target_samplesize, chr
    )
    ## load auxiliary population eQTL data
    auxiliary_population_dir = args.auxiliary_dir
    aux_cell_all_df, aux_cell_snps, aux_cell_genes = load_eQTLCatalogue(
        auxiliary_population_dir, study, chr
    )
    ## load tissue eQTL data
    tissue_population_dir = args.tissue_dir
    tissue_samplesize = args.tissue_samplesize
    if "gtex" in tissue_population_dir.lower():
        load_tissue_function = load_gtex
    elif "eqtlgen" in tissue_population_dir.lower():
        load_tissue_function = load_eQTLGen
    else:
        raise ValueError(
            f"Unsupported tissue population directory: {tissue_population_dir}"
        )
    tissue_df, tissue_snps, tissue_genes = load_tissue_function(
        tissue_population_dir, tissue_samplesize, chr
    )
    ## get common snps and genes
    common_snps = (
        set(ld_snps)
        .intersection(tissue_snps)
        .intersection(aux_cell_snps)
        .intersection(tar_cell_snps)
    )
    common_genes = (
        set(tissue_genes).intersection(aux_cell_genes).intersection(tar_cell_genes)
    )

    # filter common snps and genes
    tissue_df_common = remove_from_df(tissue_df, common_snps, common_genes)

    # find unaligned snps and flip, find mismatched snps and remove
    merged_df = tissue_df_common.merge(
        aux_cell_all_df,
        on=["RSID", "GENE"],
        how="inner",
        suffixes=("", "_AUX"),
    )
    merged_df = merged_df.merge(
        tar_cell_all_df,
        on=["RSID", "GENE"],
        how="inner",
        suffixes=("", "_TAR"),
    )
    merged_df = merged_df.rename(columns={"A1": "A1_Tissue", "A2": "A2_Tissue"})

    # remove rows with NaN in any of the key columns
    merged_df.dropna(
        subset=[
            "A1_Tissue",
            "A2_Tissue",
            "A1_AUX",
            "A2_AUX",
            "A1_TAR",
            "A2_TAR",
        ],
        inplace=True,
    )
    print(f"{len(merged_df)} common SNP-gene pairs")
    ## find and flip mismatched SNPs remove SNPs that cannot be aligned
    ### compare eQTL_Tissue and tar_cell
    condition_eQTL_Tissue_same = (merged_df.A1_TAR == merged_df.A1_Tissue) & (
        merged_df.A2_TAR == merged_df.A2_Tissue
    )
    condition_eQTL_Tissue_flip = (merged_df.A1_TAR == merged_df.A2_Tissue) & (
        merged_df.A2_TAR == merged_df.A1_Tissue
    )
    remove_eQTL_Tissue = ~(condition_eQTL_Tissue_same | condition_eQTL_Tissue_flip)
    ### compare aux_cell and tar_cell
    condition_aux_same = (merged_df.A1_TAR == merged_df.A1_AUX) & (
        merged_df.A2_TAR == merged_df.A2_AUX
    )
    condition_aux_flip = (merged_df.A1_TAR == merged_df.A2_AUX) & (
        merged_df.A2_TAR == merged_df.A1_AUX
    )
    remove_AUX = ~(condition_aux_same | condition_aux_flip)
    remove_snp = remove_eQTL_Tissue | remove_AUX
    ### flip alleles
    merged_df.loc[condition_eQTL_Tissue_flip, ["A1_Tissue", "A2_Tissue"]] = (
        merged_df.loc[condition_eQTL_Tissue_flip, ["A2_Tissue", "A1_Tissue"]].values
    )
    merged_df.loc[condition_eQTL_Tissue_flip, ["BETA", "Z"]] *= -1
    merged_df.loc[condition_aux_flip, ["A1_AUX", "A2_AUX"]] = merged_df.loc[
        condition_aux_flip, ["A2_AUX", "A1_AUX"]
    ].values
    merged_df.loc[condition_aux_flip, ["BETA_AUX", "Z_AUX"]] *= -1
    ### remove mismatched SNPs
    merged_df = merged_df[~remove_snp].reset_index(drop=True)
    print(
        f"Flip {condition_eQTL_Tissue_flip.sum()} in eQTL_Tissue and {condition_aux_flip.sum()} in AUX population, remove {remove_snp.sum()} SNP-gene pairs"
    )

    ## update common snps and genes after filtering
    common_snps = merged_df.RSID.unique()
    common_genes = merged_df.GENE.unique()
    tissue_df_common = merged_df[
        [
            "RSID",
            "GENE",
            "BETA",
            "SE",
            "Z",
            "PVAL",
            "N",
        ]
    ]
    aux_cell_common = merged_df[
        [
            "RSID",
            "GENE",
            "BETA_AUX",
            "SE_AUX",
            "Z_AUX",
            "PVAL_AUX",
            "N_AUX",
        ]
    ].rename(
        columns={
            "BETA_AUX": "BETA",
            "SE_AUX": "SE",
            "Z_AUX": "Z",
            "PVAL_AUX": "PVAL",
            "N_AUX": "N",
        }
    )
    tar_cell_common = merged_df[
        [
            "RSID",
            "GENE",
            "BETA_TAR",
            "SE_TAR",
            "Z_TAR",
            "PVAL_TAR",
            "N_TAR",
        ]
    ].rename(
        columns={
            "BETA_TAR": "BETA",
            "SE_TAR": "SE",
            "Z_TAR": "Z",
            "PVAL_TAR": "PVAL",
            "N_TAR": "N",
        }
    )
    info_df = merged_df[
        [
            "GENE",
            "RSID",
            "POS",
            "A1_TAR",
            "A2_TAR",
        ]
    ].rename(
        columns={
            "A1_TAR": "A1",
            "A2_TAR": "A2",
        }
    )

    ## save common snps and genes count in tissue_df
    chr = str(chr)
    common_gene_snp_count = tissue_df_common.groupby("GENE").size().reset_index()
    common_gene_snp_count.columns = ["GENE", "SNP_COUNT"]
    common_gene_snp_count.to_csv(
        os.path.join(save_dir, "gene_snp_count", "chr" + chr + ".csv"), index=False
    )
    record_gene_snp_count(save_dir, chr, len(common_snps), len(common_genes))
    ## save common snps and genes in tissue_df
    tissue_df_common.to_csv(
        os.path.join(save_dir, "Tissue", "chr" + chr + ".csv"),
        index=False,
    )
    ## save common snps and genes in aux_cell
    aux_cell_common.to_csv(
        os.path.join(save_dir, "AUX_" + celltype, "chr" + chr + ".csv"), index=False
    )
    ## save common snps and genes in tar_cell
    tar_cell_common.to_csv(
        os.path.join(save_dir, "TAR_" + celltype, "chr" + chr + ".csv"), index=False
    )
    ## save info_df
    info_df.to_csv(
        os.path.join(save_dir, "INFO", "chr" + chr + ".csv"),
        index=False,
    )


if __name__ == "__main__":
    args = arg_parser()
    if isinstance(args.chromosome, int):
        chr_list = [int(args.chromosome)]
    else:
        chr_list = [int(i) for i in args.chromosome]
    save_dir = os.path.join(args.save_dir_main, args.studyid)

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if not os.path.exists(os.path.join(save_dir, "gene_snp_count")):
        os.makedirs(os.path.join(save_dir, "AUX_" + args.celltype))
        os.makedirs(os.path.join(save_dir, "TAR_" + args.celltype))
        os.makedirs(os.path.join(save_dir, "Tissue"))
        os.makedirs(os.path.join(save_dir, "gene_snp_count"))
        os.makedirs(os.path.join(save_dir, "INFO"))

    for chr in chr_list:
        main(args, chr, save_dir)
# END
