# run gmm.py: Run GMM for a given study and cell type and calculate effective sample size
import sys
from pathlib import Path

SRC_DIR = Path(__file__).parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import argparse
import time
import os
import numpy as np
import pandas as pd
from numba import prange, njit

from traceCB.gmm import GMM, GMMtissue
from traceCB.ldsc import Run_cross_LDSC
from traceCB.utils import (
    z2p,
    # make_pd_shrink,
    MIN_FLOAT,
    MIN_HERITABILITY,
    P_VAL_THRED,
    MAX_CORR,
    eSNP_THRESHOLD,
)

# cmd: python src/traceCB/run_gmm.py -s QTD000081 -t Monocytes -c 1 -d /home/wjiang49/group/wjiang49/data/traceCB/EAS_GTEx >> log/test_sigmao_QTD000081_Monocytes.csv 2>&1


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--studyid",
        help="eqtlCatalogue study id",
        type=str,
    )
    parser.add_argument(
        "-t",
        "--celltype",
        help="cell type name",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        nargs="+",
        help="chromosome number",
        type=int,
    )
    parser.add_argument(
        "-d",
        "--main_dir",
        help="main directory",
        type=str,
    )
    return parser.parse_args()


def calculate_Neff(z: np.ndarray, omega: float, ld: np.ndarray) -> float:
    """
    Calculate effective sample size from z-scores, heritability omega, and LD vector.
    """
    return (np.mean(z**2) - 1) / (omega + MIN_FLOAT) / (np.mean(ld).item() + MIN_FLOAT)


def count_eSNPs(pval_list: np.ndarray) -> int:
    """
    Count the number of eSNPs based on p-values.
    """
    return np.sum(pval_list < eSNP_THRESHOLD).item()


def clip_correlation(var1: float, var2: float, cov: float) -> tuple[float, float]:
    """
    Clip the correlation coefficient to be within [-MAX_CORR, MAX_CORR].
    returns the adjusted covariance and correlation.
    """
    denominator = var1**0.5 * var2**0.5 + MIN_FLOAT
    cor = cov / denominator
    if abs(cor) > MAX_CORR:
        cor = np.sign(cor) * MAX_CORR
        cov = cor * denominator
    return cov, cor


def main(args, chr):
    studyid = args.studyid
    celltype = args.celltype
    main_dir = args.main_dir
    ## Set up directories
    study_dir = os.path.join(main_dir, studyid)
    save_dir = os.path.join(study_dir, "GMM", f"chr{chr}")
    os.makedirs(save_dir, exist_ok=True)
    ### ldsc directories
    ldsc_gene_path = os.path.join(study_dir, "LDSC", "LDSC_gene")
    tarLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr}_pop1.gz")
    auxLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr}_pop2.gz")
    xLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr}_te.gz")
    ### eQTL directories
    tarQTLfile = os.path.join(study_dir, f"TAR_{celltype}", f"chr{chr}.csv")
    auxQTLfile = os.path.join(study_dir, f"AUX_{celltype}", f"chr{chr}.csv")
    tissueQTLfile = os.path.join(study_dir, "Tissue", f"chr{chr}.csv")
    celltype_proportion_file = os.path.join(study_dir, "celltype_proportion.csv")
    ## Load data
    tarLDdf = pd.read_csv(
        tarLDfile, header=0, sep="\t", compression="gzip"
    )  # CHR	SNP	BP	base	ENSG00000015475 ...
    auxLDdf = pd.read_csv(
        auxLDfile, header=0, sep="\t", compression="gzip"
    )  # CHR	SNP	BP	base	ENSG00000015475 ...
    xLDdf = pd.read_csv(
        xLDfile, header=0, sep="\t", compression="gzip"
    )  # CHR	SNP	BP	base	ENSG00000015475 ...
    tarQTLdf = pd.read_csv(tarQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
    auxQTLdf = pd.read_csv(auxQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
    tissueQTLdf = pd.read_csv(tissueQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
    celltype_proportion_df = pd.read_csv(
        celltype_proportion_file,
        header=0,
        dtype={"Cell_type": str, "Proportion": float},
    )  # Cell_type,Proportion\nGranulocyte,53.120311443234144
    celltype_proportion_percentage = celltype_proportion_df.loc[
        celltype_proportion_df["Cell_type"] == celltype, "Proportion"
    ].iloc[  # type: ignore
        0
    ]
    celltype_proportion = celltype_proportion_percentage / 100.0

    ## run gmm for each gene
    gene_list = tarLDdf.columns[4:].tolist()
    summary_df = pd.DataFrame(
        data=np.nan,
        index=gene_list,
        columns=[
            "NSNP",
            "H1SQ",  # Per SNP heritability for target population
            "H1SQSE",
            "H2SQ",
            "H2SQSE",
            "COR",  # Genetic correlation between target and auxiliary cell types
            "COV_PVAL",  # Genetic covariance between target and auxiliary cell types
            "RUN_GMM",  # Whether run cross-population GMM
            "RUN_GMM_TISSUE",  # Whether run GMM with tissue
            "TAR_SNEFF",  # effective sample size, S is for summary statistics
            "TAR_CNEFF",  # effective sample size, C is for using cross-population GMM
            "TAR_TNEFF",  # effective sample size, T is for using GMM with tissue
            "TAR_SeSNP",  # number of eSNP
            "TAR_CeSNP",
            "TAR_TeSNP",
            "AUX_SNEFF",
            "AUX_CNEFF",
            "AUX_TNEFF",
            "AUX_SeSNP",
            "AUX_CeSNP",
            "AUX_TeSNP",
            "TISSUE_SNEFF",
            "TISSUE_SeSNP",
        ],
    )
    summary_df.index.name = "GENE"
    summary_df = summary_df.astype({"RUN_GMM": "boolean", "RUN_GMM_TISSUE": "boolean"})
    summary_df.RUN_GMM = False  # init with False
    summary_df.RUN_GMM_TISSUE = False
    # for gene in gene_list:
    # breakpoint()
    for gene_idx in prange(len(gene_list)):
        gene = gene_list[gene_idx]
        ## collect data for the gene
        tarQTLgene_df = (
            tarQTLdf.loc[tarQTLdf["GENE"] == gene, :]
            .drop_duplicates(subset="RSID", keep="first")
            .drop(columns=["GENE"])
        )
        auxQTLdfgene_df = (
            auxQTLdf.loc[auxQTLdf["GENE"] == gene, :]
            .drop_duplicates(subset="RSID", keep="first")
            .drop(columns=["GENE"])
        )
        tissueQTLdfgene_df = (
            tissueQTLdf.loc[tissueQTLdf["GENE"] == gene, :]
            .drop_duplicates(subset="RSID", keep="first")
            .drop(columns=["GENE"])
        )
        tarLDdfgene_df = tarLDdf.loc[:, ["SNP", gene]].rename(
            columns={gene: "LD_tar", "SNP": "RSID"}
        )
        auxLDdfgene_df = auxLDdf.loc[:, ["SNP", gene]].rename(
            columns={gene: "LD_aux", "SNP": "RSID"}
        )
        xLDdfgene_df = xLDdf.loc[:, ["SNP", gene]].rename(
            columns={gene: "LD_x", "SNP": "RSID"}
        )
        rename_cols = {
            col: f"{col}_tissue" for col in tissueQTLdfgene_df.columns if col != "RSID"
        }
        gene_df = (
            tarQTLgene_df.merge(auxQTLdfgene_df, on="RSID", suffixes=("_tar", "_aux"))
            .merge(tissueQTLdfgene_df.rename(columns=rename_cols), on="RSID")
            .merge(tarLDdfgene_df, on="RSID")
            .merge(auxLDdfgene_df, on="RSID")
            .merge(xLDdfgene_df, on="RSID")
            .set_index("RSID")
            .dropna(axis=0, how="any")
        )  # RSID,GENE,BETA_tar,SE_tar,Z_tar,PVAL_tar,N_tar,BETA_aux,SE_aux,Z_aux,PVAL_aux,N_aux,BETA_tissue,SE_tissue,Z_tissue,PVAL_tissue,N_tissue,tarLD,auxLD,xLD
        nsnp = gene_df.shape[0]
        summary_df.loc[gene, "NSNP"] = nsnp
        if nsnp < 100:  # skip genes with less than 100 SNPs
            continue

        # run TAR + AUX to get traceC
        OmegaC, OmegaC_se = Run_cross_LDSC(
            (gene_df.BETA_tar / gene_df.SE_tar).to_numpy(),
            gene_df.N_tar.to_numpy(),
            gene_df.LD_tar.to_numpy(),
            (gene_df.BETA_aux / gene_df.SE_aux).to_numpy(),
            gene_df.N_aux.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            gene_df.LD_x.to_numpy(),
            np.array([1.0, 1.0, 0.0]),
        )
        OmegaC_p = z2p(OmegaC / OmegaC_se)
        summary_df.loc[gene, "H1SQ"] = OmegaC[0, 0]
        summary_df.loc[gene, "H1SQSE"] = OmegaC_se[0, 0]
        summary_df.loc[gene, "H2SQ"] = OmegaC[1, 1]
        summary_df.loc[gene, "H2SQSE"] = OmegaC_se[1, 1]
        summary_df.loc[gene, "COR"] = OmegaC[0, 1] / (
            (OmegaC[0, 0] * OmegaC[1, 1]) ** 0.5 + MIN_FLOAT
        )
        summary_df.loc[gene, "COV_PVAL"] = OmegaC_p[0, 1]
        summary_df.loc[gene, "TAR_SNEFF"] = calculate_Neff(
            gene_df.Z_tar.to_numpy(),
            OmegaC[0, 0],
            gene_df.LD_tar.to_numpy(),
        )
        summary_df.loc[gene, "AUX_SNEFF"] = calculate_Neff(
            gene_df.Z_aux.to_numpy(),
            OmegaC[1, 1],
            gene_df.LD_aux.to_numpy(),
        )
        summary_df.loc[gene, "TISSUE_SNEFF"] = calculate_Neff(
            gene_df.Z_tissue.to_numpy(),
            OmegaC[1, 1],
            gene_df.LD_aux.to_numpy(),
        )
        summary_df.loc[gene, "TAR_SeSNP"] = count_eSNPs(gene_df.PVAL_tar.to_numpy())
        summary_df.loc[gene, "AUX_SeSNP"] = count_eSNPs(gene_df.PVAL_aux.to_numpy())
        summary_df.loc[gene, "TISSUE_SeSNP"] = count_eSNPs(
            gene_df.PVAL_tissue.to_numpy()
        )
        # initialize gmm(traceC) results
        gmm_b_tar = np.full((nsnp), np.nan)
        gmm_se_tar = np.full((nsnp), np.nan)
        gmm_z_tar = np.full((nsnp), np.nan)
        gmm_p_tar = np.full((nsnp), np.nan)

        gmm_b_aux = np.full((nsnp), np.nan)
        gmm_se_aux = np.full((nsnp), np.nan)
        gmm_z_aux = np.full((nsnp), np.nan)
        gmm_p_aux = np.full((nsnp), np.nan)

        gmm_b_tar_cb = np.full((nsnp), np.nan)
        gmm_se_tar_cb = np.full((nsnp), np.nan)
        gmm_z_tar_cb = np.full((nsnp), np.nan)
        gmm_p_tar_cb = np.full((nsnp), np.nan)

        gmm_b_aux_cb = np.full((nsnp), np.nan)
        gmm_se_aux_cb = np.full((nsnp), np.nan)
        gmm_z_aux_cb = np.full((nsnp), np.nan)
        gmm_p_aux_cb = np.full((nsnp), np.nan)
        # judge whether to run traceC
        if np.all(OmegaC_p < P_VAL_THRED):
            summary_df.loc[gene, "RUN_GMM"] = True
            for i in prange(nsnp):
                (gmm_b_tar[i], gmm_se_tar[i], gmm_b_aux[i], gmm_se_aux[i]) = GMM(
                    OmegaC,
                    np.eye(2),
                    gene_df.BETA_tar.iloc[i].item(),
                    gene_df.SE_tar.iloc[i].item(),
                    gene_df.LD_tar.iloc[i].item(),
                    gene_df.BETA_aux.iloc[i].item(),
                    gene_df.SE_aux.iloc[i].item(),
                    gene_df.LD_aux.iloc[i].item(),
                    gene_df.LD_x.iloc[i].item(),
                )
            gmm_z_tar = gmm_b_tar / gmm_se_tar
            gmm_z_aux = gmm_b_aux / gmm_se_aux
            gmm_p_tar = z2p(gmm_z_tar)
            gmm_p_aux = z2p(gmm_z_aux)
            summary_df.loc[gene, "TAR_CNEFF"] = calculate_Neff(
                gmm_z_tar,
                OmegaC[0, 0],
                gene_df.LD_tar.to_numpy(),
            )
            summary_df.loc[gene, "AUX_CNEFF"] = calculate_Neff(
                gmm_z_aux,
                OmegaC[1, 1],
                gene_df.LD_aux.to_numpy(),
            )
            summary_df.loc[gene, "TAR_CeSNP"] = count_eSNPs(gmm_p_tar)
            summary_df.loc[gene, "AUX_CeSNP"] = count_eSNPs(gmm_p_aux)

        # run AUX + Tissue to get updated AUX
        Omega, Omega_se = Run_cross_LDSC(
            (gene_df.BETA_aux / gene_df.SE_aux).to_numpy(),
            gene_df.N_aux.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            (gene_df.BETA_tissue / gene_df.SE_tissue).to_numpy(),
            gene_df.N_tissue.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            np.array([1.0, 1.0, 0.0]),
        )
        Omega_p = z2p(Omega / Omega_se)
        if np.all(Omega_p < P_VAL_THRED):
            # update AUX summary statistics with GMMtissue
            gmm_b_aux_t = np.full((nsnp), np.nan)
            gmm_se_aux_t = np.full((nsnp), np.nan)
            OmegaCB = np.array(
                [
                    [1, 0.0, 0.0],
                    [0.0, Omega[0, 0], Omega[0, 1]],
                    [0.0, Omega[0, 1], Omega[1, 1]],
                ]
            )
            for i in prange(nsnp):
                (
                    _,
                    _,
                    gmm_b_aux_t[i],
                    gmm_se_aux_t[i],
                ) = GMMtissue(
                    OmegaCB,
                    np.eye(3),
                    gene_df.BETA_tar.iloc[i].item(),
                    gene_df.SE_tar.iloc[i].item(),
                    gene_df.LD_tar.iloc[i].item(),
                    gene_df.BETA_aux.iloc[i].item(),
                    gene_df.SE_aux.iloc[i].item(),
                    gene_df.LD_aux.iloc[i].item(),
                    gene_df.LD_x.iloc[i].item(),
                    gene_df.BETA_tissue.iloc[i].item(),
                    gene_df.SE_tissue.iloc[i].item(),
                    celltype_proportion,
                )

            ## run TAR + updated AUX to get final TAR
            Omega, Omega_se = Run_cross_LDSC(
                (gene_df.BETA_tar / gene_df.SE_tar).to_numpy(),
                gene_df.N_tar.to_numpy(),
                gene_df.LD_tar.to_numpy(),
                gmm_b_aux_t / gmm_se_aux_t,
                gene_df.N_aux.to_numpy(),
                gene_df.LD_aux.to_numpy(),
                gene_df.LD_x.to_numpy(),
                np.array([1.0, 1.0, 0.0]),
            )
            Omega_p = z2p(Omega / Omega_se)
            if np.all(Omega_p < P_VAL_THRED):
                summary_df.loc[gene, "RUN_GMM_TISSUE"] = True
                for i in prange(nsnp):
                    (
                        gmm_b_tar_cb[i],
                        gmm_se_tar_cb[i],
                        gmm_b_aux_cb[i],
                        gmm_se_aux_cb[i],
                    ) = GMM(
                        Omega,
                        np.eye(2),
                        gene_df.BETA_tar.iloc[i].item(),
                        gene_df.SE_tar.iloc[i].item(),
                        gene_df.LD_tar.iloc[i].item(),
                        gmm_b_aux_t[i],
                        gmm_se_aux_t[i],
                        gene_df.LD_aux.iloc[i].item(),
                        gene_df.LD_x.iloc[i].item(),
                    )
                gmm_z_tar_cb = gmm_b_tar_cb / gmm_se_tar_cb
                gmm_p_tar_cb = z2p(gmm_z_tar_cb)
                gmm_z_aux_cb = gmm_b_aux_cb / gmm_se_aux_cb
                gmm_p_aux_cb = z2p(gmm_z_aux_cb)
                summary_df.loc[gene, "TAR_TNEFF"] = calculate_Neff(
                    gmm_z_tar_cb,
                    Omega[0, 0],
                    gene_df.LD_tar.to_numpy(),
                )
                summary_df.loc[gene, "TAR_TeSNP"] = count_eSNPs(gmm_p_tar_cb)
                summary_df.loc[gene, "AUX_TNEFF"] = calculate_Neff(
                    gmm_z_aux_cb,
                    Omega[1, 1],
                    gene_df.LD_aux.to_numpy(),
                )
                summary_df.loc[gene, "AUX_TeSNP"] = count_eSNPs(gmm_p_aux_cb)

        if summary_df.loc[gene, "RUN_GMM"] or summary_df.loc[gene, "RUN_GMM_TISSUE"]:
            # save gmm results for the gene
            gmm_results = pd.DataFrame(
                {
                    "RSID": gene_df.index,
                    "TAR_SBETA": gene_df.BETA_tar,
                    "TAR_CBETA": gmm_b_tar,
                    "TAR_TBETA": gmm_b_tar_cb,
                    "TAR_SSE": gene_df.SE_tar,
                    "TAR_CSE": gmm_se_tar,
                    "TAR_TSE": gmm_se_tar_cb,
                    "TAR_SZ": gene_df.Z_tar,
                    "TAR_CZ": gmm_z_tar,
                    "TAR_TZ": gmm_z_tar_cb,
                    "TAR_SPVAL": gene_df.PVAL_tar,
                    "TAR_CPVAL": gmm_p_tar,
                    "TAR_TPVAL": gmm_p_tar_cb,
                    "AUX_SBETA": gene_df.BETA_aux,
                    "AUX_CBETA": gmm_b_aux,
                    "AUX_TBETA": gmm_b_aux_cb,
                    "AUX_SSE": gene_df.SE_aux,
                    "AUX_CSE": gmm_se_aux,
                    "AUX_TSE": gmm_se_aux_cb,
                    "AUX_SZ": gene_df.Z_aux,
                    "AUX_CZ": gmm_z_aux,
                    "AUX_TZ": gmm_z_aux_cb,
                    "AUX_SPVAL": gene_df.PVAL_aux,
                    "AUX_CPVAL": gmm_p_aux,
                    "AUX_TPVAL": gmm_p_aux_cb,
                    "TISSUE_BETA": gene_df.BETA_tissue,
                    "TISSUE_SE": gene_df.SE_tissue,
                    "TISSUE_Z": gene_df.Z_tissue,
                    "TISSUE_PVAL": gene_df.PVAL_tissue,
                }
            )

            gmm_results.to_csv(os.path.join(save_dir, f"{gene}.csv"), index=False)

    ## finish gene loop, save summary df
    summary_df.to_csv(os.path.join(save_dir, "summary.csv"))
    return


if __name__ == "__main__":
    args = parse_args()
    if isinstance(args.chromosome, int):
        chr_list = [int(args.chromosome)]
    else:
        chr_list = [int(i) for i in args.chromosome]
    ## Set up directories
    save_dir = os.path.join(args.main_dir, args.studyid, "GMM")
    os.makedirs(save_dir, exist_ok=True)

    start_time = time.time()
    print(
        f"Running GMM for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(start_time)}"
    )
    for chr in chr_list:
        main(args, chr)  # run GMM for each chromosome
    end_time = time.time()
    print(
        f"GMM finished for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(end_time)} which took {end_time - start_time:.2f} seconds"
    )
