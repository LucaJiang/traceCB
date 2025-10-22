# run gmm.py: Run GMM for a given study and cell type and calculate effective sample size
import argparse
import time
import os
import numpy as np
import pandas as pd
from numba import prange

from ldsc import Run_single_LDSC, Run_cross_LDSC
from gmm import GMM, GMMtissue
from utils import make_pd_shrink, z2p, P_VAL_THRED, MIN_FLOAT, eSNP_THRESHOLD

# cmd: python -m pdb src/run_gmm.py -s QTD000066 -t CD8+T_cells -c 22 -d /gpfs1/scratch/wjiang49/xpmm/EAS_GTEx/


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
    return ((np.mean(z**2) - 1) / omega / np.mean(ld)).item()


def count_eSNPs(pval_list: np.ndarray) -> int:
    """
    Count the number of eSNPs based on p-values.
    """
    return np.sum(pval_list < eSNP_THRESHOLD).item()


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
            "COV",
            "COV_PVAL",
            "TAR_SNEFF",
            "TAR_CNEFF",
            "TAR_TNEFF",
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
    for gene in gene_list:
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
        if nsnp < 10:
            continue  # skip genes with less than 10 SNPs

        ## run ldsc to get Omega
        Omega, Omega_se = Run_cross_LDSC(
            (gene_df.BETA_tar / gene_df.SE_tar).to_numpy(),
            gene_df.N_tar.to_numpy(),
            gene_df.LD_tar.to_numpy(),
            (gene_df.BETA_aux / gene_df.SE_aux).to_numpy(),
            gene_df.N_aux.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            gene_df.LD_x.to_numpy(),
            np.array([1.0, 1.0, 0.0]),
        )
        Omega_p = z2p(Omega / Omega_se)
        ## get Sigma_o
        propt_weighted_omega, _ = Run_single_LDSC(
            gene_df.BETA_tissue.to_numpy() / gene_df.SE_tissue.to_numpy(),
            gene_df.N_tissue.to_numpy(),
            gene_df.LD_aux.to_numpy(),
            1.0,
        )
        weighted_omega_o = propt_weighted_omega - celltype_proportion**2 * Omega[1, 1]
        Sigma_o = weighted_omega_o * gene_df.LD_aux.to_numpy()
        summary_df.loc[gene, "H1SQ"] = Omega[0, 0]
        summary_df.loc[gene, "H1SQSE"] = Omega_se[0, 0]
        summary_df.loc[gene, "H2SQ"] = Omega[1, 1]
        summary_df.loc[gene, "H2SQSE"] = Omega_se[1, 1]
        summary_df.loc[gene, "COV"] = Omega[0, 1]
        summary_df.loc[gene, "COV_PVAL"] = Omega_p[0, 1]
        ### get neff of summary statistics
        summary_df.loc[gene, "TAR_SNEFF"] = calculate_Neff(
            gene_df.Z_tar.to_numpy(),
            Omega[0, 0],
            gene_df.LD_tar.to_numpy(),
        )
        summary_df.loc[gene, "AUX_SNEFF"] = calculate_Neff(
            gene_df.Z_aux.to_numpy(),
            Omega[1, 1],
            gene_df.LD_aux.to_numpy(),
        )
        summary_df.loc[gene, "TISSUE_SNEFF"] = calculate_Neff(
            gene_df.Z_tissue.to_numpy(),
            Omega[1, 1],
            gene_df.LD_aux.to_numpy(),
        )
        summary_df.loc[gene, "TAR_SeSNP"] = count_eSNPs(gene_df.PVAL_tar.to_numpy())
        summary_df.loc[gene, "AUX_SeSNP"] = count_eSNPs(gene_df.PVAL_aux.to_numpy())
        summary_df.loc[gene, "TISSUE_SeSNP"] = count_eSNPs(
            gene_df.PVAL_tissue.to_numpy()
        )
        ### threshold Omega
        for i in range(Omega.shape[0]):
            for j in range(Omega.shape[1]):
                if Omega_p[i, j] > P_VAL_THRED:
                    Omega[i, j] = MIN_FLOAT
        if np.any(Omega <= MIN_FLOAT):
            continue  # skip genes with unsignificant heritability or correlation
        Omega = make_pd_shrink(Omega, shrink=0.9)

        ## run GMM for each SNP in the gene
        gmm_b_tar = np.full((nsnp, 2), np.nan)  # cross, cross+tissue
        gmm_b_aux = np.full((nsnp, 2), np.nan)
        gmm_se_tar = np.full((nsnp, 2), np.nan)
        gmm_se_aux = np.full((nsnp, 2), np.nan)
        for i in prange(nsnp):
            try:
                gmm_b_tar[i, 0], gmm_se_tar[i, 0], gmm_b_aux[i, 0], gmm_se_aux[i, 0] = (
                    GMM(
                        Omega,
                        np.eye(2),
                        gene_df.BETA_tar.iloc[i].item(),
                        gene_df.SE_tar.iloc[i].item(),
                        gene_df.LD_tar.iloc[i].item(),
                        gene_df.BETA_aux.iloc[i].item(),
                        gene_df.SE_aux.iloc[i].item(),
                        gene_df.LD_aux.iloc[i].item(),
                        gene_df.LD_x.iloc[i].item(),
                    )
                )
                gmm_b_tar[i, 1], gmm_se_tar[i, 1], gmm_b_aux[i, 1], gmm_se_aux[i, 1] = (
                    GMMtissue(
                        Omega,
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
                        Sigma_o[i],
                        celltype_proportion,
                    )
                )
            except Exception as e:
                print(
                    f"Error in GMM for gene {gene}, SNP {gene_df.index[i]}: {e}. \
h1sq: {Omega[0, 0]}, h2sq: {Omega[1, 1]}, \
Cor: {Omega[0, 1]}, tar_beta: {gene_df.BETA_tar.iloc[i]}, \
aux_beta: {gene_df.BETA_aux.iloc[i]}, tissue_beta: \
{gene_df.BETA_tissue.iloc[i]}, tar_se: {gene_df.SE_tar.iloc[i]}, \
aux_se: {gene_df.SE_aux.iloc[i]}, tissue_se: {gene_df.SE_tissue.iloc[i]},\
tar_ld: {gene_df.LD_tar.iloc[i]}, aux_ld: {gene_df.LD_aux.iloc[i]}, \
x_ld: {gene_df.LD_x.iloc[i]}"
                )
                gmm_b_tar[i, :] = gene_df.BETA_tar.iloc[i]
                gmm_b_aux[i, :] = gene_df.BETA_aux.iloc[i]
                gmm_se_tar[i, :] = gene_df.SE_tar.iloc[i]
                gmm_se_aux[i, :] = gene_df.SE_aux.iloc[i]
        gmm_z_tar = gmm_b_tar / gmm_se_tar
        gmm_z_aux = gmm_b_aux / gmm_se_aux
        gmm_p_tar = z2p(gmm_z_tar)
        gmm_p_aux = z2p(gmm_z_aux)
        ## save gmm results
        gmm_results = pd.DataFrame(
            data=np.hstack(
                [
                    gene_df.index.to_numpy().reshape(-1, 1),
                    gene_df.BETA_tar.to_numpy().reshape(-1, 1),
                    gmm_b_tar,
                    gene_df.SE_tar.to_numpy().reshape(-1, 1),
                    gmm_se_tar,
                    gene_df.Z_tar.to_numpy().reshape(-1, 1),
                    gmm_z_tar,
                    gene_df.PVAL_tar.to_numpy().reshape(-1, 1),
                    gmm_p_tar,
                    gene_df.BETA_aux.to_numpy().reshape(-1, 1),
                    gmm_b_aux,
                    gene_df.SE_aux.to_numpy().reshape(-1, 1),
                    gmm_se_aux,
                    gene_df.Z_aux.to_numpy().reshape(-1, 1),
                    gmm_z_aux,
                    gene_df.PVAL_aux.to_numpy().reshape(-1, 1),
                    gmm_p_aux,
                    gene_df.BETA_tissue.to_numpy().reshape(-1, 1),
                    gene_df.SE_tissue.to_numpy().reshape(-1, 1),
                    gene_df.Z_tissue.to_numpy().reshape(-1, 1),
                    gene_df.PVAL_tissue.to_numpy().reshape(-1, 1),
                ]
            ),
            columns=[
                "RSID",
                "TAR_SBETA",
                "TAR_CBETA",
                "TAR_TBETA",
                "TAR_SSE",
                "TAR_CSE",
                "TAR_TSE",
                "TAR_SZ",
                "TAR_CZ",
                "TAR_TZ",
                "TAR_SPVAL",
                "TAR_CPVAL",
                "TAR_TPVAL",
                "AUX_SBETA",
                "AUX_CBETA",
                "AUX_TBETA",
                "AUX_SSE",
                "AUX_CSE",
                "AUX_TSE",
                "AUX_SZ",
                "AUX_CZ",
                "AUX_TZ",
                "AUX_SPVAL",
                "AUX_CPVAL",
                "AUX_TPVAL",
                "TISSUE_BETA",
                "TISSUE_SE",
                "TISSUE_Z",
                "TISSUE_PVAL",
            ],
        )
        gmm_results.to_csv(os.path.join(save_dir, f"{gene}.csv"), index=False)
        ## calculate effective sample size and count eSNPs
        summary_df.loc[gene, "TAR_CNEFF"] = calculate_Neff(
            gmm_z_tar[:, 0],
            Omega[0, 0],
            gene_df.LD_tar.to_numpy().reshape(-1, 1),
        )
        summary_df.loc[gene, "TAR_TNEFF"] = calculate_Neff(
            gmm_z_tar[:, 1],
            Omega[0, 0],
            gene_df.LD_tar.to_numpy().reshape(-1, 1),
        )
        summary_df.loc[gene, "AUX_CNEFF"] = calculate_Neff(
            gmm_z_aux[:, 0],
            Omega[1, 1],
            gene_df.LD_aux.to_numpy().reshape(-1, 1),
        )
        summary_df.loc[gene, "AUX_TNEFF"] = calculate_Neff(
            gmm_z_aux[:, 1],
            Omega[1, 1],
            gene_df.LD_aux.to_numpy().reshape(-1, 1),
        )
        summary_df.loc[gene, "TAR_CeSNP"] = count_eSNPs(gmm_p_tar[:, 0])
        summary_df.loc[gene, "TAR_TeSNP"] = count_eSNPs(gmm_p_tar[:, 1])
        summary_df.loc[gene, "AUX_CeSNP"] = count_eSNPs(gmm_p_aux[:, 0])
        summary_df.loc[gene, "AUX_TeSNP"] = count_eSNPs(gmm_p_aux[:, 1])
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
