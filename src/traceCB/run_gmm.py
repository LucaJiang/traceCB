import os
import time
import argparse
import pyarrow
import numpy as np
import pandas as pd
from numba import prange, njit

from traceCB.gmm import GMM, GMMtissue
from traceCB.ldsc import Run_Cross_LDSC
from traceCB.utils import (
    z2p,
    MIN_FLOAT,
    MAX_CORR,
    eSNP_THRESHOLD,
)

# ------------------------------------------------------------------------------
# Core Calculation Kernel (Numba Optimized)
# ------------------------------------------------------------------------------


@njit(parallel=True)
def run_gmm_kernel(
    nsnp,
    run_tissue,
    Omega,
    pi2_omega_sum,
    celltype_proportion,
    beta_tar,
    se_tar,
    ld_tar,
    beta_aux,
    se_aux,
    ld_aux,
    ld_x,
    beta_tissue,
    se_tissue,
):
    """
    JIT-compiled kernel to run GMM in parallel for all SNPs of a single gene.
    This replaces the slow python loop with optimized machine code.
    """
    # Pre-allocate result arrays
    gmm_b_tar = np.full((nsnp, 2), np.nan)
    gmm_se_tar = np.full((nsnp, 2), np.nan)
    gmm_b_aux = np.full((nsnp, 2), np.nan)
    gmm_se_aux = np.full((nsnp, 2), np.nan)

    # Pre-allocate constant matrices to avoid re-creation inside loop
    eye2 = np.eye(2)
    eye3 = np.eye(3)

    for i in prange(nsnp):
        # Model 0: Cross population only
        (
            gmm_b_tar[i, 0],
            gmm_se_tar[i, 0],
            gmm_b_aux[i, 0],
            gmm_se_aux[i, 0],
        ) = GMM(
            Omega,
            eye2,
            beta_tar[i],
            se_tar[i],
            ld_tar[i],
            beta_aux[i],
            se_aux[i],
            ld_aux[i],
            ld_x[i],
        )

        if run_tissue:
            # Model 1: Cross population + Tissue
            (
                gmm_b_tar[i, 1],
                gmm_se_tar[i, 1],
                gmm_b_aux[i, 1],
                gmm_se_aux[i, 1],
            ) = GMMtissue(
                Omega,
                eye3,
                beta_tar[i],
                se_tar[i],
                ld_tar[i],
                beta_aux[i],
                se_aux[i],
                ld_aux[i],
                ld_x[i],
                beta_tissue[i],
                se_tissue[i],
                pi2_omega_sum,
                celltype_proportion,
            )
        else:
            # If not running tissue, copy Model 0 results
            gmm_b_tar[i, 1] = gmm_b_tar[i, 0]
            gmm_se_tar[i, 1] = gmm_se_tar[i, 0]
            gmm_b_aux[i, 1] = gmm_b_aux[i, 0]
            gmm_se_aux[i, 1] = gmm_se_aux[i, 0]

    return gmm_b_tar, gmm_se_tar, gmm_b_aux, gmm_se_aux


# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--studyid", help="eqtlCatalogue study id", type=str)
    parser.add_argument("-t", "--celltype", help="cell type name", type=str)
    parser.add_argument(
        "-c", "--chromosome", nargs="+", help="chromosome number", type=int
    )
    parser.add_argument("-d", "--main_dir", help="main directory", type=str)
    return parser.parse_args()


def calculate_Neff(z: np.ndarray, omega: float, ld: np.ndarray) -> float:
    return (np.mean(z**2) - 1) / (omega + MIN_FLOAT) / (np.mean(ld).item() + MIN_FLOAT)


def count_eSNPs(pval_list: np.ndarray) -> int:
    return np.sum(pval_list < eSNP_THRESHOLD).item()


def clip_correlation(var1: float, var2: float, cov: float) -> tuple[float, float]:
    denominator = var1**0.5 * var2**0.5 + MIN_FLOAT
    cor = cov / denominator
    if abs(cor) > MAX_CORR:
        cor = np.sign(cor) * MAX_CORR
        cov = cor * denominator
    return cov, cor


def load_and_group_qtl(filepath: str, suffix: str = "") -> dict:
    """
    Reads a CSV file, drops duplicates, renames columns with suffix,
    and groups by 'GENE' into a dictionary for O(1) access.
    """
    df = pd.read_csv(filepath, header=0)

    # Rename columns except keys
    if suffix:
        rename_map = {
            col: f"{col}{suffix}" for col in df.columns if col not in ["RSID", "GENE"]
        }
        df.rename(columns=rename_map, inplace=True)

    # Drop duplicates
    df = df.drop_duplicates(subset=["GENE", "RSID"], keep="first")

    # Group by Gene and convert to dictionary
    # structure: {'gene_name': DataFrame}
    return dict(tuple(df.groupby("GENE")))


# ------------------------------------------------------------------------------
# Main Logic
# ------------------------------------------------------------------------------


def main(args, chr_num):
    studyid = args.studyid
    celltype = args.celltype
    main_dir = args.main_dir

    ## Set up directories
    study_dir = os.path.join(main_dir, studyid)
    save_dir = os.path.join(study_dir, "GMM", f"chr{chr_num}")
    os.makedirs(save_dir, exist_ok=True)

    ### ldsc directories
    ldsc_gene_path = os.path.join(study_dir, "LDSC", "LDSC_gene")
    tarLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr_num}_pop1.gz")
    auxLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr_num}_pop2.gz")
    xLDfile = os.path.join(ldsc_gene_path, f"TAR_AUX_std_chr{chr_num}_te.gz")

    ### eQTL directories
    tarQTLfile = os.path.join(study_dir, f"TAR_{celltype}", f"chr{chr_num}.csv")
    auxQTLfile = os.path.join(study_dir, f"AUX_{celltype}", f"chr{chr_num}.csv")
    tissueQTLfile = os.path.join(study_dir, "Tissue", f"chr{chr_num}.csv")
    celltype_proportion_file = os.path.join(study_dir, "celltype_proportion.csv")

    print(f"Loading data for chr{chr_num}...")

    # 1. Load Proportions
    celltype_proportion_df = pd.read_csv(
        celltype_proportion_file,
        header=0,
        dtype={"Cell_type": str, "Proportion": float},
    )
    # Check if cell type exists
    if celltype not in celltype_proportion_df["Cell_type"].values:
        print(f"Warning: Cell type {celltype} not found in proportions. Skipping.")
        return

    celltype_proportion_percentage = celltype_proportion_df.loc[
        celltype_proportion_df["Cell_type"] == celltype, "Proportion"
    ].iloc[0]
    celltype_proportion = celltype_proportion_percentage / 100.0

    other_celltype_proportions = (
        celltype_proportion_df.loc[
            celltype_proportion_df["Cell_type"] != celltype, "Proportion"
        ]
        / 100.0
    ).to_numpy()

    pii_pij_sum = sum(
        [
            i * j
            for i in other_celltype_proportions
            for j in other_celltype_proportions
            if i != j
        ]
    )
    pi2_omega_sum_const = 2 * celltype_proportion + pii_pij_sum / (
        1 - celltype_proportion
    )

    # 2. Load LD files (Wide format: rows=SNPs, cols=Genes)
    tarLDdf = pd.read_csv(tarLDfile, sep="\t", compression="gzip")
    auxLDdf = pd.read_csv(auxLDfile, sep="\t", compression="gzip")
    xLDdf = pd.read_csv(xLDfile, sep="\t", compression="gzip")
    gene_list = tarLDdf.columns[4:].tolist()

    # 3. Load and Group QTL files
    tar_qtl_dict = load_and_group_qtl(tarQTLfile, suffix="_tar")
    aux_qtl_dict = load_and_group_qtl(auxQTLfile, suffix="_aux")

    # Tissue file
    tissue_df_raw = pd.read_csv(tissueQTLfile, header=0)
    tissue_rename_cols = {
        col: f"{col}_tissue"
        for col in tissue_df_raw.columns
        if col not in ["RSID", "GENE"]
    }
    tissue_df_raw.rename(columns=tissue_rename_cols, inplace=True)
    tissue_df_raw = tissue_df_raw.drop_duplicates(subset=["GENE", "RSID"], keep="first")
    tissue_qtl_dict = dict(tuple(tissue_df_raw.groupby("GENE")))

    # 4. Prepare Summary Data Structure
    summary_data = []  # Collect results in a list first (faster than appending to DF)
    summary_columns = [
        "NSNP",
        "H1SQ",
        "H1SQSE",
        "H2SQ",
        "H2SQSE",
        "COR_X_ORI",
        "COV_PVAL",
        "COR_X",
        "SIGMAO",
        "RUN_GMM",
        "TAR_SNEFF",
        "TAR_CNEFF",
        "TAR_TNEFF",
        "TAR_SeSNP",
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
    ]

    # print(f"Processing {len(gene_list)} genes on chr{chr_num}...")

    # 5. Main Loop
    for gene in gene_list:
        # Initialize summary row
        gene_res = {col: np.nan for col in summary_columns}
        # Check if gene exists in all QTL dictionaries
        if (
            (gene not in tar_qtl_dict)
            or (gene not in aux_qtl_dict)
            or (gene not in tissue_qtl_dict)
        ):
            continue

        # Retrieve DataFrames from dictionary (O(1) operation)
        tarQTLgene_df = tar_qtl_dict[gene].drop(columns=["GENE"])
        auxQTLgene_df = aux_qtl_dict[gene].drop(columns=["GENE"])
        tissueQTLgene_df = tissue_qtl_dict[gene].drop(columns=["GENE"])

        # Prepare LD columns for merge
        # Original logic: rename {gene} -> LD_tar/aux/x
        tarLDgene_df = tarLDdf[["SNP", gene]].rename(
            columns={gene: "LD_tar", "SNP": "RSID"}
        )
        auxLDgene_df = auxLDdf[["SNP", gene]].rename(
            columns={gene: "LD_aux", "SNP": "RSID"}
        )
        xLDgene_df = xLDdf[["SNP", gene]].rename(columns={gene: "LD_x", "SNP": "RSID"})

        # Merge all data
        gene_df = (
            tarQTLgene_df.merge(auxQTLgene_df, on="RSID")
            .merge(tissueQTLgene_df, on="RSID")
            .merge(tarLDgene_df, on="RSID")
            .merge(auxLDgene_df, on="RSID")
            .merge(xLDgene_df, on="RSID")
        )
        gene_df.dropna(axis=0, how="any", inplace=True)
        nsnp = gene_df.shape[0]
        gene_res["NSNP"] = nsnp

        if nsnp < 100:
            summary_data.append((gene, gene_res))
            continue

        # ---------------------------------------------------------
        # Convert columns to Numpy Arrays for Vectorized operations
        # ---------------------------------------------------------
        beta_tar = gene_df.BETA_tar.to_numpy()
        se_tar = gene_df.SE_tar.to_numpy()
        n_tar = gene_df.N_tar.to_numpy()
        ld_tar = gene_df.LD_tar.to_numpy()
        z_tar = gene_df.Z_tar.to_numpy()
        pval_tar = gene_df.PVAL_tar.to_numpy()

        beta_aux = gene_df.BETA_aux.to_numpy()
        se_aux = gene_df.SE_aux.to_numpy()
        n_aux = gene_df.N_aux.to_numpy()
        ld_aux = gene_df.LD_aux.to_numpy()
        z_aux = gene_df.Z_aux.to_numpy()
        pval_aux = gene_df.PVAL_aux.to_numpy()

        beta_tissue = gene_df.BETA_tissue.to_numpy()
        se_tissue = gene_df.SE_tissue.to_numpy()
        n_tissue = gene_df.N_tissue.to_numpy()
        z_tissue = gene_df.Z_tissue.to_numpy()
        pval_tissue = gene_df.PVAL_tissue.to_numpy()

        ld_x = gene_df.LD_x.to_numpy()

        # Run LDSC
        Omega, Omega_se = Run_Cross_LDSC(
            beta_tar / se_tar,
            n_tar,
            ld_tar,
            beta_aux / se_aux,
            n_aux,
            ld_aux,
            ld_x,
            np.array([1.0, 1.0, 0.0]),
        )

        cor_original = Omega[0, 1] / (
            Omega[0, 0] ** 0.5 * Omega[1, 1] ** 0.5 + MIN_FLOAT
        )
        gene_res["COR_X_ORI"] = cor_original

        Omega[0, 1], cor_clipped = clip_correlation(
            Omega[0, 0], Omega[1, 1], Omega[0, 1]
        )
        Omega[1, 0] = Omega[0, 1]
        Omega_p = z2p(Omega / Omega_se)

        gene_res["H1SQ"] = Omega[0, 0]
        gene_res["H1SQSE"] = Omega_se[0, 0]
        gene_res["H2SQ"] = Omega[1, 1]
        gene_res["H2SQSE"] = Omega_se[1, 1]
        gene_res["COV_PVAL"] = Omega_p[0, 1]

        # Neff Summary Stats
        gene_res["TAR_SNEFF"] = calculate_Neff(z_tar, Omega[0, 0], ld_tar)
        gene_res["AUX_SNEFF"] = calculate_Neff(z_aux, Omega[1, 1], ld_aux)
        gene_res["TISSUE_SNEFF"] = calculate_Neff(z_tissue, Omega[1, 1], ld_aux)

        gene_res["TAR_SeSNP"] = count_eSNPs(pval_tar)
        gene_res["AUX_SeSNP"] = count_eSNPs(pval_aux)
        gene_res["TISSUE_SeSNP"] = count_eSNPs(pval_tissue)

        # Threshold check
        if np.any(Omega_p >= 0.1):
            gene_res["RUN_GMM"] = False
            summary_data.append((gene, gene_res))
            continue

        cor_x = Omega[0, 1] / np.sqrt(Omega[0, 0] * Omega[1, 1] + MIN_FLOAT)
        gene_res["COR_X"] = cor_x
        gene_res["RUN_GMM"] = True

        # Calculate Sigma_o
        pi2_omega_sum = 0.0
        run_gmm_tissue = True  # set default to run GMM with tissue

        aux_Omega_matrix, aux_Omega_matrix_se = Run_Cross_LDSC(
            beta_aux / se_aux,
            n_aux,
            ld_aux,
            beta_tissue / se_tissue,
            n_tissue,
            ld_aux,
            ld_aux,
            np.array([1.0, 1.0, 0.0]),
        )

        aux_Omega_matrix[0, 1], _ = clip_correlation(
            aux_Omega_matrix[0, 0], aux_Omega_matrix[1, 1], aux_Omega_matrix[0, 1]
        )
        aux_Omega_matrix[1, 0] = aux_Omega_matrix[0, 1]

        if np.all(z2p(aux_Omega_matrix / aux_Omega_matrix_se) < 0.10):
            pi2_omega_sum = (
                aux_Omega_matrix[1, 1]
                - celltype_proportion**2 * Omega[1, 1]
                - pi2_omega_sum_const
                * np.maximum(
                    aux_Omega_matrix[0, 1] - celltype_proportion * Omega[1, 1], 0
                )
            )
            gene_res["SIGMAO"] = pi2_omega_sum
            if pi2_omega_sum < 0:
                run_gmm_tissue = False
        else:
            run_gmm_tissue = False

        # ------------------------------------------------------------------
        # Run GMM (Parallelized Kernel)
        # ------------------------------------------------------------------
        # Pass numpy arrays to Numba function
        gmm_b_tar, gmm_se_tar, gmm_b_aux, gmm_se_aux = run_gmm_kernel(
            nsnp,
            run_gmm_tissue,
            Omega,
            pi2_omega_sum,
            celltype_proportion,
            beta_tar,
            se_tar,
            ld_tar,
            beta_aux,
            se_aux,
            ld_aux,
            ld_x,
            beta_tissue,
            se_tissue,
        )

        # Calculate Z and P
        gmm_z_tar = gmm_b_tar / gmm_se_tar
        gmm_z_aux = gmm_b_aux / gmm_se_aux
        gmm_p_tar = z2p(gmm_z_tar)
        gmm_p_aux = z2p(gmm_z_aux)

        # ------------------------------------------------------------------
        # Save GMM Results (Parquet format)
        # ------------------------------------------------------------------
        res_data = {
            "RSID": gene_df["RSID"].values,  # Use .values to ensure array
            "TAR_SBETA": beta_tar,
            "TAR_CBETA": gmm_b_tar[:, 0],
            "TAR_TBETA": gmm_b_tar[:, 1],
            "TAR_SSE": se_tar,
            "TAR_CSE": gmm_se_tar[:, 0],
            "TAR_TSE": gmm_se_tar[:, 1],
            "TAR_SZ": z_tar,
            "TAR_CZ": gmm_z_tar[:, 0],
            "TAR_TZ": gmm_z_tar[:, 1],
            "TAR_SPVAL": pval_tar,
            "TAR_CPVAL": gmm_p_tar[:, 0],
            "TAR_TPVAL": gmm_p_tar[:, 1],
            "AUX_SBETA": beta_aux,
            "AUX_CBETA": gmm_b_aux[:, 0],
            "AUX_TBETA": gmm_b_aux[:, 1],
            "AUX_SSE": se_aux,
            "AUX_CSE": gmm_se_aux[:, 0],
            "AUX_TSE": gmm_se_aux[:, 1],
            "AUX_SZ": z_aux,
            "AUX_CZ": gmm_z_aux[:, 0],
            "AUX_TZ": gmm_z_aux[:, 1],
            "AUX_SPVAL": pval_aux,
            "AUX_CPVAL": gmm_p_aux[:, 0],
            "AUX_TPVAL": gmm_p_aux[:, 1],
            "TISSUE_BETA": beta_tissue,
            "TISSUE_SE": se_tissue,
            "TISSUE_Z": z_tissue,
            "TISSUE_PVAL": pval_tissue,
        }

        gmm_results = pd.DataFrame(res_data)
        gmm_results.to_parquet(
            os.path.join(save_dir, f"{gene}.parquet"), index=False, engine="pyarrow"
        )

        # Calculate Post-GMM Neff
        gene_res["TAR_CNEFF"] = calculate_Neff(gmm_z_tar[:, 0], Omega[0, 0], ld_tar)
        gene_res["TAR_TNEFF"] = calculate_Neff(gmm_z_tar[:, 1], Omega[0, 0], ld_tar)
        gene_res["AUX_CNEFF"] = calculate_Neff(gmm_z_aux[:, 0], Omega[1, 1], ld_aux)
        gene_res["AUX_TNEFF"] = calculate_Neff(gmm_z_aux[:, 1], Omega[1, 1], ld_aux)

        gene_res["TAR_CeSNP"] = count_eSNPs(gmm_p_tar[:, 0])
        gene_res["TAR_TeSNP"] = count_eSNPs(gmm_p_tar[:, 1])
        gene_res["AUX_CeSNP"] = count_eSNPs(gmm_p_aux[:, 0])
        gene_res["AUX_TeSNP"] = count_eSNPs(gmm_p_aux[:, 1])

        # Add to results
        summary_data.append((gene, gene_res))

    # End Gene Loop

    # Save Summary
    if summary_data:
        genes, data = zip(*summary_data)
        summary_df = pd.DataFrame(list(data), index=genes)
        summary_df.index.name = "GENE"
        summary_df["RUN_GMM"] = summary_df["RUN_GMM"].astype("boolean")
        summary_df.to_csv(os.path.join(save_dir, "summary.csv"))
    else:
        # Handle case where no genes were processed
        pd.DataFrame(columns=summary_columns).to_csv(
            os.path.join(save_dir, "summary.csv")
        )

    return


if __name__ == "__main__":
    args = parse_args()
    if isinstance(args.chromosome, int):
        chr_list = [int(args.chromosome)]
    elif args.chromosome is None:
        chr_list = []  # Or handle error
    else:
        chr_list = [int(i) for i in args.chromosome]

    start_time = time.time()
    print(
        f"Running GMM (Optimized) for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(start_time)}"
    )

    for chr_num in chr_list:
        main(args, chr_num)

    end_time = time.time()
    print(
        f"GMM finished for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(end_time)} which took {end_time - start_time:.2f} seconds"
    )
