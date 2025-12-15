# run_re2.py: Run RE2 meta-analysis for a given study and cell type.
import argparse
import time
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

z2p = lambda z: 2 * norm.sf(abs(z))
eSNP_THRESHOLD = 1e-5


def re2_meta(betas, ses, converge_tol=1e-6, max_iter=100):
    """
    Calculate RE2 meta-analysis for multiple populations using random-effects model with iterative tau**2 estimation

    :param betas: array of effect sizes (beta coefficients) for each study (n_studies,)
    :param ses: array of standard errors for each study (n_studies,)
    :param converge_tol: convergence tolerance for tau**2 and beta estimates
    :param max_iter: maximum number of iterations
    :return: tuple (beta: combined effect size, se: standard error of combined effect)
    """
    # \widehat{\mu}_{(n+1)} & =\frac{\sum \frac{X_i}{V_i+\hat{\tau}_{(n)}^2}}{\sum \frac{1}{V_i+\hat{\tau}_{(n)}^2}} \\
    # \hat{\tau}_{(n+1)}^2 & =\frac{\sum \frac{\left(X_i-\hat{\mu}_{(n+1)}\right)^2-V_i}{\left(V_i+\hat{\tau}_{(n)}^2\right)^2}}{\sum \frac{1}{\left(V_i+\hat{\tau}_{(n)}^2\right)^2}}
    variances = ses**2

    # Step 1: Initialize
    tau2 = 0.0
    weights = 1.0 / variances
    sum_weights = np.sum(weights)
    combined_beta = np.sum(weights * betas) / sum_weights

    # Initialize variables for iteration
    converged = False

    # Step 2: Iterate to estimate tau
    for i in range(max_iter):
        # Calculate weights for the current tau
        denom = variances + tau2
        weights = 1.0 / denom
        sum_weights = np.sum(weights)

        # Update combined beta using Han & Eskin formula
        new_combined_beta = np.sum(betas * weights) / sum_weights

        # Calculate new tau using Han & Eskin formula
        num_numerator = (betas - new_combined_beta) ** 2 - variances
        num_denom = denom**2
        numerator = np.sum(num_numerator / num_denom)
        denominator = np.sum(1.0 / num_denom)
        new_tau2 = numerator / denominator

        if new_tau2 < 0:  # If tau2 is negative, switch to fixed effects model
            tau2 = 0.0
            weights = 1.0 / variances
            sum_weights = np.sum(weights)
            combined_beta = np.sum(betas * weights) / sum_weights
            converged = True
            break

        beta_diff = abs(new_combined_beta - combined_beta)
        tau_diff = abs(new_tau2 - tau2)

        combined_beta = new_combined_beta
        tau2 = new_tau2

        if max(beta_diff, tau_diff) < converge_tol:
            converged = True
            break

    # Step 3: Calculate final combined beta and standard error
    if tau2 > 0:  # RE2 model
        weights = 1.0 / (variances + tau2)
    else:  # FE model
        weights = 1.0 / variances

    sum_weights = np.sum(weights)
    combined_beta = np.sum(betas * weights) / sum_weights
    combined_se = np.sqrt(1.0 / sum_weights)

    return combined_beta, combined_se


def parse_args():
    parser = argparse.ArgumentParser(description="Run RE2 meta-analysis.")
    parser.add_argument(
        "-s", "--studyid", help="eqtlCatalogue study id", type=str, required=True
    )
    parser.add_argument(
        "-t", "--celltype", help="cell type name", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        nargs="+",
        help="chromosome number",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-d", "--main_dir", help="main directory", type=str, required=True
    )
    return parser.parse_args()


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
    save_dir = os.path.join(study_dir, "RE2", f"chr{chr}")
    os.makedirs(save_dir, exist_ok=True)

    ### eQTL directories
    tarQTLfile = os.path.join(study_dir, f"TAR_{celltype}", f"chr{chr}.csv")
    auxQTLfile = os.path.join(study_dir, f"AUX_{celltype}", f"chr{chr}.csv")
    tissueQTLfile = os.path.join(study_dir, "Tissue", f"chr{chr}.csv")

    ## Load data
    try:
        tarQTLdf = pd.read_csv(tarQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
        auxQTLdf = pd.read_csv(auxQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
        tissueQTLdf = pd.read_csv(tissueQTLfile, header=0)  # RSID,GENE,BETA,SE,Z,PVAL,N
    except FileNotFoundError as e:
        print(f"Error loading data for chr{chr}: {e}")
        return

    ## run re2 for each gene
    # Use union of all genes present in the files
    gene_list = pd.concat(
        [tarQTLdf["GENE"], auxQTLdf["GENE"], tissueQTLdf["GENE"]]
    ).unique()

    summary_df = pd.DataFrame(
        data=np.nan,
        index=gene_list,
        columns=[
            "NSNP",
            "TAR_SeSNP",  # number of eSNP from summary stats
            "TAR_CeSNP",  # number of eSNP from cross-study RE2
            "TAR_TeSNP",  # number of eSNP from cross+tissue RE2
            "AUX_SeSNP",
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

        rename_cols = {
            col: f"{col}_tissue" for col in tissueQTLdfgene_df.columns if col != "RSID"
        }
        gene_df = (
            tarQTLgene_df.merge(auxQTLdfgene_df, on="RSID", suffixes=("_tar", "_aux"))
            .merge(tissueQTLdfgene_df.rename(columns=rename_cols), on="RSID")
            .set_index("RSID")
            .dropna(axis=0, how="any")
        )
        nsnp = gene_df.shape[0]
        summary_df.loc[gene, "NSNP"] = nsnp
        if nsnp == 0:
            continue

        ## run RE2 for each SNP in the gene
        re2_c_b = np.full(nsnp, np.nan)  # cross-study beta
        re2_c_se = np.full(nsnp, np.nan)  # cross-study se
        re2_t_b = np.full(nsnp, np.nan)  # cross+tissue beta
        re2_t_se = np.full(nsnp, np.nan)  # cross+tissue se

        for i in range(nsnp):
            # Analysis 1: Cross-study (tar + aux)
            betas_c = np.array([gene_df.BETA_tar.iloc[i], gene_df.BETA_aux.iloc[i]])
            ses_c = np.array([gene_df.SE_tar.iloc[i], gene_df.SE_aux.iloc[i]])
            re2_c_b[i], re2_c_se[i] = re2_meta(betas_c, ses_c)

            # Analysis 2: Cross-study + Tissue (tar + aux + tissue)
            betas_t = np.array(
                [
                    gene_df.BETA_tar.iloc[i],
                    gene_df.BETA_aux.iloc[i],
                    gene_df.BETA_tissue.iloc[i],
                ]
            )
            ses_t = np.array(
                [
                    gene_df.SE_tar.iloc[i],
                    gene_df.SE_aux.iloc[i],
                    gene_df.SE_tissue.iloc[i],
                ]
            )
            re2_t_b[i], re2_t_se[i] = re2_meta(betas_t, ses_t)

        re2_c_z = re2_c_b / re2_c_se
        re2_t_z = re2_t_b / re2_t_se
        re2_c_p = z2p(re2_c_z)
        re2_t_p = z2p(re2_t_z)

        ## save re2 results
        re2_results = pd.DataFrame(
            {
                "RSID": gene_df.index,
                "TAR_SBETA": gene_df.BETA_tar.to_numpy(),
                "TAR_SSE": gene_df.SE_tar.to_numpy(),
                "TAR_SPVAL": gene_df.PVAL_tar.to_numpy(),
                "TAR_CBETA": re2_c_b,
                "TAR_CSE": re2_c_se,
                "TAR_CZ": re2_c_z,
                "TAR_CPVAL": re2_c_p,
                "TAR_TBETA": re2_t_b,
                "TAR_TSE": re2_t_se,
                "TAR_TZ": re2_t_z,
                "TAR_TPVAL": re2_t_p,
                "AUX_BETA": gene_df.BETA_aux.to_numpy(),
                "AUX_SE": gene_df.SE_aux.to_numpy(),
                "AUX_PVAL": gene_df.PVAL_aux.to_numpy(),
                "TISSUE_BETA": gene_df.BETA_tissue.to_numpy(),
                "TISSUE_SE": gene_df.SE_tissue.to_numpy(),
                "TISSUE_PVAL": gene_df.PVAL_tissue.to_numpy(),
            }
        )
        re2_results.to_csv(os.path.join(save_dir, f"{gene}.csv"), index=False)

        ## count eSNPs
        summary_df.loc[gene, "TAR_SeSNP"] = count_eSNPs(gene_df.PVAL_tar.to_numpy())
        summary_df.loc[gene, "AUX_SeSNP"] = count_eSNPs(gene_df.PVAL_aux.to_numpy())
        summary_df.loc[gene, "TISSUE_SeSNP"] = count_eSNPs(
            gene_df.PVAL_tissue.to_numpy()
        )
        summary_df.loc[gene, "TAR_CeSNP"] = count_eSNPs(re2_c_p)
        summary_df.loc[gene, "TAR_TeSNP"] = count_eSNPs(re2_t_p)

    ## finish gene loop, save summary df
    summary_df.to_csv(os.path.join(save_dir, "summary.csv"))
    return


if __name__ == "__main__":
    args = parse_args()
    if isinstance(args.chromosome, int):
        chr_list = [args.chromosome]
    else:
        chr_list = [int(i) for i in args.chromosome]

    start_time = time.time()
    print(
        f"Running RE2 for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(start_time)}"
    )
    for chr in chr_list:
        main(args, chr)  # run RE2 for each chromosome
    end_time = time.time()
    print(
        f"RE2 finished for study {args.studyid} in cell type {args.celltype} on chromosomes {chr_list} at {time.ctime(end_time)} which took {end_time - start_time:.2f} seconds"
    )
