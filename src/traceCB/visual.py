# visualize eSNP eGene NEFF
# cmd: python -m pdb src/visual.py -d /gpfs1/scratch/wjiang49/xpmm/EAS_GTEx/ -s QTD000067 -t CD4+T_cells -tpn 100 -apn 400 -tin 1000
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import argparse
import os
import glob
import upsetplot
from utils import P_VAL_THRED_Z

MAX_NEFF_RATIO = 10.0
method_list = ["Original", "GMM", "GMM+"]

# warning comes from upsetplot, ignore it


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--main_dir",
        help="main directory",
        type=str,
    )
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
        "-tpn",
        "--target_samplesize",
        help="target population sample size",
        type=int,
    )
    parser.add_argument(
        "-apn",
        "--auxiliary_samplesize",
        help="auxiliary population sample size",
        type=int,
    )
    parser.add_argument(
        "-tin",
        "--tissue_samplesize",
        help="tissue population sample size",
        type=int,
    )

    return parser.parse_args()


def load_summary(study_dir):
    """
    Load summary data from GMM results.
    returns summary_sign_df: DataFrame with significant correlations
    returns summary_df: all results DataFrame
    """
    #     <save_path_main>/EAS_GTEx/QTD@/GMM/chr@/summary.csv
    # GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
    summary_dirs = glob.glob(os.path.join(study_dir, "GMM", "chr*", "summary.csv"))
    summary_df = pd.DataFrame()
    for summary_file in summary_dirs:
        df = pd.read_csv(summary_file)
        df["CHR"] = os.path.basename(os.path.dirname(summary_file)).replace("chr", "")
        summary_df = (
            df if summary_df.empty else pd.concat([summary_df, df], ignore_index=True)
        )
    summary_df.loc[:, "COR"] = (
        summary_df.COV / summary_df.H1SQ.apply(np.sqrt) / summary_df.H2SQ.apply(np.sqrt)
    )  # convert correlation to genetic correlation
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
    return summary_sign_df, summary_df


def plot_neff_violin(
    summary_df,
    target_samplesize,
    auxiliary_samplesize,
    tissue_samplesize,
    title,
    save_path,
):
    # NOTE: filter significant h2 genes and clip
    tar_df = summary_df.loc[
        summary_df.H1SQ / summary_df.H1SQSE
        > P_VAL_THRED_Z,  # filter significant h2 genes
        ["TAR_SNEFF", "TAR_CNEFF", "TAR_TNEFF"],
    ]
    tmp = summary_df.H2SQ / summary_df.H2SQSE > P_VAL_THRED_Z
    aux_df = summary_df.loc[
        tmp, ["AUX_SNEFF", "AUX_CNEFF", "AUX_TNEFF"]  # filter significant h2 genes
    ]
    tissue_df = summary_df.loc[tmp, "TISSUE_SNEFF"].to_frame()

    # clip effective sample sizes to make plot looks better
    tar_df.TAR_SNEFF = tar_df.TAR_SNEFF.clip(
        upper=MAX_NEFF_RATIO * target_samplesize, lower=0
    )
    tar_df.TAR_CNEFF = tar_df.TAR_CNEFF.clip(
        upper=MAX_NEFF_RATIO * auxiliary_samplesize, lower=0
    )
    tar_df.TAR_TNEFF = tar_df.TAR_TNEFF.clip(
        upper=MAX_NEFF_RATIO * tissue_samplesize, lower=0
    )
    aux_df.AUX_SNEFF = aux_df.AUX_SNEFF.clip(
        upper=MAX_NEFF_RATIO * auxiliary_samplesize, lower=0
    )
    aux_df.AUX_CNEFF = aux_df.AUX_CNEFF.clip(
        upper=MAX_NEFF_RATIO * auxiliary_samplesize, lower=0
    )
    aux_df.AUX_TNEFF = aux_df.AUX_TNEFF.clip(
        upper=MAX_NEFF_RATIO * tissue_samplesize, lower=0
    )
    tissue_df.TISSUE_SNEFF = tissue_df.TISSUE_SNEFF.clip(
        upper=MAX_NEFF_RATIO * tissue_samplesize, lower=0
    )

    # rename columns and melt DataFrames for plotting
    tar_df.columns = method_list
    aux_df.columns = method_list
    tissue_df.columns = ["Original"]
    tar_melt_df = tar_df.melt(var_name="Method", value_name="NEFF")
    aux_melt_df = aux_df.melt(var_name="Method", value_name="NEFF")
    tissue_melt_df = tissue_df.melt(var_name="Method", value_name="NEFF")

    # begin plotting
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.42])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])

    # target population violin plot
    sns.violinplot(
        x="Method",
        y="NEFF",
        data=tar_melt_df,
        ax=ax0,
        inner="box",
        cut=0,
    )
    ax0.axhline(
        target_samplesize, ls="-", color="gray", label=f"N_study={target_samplesize}"
    )  # add sample size line
    ax0.legend(loc="upper left")
    ax0.set_title("Target Population")
    ax0.set_xlabel("")
    ax0.set_ylabel("Effective Sample Size")
    # auxiliary population violin plot
    sns.violinplot(
        x="Method",
        y="NEFF",
        data=aux_melt_df,
        ax=ax1,
        inner="box",
        cut=0,
    )
    ax1.axhline(
        auxiliary_samplesize,
        ls="-",
        color="gray",
        label=f"N_study={auxiliary_samplesize}",
    )  # add sample size line
    ax1.legend(loc="upper left")
    ax1.set_title("Auxiliary Population")
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    # tissue population violin plot
    sns.violinplot(
        x="Method",
        y="NEFF",
        data=tissue_melt_df,
        ax=ax2,
        inner="box",
        cut=0,
    )
    ax2.axhline(
        tissue_samplesize, ls="-", color="gray", label=f"N_study={tissue_samplesize}"
    )  # add sample size line
    ax2.legend(loc="upper left")
    ax2.set_title("Auxiliary Tissue")
    ax2.set_xlabel("")
    ax2.set_ylabel("")

    fig.suptitle(f"Effective Sample Size of {title}")
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"Neff_violin_{title}.png"))
    plt.close()


# Define correlation groups
cor_group = [
    "0-0.1",
    "0.1-0.2",
    "0.2-0.3",
    "0.3-0.4",
    "0.4-0.5",
    "0.5-0.6",
    "0.6-0.7",
    "0.7-0.8",
    "0.8-0.9",
    "0.9-1",
]


def plot_cor_box(summary_df, title, save_path):
    summary_df = summary_df.copy()
    # calculate cross/sumstat and cross_tissue/sumstat
    summary_df.loc[:, method_list[1] + "/" + method_list[0]] = (
        summary_df.TAR_CNEFF / summary_df.TAR_SNEFF
    )
    summary_df.loc[:, method_list[2] + "/" + method_list[0]] = (
        summary_df.TAR_TNEFF / summary_df.TAR_SNEFF
    )

    # Bin the COR values into groups
    summary_df.loc[:, "COR_group"] = pd.cut(
        summary_df.COR,
        bins=[-0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        labels=cor_group,
    )

    # Ensure all categories are present
    summary_df.COR_group = summary_df.COR_group.cat.add_categories(
        [cat for cat in cor_group if cat not in summary_df.COR_group.cat.categories]
    ).cat.reorder_categories(cor_group, ordered=True)
    # Reorder categories
    summary_df.COR_group = summary_df.COR_group.cat.reorder_categories(
        cor_group, ordered=True
    )
    # Convert data to long format for plotting
    melt_df = summary_df.melt(
        id_vars=["COR_group"],
        value_vars=[
            method_list[1] + "/" + method_list[0],
            method_list[2] + "/" + method_list[0],
        ],
        var_name="Method",
        value_name="Value",
    )

    # begin plotting
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True, gridspec_kw={"height_ratios": [1, 3]}
    )

    # bar plot for number of genes in each COR group
    cor_group_counts = summary_df.COR_group.value_counts().sort_index()
    sns.barplot(
        x=cor_group_counts.index, y=cor_group_counts.values, color="skyblue", ax=ax1
    )
    ax1.bar_label(ax1.containers[0], fontsize=10)
    ax1.set_ylabel("Number of genes")
    ax1.set_xlabel("")
    ax1.set_ylim(
        0, ax1.get_ylim()[1] + 150
    )  # add some space above the bars for annotations

    # box plot for effective sample size ratios
    sns.boxplot(
        x="COR_group", y="Value", hue="Method", data=melt_df, fliersize=0, ax=ax2
    )
    ax2.set_ylim([0, 5])
    ax2.set_ylabel("Ratio of effective sample size")
    ax2.set_xlabel("COR Group")
    ax2.legend(title="Ratio", loc="upper left")

    fig.suptitle(f"Effective sample size ratio by different correlation in {title}")
    # Rotate x-ticks for better readability
    for ax in [ax1, ax2]:
        plt.sca(ax)
        plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"cor_box_{title}.png"))
    plt.close()
    return summary_df


def plot_egene_cor_bar(summary_df, save_path, title):
    melted_df = pd.melt(
        summary_df,
        id_vars=["COR_group"],
        value_vars=["TAR_SeSNP", "TAR_CeSNP", "TAR_TeSNP"],
        var_name="Method",
        value_name="eSNP_count",
    )

    # filter out egenes with zero eSNP counts
    melted_df = melted_df[melted_df["eSNP_count"] > 0]
    # rename values to method_list
    melted_df["Method"] = melted_df["Method"].replace(
        {
            "TAR_SeSNP": method_list[0],
            "TAR_CeSNP": method_list[1],
            "TAR_TeSNP": method_list[2],
        }
    )

    # count the number of eSNPs in each COR group for each method
    plot_df = (
        melted_df.groupby(["COR_group", "Method"], observed=False)
        .size()
        .reset_index(name="Count")
    )

    # count total eSNPs in each COR group
    cor_group_counts = summary_df["COR_group"].value_counts().sort_index().reset_index()
    cor_group_counts.columns = ["COR_group", "Total"]
    plot_df = plot_df.merge(cor_group_counts, on="COR_group", how="left")
    plot_df.loc[:, "Proportion"] = plot_df.Count / plot_df.Total

    # begin plotting
    fig, ax1 = plt.subplots(figsize=(12, 6))
    sns.set_style("ticks")
    sns.barplot(
        x="COR_group",
        y="Count",
        hue="Method",
        data=plot_df,
        palette="viridis",
        hue_order=method_list,
        ax=ax1,
    )

    ax1.set_title(f"eGene Detection and Proportion by COR Group of {title}")
    ax1.set_xlabel("Genetic Correlation (COR) Group")
    ax1.set_ylabel("Number of eGenes")
    ax1.tick_params(axis="x", rotation=45)

    # second y-axis for proportion
    ax2 = ax1.twinx()
    sns.lineplot(
        x="COR_group",
        y="Proportion",
        hue="Method",
        data=plot_df,
        marker="o",
        ax=ax2,
        palette="viridis",
        hue_order=method_list,
        legend=False,  # 避免重复图例
    )

    ax2.set_ylabel("Proportion of eGenes")
    ax2.set_ylim(0, 1)

    # Add legend for the first axis
    handles1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(handles1, labels1, title="Method", loc="upper right")

    plt.tight_layout()
    plt.savefig(os.path.join(save_path, f"egene_cor_bar_proportion_{title}.png"))
    plt.close()


def format_eGene_df(df):
    gene_indicators = pd.DataFrame(
        {
            method_list[0]: df["TAR_SeSNP"] > 0,
            method_list[1]: df["TAR_CeSNP"] > 0,
            method_list[2]: df["TAR_TeSNP"] > 0,
        }
    )
    result = gene_indicators.value_counts().sort_index()
    # remove rows with all False
    result = result.drop((False, False, False))
    return result.to_frame(name="value")


def plot_upset(df, title, save_file):
    plt.figure(figsize=(10, 6))
    upsetplot.plot(df["value"], sort_categories_by="-input", fig=plt.gcf())
    plt.suptitle(title)
    plt.savefig(save_file)
    plt.close()


if __name__ == "__main__":
    args = parse_args()
    main_dir = args.main_dir
    studyid = args.studyid
    celltype = args.celltype
    target_samplesize = args.target_samplesize
    auxiliary_samplesize = args.auxiliary_samplesize
    tissue_samplesize = args.tissue_samplesize

    study_dir = os.path.join(main_dir, studyid)
    save_dir = os.path.join(study_dir, "results", "visualization")
    os.makedirs(save_dir, exist_ok=True)

    summary_sign_df, summary_df = load_summary(study_dir)
    title = f"{studyid} {celltype}"
    plot_neff_violin(
        summary_df,
        target_samplesize,
        auxiliary_samplesize,
        tissue_samplesize,
        title,
        save_dir,
    )
    plot_neff_violin(
        summary_sign_df,  # significant correlations only
        target_samplesize,
        auxiliary_samplesize,
        tissue_samplesize,
        title + " significant cor",
        save_dir,
    )

    summary_sign_df_annotated = plot_cor_box(summary_sign_df, title, save_dir)
    plot_egene_cor_bar(summary_sign_df_annotated, save_dir, title)

    egene_df = format_eGene_df(summary_sign_df)
    plot_upset(
        egene_df,
        f"eGene of {title}",
        os.path.join(save_dir, f"egene_{title}.png"),
    )
    # # Filter out total h2 <= 0.001
    # summary_sign_df = summary_sign_df.copy()
    # summary_sign_df.loc[:, "H1SQ_TOTAL"] = summary_sign_df.H1SQ * summary_sign_df.NSNP
    # summary_sign_df.loc[:, "H2SQ_TOTAL"] = summary_sign_df.H2SQ * summary_sign_df.NSNP
    # summary_sign_high_h2_df = summary_sign_df[
    #     (summary_sign_df.H1SQ_TOTAL > 0.001) & (summary_sign_df.H2SQ_TOTAL > 0.001)
    # ]
    # summary_sign_high_h2_df_annotated = plot_cor_box(
    #     summary_sign_high_h2_df, title, save_dir
    # )
    # plot_egene_cor_bar(summary_sign_high_h2_df_annotated, save_dir, title)
