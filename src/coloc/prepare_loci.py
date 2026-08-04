"""Prepare GWAS loci and lead variants for colocalization analysis."""

# Find the leading SNP in each GWAS cytoband.
import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd


def argparser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwas", type=str, required=True, help="GWAS summary statistics.")
    parser.add_argument(
        "--cytobands",
        type=str,
        required=True,
        help="Genome-build-matched cytoband table.",
    )
    parser.add_argument(
        "--p-value-threshold",
        type=float,
        default=1e-6,
        help="Significance threshold used to select lead variants.",
    )
    parser.add_argument(
        "--gwas-format",
        choices=("auto", "standard", "bcx", "replication"),
        default="auto",
        help=(
            "Input schema. 'standard' expects "
            "RSID/PVAL/CHR/POS/REF/ALT/BETA/SE; "
            "'auto' recognizes BCX and replication filenames and otherwise "
            "uses the standard schema."
        ),
    )
    parser.add_argument(
        "--ldlink-script",
        type=str,
        default=str(Path(__file__).with_name("query_ldlink.R")),
    )
    parser.add_argument(
        "--bedtools-closest",
        type=str,
        default="closestBed",
    )
    parser.add_argument(
        "--genes",
        type=str,
        required=True,
        help="Protein-coding gene intervals in BED format.",
    )
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--output-prefix", type=str, required=True)
    return parser.parse_args()


def read_format_GWAS_standard(
    gwas_sumstats_path, p_val_threshold, output_dir, output_prefix
):
    # output: gwas_sig_df: significant SNPs in GWAS
    # also format GWAS data for coloc
    # 'VARIANT', 'CHR', 'POS', 'REF', 'ALT', 'RSID', 'BETA', 'SE', 'Z',
    #    'PVAL', 'MAF', 'N'
    gwas_df = pd.read_csv(gwas_sumstats_path)
    gwas_df = gwas_df.rename(
        columns={
            "RSID": "SNP",
            "PVAL": "PVALUE",
        }
    )
    # save GWAS data for coloc
    gwas_df.to_csv(
        os.path.join(output_dir, f"{output_prefix}_GWAS.csv"),
        sep=",",
        header=True,
        index=False,
    )
    # select significant SNPs
    gwas_sig_df_tmp = gwas_df[gwas_df["PVALUE"] < p_val_threshold].copy()
    print(
        f"Using p-value threshold of {p_val_threshold} and {len(gwas_sig_df_tmp)} significant SNPs"
    )
    gwas_sig_df = gwas_sig_df_tmp[["SNP", "PVALUE", "CHR", "POS", "REF", "ALT"]]
    return gwas_sig_df


def read_format_GWAS_bcx(
    gwas_sumstats_path, p_val_threshold, output_dir, output_prefix
):
    # output: gwas_sig_df: significant SNPs in GWAS
    # also format GWAS data for coloc
    #     rs_number	reference_allele	other_allele	eaf	beta	se	beta_95L	beta_95U	z	p-value	_-log10_p-value	q_statistic	q_p-value	i2	n_studies	n	effects	rsid
    # 10:10000018_A_G	A	G	0.654829	-0.006301	0.005041	-0.016181	0.003579	-1.249934	0.211303	0.675094	0.123154	0.725639	0.0	2	89266	--	rs6602381
    gwas_df = pd.read_csv(gwas_sumstats_path, compression="gzip", sep="\t")
    if "z" not in gwas_df.columns:
        gwas_df["z"] = gwas_df["beta"] / gwas_df["se"]
    # gwas_df = gwas_df.loc[gwas_df["eaf"] > 0.01, :]
    gwas_df4coloc = gwas_df[["rsid", "z"]]
    gwas_df4coloc = gwas_df4coloc.rename(columns={"rsid": "SNP", "z": "Z"})
    gwas_df4coloc = gwas_df4coloc.assign(
        CHR=gwas_df["rs_number"].str.split(":").str[0].astype(int),
        POS=gwas_df["rs_number"]
        .str.split(":")
        .str[1]
        .str.split("_")
        .str[0]
        .astype(int),
        MAF=gwas_df["eaf"],
        BETA=gwas_df["beta"],
        SE=gwas_df["se"],
    )
    gwas_df4coloc.to_csv(
        os.path.join(output_dir, f"{output_prefix}_GWAS.csv"),
        sep=",",
        header=True,
        index=False,
    )
    gwas_sig_df_tmp = gwas_df[gwas_df["p-value"] < p_val_threshold].copy()
    if len(gwas_sig_df_tmp) >= 5e4:
        p_val_threshold = 5e-8
        gwas_sig_df_tmp = gwas_df[gwas_df["p-value"] < p_val_threshold].copy()
    print(
        f"Using p-value threshold of {p_val_threshold} and {len(gwas_sig_df_tmp)} significant SNPs"
    )

    gwas_sig_df_tmp.loc[:, "chr"] = (
        gwas_sig_df_tmp["rs_number"].str.split(":").str[0]
    ).astype(int)
    gwas_sig_df_tmp.loc[:, "pos"] = (
        gwas_sig_df_tmp["rs_number"].str.split(":").str[1].str.split("_").str[0]
    )
    gwas_sig_df = gwas_sig_df_tmp[
        ["rsid", "p-value", "chr", "pos", "reference_allele", "other_allele"]
    ]
    gwas_sig_df.columns = ["SNP", "PVALUE", "CHR", "POS", "REF", "ALT"]
    gwas_sig_df.loc[:, "POS"] = gwas_sig_df["POS"].astype(int)
    return gwas_sig_df


def read_replicate_eQTL(
    eQTL_path, p_val_threshold, output_dir, output_prefix
):
    #     chr,pos,variant_id,ref,alt,gene,gene_name,beta,se,pval,pip_susie,rsid
    # 1,100000012,1_100000012_G_T,G,T,ENSG00000162688,AGL,0.0505631,0.0160197,0.00164814,0.0011359478923747,rs10875231
    eQTL_df = pd.read_csv(eQTL_path, sep=",")
    eQTL_df.columns = eQTL_df.columns.str.upper()
    eQTL_df = eQTL_df.rename(columns={"RSID": "SNP", "PVAL": "PVALUE"})
    eQTL_df.to_csv(
        os.path.join(output_dir, f"{output_prefix}_GWAS.csv"), index=False
    )
    eQTL_sig_df = eQTL_df[eQTL_df["PVALUE"] < p_val_threshold].copy()
    print(
        f"Using p-value threshold of {p_val_threshold} and {len(eQTL_sig_df)} significant eQTLs"
    )
    eQTL_sig_df = eQTL_sig_df[["SNP", "PVALUE", "CHR", "POS", "REF", "ALT"]]
    return eQTL_sig_df


def load_cytoband(cytoband_path):
    # output: cytoband_df: cytoband info
    # chr1	0	2300000	p36.33	gneg
    cytoband_df = pd.read_csv(cytoband_path, sep="\t", header=None)
    cytoband_df.columns = ["chr", "start", "end", "band", "stain"]
    # exclude chrX and chrY
    cytoband_df = cytoband_df[cytoband_df["chr"].str.contains("chr[0-9]+")]
    cytoband_df.loc[:, "chr"] = cytoband_df["chr"].str.replace("chr", "").astype(int)
    return cytoband_df


def find_band(gwas_sig_df, cytoband_df):
    # output: gwas_sig_df: significant SNPs in GWAS with band info
    gwas_sig_df = gwas_sig_df.copy()
    gwas_sig_df.loc[:, "band"] = ""
    for idx, row in gwas_sig_df.iterrows():
        chr = row["CHR"]
        pos = row["POS"]
        band = cytoband_df[
            (cytoband_df["chr"] == chr)
            & (cytoband_df["start"] < pos)
            & (cytoband_df["end"] > pos)
        ]
        if band.shape[0] > 0:
            gwas_sig_df.loc[idx, "band"] = band.iloc[-1, :]["band"]
    return gwas_sig_df


def find_lead_variants(gwas_sig_df, save_path, save_prefix):
    # output: leadingSNP_df: leading SNP in each band of GWAS
    gwas_sig_df = gwas_sig_df.sort_values(by="PVALUE")
    leadingSNP_df = gwas_sig_df.drop_duplicates(subset=["CHR", "band"], keep="first")
    leadingSNP_df.to_csv(
        f"{save_path}/{save_prefix}_lead_variants.csv",
        sep=",",
        header=True,
        index=False,
    )
    print(f"Found {len(leadingSNP_df)} bands with leading SNPs")
    return leadingSNP_df


def run_ldlink_query(rcode_path, leadingSNP_path, save_prefix, save_path):
    subprocess.run(
        [
            "Rscript",
            rcode_path,
            leadingSNP_path,
            os.path.join(save_path, f"{save_prefix}_ldlink.tsv"),
        ],
        check=True,
    )


def create_index_snp_bed(ldlink_path, save_prefix, save_path):
    ldlink_loci_df = pd.read_csv(ldlink_path, sep="\t", header=0)
    positions = ldlink_loci_df["GWAS_snp_pos"].str.split(":").str[1].astype(int)
    index_snp_bed_df = pd.DataFrame(
        {
            "chr": ldlink_loci_df["chrom"],
            "start": positions - 1,
            "end": positions,
            "rsid": ldlink_loci_df["GWAS_snp"],
            "score": 0,
        }
    )
    index_snp_bed_df = index_snp_bed_df.sort_values(by=["chr", "start"])
    index_snp_bed_df.to_csv(
        os.path.join(save_path, f"{save_prefix}_GWAS_index_snps.bed"),
        sep="\t",
        header=False,
        index=False,
    )


def run_closestBed(
    closestBed_path, gene_bed_path, index_snp_bed_path, prefix, save_path
):
    with open(
        os.path.join(save_path, f"{prefix}.closest.protein_coding.bed"), "w"
    ) as f:
        subprocess.run(
            [
                closestBed_path,
                "-d",
                "-wa",
                "-a",
                index_snp_bed_path,
                "-b",
                gene_bed_path,
            ],
            stdout=f,
            check=True,
        )


def annot_closestBed(cloestBed_path, ldlink_path, save_path, prefix):
    # GWAS_snp	GWAS_snp_pos	locus_name	chrom	start	end	locus_name_gene	ensembl top_pval
    cloestBed_df = pd.read_csv(cloestBed_path, sep="\t", header=None)
    cloestBed_df.columns = [
        "chr_snp",
        "start_snp",
        "end_snp",
        "GWAS_snp",
        "beta",
        "chr_gene",
        "start_gene",
        "end_gene",
        "gene",
        "distance",
    ]
    cloestBed_df["GWAS_snp_pos"] = (
        cloestBed_df["chr_snp"].astype(str)
        + ":"
        + cloestBed_df["end_snp"].astype(str)
    )
    cloestBed_df = cloestBed_df[["GWAS_snp", "GWAS_snp_pos", "gene", "distance"]]
    cloestBed_df[["symbol", "ensembl"]] = cloestBed_df["gene"].str.split(
        "_", expand=True
    )
    result_df = cloestBed_df[["GWAS_snp_pos", "ensembl"]]

    ldlink_df = pd.read_csv(ldlink_path, sep="\t", header=0)
    # GWAS_snp	GWAS_snp_pos	locus_name	chrom	start	end
    ldlink_df = ldlink_df.merge(result_df, on="GWAS_snp_pos", how="inner")
    ldlink_df.to_csv(
        os.path.join(save_path, f"{prefix}_loci.csv"),
        sep=",",
        header=True,
        index=False,
    )


if __name__ == "__main__":
    args = argparser()
    gwas_sumstats_path = args.gwas
    cytoband_path = args.cytobands
    p_val_threshold = args.p_value_threshold
    save_path = args.output_dir
    save_prefix = args.output_prefix
    ldlinkr_src = args.ldlink_script
    closestBed_path = args.bedtools_closest
    gencode_annotation_path = args.genes

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    try:
        gwas_format = args.gwas_format
        if gwas_format == "auto":
            input_name = Path(gwas_sumstats_path).name.lower()
            if "bcx" in input_name:
                gwas_format = "bcx"
            elif "hum" in input_name:
                gwas_format = "replication"
            else:
                gwas_format = "standard"

        if gwas_format == "standard":
            gwas_sig_df = read_format_GWAS_standard(
                gwas_sumstats_path,
                p_val_threshold,
                save_path,
                save_prefix,
            )
        elif gwas_format == "bcx":
            gwas_sig_df = read_format_GWAS_bcx(
                gwas_sumstats_path,
                p_val_threshold,
                save_path,
                save_prefix,
            )
        else:
            gwas_sig_df = read_replicate_eQTL(
                gwas_sumstats_path,
                p_val_threshold,
                save_path,
                save_prefix,
            )
        cytoband_df = load_cytoband(cytoband_path)
        gwas_sig_df = find_band(gwas_sig_df, cytoband_df)
        find_lead_variants(gwas_sig_df, save_path, save_prefix)
        run_ldlink_query(
            ldlinkr_src,
            os.path.join(save_path, f"{save_prefix}_lead_variants.csv"),
            save_prefix,
            save_path,
        )
        create_index_snp_bed(
            os.path.join(save_path, f"{save_prefix}_ldlink.tsv"),
            save_prefix,
            save_path,
        )
        run_closestBed(
            closestBed_path,
            gencode_annotation_path,
            os.path.join(save_path, f"{save_prefix}_GWAS_index_snps.bed"),
            save_prefix,
            save_path,
        )
        annot_closestBed(
            os.path.join(save_path, f"{save_prefix}.closest.protein_coding.bed"),
            os.path.join(save_path, f"{save_prefix}_ldlink.tsv"),
            save_path,
            save_prefix,
        )

        print("For GWAS data at", gwas_sumstats_path)
        print(
            "Created GWAS file for coloc at",
            os.path.join(save_path, f"{save_prefix}_GWAS.csv"),
        )
        print(
            "Created lead-variant file at", f"{save_path}/{save_prefix}_lead_variants.csv"
        )
        print(
            "Created closest protein coding file at",
            os.path.join(save_path, f"{save_prefix}.closest.protein_coding.bed"),
        )
        print(
            "Created final loci file for coloc at",
            os.path.join(save_path, f"{save_prefix}_loci.csv"),
        )
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        raise
    except pd.errors.EmptyDataError:
        print("Error: Input file is empty")
        raise
    except subprocess.CalledProcessError as e:
        print(f"Error running external process: {e}")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise
