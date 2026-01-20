# find leading SNP in each band of gwas
import numpy as np
import pandas as pd
import os
import subprocess
import argparse

# # for bcx trait:
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_LYM_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_lym --save_prefix bcx_lym > /data/coloc/bcx_lym/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_NEU_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_neu --save_prefix bcx_neu > /data/coloc/bcx_neu/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_MON_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_mon --save_prefix bcx_mon > /data/coloc/bcx_mon/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_WBC_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_wbc --save_prefix bcx_wbc > /data/coloc/bcx_wbc/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_MCHC_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_mchc --save_prefix bcx_mchc > /data/coloc/bcx_mchc/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_RBC_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_rbc --save_prefix bcx_rbc > /data/coloc/bcx_rbc/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_MCV_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_mcv --save_prefix bcx_mcv > /data/coloc/bcx_mcv/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_MPV_EAS_UKBB_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_mpv --save_prefix bcx_mpv > /data/coloc/bcx_mpv/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_HGB_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_hgb --save_prefix bcx_hgb > /data/coloc/bcx_hgb/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_HCT_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_hct --save_prefix bcx_hct > /data/coloc/bcx_hct/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_EOS_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_eos --save_prefix bcx_eos > /data/coloc/bcx_eos/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_BAS_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_bas --save_prefix bcx_bas > /data/coloc/bcx_bas/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_PLT_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_plt --save_prefix bcx_plt > /data/coloc/bcx_plt/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BCX/BCX2_MCH_EAS_GWAMA_mapped.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bcx_mch --save_prefix bcx_mch > /data/coloc/bcx_mch/find_leadingSNP.log

# # for bbj trait:
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BBJ_GWAS/formatted/RA.csv.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bbj_ra --save_prefix bbj_ra > /data/coloc/bbj_ra/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/BBJ_GWAS/formatted/Asthma.csv.gz --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bbj_asthma --save_prefix bbj_asthma > /data/coloc/bbj_asthma/find_leadingSNP.log
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/coloc/mama_atopy/Atopy_bbj.csv --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 1e-6 --save_path /data/coloc/bbj_atopy --save_prefix bbj_atopy > /data/coloc/bbj_atopy/find_leadingSNP.log

# replicate eQTL:
# python src/coloc/find_leadingSNP.py --gwas_sumstats_path /data/hum0343.v3.qtl.v1/hum0304.csv --cytoband_path /data/coloc/hg19_cytoBand.txt --p_val_threshold 5e-8 --save_path /data/coloc/replicate --save_prefix replicate > /data/coloc/replicate/find_leadingSNP.log


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas_sumstats_path", type=str, required=True)
    parser.add_argument(
        "--cytoband_path",
        type=str,
        default="/data/coloc/hg19_cytoBand.txt",
    )
    parser.add_argument(
        "--p_val_threshold",
        "-p",
        type=float,
        default=1e-6,
        help="p-value threshold for significant SNPs. Should be 5e-8 for BCX and 1e-6 for BBJ",
    )
    parser.add_argument(
        "--ldlinkr_src",
        type=str,
        default="traceCB/src/coloc/run_ldlinkr.r",
    )
    parser.add_argument(
        "--closestBed_path",
        type=str,
        default="/src/coloc/bedtools2/bin/closestBed",
    )
    parser.add_argument(
        "--gencode_annotation_path",
        type=str,
        default="/data/coloc/gencode.v39lift37.annotation.protein_coding.1_22.bed",
    )
    parser.add_argument("--save_path", type=str, required=True)
    parser.add_argument("--save_prefix", type=str, required=True)
    return parser.parse_args()


def read_format_GWAS_bbj(gwas_sumstats_path, p_val_threshold):
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
        f"{save_path}/{save_prefix}_GWAS.csv",
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


def read_format_GWAS_bcx(gwas_sumstats_path, p_val_threshold):
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
        f"{save_path}/{save_prefix}_GWAS.csv",
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


def read_replicate_eQTL(eQTL_path, p_val_threshold):
    #     chr,pos,variant_id,ref,alt,gene,gene_name,beta,se,pval,pip_susie,rsid
    # 1,100000012,1_100000012_G_T,G,T,ENSG00000162688,AGL,0.0505631,0.0160197,0.00164814,0.0011359478923747,rs10875231
    eQTL_df = pd.read_csv(eQTL_path, sep=",")
    eQTL_df.columns = eQTL_df.columns.str.upper()
    eQTL_df = eQTL_df.rename(columns={"RSID": "SNP", "PVAL": "PVALUE"})
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


def find_leadingSNP(gwas_sig_df, save_path, save_prefix):
    # output: leadingSNP_df: leading SNP in each band of GWAS
    gwas_sig_df = gwas_sig_df.sort_values(by="PVALUE")
    leadingSNP_df = gwas_sig_df.drop_duplicates(subset=["CHR", "band"], keep="first")
    leadingSNP_df.to_csv(
        f"{save_path}/{save_prefix}_leadingSNP.csv",
        sep=",",
        header=True,
        index=False,
    )
    print(f"Found {len(leadingSNP_df)} bands with leading SNPs")
    return leadingSNP_df


def run_ldlinkr(rcode_path, leadingSNP_path, save_prefix, save_path):
    subprocess.run(
        [
            "Rscript",
            rcode_path,
            leadingSNP_path,
            os.path.join(save_path, f"{save_prefix}_ldlink.csv"),
        ]
    )


def create_index_snp_bed(ldlink_path, save_prefix, save_path):
    ldlink_loci_df = pd.read_csv(ldlink_path, sep="\t", header=0)
    # sort by chr and start
    ldlink_loci_df.loc[:, "chr"] = ldlink_loci_df["chrom"].str.replace("chr", "")
    ldlink_loci_df.head()

    index_snp_bed_df = pd.DataFrame(
        {
            "chr": ldlink_loci_df["chrom"],
            "pos": ldlink_loci_df["GWAS_snp_pos"].str.split(":").str[1],
            "pos_": ldlink_loci_df["GWAS_snp_pos"].str.split(":").str[1],
            "rsid": ldlink_loci_df["GWAS_snp"],
            "score": 0,
        }
    )
    index_snp_bed_df["pos"] = index_snp_bed_df["pos"].astype(int)
    index_snp_bed_df = index_snp_bed_df.sort_values(by=["chr", "pos"])
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
        + cloestBed_df["start_snp"].astype(str)
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
        sep="\t",
        header=True,
        index=False,
    )


if __name__ == "__main__":
    args = argparser()
    gwas_sumstats_path = args.gwas_sumstats_path
    cytoband_path = args.cytoband_path
    p_val_threshold = args.p_val_threshold
    save_path = args.save_path
    save_prefix = args.save_prefix
    ldlinkr_src = args.ldlinkr_src
    closestBed_path = args.closestBed_path
    gencode_annotation_path = args.gencode_annotation_path

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    try:
        if "BBJ" in gwas_sumstats_path or "bbj" in gwas_sumstats_path:
            gwas_sig_df = read_format_GWAS_bbj(gwas_sumstats_path, p_val_threshold)
        elif "BCX" in gwas_sumstats_path:
            gwas_sig_df = read_format_GWAS_bcx(gwas_sumstats_path, p_val_threshold)
        elif "eur" in gwas_sumstats_path:
            gwas_sig_df = read_format_GWAS_eur(gwas_sumstats_path, p_val_threshold)
        elif "hum" in gwas_sumstats_path:
            gwas_sig_df = read_replicate_eQTL(gwas_sumstats_path, p_val_threshold)
        else:
            raise ValueError("Unknown GWAS data source")
        cytoband_df = load_cytoband(cytoband_path)
        gwas_sig_df = find_band(gwas_sig_df, cytoband_df)
        leadingSNP_df = find_leadingSNP(gwas_sig_df, save_path, save_prefix)
        run_ldlinkr(
            ldlinkr_src,
            os.path.join(save_path, f"{save_prefix}_leadingSNP.csv"),
            save_prefix,
            save_path,
        )
        create_index_snp_bed(
            os.path.join(save_path, f"{save_prefix}_ldlink.csv"),
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
            os.path.join(save_path, f"{save_prefix}_ldlink.csv"),
            save_path,
            save_prefix,
        )

        print("For GWAS data at", gwas_sumstats_path)
        print(
            "Created GWAS file for coloc at",
            os.path.join(save_path, f"{save_prefix}_GWAS.csv"),
        )
        print(
            "Created leading SNP file at", f"{save_path}/{save_prefix}_leadingSNP.csv"
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
    except pd.errors.EmptyDataError:
        print("Error: Input file is empty")
    except subprocess.CalledProcessError as e:
        print(f"Error running external process: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# TODO: create protein coding bed file from gencode annotation
# which use the following code:
# if(!file.exists('data/gwas/ad/closest.protein.coding.bed')){
#   gtf_b37 <- rtracklayer::import('data/gencode/gencode.v39lift37.annotation.gtf.gz') %>%
#     as.data.frame() %>%
#     as_tibble()
# }
# if(!file.exists('data/gencode/gencode.v39lift37.annotation.protein_coding.1_22.bed')){
#   protein_coding <- filter(gtf_b37,type=='gene',gene_type=='protein_coding') %>%
#     mutate(gene_label=paste0(gene_name,'_',gsub('\\..+','',gene_id))) %>%
#     dplyr::select(seqnames,start,end,gene_label) %>%
#     filter(seqnames %in% paste0('chr',1:22)) %>%
#     mutate(seqnames=as.character(seqnames)) %>%
#     arrange(seqnames,start) %>%
#     write_tsv('data/gencode/gencode.v39lift37.annotation.protein_coding.1_22.bed',col_names = FALSE)
# }
