#!/bin/bash

GTEx_main_path="/gpfs1/scratch/ResearchGroups/bios_mingxcai/data/GTEx"
GTEx_raw_file="${GTEx_main_path}/GTEx_Analysis_v8_QTLs-GTEx_Analysis_v8_eQTL_all_associations-Whole_Blood.allpairs.txt.gz"
# gene_id	variant_id	tss_distance	ma_samples	ma_count	maf	pval_nominal	slope	slope_se
# ENSG00000227232.5	chr1_13550_G_A_b38	-16003	19	19	0.0141791	0.734151	0.0587242	0.172837
# ENSG00000227232.5	chr1_14671_G_C_b38	-14882	17	17	0.0126866	0.876478	-0.0282343	0.181569
# ENSG00000227232.5	chr1_14677_G_A_b38	-14876	69	69	0.0514925	0.0487147	-0.185212	0.0937735
GTEx_lookup_file="${GTEx_main_path}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz"
# variant_id	chr	variant_pos	ref	alt	num_alt_per_site	rs_id_dbSNP151_GRCh38p7	variant_id_b37
# chr1_13526_C_T_b38	chr1	13526	C	T	1	rs1209314672	1_13526_C_T_b37
# chr1_13550_G_A_b38	chr1	13550	G	A	1	rs554008981	1_13550_G_A_b37
# chr17_13496299_A_T_b38	chr17	13496299	A	T	1	rs61741326,rs62057033	17_13399616_A_T_b37
# chr1_14451_CTCT_C_b38	chr1	14451	CTCT	C	1	.	1_14451_CTCT_C_b37

OUT_PATH="${GTEx_main_path}/GTEx_Whole_Blood_by_chr"
mkdir -p $OUT_PATH
# targrt structure: $OUT_PATH/chr1.csv

# Map according to variant_id, format to following columns:
# GENE  RSID    CHR  POS    TSS_DISTANCE    A1  A2  MAF PVAL    BETA    SE
# ENSG00000227232    rs1209314672    1    13550    -16003    T    C    0.0141791    0.734151    0.0587242    0.172837

header="GENE,RSID,CHR,POS,TSS_DISTANCE,A1,A2,MAF,PVAL,BETA,SE"

# Create temporary files
temp_dir=$(mktemp -d)
echo "Using temporary directory: $temp_dir"

# Check input files exist
[ ! -f "$GTEx_raw_file" ] && echo "Error: $GTEx_raw_file not found" && exit 1
[ ! -f "$GTEx_lookup_file" ] && echo "Error: $GTEx_lookup_file not found" && exit 1

# Create lookup dictionary
echo "Creating lookup dictionary..."
zcat "$GTEx_lookup_file" | tail -n +2 |
    awk 'length($4)==1 && length($5)==1 && $7 !~ /,/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7}' |
    sort -k1,1 >"$temp_dir/lookup_sorted.txt"
# variant_id    chr	variant_pos    ref	alt	rs_id_dbSNP151_GRCh38p7
# head -n 3 "$temp_dir/lookup_sorted.txt"
# chr10_10000019_G_C_b38	chr10	10000019	G	C	rs571911750
# chr10_100000222_G_A_b38	chr10	100000222	G	A	rs138880521
# chr10_100000235_C_T_b38	chr10	100000235	C	T	rs11596870

for chr in {1..22}; do
(
    echo "Processing chromosome $chr"
    output_file="$OUT_PATH/chr${chr}.csv"

    # Write header
    echo "GENE,RSID,CHR,POS,TSS_DISTANCE,A1,A2,MAF,PVAL,BETA,SE" >"$output_file"

    # Process raw file for current chromosome
    zcat "$GTEx_raw_file" | tail -n +2 |
        awk -v chr="$chr" '
    BEGIN {FS="\t"; OFS="\t"}
    $2 ~ "^chr"chr"[^0-9]" || $2 ~ "^chr"chr"$" {
        gene=$1
        sub(/\.[0-9]+$/, "", gene)
        variant=$2
        print variant, gene, $3, $6, $7, $8, $9
    }' >"$temp_dir/raw_chr${chr}.txt"

    # Print first few lines of raw file for verification
    # variant_id    gene_id    tss_distance    maf    pval_nominal    slope    slope_se
    # head -n 3 "$temp_dir/raw_chr${chr}.txt"

    # Sort raw file
    sort -k1,1 "$temp_dir/raw_chr${chr}.txt" >"$temp_dir/raw_chr${chr}_sorted.txt"

    # Join files and format output
    join -t $'\t' "$temp_dir/raw_chr${chr}_sorted.txt" "$temp_dir/lookup_sorted.txt" |
        awk -F'\t' -v OFS=',' '
    {
        gsub(/^chr/, "", $8)  # Remove chr prefix from chromosome number
        print $2,         # GENE
              $12,         # RSID
              $8,         # CHR
              $9,        # POS
              int($3),    # TSS_DISTANCE
              $10,        # A1
              $11,        # A2
              $4,         # MAF
              $5,         # PVAL
              $6,         # BETA
              $7         # SE
    }' >"$temp_dir/chr${chr}_formatted.txt"

    # head -n 3 "$temp_dir/chr${chr}_formatted.txt"
    # GENE  RSID    CHR  POS    TSS_DISTANCE    A1  A2  MAF PVAL    BETA    SE
    # ENSG00000227232    rs1209314672    1    13550    -16003    C    T    0.0141791    0.734151    0.0587242    0.172837

    # Sort and append to output file
    sort -t',' -k1,1 -k4,4n "$temp_dir/chr${chr}_formatted.txt" >>"$output_file"
) &
done

wait

# Clean up
rm -rf "$temp_dir"
echo "All chromosomes processed"