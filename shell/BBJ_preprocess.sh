#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J BBJ_preprocess
#SBATCH -o log/BBJ_preprocess.log
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH -t 10:00:00
#SBATCH --array=0-4

BBJ_path_main="/gpfs1/scratch/ResearchGroups/bios_mingxcai/data/BBJ_eQTL"
OUT_path="${BBJ_path_main}/by_celltype_chr"
cell_type_list=("B_cells" "CD4+T_cells" "CD8+T_cells" "Monocytes" "NK_cells")
# target path: ~/ResearchGroup/data/BBJ_eQTL/cell_type/chr1.csv

if [ ! -d "$OUT_path" ]; then
    mkdir -p "$OUT_path"
fi
# structure of eQTL files:
# tar -tzf eQTL_B_cells.tar.gz | head
# B_cells/
# B_cells/chr1_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz
# tar -xOzf eQTL_B_cells.tar.gz B_cells/chr1_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz | zcat | head
# SNP	POS	REF	ALT	gene	beta	t-stat	p-value
# rs11090253	100017147	TC	T	ENSG00000036054.8	-0.368564601516914	-1.22853496743858	0.222073259733985
# chr3:100017147:D	100017147	TC	T	ENSG00000081148.11	0.101998059353839	0.337710406631037	0.736275455591437
# chr3:100017147:D	100017147	TC	T	ENSG00000114021.7	0.269773020022972	0.896160599960745	0.372277126551398

cell_type=${cell_type_list[${SLURM_ARRAY_TASK_ID}]}
if [ ! -d "$OUT_path/${cell_type}" ]; then
    mkdir -p "$OUT_path/${cell_type}"
fi
cell_type_path="$BBJ_path_main/eQTL_${cell_type}.tar.gz"
header="CHR,RSID,POS,A2,A1,GENE,BETA,Z,PVAL"

for chr in {1..22}; do
    echo "Processing chromosome $chr for ${cell_type}"
    output_file="$OUT_path/${cell_type}/chr${chr}.csv"
    # Write header to output file
    echo "$header" >"$output_file"
    # Process the chromosome data and append to output file
    # Input columns: SNP POS REF ALT gene beta t-stat p-value
    # Output columns: CHR,RSID,POS,A2,A1,GENE,BETA,Z,PVAL
    tar -xOzf "$cell_type_path" "${cell_type}/chr${chr}_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz" |
        zcat |
        awk -v OFS=',' -v chr=$chr 'NR>1 {
        sub(/\.[0-9]+$/, "", $5);  # Remove version number from gene ID
        if (substr($1, 1, 2) == "rs") {
            print chr,    # CHR
                $1,       # RSID
                $2,       # POS
                $3,       # A2 (REF)
                $4,       # A1 (ALT)
                $5,       # Gene
                $6,       # BETA
                $7,       # Z (t-stat)
                $8        # PVAL (p-value)
        }
    }' >>"$output_file"

    echo "Finished processing chromosome $chr for ${cell_type}"
done

echo "All chromosomes processed for ${cell_type}"
