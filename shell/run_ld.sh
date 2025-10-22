#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J ld
#SBATCH -o /dev/null
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH -t 40:00:00
#SBATCH --array=0-9

source ./shell/setting.sh
source activate $python_env

sleep $(((RANDOM % 10) + 1))
i=$SLURM_ARRAY_TASK_ID
QTDid=${QTDids[$i]}
Celltype=${Celltypes[$i]}
annotation_path="$save_path_main/${QTDid}/LDSC/LD_annotation" # chr.print_snps.txt
ld_save_path="$save_path_main/${QTDid}/LDSC"

## create annotation LD files
python $src_path/create_ld_annot.py -s ${QTDid} -sd $save_path_main -tpld $tar_LDSC_path -apld $aux_LDSC_path >>${log_path}/create_ld_annot.log 2>&1

mkdir -p "${ld_save_path}/TAR" "${ld_save_path}/AUX"
## target population LD files
for chr in ${chrs[@]}; do
    target_file=$(find "$tar_LDSC_path" -name "1000G.*.QC.maf.${chr}.bim" -print -quit)
    base_filename=$(basename "$target_file" .bim)
    ${plink_path} \
        --bfile "${tar_LDSC_path}/${base_filename}" \
        --extract "${annotation_path}/${chr}.print_snps.txt" \
        --keep-allele-order \
        --make-bed \
        --out "${ld_save_path}/TAR/1000G.TAR.QC.maf.${chr}"
done

## auxiliary population LD files
for chr in ${chrs[@]}; do
    aux_file=$(find "$aux_LDSC_path" -name "1000G.*.QC.maf.${chr}.bim" -print -quit)
    base_filename=$(basename "$aux_file" .bim)
    ${plink_path} \
        --bfile "${aux_LDSC_path}/${base_filename}" \
        --extract "${annotation_path}/${chr}.print_snps.txt" \
        --keep-allele-order \
        --make-bed \
        --out "${ld_save_path}/AUX/1000G.AUX.QC.maf.${chr}"
done

### run s-ldxr
ldsc_save_path="${save_path_main}/${QTDid}/LDSC/LDSC_gene"
mkdir -p ${ldsc_save_path}
for chr in ${chrs[@]}; do
    python ${sldxr_path}/s-ldxr.py \
        --bfile "${save_path_main}/${QTDid}/LDSC/TAR/1000G.TAR.QC.maf.${chr}" "${save_path_main}/${QTDid}/LDSC/AUX/1000G.AUX.QC.maf.${chr}" \
        --print-snps ${annotation_path}/${chr}.print_snps.txt \
        --annot ${annotation_path}/${chr}.annot.gz \
        --ld-wind-cm 1.0 \
        --score standardized \
        --out ${ldsc_save_path}/TAR_AUX_std_chr${chr} \
        >>${log_path}/ldscore.log 2>&1
done
