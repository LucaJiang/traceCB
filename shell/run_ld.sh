#!/bin/bash
source ./shell/setting.sh
conda activate $python_env

# 最大并行任务数
MAX_JOBS=10

for i in ${!QTDids[@]}; do
(
    sleep $(((RANDOM % 10) + 1))

    QTDid=${QTDids[$i]}
    Celltype=${Celltypes[$i]}
    annotation_path="$save_path_main/${QTDid}/LDSC/LD_annotation" # chr.print_snps.txt
    ld_save_path="$save_path_main/${QTDid}/LDSC"

    ## create annotation LD files
    python $src_path/preprocess/create_ld_annot.py -s ${QTDid} -sd $save_path_main -tpld $tar_LDSC_path -apld $aux_LDSC_path >>${log_path}/create_ld_annot.log 2>&1

    mkdir -p "${ld_save_path}/TAR" "${ld_save_path}/AUX"
    
    for chr in ${chrs[@]}; do
        while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
            wait -n
        done

    (
        ## target population LD files
        target_file=$(find "$tar_LDSC_path" -name "1000G.*.QC.maf.${chr}.bim" -print -quit)
        base_filename=$(basename "$target_file" .bim)
        ${plink_path} \
            --bfile "${tar_LDSC_path}/${base_filename}" \
            --extract "${annotation_path}/${chr}.print_snps.txt" \
            --keep-allele-order \
            --make-bed \
            --out "${ld_save_path}/TAR/1000G.TAR.QC.maf.${chr}"

        ## auxiliary population LD files
        aux_file=$(find "$aux_LDSC_path" -name "1000G.*.QC.maf.${chr}.bim" -print -quit)
        base_filename=$(basename "$aux_file" .bim)
        ${plink_path} \
            --bfile "${aux_LDSC_path}/${base_filename}" \
            --extract "${annotation_path}/${chr}.print_snps.txt" \
            --keep-allele-order \
            --make-bed \
            --out "${ld_save_path}/AUX/1000G.AUX.QC.maf.${chr}"

        ### run s-ldxr
        ldsc_save_path="${save_path_main}/${QTDid}/LDSC/LDSC_gene"
        mkdir -p ${ldsc_save_path}
        python ${sldxr_path}/s-ldxr.py \
            --bfile "${save_path_main}/${QTDid}/LDSC/TAR/1000G.TAR.QC.maf.${chr}" "${save_path_main}/${QTDid}/LDSC/AUX/1000G.AUX.QC.maf.${chr}" \
            --print-snps ${annotation_path}/${chr}.print_snps.txt \
            --annot ${annotation_path}/${chr}.annot.gz \
            --ld-wind-cm 1.0 \
            --score standardized \
            --out ${ldsc_save_path}/TAR_AUX_std_chr${chr} \
            >>${log_path}/ldscore.log 2>&1
    ) &
    done
    wait

) &
done

wait

echo "LD score calculation for ${target_population} using ${use_tissue} data completed!" >>${log_path}/ldscore.log
