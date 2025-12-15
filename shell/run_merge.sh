#!/bin/bash

source ./shell/setting.sh
conda activate $python_env

for i in ${!QTDids[@]}; do
(
    sleep $(((RANDOM % 10) + 1))
    QTDid=${QTDids[$i]}
    Celltype=${Celltypes[$i]}
    for chr in ${chrs[@]}; do
        echo "Processing QTDid: $QTDid, Celltype: $Celltype, Chr: $chr" >>${log_path}/merge.log
        python $src_path/preprocess/merge.py -s $QTDid -t $Celltype -c $chr \
            -tpd $target_population_data_path \
            -tpn ${target_population_samplesize[$i]} \
            -tid $tissue_data_path -tin ${tissue_samplesize[$i]} \
            -apd $eQTLCatalog_path \
            -tpld $tar_LDSC_path -apld $aux_LDSC_path -sd $save_path_main >>${log_path}/merge.log 2>&1
    done

    # copy Celltype proportion
    cp ${celltype_proportion_path} ${save_path_main}/${QTDid}
) &
done

wait
echo "All merging tasks for ${target_population} using ${use_tissue} data completed!" >>${log_path}/merge.log
