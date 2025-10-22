#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J eQTL_merge
#SBATCH -o /dev/null
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH -t 10:00:00
#SBATCH --array=0-9

source ./shell/setting.sh
source activate $python_env

sleep $(((RANDOM % 10) + 1))
i=$SLURM_ARRAY_TASK_ID
QTDid=${QTDids[$i]}
Celltype=${Celltypes[$i]}
for chr in ${chrs[@]}; do
    echo "Processing QTDid: $QTDid, Celltype: $Celltype, Chr: $chr" >>${log_path}/merge.log
    python $src_path/MergeChr.py -s $QTDid -t $Celltype -c $chr \
        -tpd $target_population_data_path \
        -tpn ${target_population_samplesize[$i]} \
        -tid $tissue_data_path -tin ${tissue_samplesize[$i]} \
        -apd $eQTLCatalog_path \
        -tpld $tar_LDSC_path -apld $aux_LDSC_path -sd $save_path_main >>${log_path}/merge.log 2>&1
done

# copy Celltype proportion
cp ${celltype_proportion_path} ${save_path_main}/${QTDid}
# End of file
