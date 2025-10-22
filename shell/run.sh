#!/bin/bash
#SBATCH -A pa_bios_department
#SBATCH -p special_bios
#SBATCH -J gmm
#SBATCH -o /dev/null
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH -t 10:00:00
#SBATCH --array=0-9

source ./shell/setting.sh
source activate $python_env

sleep $(((RANDOM % 10) + 1))
i=$SLURM_ARRAY_TASK_ID
# visualization
python -m pdb $src_path/visual.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -d $save_path_main -tin ${tissue_samplesize[$i]} -tpn ${target_population_samplesize[$i]} -apn ${aux_population_samplesize[$i]} >>${log_path}/visual.log 2>&1
