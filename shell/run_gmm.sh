
# use screen to run this script, parallel run 3 jobs
source ./shell/setting.sh

for i in $(seq 0 9); do
    echo "Running GMM for ${QTDids[$i]} - ${Celltypes[$i]} of $target_population using $use_tissue data..." >>${log_path}/run_gmm.log 2>&1
    python $src_path/traceCB/run_gmm.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -c ${chrs[@]} -d $save_path_main >>${log_path}/run_gmm.log 2>&1
    # visualization
    python $src_path/traceCB/visual.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -d $save_path_main -tin ${tissue_samplesize[$i]} -tpn ${target_population_samplesize[$i]} -apn ${aux_population_samplesize[$i]} >>${log_path}/visual.log 2>&1
done
