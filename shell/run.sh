# use screen to run this script, parallel run 2 jobs
source ./shell/setting.sh

# 最大并行任务数
MAX_JOBS=2

# 函数：运行单个任务
run_task() {
    local i=$1
    echo "Running GMM for ${QTDids[$i]} - ${Celltypes[$i]} of $target_population using $use_tissue data..." >>${log_path}/run_gmm.log 2>&1
    python $src_path/traceCB/run_gmm.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -c ${chrs[@]} -d $save_path_main >>${log_path}/run_gmm_${QTDids[$i]}.log 2>&1
    # visualization
    python $src_path/traceCB/visual.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -d $save_path_main -tin ${tissue_samplesize[$i]} -tpn ${target_population_samplesize[$i]} -apn ${aux_population_samplesize[$i]} >>${log_path}/visual_${QTDids[$i]}.log 2>&1
}

# 并行运行
for i in $(seq 4 9); do
    # 等待直到运行的任务数小于 MAX_JOBS
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 10
    done
    
    # 在后台运行任务
    run_task $i &
done

# 等待所有后台任务完成
wait
echo "All tasks completed!" >>${log_path}/run_gmm.log 2>&1