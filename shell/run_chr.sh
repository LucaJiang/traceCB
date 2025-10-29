# use screen to run this script, parallel run 2 jobs
source ./shell/setting.sh

# 最大并行任务数
MAX_JOBS=12
i=6
# 函数：运行单个任务
run_task() {
    local j=$1
    echo "Running GMM for chr${j} ${QTDids[$i]} - ${Celltypes[$i]} of $target_population using $use_tissue data..." >>${log_path}/run_gmm_test_sign_sigma.log 2>&1
    python $src_path/traceCB/run_gmm.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -c ${j} -d $save_path_main
    # visualization
}

# 并行运行
for j in $(seq 1 22); do
    # 等待直到运行的任务数小于 MAX_JOBS
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 10
    done
    
    # 在后台运行任务
    run_task $j &
done

# 等待所有后台任务完成
wait

# python $src_path/traceCB/visual.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -d $save_path_main -tin ${tissue_samplesize[$i]} -tpn ${target_population_samplesize[$i]} -apn ${aux_population_samplesize[$i]}
echo "All tasks completed!"