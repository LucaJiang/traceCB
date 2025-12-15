# use screen to run this script, parallel run 2 jobs
source ./shell/setting.sh
conda activate $python_env

# 最大并行任务数
MAX_JOBS=10
i=6
log_file="${log_path}/run_gmm_chr.log"
# ensure log dir
mkdir -p "$log_path"
TIME_LOG="${log_path}/run_gmm_time.log"
: > "$TIME_LOG"  # truncate or create

# helper to format seconds as HH:MM:SS
format_time() {
    local T=$1
    printf '%02d:%02d:%02d' $((T/3600)) $(((T%3600)/60)) $((T%60))
}

script_start=$(date +%s)

# 函数：运行单个任务
run_task() {
        local j=$1
        echo "Starting GMM for chr${j} ${QTDids[$i]} - ${Celltypes[$i]} of $target_population using $use_tissue data..."
        task_start=$(date +%s)

        python $src_path/traceCB/run_double_gmm.py -s ${QTDids[$i]} -t ${Celltypes[$i]} -c ${j} -d $save_path_main >>"$log_file" 2>&1

        task_end=$(date +%s)
        task_elapsed=$((task_end - task_start))
        task_elapsed_fmt=$(format_time $task_elapsed)
        msg="Chr${j} finished. Elapsed: ${task_elapsed_fmt} (${task_elapsed}s)"
        echo "$msg"
        echo "$(date '+%F %T') - $msg" >>"$TIME_LOG"
}

# 并行运行
for j in $(seq 1 22); do
        # 等待直到运行的任务数小于 MAX_JOBS
        while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
                sleep 10
        done
        
        # 在后台运行任务
        run_task $j &
done

# 等待所有后台任务完成
wait

python src/visual/f3egene.py >>"$log_file" 2>&1
python src/visual/f3violin.py >>"$log_file" 2>&1

script_end=$(date +%s)
total_elapsed=$((script_end - script_start))
total_elapsed_fmt=$(format_time $total_elapsed)
echo "All tasks completed! Total elapsed: ${total_elapsed_fmt} (${total_elapsed}s)"
echo "$(date '+%F %T') - All tasks completed! Total elapsed: ${total_elapsed_fmt} (${total_elapsed}s)" >>"$TIME_LOG"
