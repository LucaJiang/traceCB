#!/bin/bash
# filepath: /home/wjiang49/traceCB/shell/run_gmm.sh

# This script processes all chromosomes from all studies in a global parallel pool.

# Load settings and activate conda environment
source ./shell/setting.sh
conda activate "$python_env"

# --- Parallel Execution Settings ---
# Maximum number of parallel jobs across ALL studies
MAX_JOBS=256
# Log file setup
log_dir="$log_path"
time_log_file="${log_dir}/run_gmm_timing.log"
# Ensure the log directory exists
mkdir -p "$log_dir"
# Truncate or create the timing log file
: > "$time_log_file"

# --- Helper Functions ---
# Formats seconds into HH:MM:SS
format_time() {
    local T=$1
    printf '%02d:%02d:%02d' $((T/3600)) $(((T%3600)/60)) $((T%60))
}

# --- Task Definition ---
# This function runs the GMM analysis for a single study and a single chromosome.
# Arguments: $1=study_index, $2=chromosome_number
run_gmm_task() {
    local i=$1
    local j=$2
    local study_id="${QTDids[$i]}"
    local cell_type="${Celltypes[$i]}"
    
    local task_start_time=$(date +%s)
    echo "Starting GMM for study '${study_id}', chr${j}..."
    
    # Run GMM for one chromosome. Redirect output to a per-study log file.
    python "$src_path/traceCB/run_gmm.py" \
        -s "$study_id" \
        -t "$cell_type" \
        -c "$j" \
        -d "$save_path_main" >> "${log_dir}/run_gmm_study.log" 2>&1
    
    local task_end_time=$(date +%s)
    local elapsed=$((task_end_time - task_start_time))
    echo "Finished GMM for study '${study_id}', chr${j}. Elapsed: $(format_time $elapsed)"
}

# --- Main Execution Logic ---
script_start_time=$(date +%s)
echo "Script started. Launching all chromosome jobs from all studies into a global pool (max ${MAX_JOBS} jobs)."

num_studies=${#QTDids[@]}

# --- Part 1: Launch all GMM tasks in parallel ---
# Use nested loops to generate all tasks, but control them with a single job limit.
for i in $(seq 0 $((num_studies - 1))); do
    for j in "${chrs[@]}"; do
        # Wait until the number of running jobs is less than the maximum limit
        while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
            sleep 5
        done
        
        # Launch the GMM task for the current (study, chromosome) pair in the background
        run_gmm_task "$i" "$j" &
    done
done

# Wait for ALL background GMM jobs from ALL studies to complete
echo "All GMM tasks have been launched. Waiting for completion..."
wait
gmm_end_time=$(date +%s)
gmm_elapsed_fmt=$(format_time $((gmm_end_time - script_start_time)))
echo "All GMM tasks completed. Total GMM time: ${gmm_elapsed_fmt}"
echo "$(date '+%F %T') - All GMM tasks completed in ${gmm_elapsed_fmt}" >> "$time_log_file"


# # --- Part 2: Run post-processing scripts sequentially ---
# echo "-----------------------------------------------------------------"
# # Run final summary scripts
# echo "Running final summary scripts..."
# echo "$(date '+%F %T') with tissue=${use_tissue}, target_population=${target_population}" >> "${log_dir}/final_f3egene.log"
# python "$src_path/visual/f3egene.py" >> "${log_dir}/final_f3egene.log" 2>&1
# echo "$(date '+%F %T') with tissue=${use_tissue}, target_population=${target_population}" >> "${log_dir}/final_f3violin.log"
# python "$src_path/visual/f3violin.py" >> "${log_dir}/final_f3violin.log" 2>&1

# --- Final Summary ---
script_end_time=$(date +%s)
total_elapsed=$(format_time $((script_end_time - script_start_time)))
final_msg="All tasks completed! Total script time: ${total_elapsed}"
echo "================================================================="
echo "$final_msg"
echo "$(date '+%F %T') - $final_msg" >> "$time_log_file"