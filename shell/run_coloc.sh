#!/bin/bash
# filepath: /home/wjiang49/traceCB/shell/run_coloc.sh

# Load settings and activate conda environment
source ./shell/setting.sh
conda activate "$r_env"

# Using save_path_main to keep results with other traceCB outputs
coloc_result_save_path="${save_path_main}/coloc/"
if [ ! -d "$coloc_result_save_path" ]; then
    mkdir -p "$coloc_result_save_path"
fi

# Generate eQTL paths
# QTDids and Celltypes are loaded from setting.sh
eqtl_paths=()
eqtl_save_names=()
for i in "${!QTDids[@]}"; do
    # Assuming GMM results are stored in study-specific folders under save_path_main
    eqtl_paths+=("${save_path_main}/${QTDids[$i]}/GMM")
    # Construct save name
    eqtl_save_names+=("${use_tissue}_${QTDids[$i]}_${Celltypes[$i]}")
done

# GWAS and loci paths
gwas_paths=()
loci_paths=()

# Helper function to add paths
add_trait_paths() {
    local prefix=$1
    local trait=$2
    gwas_paths+=("${coloc_data_path}/${prefix}/${prefix}_${trait}_GWAS.csv")
    loci_paths+=("${coloc_data_path}/${prefix}/${prefix}_${trait}_loci.csv")
}

# EAS traits (bcx/bbj)
for trait in "mon" "lym" "neu" "ra" "asthma" "atopy"; do
    if [[ $trait =~ ^(mon|lym|neu)$ ]]; then
        prefix="bcx"
    else
        prefix="bbj"
    fi
    add_trait_paths "$prefix" "$trait"
done

# BCX traits
for trait in "wbc" "mchc" "rbc" "mcv" "mpv" "plt" "hct" "hgb" "mch" "eos" "bas"; do
    prefix="bcx"
    add_trait_paths "$prefix" "$trait"
done


# --- Parallel Execution Settings ---
MAX_JOBS=64

# --- Task Definition ---
run_coloc_task() {
    local i=$1
    local j=$2
    
    local eqtl_path="${eqtl_paths[$j]}"
    local eqtl_save_name="${eqtl_save_names[$j]}"
    local gwas_path="${gwas_paths[$i]}"
    local loci_path="${loci_paths[$i]}"

    local gwas_base=$(basename "$gwas_path" _GWAS.csv)
    local output_file="${coloc_result_save_path}${gwas_base}_${eqtl_save_name}_coloc.csv"
    local log_file="${coloc_result_save_path}${gwas_base}_${eqtl_save_name}.log"

    echo "Running coloc for ${gwas_base} and ${eqtl_save_name}..."
    Rscript "${src_path}/coloc/run_coloc.r" "$loci_path" "$gwas_path" "$eqtl_path" "$output_file" >>"$log_file" 2>&1
    echo "Finished coloc for ${gwas_base} and ${eqtl_save_name}"
}

# --- Main Execution Logic ---
echo "Launching coloc tasks in parallel (max ${MAX_JOBS} jobs)..."
# Run coloc for each pair of GWAS and eQTL
for i in "${!gwas_paths[@]}"; do
    # Inner loop over eQTLs
    for j in "${!eqtl_paths[@]}"; do
        # Wait until the number of running jobs is less than the maximum limit
        while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
            sleep 5
        done
        
        run_coloc_task "$i" "$j" &
    done
done

wait
echo "All coloc tasks completed."

# Summary results
conda activate "$python_env"
python "${src_path}/coloc/summary_results.py"
