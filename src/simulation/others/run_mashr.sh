#!/usr/bin/env bash
set -euo pipefail

# Generate mashr-ready simulations for selected nt_n2_propt settings,
# run mashr with two input condition sets, and draw true-sign diagnostics.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

BASE_PATH="${BASE_PATH:-bench/result_mashr}"
RUNNAME="${RUNNAME:-nt_n2_propt_mashr}"
MASHR_COV_METHOD="${MASHR_COV_METHOD:-canonical_pca}"
FORCE="${FORCE:-0}"

run_simulation() {
    local args=(
        src/simulation/others/simulation_mashr.py
        --out_dir "${BASE_PATH}"
        --runname "${RUNNAME}"
        --n2 400
        --nt 5000
        --propt 0.4 0.8
        --nrep 100
    )
    if [[ "${FORCE}" == "1" ]]; then
        args+=(--force)
    fi
    python3 "${args[@]}"
}

run_mashr_condition_set() {
    local condition_set="$1"
    local output_prefix="$2"
    local args=(
        src/simulation/others/run_mashr_benchmark.R
        --base_path "${BASE_PATH}"
        --runname "${RUNNAME}"
        --cov_method "${MASHR_COV_METHOD}"
        --condition_set "${condition_set}"
        --output_prefix "${output_prefix}"
        --lfsr_threshold 0.05
        --pvalue_threshold 0.05
    )
    if [[ "${FORCE}" == "1" ]]; then
        args+=(--force)
    fi
    Rscript "${args[@]}"
}

run_simulation
run_mashr_condition_set sc mashr_sc
run_mashr_condition_set sc_bulk mashr_sc_bulk
