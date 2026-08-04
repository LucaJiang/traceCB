#!/usr/bin/env bash
set -euo pipefail

# Generate mashr-ready simulations for selected nt_n2_propt settings,
# run mashr with two input condition sets, and draw true-sign diagnostics.
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data BASE_PATH=/path/to/mashr_result \
#     NREP=100 NSNP=2000 bash src/simulation/experiments/run_mashr.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
setup_simulation_env

BASE_PATH="${BASE_PATH:-${OUT_DIR:-bench/result_mashr}}"
RUNNAME="${RUNNAME:-nt_n2_propt_mashr}"
MASHR_COV_METHOD="${MASHR_COV_METHOD:-canonical_pca}"
MASHR_SEED="${MASHR_SEED:-20260604}"
FORCE="${FORCE:-0}"

run_simulation() {
    local args=(
        src/simulation/experiments/simulate_mashr.py
        --pop1_geno "${POP1_GENO}"
        --pop2_geno "${POP2_GENO}"
        --out_dir "${BASE_PATH}"
        --runname "${RUNNAME}"
        --n2 400
        --nt 5000
        --nsnp "${NSNP}"
        --propt 0.4 0.8
        --nrep "${NREP}"
        --seed "${MASHR_SEED}"
    )
    if [[ "${FORCE}" == "1" ]]; then
        args+=(--force)
    fi
    run_cmd "${PYTHON}" "${args[@]}"
}

run_mashr_condition_set() {
    local condition_set="$1"
    local output_prefix="$2"
    local args=(
        src/simulation/experiments/benchmark_mashr.R
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
    run_cmd Rscript "${args[@]}"
}

run_simulation
run_mashr_condition_set sc mashr_sc
run_mashr_condition_set sc_bulk mashr_sc_bulk
if [[ "${RUN_VISUALS}" == "1" ]]; then
    run_cmd "${PYTHON}" src/simulation/experiments/plot_mashr_fsr.py \
        --base_path "${BASE_PATH}" \
        --runname "${RUNNAME}" \
        --img_dir "${IMG_DIR}"
fi
