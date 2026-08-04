#!/usr/bin/env bash
set -euo pipefail

# Run masked-omega comparisons for pop1 target.
#
# Usage:
#   bash src/simulation/experiments/run_masked_omega.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 bash src/simulation/experiments/run_masked_omega.sh
#
# This curated entry point uses estimated omega. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result/masked_omega_compare}"
MASKED_OMEGA_SEED="${MASKED_OMEGA_SEED:-20260525}"

ALPHA_GC=0
ALPHA_RUNNAME="alpha_masked_omega_trueomega_gc0_h2sq_pcausal_propt"

masked_runname_for_omega() {
    local omega_kind="$1"
    local base_runname="$2"
    if [[ "${omega_kind}" == "true" ]]; then
        printf "%s\n" "${base_runname}"
    elif [[ "${base_runname}" == *trueomega* ]]; then
        printf "%s\n" "${base_runname/trueomega/estomega}"
    else
        printf "%s_estomega\n" "${base_runname}"
    fi
}

run_masked_for_omega() {
    local omega_kind="$1"
    local runname
    runname="$(masked_runname_for_omega "${omega_kind}" "${ALPHA_RUNNAME}")"
    local sim_args=(
        "${PYTHON}" src/simulation/experiments/simulate_masked_omega.py
        --pop1_geno "${POP1_GENO}"
        --pop2_geno "${POP2_GENO}"
        --runname "${runname}"
        --h1sq 0
        --h2sq 0.1 0.2
        --gc "${ALPHA_GC}"
        --n1 100
        --n2 400
        --nt 5000
        --nsnp "${NSNP}"
        --propt 0.01 0.2 0.4 0.6 0.8
        --pcausal 0.005 0.01 0.02
        --out_dir "${OUT_DIR}"
        --nrep "${NREP}"
        --seed "${MASKED_OMEGA_SEED}"
    )
    if [[ "${omega_kind}" == "estimate" ]]; then
        sim_args+=(--estimate_omega)
    fi
    run_cmd "${sim_args[@]}"
}

plot_masked_for_omega() {
    local omega_kind="$1"
    if [[ "${RUN_VISUALS}" != "1" ]]; then
        return
    fi

    local visual_args=(
        "${PYTHON}" src/simulation/experiments/plot_masked_omega.py
        --base_path "${OUT_DIR}"
        --img_dir "${IMG_DIR}"
        --runname
        "$(masked_runname_for_omega "${omega_kind}" "${ALPHA_RUNNAME}")"
    )

    visual_args+=(
        --save_prefix "alpha_masked_omega_$(omega_label "${omega_kind}")_gc_h2sq_pcausal_propt"
        --metric alpha
        --x propt
        --row h2sq
        --col pcausal
        --alpha_ymax 0.48
    )
    run_cmd "${visual_args[@]}"
}


run_masked_for_omega estimate
plot_masked_for_omega estimate
