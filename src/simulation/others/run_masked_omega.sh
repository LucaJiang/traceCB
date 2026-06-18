#!/usr/bin/env bash
set -euo pipefail

# Run masked-omega comparisons for pop1 target.
#
# Usage:
#   bash src/simulation/others/run_masked_omega.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/others/run_masked_omega.sh
#
# OMEGA_MODE can be both, estimate, or true. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../run_common.sh
source "${SCRIPT_DIR}/../run_common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result/masked_omega_compare}"
MASKED_OMEGA_SEED="${MASKED_OMEGA_SEED:-20260525}"

ALPHA_GCS="${ALPHA_GCS:-0.3 0.7}"
ALPHA_RUNNAME_GC03="${ALPHA_RUNNAME_GC03:-alpha_masked_omega_trueomega_gc0.3_h2sq_pcausal_propt}"
ALPHA_RUNNAME_GC07="${ALPHA_RUNNAME_GC07:-alpha_masked_omega_trueomega_gc0.7_h2sq_pcausal_propt}"

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

base_alpha_runname_for_gc() {
    local gc="$1"
    if [[ "${gc}" == "0.3" ]]; then
        printf "%s\n" "${ALPHA_RUNNAME_GC03}"
    elif [[ "${gc}" == "0.7" ]]; then
        printf "%s\n" "${ALPHA_RUNNAME_GC07}"
    else
        echo "Unexpected alpha GC value: ${gc}" >&2
        return 1
    fi
}

run_masked_for_omega() {
    local omega_kind="$1"

    for GC in ${ALPHA_GCS}; do
        local base_runname
        local runname
        base_runname="$(base_alpha_runname_for_gc "${GC}")"
        runname="$(masked_runname_for_omega "${omega_kind}" "${base_runname}")"
        local sim_args=(
            "${PYTHON}" src/simulation/others/simulation_masked_omega.py
            --pop1_geno "${POP1_GENO}"
            --pop2_geno "${POP2_GENO}"
            --runname "${runname}"
            --h1sq 0
            --h2sq 0.1 0.2
            --gc "${GC}"
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
    done
}

plot_masked_for_omega() {
    local omega_kind="$1"
    if [[ "${RUN_VISUALS}" != "1" ]]; then
        return
    fi

    local runname_count=0
    local visual_args=(
        "${PYTHON}" src/simulation/others/visual_masked_omega.py
        --base_path "${OUT_DIR}"
        --img_dir "${IMG_DIR}"
        --runname
    )
    for GC in ${ALPHA_GCS}; do
        local base_runname
        base_runname="$(base_alpha_runname_for_gc "${GC}")"
        visual_args+=("$(masked_runname_for_omega "${omega_kind}" "${base_runname}")")
        runname_count=$((runname_count + 1))
    done
    if [[ "${runname_count}" == "0" ]]; then
        echo "ALPHA_GCS must include at least one GC value for plotting." >&2
        return 1
    fi

    visual_args+=(
        --save_prefix "alpha_masked_omega_$(omega_label "${omega_kind}")_gc_h2sq_pcausal_propt"
        --metric alpha
        --x propt
        --row gc
        --col pcausal
        --alpha_ymax 0.48
    )
    run_cmd "${visual_args[@]}"
}


run_masked_for_omega estimate
plot_masked_for_omega estimate
