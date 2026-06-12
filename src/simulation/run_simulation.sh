#!/usr/bin/env bash
set -euo pipefail

# Curated small-window simulation.py grid for the traceCB paper.
#
# Usage:
#   bash src/simulation/run_simulation.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/run_simulation.sh
#
# OMEGA_MODE can be both, estimate, or true. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=run_common.sh
source "${SCRIPT_DIR}/run_common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result}"
SIMULATION_SEED="${SIMULATION_SEED:-20260525}"

run_simulation_grid() {
    local omega_kind="$1"
    local runname="$2"
    local metric="$3"
    local ymin="$4"
    local ymax="$5"
    shift 5

    local omega_args=()
    if [[ "${omega_kind}" == "estimate" ]]; then
        omega_args+=(--estimate_omega)
    fi

    run_cmd "${PYTHON}" src/simulation/simulation.py \
        --pop1_geno "${POP1_GENO}" \
        --pop2_geno "${POP2_GENO}" \
        --runname "${runname}" \
        "$@" \
        --nsnp "${NSNP}" \
        --out_dir "${OUT_DIR}" \
        --nrep "${NREP}" \
        --seed "${SIMULATION_SEED}" \
        "${omega_args[@]}"

    if [[ "${RUN_VISUALS}" == "1" ]]; then
        local visual_args=(
            src/simulation/visual_simulation.py
            --metric "${metric}"
            --runname "${runname}"
            --base_path "${OUT_DIR}"
            --omega "$(omega_filter_arg "${omega_kind}")"
        )
        if [[ -n "${ymin}" ]]; then
            visual_args+=(--ymin "${ymin}")
        fi
        if [[ -n "${ymax}" ]]; then
            visual_args+=(--ymax "${ymax}")
        fi
        run_cmd "${PYTHON}" "${visual_args[@]}"
    fi
}

run_all_grids_for_omega() {
    local omega_kind="$1"
    local omega_name
    omega_name="$(omega_label "${omega_kind}")"
    local nt_ymin="0.16"
    local nt_ymax="0.38"
    local h2_ymin="0.12"
    local h2_ymax="0.42"
    local n1_ymin="0.05"
    local n1_ymax="0.32"
    if [[ "${omega_kind}" == "true" ]]; then
        nt_ymin=""
        nt_ymax="0.88"
        h2_ymin=""
        h2_ymax="0.88"
        n1_ymin=""
        n1_ymax="0.88"
    fi
    echo "Running main small-window grids with ${omega_name}"

    run_simulation_grid "${omega_kind}" nt_n2_propt power "${nt_ymin}" "${nt_ymax}" \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 100 \
        --n2 100 200 400 \
        --nt 1000 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.05

    run_simulation_grid "${omega_kind}" h2sq_gc_propt power "${h2_ymin}" "${h2_ymax}" \
        --h1sq 0.1 \
        --h2sq 0.1 0.2 \
        --gc 0.01 0.5 0.9 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.05

    run_simulation_grid "${omega_kind}" n1_pcausal_propt power "${n1_ymin}" "${n1_ymax}" \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 50 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.01 0.02 0.05

    run_simulation_grid "${omega_kind}" alpha_h2sq_pcausal_propt alpha "" 0.48 \
        --h1sq 0.000000000001 \
        --h2sq 0.1 0.2 \
        --gc 0 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.01 0.02 0.05
}

while IFS= read -r omega_kind; do
    run_all_grids_for_omega "${omega_kind}"
done < <(omega_modes)
