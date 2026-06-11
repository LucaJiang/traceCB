#!/usr/bin/env bash
set -euo pipefail

# Curated traceCB^2 simulation grid.
#
# Usage:
#   bash src/simulation/others/run_tracecb2.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/others/run_tracecb2.sh
#
# OMEGA_MODE can be both, estimate, or true. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../run_common.sh
source "${SCRIPT_DIR}/../run_common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result}"
TRACECB2_SEED="${TRACECB2_SEED:-20260525}"

run_tracecb2_grid() {
    local omega_kind="$1"
    local runname="$2"
    shift 2

    local omega_args=()
    if [[ "${omega_kind}" == "estimate" ]]; then
        omega_args+=(--estimate_omega)
    fi

    run_cmd "${PYTHON}" src/simulation/others/simulation_tracecb2.py \
        --pop1_geno "${POP1_GENO}" \
        --pop2_geno "${POP2_GENO}" \
        --runname "${runname}" \
        "$@" \
        --nsnp "${NSNP}" \
        --out_dir "${OUT_DIR}" \
        --nrep "${NREP}" \
        --seed "${TRACECB2_SEED}" \
        "${omega_args[@]}"
}

plot_tracecb2_grid() {
    local omega_kind="$1"
    local runname="$2"
    local metric="$3"
    local ymin="$4"
    local ymax="$5"

    if [[ "${RUN_VISUALS}" != "1" ]]; then
        return
    fi

    local visual_args=(
        src/simulation/others/visual_tracecb2.py
        --metric "${metric}"
        --runname "${runname}"
        --omega "$(omega_filter_arg "${omega_kind}")"
        --base_path "${OUT_DIR}"
    )
    if [[ -n "${ymin}" ]]; then
        visual_args+=(--ymin "${ymin}")
    fi
    if [[ -n "${ymax}" ]]; then
        visual_args+=(--ymax "${ymax}")
    fi
    run_cmd "${PYTHON}" "${visual_args[@]}"
}

run_all_tracecb2_for_omega() {
    local omega_kind="$1"
    local omega_name
    omega_name="$(omega_label "${omega_kind}")"
    echo "Running traceCB^2 grids with ${omega_name}"

    run_tracecb2_grid "${omega_kind}" nt1_nt2_propt \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 100 \
        --n2 400 \
        --nt1 500 1000 2000 \
        --nt2 1000 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.005

    run_tracecb2_grid "${omega_kind}" n1_nt1_propt \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 100 200 \
        --n2 400 \
        --nt1 500 1000 2000 \
        --nt2 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.005

    run_tracecb2_grid "${omega_kind}" n2_nt1_propt \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 100 \
        --n2 200 400 \
        --nt1 500 1000 2000 \
        --nt2 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.005

    run_tracecb2_grid "${omega_kind}" alpha_h2sq_pcausal_propt \
        --h1sq 0.000000000001 \
        --h2sq 0.1 0.2 \
        --gc 0 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.005 0.01 0.02
}

plot_all_tracecb2_for_omega() {
    local omega_kind="$1"
    local nt_ymax="0.39"
    if [[ "${omega_kind}" == "true" ]]; then
        nt_ymax=""
    fi
    plot_tracecb2_grid "${omega_kind}" nt1_nt2_propt power "" "${nt_ymax}"
    plot_tracecb2_grid "${omega_kind}" n1_nt1_propt power "" ""
    plot_tracecb2_grid "${omega_kind}" n2_nt1_propt power "" ""
    plot_tracecb2_grid "${omega_kind}" alpha_h2sq_pcausal_propt alpha 0 0.45
}

while IFS= read -r omega_kind; do
    run_all_tracecb2_for_omega "${omega_kind}"
done < <(omega_modes)

while IFS= read -r omega_kind; do
    plot_all_tracecb2_for_omega "${omega_kind}"
done < <(omega_modes)
