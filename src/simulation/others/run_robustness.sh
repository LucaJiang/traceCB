#!/usr/bin/env bash
set -euo pipefail

# Supplementary robustness grids using simulation_robustness.py.
#
# Usage:
#   bash src/simulation/others/run_robustness.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 bash src/simulation/others/run_robustness.sh
#
# This curated entry point pins omega per robustness suite. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../run_common.sh
source "${SCRIPT_DIR}/../run_common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result}"
GMM_PROPT_SEED="${GMM_PROPT_SEED:-20260525}"

run_robustness_grid() {
    local omega_kind="$1"
    local runname="$2"
    shift 2

    local sim_args=(
        "${PYTHON}" src/simulation/others/simulation_robustness.py
        --pop1_geno "${POP1_GENO}"
        --pop2_geno "${POP2_GENO}"
        --runname "${runname}"
        "$@"
        --nsnp "${NSNP}"
        --out_dir "${OUT_DIR}"
        --nrep "${NREP}"
        --seed "${GMM_PROPT_SEED}"
    )
    if [[ "${omega_kind}" == "estimate" ]]; then
        sim_args+=(--estimate_omega)
    fi
    run_cmd "${sim_args[@]}"
}

plot_robustness_grid() {
    local omega_kind="$1"
    local metric="$2"
    local runname="$3"
    local ymin="$4"
    local ymax="$5"
    shift 5

    if [[ "${RUN_VISUALS}" != "1" ]]; then
        return
    fi

    local visual_args=(
        src/simulation/visual_simulation.py
        --metric "${metric}"
        --runname "${runname}"
        --omega "$(omega_filter_arg "${omega_kind}")"
        --base_path "${OUT_DIR}"
        --img_dir "${IMG_DIR}"
        "$@"
    )
    if [[ -n "${ymin}" ]]; then
        visual_args+=(--ymin "${ymin}")
    fi
    if [[ -n "${ymax}" ]]; then
        visual_args+=(--ymax "${ymax}")
    fi
    run_cmd "${PYTHON}" "${visual_args[@]}"
}

run_gmm_propt_robustness() { # whether robust to error in cell type proportion estimation
    local omega_kind="$1"
    local runname="alpha_pcausal_propt_gmmproptmode"

    run_robustness_grid "${omega_kind}" "${runname}" \
        --h1sq 0.000000000001 \
        --h2sq 0.1 \
        --gc 0 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --gmm_propt_mode exact \
        --pcausal 0.005

    run_robustness_grid "${omega_kind}" "${runname}" \
        --h1sq 0.000000000001 \
        --h2sq 0.1 \
        --gc 0 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --gmm_propt_mode underestimate overestimate \
        --gmm_propt_mode_scale 0.1 0.2 \
        --pcausal 0.005

    run_robustness_grid "${omega_kind}" "${runname}" \
        --h1sq 0.000000000001 \
        --h2sq 0.1 \
        --gc 0 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --gmm_propt_mode normal \
        --gmm_propt_normal_var 0.1 \
        --pcausal 0.005

    plot_robustness_grid "${omega_kind}" alpha "${runname}" "" 0.4
}

run_segmented_null_robustness() {
    local omega_kind="$1"
    # h1sq is for shared region |B|
    run_robustness_grid "${omega_kind}" alpha_a_h2sq_nullregionprop_propt_nocormax \
        --h1sq 0.1 \
        --h2sq 0.1 0.2 \
        --gc 0.3 0.5 0.9 \
        --causal_partition_mode pop2_a_shared_b \
        --null_region_prop 0.2 0.5 \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --propt 0.01 0.3 0.6 \
        --pcausal 0.005
    plot_robustness_grid "${omega_kind}" alpha_a alpha_a_h2sq_nullregionprop_propt_nocormax 0.0 0.4 \
        --row h2sq_nullregionprop \
        --col gc \
        --x propt \
        --plot_filter gc=0.3,0.5,0.9
}

echo "Running GMM proportion robustness with $(omega_label estimate)"
run_gmm_propt_robustness estimate

echo "Running segmented-null robustness with $(omega_label true)"
run_segmented_null_robustness true
