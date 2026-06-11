#!/usr/bin/env bash
set -euo pipefail

# Rebuttal simulation grids for showing when traceC and traceCB gain power.
#
# Usage:
#   bash src/simulation/others/run_power_gain.sh
#
# Useful server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/others/run_power_gain.sh
#
# OMEGA_MODE can be both, estimate, or true. RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../run_common.sh
source "${SCRIPT_DIR}/../run_common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result_power_gain}"
POWER_GAIN_SEED="${POWER_GAIN_SEED:-20260605}"

TRACEC_RUNNAME="${TRACEC_RUNNAME:-power_gain_tracec_n2_rho}"
TRACECB_RUNNAME="${TRACECB_RUNNAME:-power_gain_tracecb_propt_nt}"
GAIN_OUTPUT_PREFIX="${GAIN_OUTPUT_PREFIX:-${OUT_DIR}/img/power_gain}"

mkdir -p "${OUT_DIR}/img"

run_power_gain_grid() {
    local omega_kind="$1"
    local omega_args=()
    if [[ "${omega_kind}" == "estimate" ]]; then
        omega_args+=(--estimate_omega)
    fi

    echo "Running traceC gain grid with $(omega_label "${omega_kind}"): x=N2, group=rho"
    run_cmd "${PYTHON}" src/simulation/simulation.py \
        --pop1_geno "${POP1_GENO}" \
        --pop2_geno "${POP2_GENO}" \
        --runname "${TRACEC_RUNNAME}" \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0 0.3 0.6 0.9 \
        --n1 100 \
        --n2 100 200 400 800 \
        --nt 5000 \
        --nsnp "${NSNP}" \
        --propt 0.4 \
        --pcausal 0.005 \
        --out_dir "${OUT_DIR}" \
        --nrep "${NREP}" \
        --seed "${POWER_GAIN_SEED}" \
        "${omega_args[@]}"

    if [[ "${RUN_VISUALS}" == "1" ]]; then
        run_cmd "${PYTHON}" src/simulation/visual_simulation.py \
            --metric power \
            --runname "${TRACEC_RUNNAME}" \
            --row h2sq \
            --col gc \
            --x n2 \
            --omega "$(omega_filter_arg "${omega_kind}")" \
            --base_path "${OUT_DIR}" \
            --ymin 0 \
            --ymax 1
    fi

    echo "Running traceCB gain grid with $(omega_label "${omega_kind}"): x=propt, group=Nt"
    run_cmd "${PYTHON}" src/simulation/simulation.py \
        --pop1_geno "${POP1_GENO}" \
        --pop2_geno "${POP2_GENO}" \
        --runname "${TRACECB_RUNNAME}" \
        --h1sq 0.1 \
        --h2sq 0.1 \
        --gc 0.7 \
        --n1 100 \
        --n2 400 \
        --nt 1000 5000 10000 \
        --nsnp "${NSNP}" \
        --propt 0.01 0.3 0.6 0.9 \
        --pcausal 0.005 \
        --out_dir "${OUT_DIR}" \
        --nrep "${NREP}" \
        --seed "${POWER_GAIN_SEED}" \
        "${omega_args[@]}"

    if [[ "${RUN_VISUALS}" == "1" ]]; then
        run_cmd "${PYTHON}" src/simulation/visual_simulation.py \
            --metric power \
            --runname "${TRACECB_RUNNAME}" \
            --row h2sq \
            --col nt \
            --x propt \
            --omega "$(omega_filter_arg "${omega_kind}")" \
            --base_path "${OUT_DIR}" \
            --ymin 0 \
            --ymax 1

        local output_prefix="${GAIN_OUTPUT_PREFIX}"
        if [[ "${OMEGA_MODE}" == "both" ]]; then
            output_prefix="${GAIN_OUTPUT_PREFIX}_$(omega_label "${omega_kind}")"
        fi
        run_cmd "${PYTHON}" src/simulation/others/visual_power_gain.py \
            --base_path "${OUT_DIR}" \
            --tracec_runname "${TRACEC_RUNNAME}" \
            --tracecb_runname "${TRACECB_RUNNAME}" \
            --omega "$(omega_filter_arg "${omega_kind}")" \
            --output_prefix "${output_prefix}"
        echo "Done. Gain figure: ${output_prefix}.pdf"
    fi
}

while IFS= read -r omega_kind; do
    run_power_gain_grid "${omega_kind}"
done < <(omega_modes)
