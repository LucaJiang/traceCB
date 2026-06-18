#!/usr/bin/env bash

# Shared environment setup for simulation shell entrypoints.
# Source this file from a run script, then call setup_simulation_env.

setup_simulation_env() {
    local common_dir
    common_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    REPO_ROOT="${REPO_ROOT:-$(cd "${common_dir}/../.." && pwd)}"
    cd "${REPO_ROOT}"

    CONDA_ENV="${CONDA_ENV:-py312}"
    if [[ "${SKIP_CONDA:-0}" != "1" ]]; then
        # shellcheck source=/dev/null
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "${CONDA_ENV}"
    fi

    export PYTHONUNBUFFERED=1
    PYTHON="${PYTHON:-python3}"

    SIM_DATA_DIR="${SIM_DATA_DIR:-data/simulation}"
    POP1_GENO="${POP1_GENO:-${SIM_DATA_DIR}/EAS_n5000_chr22_loci29.npy}"
    POP2_GENO="${POP2_GENO:-${SIM_DATA_DIR}/EUR_n20000_chr22_loci29.npy}"

    NREP="${NREP:-100}"
    NSNP="${NSNP:-2000}"
    RUN_VISUALS="${RUN_VISUALS:-1}"
    OMEGA_MODE="${OMEGA_MODE:-both}"
    IMG_DIR="${IMG_DIR:-${REPO_ROOT}/bench/result/img}"

    export REPO_ROOT CONDA_ENV PYTHON
    export SIM_DATA_DIR POP1_GENO POP2_GENO
    export NREP NSNP RUN_VISUALS OMEGA_MODE IMG_DIR
}

omega_modes() {
    case "${OMEGA_MODE}" in
        both)
            printf "estimate\ntrue\n"
            ;;
        estimate|estimated|est)
            printf "estimate\n"
            ;;
        true|true_omega|trueomega)
            printf "true\n"
            ;;
        *)
            echo "Invalid OMEGA_MODE=${OMEGA_MODE}; expected both, estimate, or true." >&2
            return 1
            ;;
    esac
}

omega_flag_args() {
    local omega_kind="$1"
    if [[ "${omega_kind}" == "estimate" ]]; then
        printf "%s\n" "--estimate_omega"
    fi
}

omega_filter_arg() {
    local omega_kind="$1"
    if [[ "${omega_kind}" == "true" ]]; then
        printf "true\n"
    else
        printf "false\n"
    fi
}

omega_label() {
    local omega_kind="$1"
    if [[ "${omega_kind}" == "true" ]]; then
        printf "trueomega\n"
    else
        printf "estomega\n"
    fi
}

run_cmd() {
    printf "+ "
    printf "%q " "$@"
    printf "\n"
    "$@"
}
