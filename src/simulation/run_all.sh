#!/usr/bin/env bash
set -euo pipefail

# Simulation shell entrypoint.
#
# Usage:
#   bash src/simulation/run_all.sh main          # main paper grids
#   bash src/simulation/run_all.sh robustness    # supplementary robustness grids
#   bash src/simulation/run_all.sh tracecb2      # traceCB^2 grids
#   bash src/simulation/run_all.sh masked-omega  # masked-omega grids
#   bash src/simulation/run_all.sh mashr         # mashr benchmark grids
#   bash src/simulation/run_all.sh power-gain    # power-gain figure grids
#   bash src/simulation/run_all.sh chr22         # chromosome 22 experiment
#   bash src/simulation/run_all.sh all           # all small-window grids
#
# Common server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/run_all.sh all
#
# OMEGA_MODE can be both, estimate, or true for scripts with omega modes.
# RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
setup_simulation_env

TARGET="${1:-main}"

run_target() {
    local target="$1"
    case "${target}" in
        main)
            run_cmd bash src/simulation/run_main.sh
            ;;
        robustness)
            run_cmd bash src/simulation/experiments/run_robustness.sh
            ;;
        tracecb2)
            run_cmd bash src/simulation/experiments/run_tracecb2.sh
            ;;
        masked-omega)
            run_cmd bash src/simulation/experiments/run_masked_omega.sh
            ;;
        mashr)
            run_cmd bash src/simulation/experiments/run_mashr.sh
            ;;
        power-gain)
            run_cmd bash src/simulation/experiments/run_power_gain.sh
            ;;
        chr22)
            run_cmd bash src/simulation/chr22/run.sh
            ;;
        *)
            echo "Unknown target: ${target}" >&2
            echo "Expected: main, robustness, tracecb2, masked-omega, mashr, power-gain, chr22, or all" >&2
            return 1
            ;;
    esac
}

if [[ "${TARGET}" == "all" ]]; then
    run_target main
    run_target robustness
    run_target tracecb2
    run_target masked-omega
    run_target mashr
    run_target power-gain
else
    run_target "${TARGET}"
fi
