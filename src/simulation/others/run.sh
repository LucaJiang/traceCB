#!/usr/bin/env bash
set -euo pipefail

# Simulation shell entrypoint.
#
# Usage:
#   bash src/simulation/others/run.sh simulation      # simulation.py paper grids
#   bash src/simulation/others/run.sh robustness      # robustness supplementary grids
#   bash src/simulation/others/run.sh tracecb2        # traceCB^2 grids
#   bash src/simulation/others/run.sh masked_omega    # masked-omega grids
#   bash src/simulation/others/run.sh mashr           # mashr benchmark grids
#   bash src/simulation/others/run.sh power_gain      # power-gain figure grids
#   bash src/simulation/others/run.sh all             # run all of the above
#
# Common server overrides:
#   SIM_DATA_DIR=/path/to/simulation/data OUT_DIR=/path/to/result \
#     NREP=100 NSNP=2000 OMEGA_MODE=both bash src/simulation/others/run.sh all
#
# OMEGA_MODE can be both, estimate, or true for scripts with omega modes.
# RUN_VISUALS=0 skips plotting.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../run_common.sh
source "${SCRIPT_DIR}/../run_common.sh"
setup_simulation_env

TARGET="${1:-simulation}"

run_target() {
    local target="$1"
    case "${target}" in
        simulation)
            run_cmd bash src/simulation/run_simulation.sh
            ;;
        robustness)
            run_cmd bash src/simulation/others/run_robustness.sh
            ;;
        tracecb2)
            run_cmd bash src/simulation/others/run_tracecb2.sh
            ;;
        masked_omega)
            run_cmd bash src/simulation/others/run_masked_omega.sh
            ;;
        mashr)
            run_cmd bash src/simulation/others/run_mashr.sh
            ;;
        power_gain)
            run_cmd bash src/simulation/others/run_power_gain.sh
            ;;
        *)
            echo "Unknown target: ${target}" >&2
            echo "Expected one of: simulation, robustness, tracecb2, masked_omega, mashr, power_gain, all" >&2
            return 1
            ;;
    esac
}

if [[ "${TARGET}" == "all" ]]; then
    run_target simulation
    run_target robustness
    run_target tracecb2
    run_target masked_omega
    run_target mashr
    run_target power_gain
else
    run_target "${TARGET}"
fi
