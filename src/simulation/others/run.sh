#!/usr/bin/env bash
set -euo pipefail

# Simulation shell entrypoint.
#
# Usage:
#   bash src/simulation/others/run.sh simulation      # simulation.py paper grids
#   bash src/simulation/others/run.sh robustness      # simulation_robustness.py supplementary grids
#   bash src/simulation/others/run.sh tracecb2        # simulation_tracecb2.py grids
#   bash src/simulation/others/run.sh masked_omega    # simulation_masked_omega.py grids
#   bash src/simulation/others/run.sh mashr           # mashr benchmark grids
#   bash src/simulation/others/run.sh power_gain      # rebuttal power-gain figure
#   bash src/simulation/others/run.sh all             # run all of the above in order
#
# Each sub-script activates conda env py312 and runs with python3.

TARGET="${1:-simulation}" # Default to "simulation" if no argument is provided.

case "${TARGET}" in
    simulation)
        bash src/simulation/run_simulation.sh
        ;;
    robustness)
        bash src/simulation/others/run_robustness.sh
        ;;
    tracecb2)
        bash src/simulation/others/run_tracecb2.sh
        ;;
    masked_omega)
        bash src/simulation/others/run_masked_omega.sh
        ;;
    mashr)
        bash src/simulation/others/run_mashr.sh
        ;;
    power_gain)
        bash src/simulation/others/run_power_gain.sh
        ;;
    all)
        bash src/simulation/run_simulation.sh
        bash src/simulation/others/run_robustness.sh
        bash src/simulation/others/run_tracecb2.sh
        bash src/simulation/others/run_masked_omega.sh
        bash src/simulation/others/run_mashr.sh
        ;;
    *)
        echo "Unknown target: ${TARGET}" >&2
        echo "Expected one of: simulation, robustness, tracecb2, masked_omega, mashr, power_gain, all" >&2
        exit 1
        ;;
esac
