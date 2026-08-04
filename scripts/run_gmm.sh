#!/usr/bin/env bash
set -euo pipefail

# Run traceCB for every configured study and chromosome.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
activate_conda_env "${PYTHON_ENV}"
PYTHON_BIN="${PYTHON_BIN:-python}"
MAX_JOBS="${MAX_JOBS:-8}"
if [[ ! "${MAX_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_JOBS must be a positive integer; got ${MAX_JOBS}." >&2
    exit 2
fi
TIME_LOG="${LOG_DIR}/run_gmm_timing.log"
: > "${TIME_LOG}"
pids=()

wait_for_batch() {
    local pid status=0
    for pid in "$@"; do
        if ! wait "${pid}"; then
            status=1
        fi
    done
    return "${status}"
}

format_time() {
    local seconds="$1"
    printf '%02d:%02d:%02d' \
        $((seconds / 3600)) $(((seconds % 3600) / 60)) $((seconds % 60))
}

run_gmm_task() {
    local study_index="$1"
    local chromosome="$2"
    local study_id="${STUDY_IDS[$study_index]}"
    local cell_type="${CELL_TYPES[$study_index]}"
    local started finished elapsed

    started="$(date +%s)"
    echo "Starting ${study_id}, chr${chromosome}"
    "${PYTHON_BIN}" "${SRC_DIR}/traceCB/run_gmm.py" \
        --study "${study_id}" \
        --cell-type "${cell_type}" \
        --chromosome "${chromosome}" \
        --data-dir "${OUTPUT_DIR}" \
        >> "${LOG_DIR}/run_gmm.log" 2>&1
    finished="$(date +%s)"
    elapsed=$((finished - started))
    echo "Finished ${study_id}, chr${chromosome} in $(format_time "${elapsed}")"
}

started="$(date +%s)"
for study_index in "${!STUDY_IDS[@]}"; do
    for chromosome in "${CHROMOSOMES[@]}"; do
        if (( ${#pids[@]} >= MAX_JOBS )); then
            wait_for_batch "${pids[@]}"
            pids=()
        fi
        run_gmm_task "${study_index}" "${chromosome}" &
        pids+=("$!")
    done
done
wait_for_batch "${pids[@]}"

elapsed=$(($(date +%s) - started))
echo "All traceCB jobs completed in $(format_time "${elapsed}")." | tee -a "${TIME_LOG}"
