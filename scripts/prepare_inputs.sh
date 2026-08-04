#!/usr/bin/env bash
set -euo pipefail

# Harmonize target, auxiliary, and tissue summary statistics for every study.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
activate_conda_env "${PYTHON_ENV}"
PYTHON_BIN="${PYTHON_BIN:-python}"

pids=()
for i in "${!STUDY_IDS[@]}"; do
    (
        study_id="${STUDY_IDS[$i]}"
        cell_type="${CELL_TYPES[$i]}"
        for chromosome in "${CHROMOSOMES[@]}"; do
            echo "Preparing ${study_id} (${cell_type}), chr${chromosome}" \
                >> "${LOG_DIR}/prepare_inputs.log"
            "${PYTHON_BIN}" "${SRC_DIR}/preprocess/harmonize_inputs.py" \
                --study "${study_id}" \
                --cell-type "${cell_type}" \
                --chromosome "${chromosome}" \
                --target-data-dir "${TARGET_EQTL_DIR}" \
                --target-sample-size "${TARGET_SAMPLE_SIZES[$i]}" \
                --tissue-data-dir "${TISSUE_DIR}" \
                --tissue-sample-size "${TISSUE_SAMPLE_SIZES[$i]}" \
                --aux-data-dir "${AUX_EQTL_DIR}" \
                --target-ld-dir "${TARGET_LD_DIR}" \
                --aux-ld-dir "${AUX_LD_DIR}" \
                --output-dir "${OUTPUT_DIR}" \
                >> "${LOG_DIR}/prepare_inputs.log" 2>&1
        done

        cp "${CELL_TYPE_PROPORTION_FILE}" "${OUTPUT_DIR}/${study_id}/"
    ) &
    pids+=("$!")
done

wait_status=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        wait_status=1
    fi
done
if (( wait_status != 0 )); then
    echo "One or more input-preparation jobs failed; inspect ${LOG_DIR}/prepare_inputs.log." >&2
    exit "${wait_status}"
fi
echo "Prepared ${TARGET_POPULATION} inputs with ${TISSUE_SOURCE} tissue data."
