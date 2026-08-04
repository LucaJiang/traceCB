#!/usr/bin/env bash
set -euo pipefail

# Run colocalization for every configured GWAS/eQTL pair.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
activate_conda_env "${R_ENV}"
MAX_JOBS="${MAX_JOBS:-8}"
if [[ ! "${MAX_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_JOBS must be a positive integer; got ${MAX_JOBS}." >&2
    exit 2
fi
COLOC_OUTPUT_DIR="${OUTPUT_DIR}/coloc"
mkdir -p "${COLOC_OUTPUT_DIR}"
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

EQTL_DIRS=()
EQTL_LABELS=()
for i in "${!STUDY_IDS[@]}"; do
    EQTL_DIRS+=("${OUTPUT_DIR}/${STUDY_IDS[$i]}/GMM")
    EQTL_LABELS+=("${TISSUE_SOURCE}_${STUDY_IDS[$i]}_${CELL_TYPES[$i]}")
done

GWAS_FILES=()
LOCI_FILES=()
add_trait() {
    local cohort="$1"
    local trait="$2"
    GWAS_FILES+=("${COLOC_INPUT_DIR}/${cohort}/${cohort}_${trait}_GWAS.csv")
    LOCI_FILES+=("${COLOC_INPUT_DIR}/${cohort}/${cohort}_${trait}_loci.csv")
}

for trait in mon lym neu; do add_trait bcx "${trait}"; done
for trait in ra asthma atopy; do add_trait bbj "${trait}"; done
for trait in wbc mchc rbc mcv mpv plt hct hgb mch eos bas; do
    add_trait bcx "${trait}"
done

run_coloc_task() {
    local trait_index="$1"
    local study_index="$2"
    local gwas_file="${GWAS_FILES[$trait_index]}"
    local loci_file="${LOCI_FILES[$trait_index]}"
    local trait_name="$(basename "${gwas_file}" _GWAS.csv)"
    local label="${EQTL_LABELS[$study_index]}"

    Rscript "${SRC_DIR}/coloc/run_colocalization.R" \
        "${loci_file}" \
        "${gwas_file}" \
        "${EQTL_DIRS[$study_index]}" \
        "${COLOC_OUTPUT_DIR}/${trait_name}_${label}_coloc.csv" \
        >> "${COLOC_OUTPUT_DIR}/${trait_name}_${label}.log" 2>&1
}

for trait_index in "${!GWAS_FILES[@]}"; do
    for study_index in "${!EQTL_DIRS[@]}"; do
        if (( ${#pids[@]} >= MAX_JOBS )); then
            wait_for_batch "${pids[@]}"
            pids=()
        fi
        run_coloc_task "${trait_index}" "${study_index}" &
        pids+=("$!")
    done
done
wait_for_batch "${pids[@]}"

activate_conda_env "${PYTHON_ENV}"
PYTHON_BIN="${PYTHON_BIN:-python}"
"${PYTHON_BIN}" "${SRC_DIR}/coloc/summarize_results.py" \
    --input-dir "${COLOC_OUTPUT_DIR}"
