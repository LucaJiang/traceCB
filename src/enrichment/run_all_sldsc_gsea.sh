#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RESULT_DIR="${RESULT_DIR:-/home/group1/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea}"
PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/envs/py312/bin/python}"
LDSC_L2_MAX_JOBS="${LDSC_L2_MAX_JOBS:-72}"
LDSC_H2_MAX_JOBS="${LDSC_H2_MAX_JOBS:-72}"
OVERWRITE="${OVERWRITE:-1}"

mkdir -p "${RESULT_DIR}/logs"
LOG_FILE="${RESULT_DIR}/logs/run_all_sldsc_gsea.$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[start] $(date)"
echo "[config] RESULT_DIR=${RESULT_DIR}"
echo "[config] LDSC_L2_MAX_JOBS=${LDSC_L2_MAX_JOBS}"
echo "[config] LDSC_H2_MAX_JOBS=${LDSC_H2_MAX_JOBS}"
echo "[config] OVERWRITE=${OVERWRITE}"

PREPARE_ARGS=("--result-dir" "${RESULT_DIR}")
if [[ "${OVERWRITE}" == "1" ]]; then
  PREPARE_ARGS+=("--overwrite")
fi
"${PYTHON_BIN}" "${SCRIPT_DIR}/01_prepare_annotations.py" "${PREPARE_ARGS[@]}"

RESULT_DIR="${RESULT_DIR}" \
MAX_JOBS="${LDSC_L2_MAX_JOBS}" \
OVERWRITE="${OVERWRITE}" \
bash "${SCRIPT_DIR}/02_compute_ldscores.sh"

RESULT_DIR="${RESULT_DIR}" \
MAX_JOBS="${LDSC_H2_MAX_JOBS}" \
OVERWRITE="${OVERWRITE}" \
bash "${SCRIPT_DIR}/03_run_h2.sh"

"${PYTHON_BIN}" "${SCRIPT_DIR}/04_aggregate_results.py" --result-dir "${RESULT_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/05_visualize_results.py" --result-dir "${RESULT_DIR}"

echo "[check] LD score files: $(find "${RESULT_DIR}/annotations/ldscores" -name '*.l2.ldscore.gz' | wc -l)"
echo "[check] h2 result files: $(find "${RESULT_DIR}/results/raw" -name '*.results' | wc -l)"
echo "[done] $(date)"
echo "[log] ${LOG_FILE}"
