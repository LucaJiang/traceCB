#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RESULT_DIR="${RESULT_DIR:-${REPO_ROOT}/output/sldsc_gsea_eur_release_matched}"
PYTHON_BIN="${PYTHON_BIN:-python}"
STUDY_DIR="${STUDY_DIR:?Set STUDY_DIR to the EAS_eQTLGen study-data directory}"
GWAS_ROOT="${GWAS_ROOT:?Set GWAS_ROOT to the Pan-UK Biobank sumstats directory}"
LDSC_DIR="${LDSC_DIR:?Set LDSC_DIR to an LDSC source checkout}"
LDSC_L2_MAX_JOBS="${LDSC_L2_MAX_JOBS:-72}"
LDSC_H2_MAX_JOBS="${LDSC_H2_MAX_JOBS:-72}"
OVERWRITE="${OVERWRITE:-1}"
BFILE_PREFIX="${BFILE_PREFIX:-${RESULT_DIR}/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.}"
FRQ_PREFIX="${FRQ_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_frq/1000G.EUR.QC.}"
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_baselineLD_v2.2_exact_hm3/baselineLD.}"
WEIGHTS_LD_PREFIX="${WEIGHTS_LD_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.}"
PRINT_SNPS="${PRINT_SNPS:-${RESULT_DIR}/reference/hm3_no_MHC.list.txt}"
TRACECB_REVISION="$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo unavailable)"
LDSC_REVISION="$(git -C "${LDSC_DIR}" rev-parse HEAD 2>/dev/null || echo unavailable)"
if [[ -n "$(git -C "${REPO_ROOT}" status --porcelain 2>/dev/null || true)" ]]; then
  TRACECB_GIT_DIRTY=1
else
  TRACECB_GIT_DIRTY=0
fi

mkdir -p "${RESULT_DIR}/logs"
LOG_FILE="${RESULT_DIR}/logs/run_all_sldsc_gsea.$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[start] $(date)"
echo "[config] RESULT_DIR=${RESULT_DIR}"
echo "[config] LDSC_L2_MAX_JOBS=${LDSC_L2_MAX_JOBS}"
echo "[config] LDSC_H2_MAX_JOBS=${LDSC_H2_MAX_JOBS}"
echo "[config] OVERWRITE=${OVERWRITE}"
echo "[config] STUDY_DIR=${STUDY_DIR}"
echo "[config] GWAS_ROOT=${GWAS_ROOT}"
echo "[config] LDSC_DIR=${LDSC_DIR}"
echo "[config] TRACECB_REVISION=${TRACECB_REVISION}"
echo "[config] TRACECB_GIT_DIRTY=${TRACECB_GIT_DIRTY}"
echo "[config] LDSC_REVISION=${LDSC_REVISION}"
echo "[config] PYTHON_VERSION=$("${PYTHON_BIN}" --version 2>&1)"
echo "[config] BFILE_PREFIX=${BFILE_PREFIX}"
echo "[config] BASELINE_LD_PREFIX=${BASELINE_LD_PREFIX}"
echo "[config] WEIGHTS_LD_PREFIX=${WEIGHTS_LD_PREFIX}"
echo "[config] FRQ_PREFIX=${FRQ_PREFIX}"

RESULT_DIR="${RESULT_DIR}" BFILE_PREFIX="${BFILE_PREFIX}" \
PYTHON_BIN="${PYTHON_BIN}" \
  bash "${SCRIPT_DIR}/00_prepare_eur_reference.sh"

PREPARE_ARGS=(
  "--result-dir" "${RESULT_DIR}"
  "--study-dir" "${STUDY_DIR}"
  "--gwas-root" "${GWAS_ROOT}"
  "--bim-prefix" "${BFILE_PREFIX}"
)
if [[ "${OVERWRITE}" == "1" ]]; then
  PREPARE_ARGS+=("--overwrite")
fi
"${PYTHON_BIN}" "${SCRIPT_DIR}/01_prepare_annotations.py" "${PREPARE_ARGS[@]}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_validate_eur_stack.py" \
  --result-dir "${RESULT_DIR}" \
  --study-dir "${STUDY_DIR}" \
  --bfile-prefix "${BFILE_PREFIX}" \
  --frq-prefix "${FRQ_PREFIX}"

RESULT_DIR="${RESULT_DIR}" \
MAX_JOBS="${LDSC_L2_MAX_JOBS}" \
OVERWRITE="${OVERWRITE}" \
BFILE_PREFIX="${BFILE_PREFIX}" \
PRINT_SNPS="${PRINT_SNPS}" \
LDSC_DIR="${LDSC_DIR}" \
bash "${SCRIPT_DIR}/02_compute_ldscores.sh"

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_filter_baseline_to_regression_snps.py" \
  --source-prefix "${RESULT_DIR}/reference/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD." \
  --output-prefix "${BASELINE_LD_PREFIX}" \
  --regression-snps "${PRINT_SNPS}" \
  --custom-prefix "${RESULT_DIR}/annotations/ldscores/incremental_QTD000021_original/incremental_QTD000021_original."

"${PYTHON_BIN}" "${SCRIPT_DIR}/00_validate_custom_ldscores.py" \
  --result-dir "${RESULT_DIR}"

RESULT_DIR="${RESULT_DIR}" \
MAX_JOBS="${LDSC_H2_MAX_JOBS}" \
OVERWRITE="${OVERWRITE}" \
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX}" \
WEIGHTS_LD_PREFIX="${WEIGHTS_LD_PREFIX}" \
FRQ_PREFIX="${FRQ_PREFIX}" \
LDSC_DIR="${LDSC_DIR}" \
bash "${SCRIPT_DIR}/03_run_h2.sh"

"${PYTHON_BIN}" "${SCRIPT_DIR}/04_aggregate_results.py" --result-dir "${RESULT_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/05_visualize_results.py" --result-dir "${RESULT_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/06_write_manuscript_sections.py" --result-dir "${RESULT_DIR}"

echo "[check] LD score files: $(find "${RESULT_DIR}/annotations/ldscores" -name '*.l2.ldscore.gz' | wc -l)"
echo "[check] h2 result files: $(find "${RESULT_DIR}/results/raw" -name '*.results' | wc -l)"
echo "[done] $(date)"
echo "[log] ${LOG_FILE}"
