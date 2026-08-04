#!/usr/bin/env bash
# Step 04: compute the custom LD scores.
set -euo pipefail

if [[ "${SKIP_CONDA_ACTIVATE:-0}" == "1" ]]; then
  :
elif [[ -f /opt/anaconda3/etc/profile.d/conda.sh ]]; then
  source /opt/anaconda3/etc/profile.d/conda.sh
  conda activate ldsc
else
  source activate ldsc
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export PYTHONWARNINGS="${PYTHONWARNINGS:-ignore}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RESULT_DIR="${RESULT_DIR:-${TRACECB_ENRICHMENT_DIR:-${REPO_ROOT}/results/enrichment}}"
LDSC_DIR="${LDSC_DIR:-${REPO_ROOT}/external/ldsc}"
BFILE_PREFIX="${BFILE_PREFIX:-${RESULT_DIR}/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.}"
PRINT_SNPS="${PRINT_SNPS:-${RESULT_DIR}/reference/hm3_no_MHC.list.txt}"
MAX_JOBS="${MAX_JOBS:-8}"
OVERWRITE="${OVERWRITE:-0}"
MANIFEST="${RESULT_DIR}/metadata/annotation_manifest.tsv"
LOG_DIR="${RESULT_DIR}/logs"

mkdir -p "${LOG_DIR}"
if [[ ! -f "${MANIFEST}" ]]; then
  echo "[error] Missing annotation manifest: ${MANIFEST}" >&2
  exit 1
fi

FAILURES=0
ACTIVE_JOBS=0

if [[ ! "${MAX_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "[error] MAX_JOBS must be a positive integer: ${MAX_JOBS}" >&2
  exit 1
fi

wait_for_slot() {
  if (( ACTIVE_JOBS >= MAX_JOBS )); then
    if ! wait -n; then
      FAILURES=1
    fi
    ACTIVE_JOBS=$((ACTIVE_JOBS - 1))
  fi
}

wait_for_all() {
  while (( ACTIVE_JOBS > 0 )); do
    if ! wait -n; then
      FAILURES=1
    fi
    ACTIVE_JOBS=$((ACTIVE_JOBS - 1))
  done
}

run_one() {
  local annot_id="$1"
  local prefix="$2"
  local chrom="$3"
  local out_prefix="${prefix}${chrom}"
  local log_file="${LOG_DIR}/02_l2.${annot_id}.chr${chrom}.log"

  if [[ "${OVERWRITE}" != "1" && -f "${out_prefix}.l2.ldscore.gz" ]]; then
    echo "[skip] ${annot_id} chr${chrom}"
    return 0
  fi

  echo "[l2] ${annot_id} chr${chrom}"
  python "${LDSC_DIR}/ldsc.py" \
    --l2 \
    --bfile "${BFILE_PREFIX}${chrom}" \
    --ld-wind-cm 1 \
    --annot "${out_prefix}.annot.gz" \
    --thin-annot \
    --print-snps "${PRINT_SNPS}" \
    --out "${out_prefix}" \
    > "${log_file}" 2>&1

  if [[ ! -f "${out_prefix}.l2.ldscore.gz" ]]; then
    echo "[error] Missing output ${out_prefix}.l2.ldscore.gz" >&2
    return 1
  fi
}

python - "${MANIFEST}" <<'PY' > "${RESULT_DIR}/metadata/annotation_prefixes.tsv"
import csv
import sys

seen = set()
print("AnnotID\tAnnotPrefix")
with open(sys.argv[1], newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        key = row["AnnotID"]
        if key in seen:
            continue
        seen.add(key)
        print(f"{row['AnnotID']}\t{row['AnnotPrefix']}")
PY

{
  IFS= read -r _header
  while IFS=$'\t' read -r annot_id prefix; do
    for chrom in {1..22}; do
      wait_for_slot
      run_one "${annot_id}" "${prefix}" "${chrom}" &
      ACTIVE_JOBS=$((ACTIVE_JOBS + 1))
    done
  done
} < "${RESULT_DIR}/metadata/annotation_prefixes.tsv"

wait_for_all
if [[ "${FAILURES}" != "0" ]]; then
  echo "[error] One or more LD score jobs failed. See ${LOG_DIR}" >&2
  exit 1
fi

echo "[done] LD scores written under ${RESULT_DIR}/annotations/ldscores"
