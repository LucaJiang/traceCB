#!/usr/bin/env bash
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
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_baselineLD_v2.2_exact_hm3/baselineLD.}"
WEIGHTS_LD_PREFIX="${WEIGHTS_LD_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.}"
FRQ_PREFIX="${FRQ_PREFIX:-${RESULT_DIR}/reference/1000G_Phase3_frq/1000G.EUR.QC.}"
MAX_JOBS="${MAX_JOBS:-8}"
OVERWRITE="${OVERWRITE:-0}"

TRAIT_MANIFEST="${RESULT_DIR}/metadata/trait_manifest.tsv"
ANNOT_MANIFEST="${RESULT_DIR}/metadata/annotation_manifest.tsv"
RAW_DIR="${RESULT_DIR}/results/raw"
LOG_DIR="${RESULT_DIR}/logs"

mkdir -p "${RAW_DIR}" "${LOG_DIR}"
if [[ ! -f "${TRAIT_MANIFEST}" || ! -f "${ANNOT_MANIFEST}" ]]; then
  echo "[error] Missing manifests under ${RESULT_DIR}/metadata" >&2
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
  local trait="$1"
  local sumstats="$2"
  local annot_id="$3"
  local prefix="$4"
  local out_prefix="${RAW_DIR}/${trait}__${annot_id}"
  local log_file="${LOG_DIR}/03_h2.${trait}.${annot_id}.log"

  if [[ "${OVERWRITE}" != "1" && -f "${out_prefix}.results" ]]; then
    echo "[skip] ${trait} ${annot_id}"
    return 0
  fi
  if [[ ! -f "${sumstats}" ]]; then
    echo "[error] Missing sumstats: ${sumstats}" >&2
    return 1
  fi
  if [[ ! -f "${prefix}1.l2.ldscore.gz" ]]; then
    echo "[error] Missing custom LD score prefix: ${prefix}" >&2
    return 1
  fi

  echo "[h2] ${trait} ${annot_id}"
  python "${LDSC_DIR}/ldsc.py" \
    --h2 "${sumstats}" \
    --ref-ld-chr "${BASELINE_LD_PREFIX},${prefix}" \
    --w-ld-chr "${WEIGHTS_LD_PREFIX}" \
    --overlap-annot \
    --frqfile-chr "${FRQ_PREFIX}" \
    --out "${out_prefix}" \
    --print-coefficients \
    > "${log_file}" 2>&1

  if [[ ! -f "${out_prefix}.results" ]]; then
    echo "[error] Missing output ${out_prefix}.results" >&2
    return 1
  fi
}

python - "${TRAIT_MANIFEST}" "${ANNOT_MANIFEST}" <<'PY' > "${RESULT_DIR}/metadata/h2_jobs.tsv"
import csv
import sys

traits = []
with open(sys.argv[1], newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        traits.append((row["Trait"], row["SumstatsPath"]))

annots = []
seen = set()
with open(sys.argv[2], newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["AnnotID"] in seen:
            continue
        seen.add(row["AnnotID"])
        annots.append((row["AnnotID"], row["AnnotPrefix"]))

print("Trait\tSumstatsPath\tAnnotID\tAnnotPrefix")
for trait, sumstats in traits:
    for annot_id, prefix in annots:
        print(f"{trait}\t{sumstats}\t{annot_id}\t{prefix}")
PY

{
  IFS= read -r _header
  while IFS=$'\t' read -r trait sumstats annot_id prefix; do
    wait_for_slot
    run_one "${trait}" "${sumstats}" "${annot_id}" "${prefix}" &
    ACTIVE_JOBS=$((ACTIVE_JOBS + 1))
  done
} < "${RESULT_DIR}/metadata/h2_jobs.tsv"

wait_for_all
if [[ "${FAILURES}" != "0" ]]; then
  echo "[error] One or more h2 jobs failed. See ${LOG_DIR}" >&2
  exit 1
fi

echo "[done] h2 results written under ${RAW_DIR}"
