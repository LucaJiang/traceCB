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

RESULT_DIR="${RESULT_DIR:-/home/group1/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea}"
LDSC_DIR="${LDSC_DIR:-/home/group1/wjiang49/software/ldsc}"
BASELINE_LD_PREFIX="${BASELINE_LD_PREFIX:-/home/group1/wjiang49/data/traceCB/EAS_eQTLGen/results/disease_heritability_sldsc/reference/baselineLD_joint/baselineLD.}"
WEIGHTS_LD_PREFIX="${WEIGHTS_LD_PREFIX:-/home/group1/wjiang49/data/1000G/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.}"
FRQ_PREFIX="${FRQ_PREFIX:-/home/group1/wjiang49/data/1000G/1000G_EAS_EUR/EAS/1000G.EAS.QC.}"
MAX_JOBS="${MAX_JOBS:-72}"
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

wait_for_slot() {
  while [[ "$(jobs -r -p | wc -l)" -ge "${MAX_JOBS}" ]]; do
    if ! wait -n; then
      FAILURES=1
    fi
  done
}

wait_for_all() {
  while [[ "$(jobs -r -p | wc -l)" -gt 0 ]]; do
    if ! wait -n; then
      FAILURES=1
    fi
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

for trait, sumstats in traits:
    for annot_id, prefix in annots:
        print(f"{trait}\t{sumstats}\t{annot_id}\t{prefix}")
PY

while IFS=$'\t' read -r trait sumstats annot_id prefix; do
  wait_for_slot
  run_one "${trait}" "${sumstats}" "${annot_id}" "${prefix}" &
done < "${RESULT_DIR}/metadata/h2_jobs.tsv"

wait_for_all
if [[ "${FAILURES}" != "0" ]]; then
  echo "[error] One or more h2 jobs failed. See ${LOG_DIR}" >&2
  exit 1
fi

echo "[done] h2 results written under ${RAW_DIR}"
