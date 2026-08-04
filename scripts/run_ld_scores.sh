#!/usr/bin/env bash
set -euo pipefail

# Build gene annotations and compute cross-population LD scores.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
activate_conda_env "${PYTHON_ENV}"
PYTHON_BIN="${PYTHON_BIN:-python}"
MAX_JOBS="${MAX_JOBS:-10}"
if [[ ! "${MAX_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_JOBS must be a positive integer; got ${MAX_JOBS}." >&2
    exit 2
fi
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

# Annotations are prerequisites for every chromosome task. Build them first so
# MAX_JOBS below is a repository-wide limit rather than a per-study limit.
for i in "${!STUDY_IDS[@]}"; do
    study_id="${STUDY_IDS[$i]}"
    ld_dir="${OUTPUT_DIR}/${study_id}/LDSC"
    "${PYTHON_BIN}" "${SRC_DIR}/preprocess/build_ld_annotations.py" \
        --study "${study_id}" \
        --output-dir "${OUTPUT_DIR}" \
        --target-ld-dir "${TARGET_LD_DIR}" \
        --aux-ld-dir "${AUX_LD_DIR}" \
        >> "${LOG_DIR}/build_ld_annotations.log" 2>&1
    mkdir -p "${ld_dir}/TAR" "${ld_dir}/AUX" "${ld_dir}/LDSC_gene"
done

run_ld_task() {
    local study_index="$1"
    local chromosome="$2"
    local study_id="${STUDY_IDS[$study_index]}"
    local annotation_dir="${OUTPUT_DIR}/${study_id}/LDSC/LD_annotation"
    local ld_dir="${OUTPUT_DIR}/${study_id}/LDSC"
    local target_bim target_prefix aux_bim aux_prefix

    target_bim="$(find "${TARGET_LD_DIR}" -name "1000G.*.QC.maf.${chromosome}.bim" -print -quit)"
    aux_bim="$(find "${AUX_LD_DIR}" -name "1000G.*.QC.maf.${chromosome}.bim" -print -quit)"
    if [[ -z "${target_bim}" || -z "${aux_bim}" ]]; then
        echo "Missing LD reference BIM for chromosome ${chromosome}." >&2
        return 1
    fi
    target_prefix="$(basename "${target_bim}" .bim)"
    aux_prefix="$(basename "${aux_bim}" .bim)"

    "${PLINK_BIN}" \
        --bfile "${TARGET_LD_DIR}/${target_prefix}" \
        --extract "${annotation_dir}/${chromosome}.print_snps.txt" \
        --keep-allele-order \
        --make-bed \
        --out "${ld_dir}/TAR/1000G.TAR.QC.maf.${chromosome}"
    "${PLINK_BIN}" \
        --bfile "${AUX_LD_DIR}/${aux_prefix}" \
        --extract "${annotation_dir}/${chromosome}.print_snps.txt" \
        --keep-allele-order \
        --make-bed \
        --out "${ld_dir}/AUX/1000G.AUX.QC.maf.${chromosome}"
    "${PYTHON_BIN}" "${SLDXR_DIR}/s-ldxr.py" \
        --bfile \
            "${ld_dir}/TAR/1000G.TAR.QC.maf.${chromosome}" \
            "${ld_dir}/AUX/1000G.AUX.QC.maf.${chromosome}" \
        --print-snps "${annotation_dir}/${chromosome}.print_snps.txt" \
        --annot "${annotation_dir}/${chromosome}.annot.gz" \
        --ld-wind-cm 1.0 \
        --score standardized \
        --out "${ld_dir}/LDSC_gene/TAR_AUX_std_chr${chromosome}" \
        >> "${LOG_DIR}/ld_scores.log" 2>&1
}

for i in "${!STUDY_IDS[@]}"; do
    for chromosome in "${CHROMOSOMES[@]}"; do
        if (( ${#pids[@]} >= MAX_JOBS )); then
            wait_for_batch "${pids[@]}"
            pids=()
        fi
        run_ld_task "${i}" "${chromosome}" &
        pids+=("$!")
    done
done
wait_for_batch "${pids[@]}"
echo "Computed LD scores for ${TARGET_POPULATION} with ${TISSUE_SOURCE} tissue data."
