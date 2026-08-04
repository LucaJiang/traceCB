#!/usr/bin/env bash

# Shared configuration for the full-data workflows.
# Override any path with the corresponding environment variable; no personal
# filesystem paths are assumed.

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${TRACECB_REPO_ROOT:-$(cd "${CONFIG_DIR}/.." && pwd)}"

TISSUE_SOURCE="${TISSUE_SOURCE:-eQTLGen}"  # eQTLGen or GTEx
TARGET_POPULATION="${TARGET_POPULATION:-EAS}"  # EAS or AFR
CHROMOSOMES=({1..22})

DATA_ROOT="${TRACECB_DATA_ROOT:-${REPO_ROOT}/data}"
OUTPUT_ROOT="${TRACECB_OUTPUT_ROOT:-${REPO_ROOT}/results}"
GTEX_DIR="${GTEX_DIR:-${DATA_ROOT}/GTEx/GTEx_Whole_Blood_by_chr}"
EQTLGEN_DIR="${EQTLGEN_DIR:-${DATA_ROOT}/eQTLGen}"
CELL_TYPE_PROPORTION_FILE="${CELL_TYPE_PROPORTION_FILE:-${DATA_ROOT}/GTEx/celltype_proportion.csv}"
EQTL_CATALOGUE_DIR="${EQTL_CATALOGUE_DIR:-${DATA_ROOT}/eQTLCatalogue/by_celltype_chr}"
BBJ_DIR="${BBJ_DIR:-${DATA_ROOT}/BBJ_eQTL/by_celltype_chr}"
AFR_DIR="${AFR_DIR:-${DATA_ROOT}/popcell/AFB_NS}"
LD_REFERENCE_DIR="${LD_REFERENCE_DIR:-${DATA_ROOT}/1000G}"
AUX_LD_DIR="${AUX_LD_DIR:-${LD_REFERENCE_DIR}/1000G_EUR}"
COLOC_INPUT_DIR="${COLOC_INPUT_DIR:-${DATA_ROOT}/coloc}"

SRC_DIR="${REPO_ROOT}/src"
LOG_DIR="${TRACECB_LOG_DIR:-${OUTPUT_ROOT}/logs}"
PYTHON_ENV="${PYTHON_ENV:-py312}"
R_ENV="${R_ENV:-r4}"
PLINK_BIN="${PLINK_BIN:-plink}"
SLDXR_DIR="${SLDXR_DIR:-${REPO_ROOT}/external/s-ldxr}"

STUDY_IDS=(
    "QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069"
    "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372"
)
CELL_TYPES=(
    "Monocytes" "CD4+T_cells" "CD8+T_cells" "CD4+T_cells" "Monocytes"
    "B_cells" "Monocytes" "NK_cells" "CD4+T_cells" "CD8+T_cells"
)
GTEX_SAMPLE_SIZE=670
EQTLGEN_SAMPLE_SIZE=30000
BBJ_SAMPLE_SIZES=(105 103 103 103 105 104 105 104 103 103)
EQTL_CATALOGUE_SAMPLE_SIZES=(191 167 277 290 286 262 420 247 280 269)
AFR_SAMPLE_SIZES=(80 80 80 80 80 80 80 80 80 80)

AUX_SAMPLE_SIZES=("${EQTL_CATALOGUE_SAMPLE_SIZES[@]}")
AUX_EQTL_DIR="${EQTL_CATALOGUE_DIR}"

case "${TARGET_POPULATION}" in
    EAS)
        TARGET_EQTL_DIR="${BBJ_DIR}"
        TARGET_SAMPLE_SIZES=("${BBJ_SAMPLE_SIZES[@]}")
        TARGET_LD_DIR="${LD_REFERENCE_DIR}/1000G_EAS"
        ;;
    AFR)
        TARGET_EQTL_DIR="${AFR_DIR}"
        TARGET_SAMPLE_SIZES=("${AFR_SAMPLE_SIZES[@]}")
        TARGET_LD_DIR="${LD_REFERENCE_DIR}/1000G_AFR"
        ;;
    *)
        echo "TARGET_POPULATION must be EAS or AFR; got ${TARGET_POPULATION}." >&2
        return 2 2>/dev/null || exit 2
        ;;
esac

case "${TISSUE_SOURCE}" in
    GTEx)
        TISSUE_DIR="${GTEX_DIR}"
        TISSUE_SAMPLE_SIZES=()
        for _ in "${STUDY_IDS[@]}"; do
            TISSUE_SAMPLE_SIZES+=("${GTEX_SAMPLE_SIZE}")
        done
        ;;
    eQTLGen)
        TISSUE_DIR="${EQTLGEN_DIR}"
        TISSUE_SAMPLE_SIZES=()
        for _ in "${STUDY_IDS[@]}"; do
            TISSUE_SAMPLE_SIZES+=("${EQTLGEN_SAMPLE_SIZE}")
        done
        ;;
    *)
        echo "TISSUE_SOURCE must be eQTLGen or GTEx; got ${TISSUE_SOURCE}." >&2
        return 2 2>/dev/null || exit 2
        ;;
esac

OUTPUT_DIR="${OUTPUT_ROOT}/${TARGET_POPULATION}_${TISSUE_SOURCE}"

activate_conda_env() {
    local env_name="$1"
    # shellcheck source=/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}"
