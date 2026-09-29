#!/usr/bin/env bash

# Shared configuration for the full-data workflows and all manuscript figures.
# Review the user settings below before running. Derived paths are grouped in
# the defaults section and can be overridden with their environment variables.
# Apply the configuration to the current shell: source scripts/config.sh

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Run 'source scripts/config.sh' from the repository root to configure the current shell." >&2
    exit 2
fi

# 1. User settings: machine-specific roots and Conda environments.
# These values match this machine; update them when using another filesystem.
DATA_ROOT="${TRACECB_DATA_ROOT:-/home/wjiang49/group/wjiang49/data}"
export TRACECB_DATA_ROOT="${DATA_ROOT}"
export TRACECB_SOFTWARE_ROOT="${TRACECB_SOFTWARE_ROOT:-/home/wjiang49/group/wjiang49/software}"
PYTHON_ENV="${PYTHON_ENV:-py312}"
R_ENV="${R_ENV:-r4}"

# 2. Defaults: workflow options and paths derived from the roots above.
# Usually no edits are needed if the data/software directory layout matches.
# Override individual paths below if your files use a different layout.
export TISSUE_SOURCE="${TISSUE_SOURCE:-eQTLGen}"  # eQTLGen or GTEx
export TARGET_POPULATION="${TARGET_POPULATION:-EAS}"  # EAS or AFR; also selects figure bins
CHROMOSOMES=({1..22})

# Repository, output, and software paths.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${TRACECB_REPO_ROOT:-$(cd "${CONFIG_DIR}/.." && pwd)}"
export TRACECB_REPO_ROOT="${REPO_ROOT}"
SRC_DIR="${REPO_ROOT}/src"
OUTPUT_ROOT="${TRACECB_OUTPUT_ROOT:-${REPO_ROOT}/results}"
export TRACECB_OUTPUT_ROOT="${OUTPUT_ROOT}"
OUTPUT_DIR="${OUTPUT_ROOT}/${TARGET_POPULATION}_${TISSUE_SOURCE}"
LOG_DIR="${TRACECB_LOG_DIR:-${OUTPUT_ROOT}/logs}"
PLINK_BIN="${PLINK_BIN:-${TRACECB_SOFTWARE_ROOT}/plink}"
PLINK2="${PLINK2:-${TRACECB_SOFTWARE_ROOT}/plink2}"
SLDXR_DIR="${SLDXR_DIR:-${TRACECB_SOFTWARE_ROOT}/s-ldxr-master}"

# Full-data workflow inputs.
# Existing study results are stored separately from newly generated outputs.
export TRACECB_STUDY_ROOT="${TRACECB_STUDY_ROOT:-${DATA_ROOT}/traceCB}"
GTEX_SOURCE_DIR="${GTEX_SOURCE_DIR:-${DATA_ROOT}/GTEx}"
GTEX_DIR="${GTEX_DIR:-${GTEX_SOURCE_DIR}/GTEx_Whole_Blood_by_chr}"
EQTLGEN_DIR="${EQTLGEN_DIR:-${DATA_ROOT}/eQTLGen}"
CELL_TYPE_PROPORTION_FILE="${CELL_TYPE_PROPORTION_FILE:-${GTEX_SOURCE_DIR}/celltype_proportion.csv}"
EQTL_CATALOGUE_DIR="${EQTL_CATALOGUE_DIR:-${DATA_ROOT}/eQTLCatalogue/by_celltype_chr}"
BBJ_DIR="${BBJ_DIR:-${DATA_ROOT}/BBJ_eQTL/by_celltype_chr}"
AFR_DIR="${AFR_DIR:-${DATA_ROOT}/popcell/AFB_NS}"
LD_REFERENCE_DIR="${LD_REFERENCE_DIR:-${DATA_ROOT}/1000G}"
AUX_LD_DIR="${AUX_LD_DIR:-${LD_REFERENCE_DIR}/1000G_EUR}"
COLOC_INPUT_DIR="${COLOC_INPUT_DIR:-${TRACECB_COLOC_INPUT_DIR:-${TRACECB_STUDY_ROOT}/coloc}}"

# Figure inputs: external data under DATA_ROOT / TRACECB_STUDY_ROOT.
# STUDY_DIR must contain QTD*/GMM/chr*/summary.csv; GENE_ANNOTATION is a GTF file.
export TRACECB_STUDY_DIR="${TRACECB_STUDY_DIR:-${TRACECB_STUDY_ROOT}/${TARGET_POPULATION}_${TISSUE_SOURCE}}"
export TRACECB_AFR_STUDY_DIR="${TRACECB_AFR_STUDY_DIR:-${TRACECB_STUDY_ROOT}/AFR_${TISSUE_SOURCE}}"
export TRACECB_GTEX_GENE_ANNOTATION="${TRACECB_GTEX_GENE_ANNOTATION:-${GTEX_SOURCE_DIR}/gencode.v26.GRCh38.genes.gtf}"
export TRACECB_ONEK1K_FILE="${TRACECB_ONEK1K_FILE:-${TRACECB_STUDY_ROOT}/onek1k_supp/onek1k_esnp.csv}"
export TRACECB_GTEX_LOOKUP="${TRACECB_GTEX_LOOKUP:-${GTEX_LOOKUP_FILE:-${GTEX_SOURCE_DIR}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz}}"
GTEX_LOOKUP_FILE="${GTEX_LOOKUP_FILE:-${TRACECB_GTEX_LOOKUP}}"
export TRACECB_OASIS_DIR="${TRACECB_OASIS_DIR:-${DATA_ROOT}/hum0197/eQTL_summary_statistics}"
export TRACECB_CIMA_DIR="${TRACECB_CIMA_DIR:-${DATA_ROOT}/CIMA}"
export TRACECB_CIMA_LEAD_EQTL="${TRACECB_CIMA_LEAD_EQTL:-${TRACECB_CIMA_DIR}/xQTL/CIMA_Lead_cis-xQTL.csv}"
export TRACECB_CIMA_CELL_TYPES="${TRACECB_CIMA_CELL_TYPES:-${TRACECB_CIMA_DIR}/Cell_Atlas/CIMA_Cell_Type_Level_and_Marker.xlsx}"
export TRACECB_REPLICATION_EGENES="${TRACECB_REPLICATION_EGENES:-${DATA_ROOT}/hum0343/hum0343_eGene.csv}"
export TRACECB_REPLICATION_ESNPS="${TRACECB_REPLICATION_ESNPS:-${DATA_ROOT}/hum0343/hum0343_eSNP.csv}"
export TRACECB_CELL_PROPORTIONS="${TRACECB_CELL_PROPORTIONS:-${TRACECB_STUDY_ROOT}/cell_type_proportion/ind_celltype_proportion.csv}"
export TRACECB_COLOC_INPUT_DIR="${TRACECB_COLOC_INPUT_DIR:-${COLOC_INPUT_DIR}}"
export TRACECB_COLOC_GENES="${TRACECB_COLOC_GENES:-${TRACECB_COLOC_INPUT_DIR}/bcx/bcx_mon.closest.protein_coding.bed}"
export TRACECB_LOCUS_GWAS="${TRACECB_LOCUS_GWAS:-${TRACECB_COLOC_INPUT_DIR}/bcx/bcx_mon_GWAS.csv}"
# Locus plots also require exported per-gene CSVs and these external bigWigs.
export TRACECB_LOCUS_EQTL_DIR="${TRACECB_LOCUS_EQTL_DIR:-${TRACECB_STUDY_DIR}}"
export TRACECB_LOCUS_TRACK_DIR="${TRACECB_LOCUS_TRACK_DIR:-${DATA_ROOT}/locuszoom}"

# Figure outputs, repository metadata, and generated result tables.
# Keep each population/bulk-source configuration's figures and composites together.
export TRACECB_FIGURE_DIR="${TRACECB_FIGURE_DIR:-${OUTPUT_ROOT}/figures/${TARGET_POPULATION}_${TISSUE_SOURCE}}"
export TRACECB_AFR_FIGURE_DIR="${TRACECB_AFR_FIGURE_DIR:-${TRACECB_FIGURE_DIR}/afr}"
export TRACECB_ESNP_REPLICATION_DIR="${TRACECB_ESNP_REPLICATION_DIR:-${TRACECB_FIGURE_DIR}/esnp_replication}"
export TRACECB_COLOC_FIGURE_DIR="${TRACECB_COLOC_FIGURE_DIR:-${TRACECB_FIGURE_DIR}/coloc}"
export TRACECB_FIGURE_METADATA="${TRACECB_FIGURE_METADATA:-${SRC_DIR}/figures/metadata.json}"
export TRACECB_TIMING_FILE="${TRACECB_TIMING_FILE:-${REPO_ROOT}/tmp/timing/summary_timing.csv}"
# Colocalization result tables must be generated by run_colocalization.sh.
export TRACECB_COLOC_DIR="${TRACECB_COLOC_DIR:-${OUTPUT_DIR}/coloc}"
export TRACECB_COLOC_REPLICATION="${TRACECB_COLOC_REPLICATION:-${TRACECB_COLOC_DIR}/replication.csv}"

# Study metadata and paths selected by population / tissue source.
STUDY_IDS=(
    "QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069"
    "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372"
)
CELL_TYPES=(
    "Monocytes" "CD4+T_cells" "CD8+T_cells" "CD4+T_cells" "Monocytes"
    "B_cells" "Monocytes" "NK_cells" "CD4+T_cells" "CD8+T_cells"
)
GTEX_SAMPLE_SIZE=670
EQTLGEN_SAMPLE_SIZE=30000  # Already present in eQTLGen files; used for plotting.
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

# Runtime setup.
# Make the source-only figures package importable from the active Python
# environment. Re-sourcing config.sh must not duplicate the source path.
case ":${PYTHONPATH:-}:" in
    *":${SRC_DIR}:"*) ;;
    *) export PYTHONPATH="${SRC_DIR}${PYTHONPATH:+:${PYTHONPATH}}" ;;
esac
export PYTHONPATH
export MPLBACKEND="${MPLBACKEND:-Agg}"

activate_conda_env() {
    local env_name="$1"
    # shellcheck source=/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
}

mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}" "${TRACECB_FIGURE_DIR}"
