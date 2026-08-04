#!/usr/bin/env bash
set -euo pipefail

# Split one BBJ cell-type archive into chromosome-level CSV files.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"

CELL_TYPE="${1:-${CELL_TYPE:-}}"
if [[ -z "${CELL_TYPE}" ]]; then
    echo "Usage: bash scripts/preprocess_bbj.sh <cell-type>" >&2
    echo "Cell types: B_cells, CD4+T_cells, CD8+T_cells, Monocytes, NK_cells" >&2
    exit 2
fi

BBJ_SOURCE_DIR="${BBJ_SOURCE_DIR:-${DATA_ROOT}/BBJ_eQTL}"
BBJ_OUTPUT_DIR="${BBJ_OUTPUT_DIR:-${BBJ_SOURCE_DIR}/by_celltype_chr}"
ARCHIVE="${BBJ_SOURCE_DIR}/eQTL_${CELL_TYPE}.tar.gz"
CELL_OUTPUT_DIR="${BBJ_OUTPUT_DIR}/${CELL_TYPE}"
mkdir -p "${CELL_OUTPUT_DIR}"

if [[ ! -f "${ARCHIVE}" ]]; then
    echo "Missing BBJ archive: ${ARCHIVE}" >&2
    exit 2
fi

for chromosome in {1..22}; do
    output="${CELL_OUTPUT_DIR}/chr${chromosome}.csv"
    printf 'CHR,RSID,POS,A2,A1,GENE,BETA,Z,PVAL\n' > "${output}"
    tar -xOzf "${ARCHIVE}" \
        "${CELL_TYPE}/chr${chromosome}_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz" |
        gzip -dc |
        awk -v OFS=',' -v chromosome="${chromosome}" 'NR > 1 {
            sub(/\.[0-9]+$/, "", $5)
            if (substr($1, 1, 2) == "rs") {
                print chromosome, $1, $2, $3, $4, $5, $6, $7, $8
            }
        }' >> "${output}"
    echo "Prepared ${CELL_TYPE}, chr${chromosome}"
done
