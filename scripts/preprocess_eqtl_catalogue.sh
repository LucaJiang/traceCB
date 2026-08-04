#!/usr/bin/env bash
set -euo pipefail

# Format eQTL Catalogue downloads as chromosome-level traceCB inputs.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${SCRIPT_DIR}/config.sh"
activate_conda_env "${PYTHON_ENV}"

python "${SRC_DIR}/preprocess/format_eqtl_catalogue.py" \
    --input-dir "${EQTL_CATALOGUE_SOURCE_DIR:-${DATA_ROOT}/eQTLCatalogue}" \
    --output-dir "${EQTL_CATALOGUE_DIR}"
