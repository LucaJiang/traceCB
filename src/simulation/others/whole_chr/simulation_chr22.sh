#!/usr/bin/env bash
set -euo pipefail

# Curated whole-chromosome chr22 mixture-architecture simulation and figures.
#
# Usage:
#   bash src/simulation/others/whole_chr/simulation_chr22.sh
#
# Run from the repository root. The script activates conda env py312.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "${SCRIPT_DIR}/../../../.." && pwd)}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

OUT_DIR="${OUT_DIR:-bench/result/chr22_eqtl_mixture}"
IMG_DIR="${IMG_DIR:-${REPO_ROOT}/bench/result/img}"
NREP="${NREP:-1}"
MAX_GENES="${MAX_GENES:-0}"
MAX_SNPS_PER_GENE="${MAX_SNPS_PER_GENE:-0}"
ERROR_UNIT="${ERROR_UNIT:-gene}"

python3 src/simulation/others/whole_chr/simulation_chr22.py \
    --run_prefix mixture_chr22 \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.1 0.2 \
    --gc 0.7 \
    --n1 100 \
    --n2 100 400 \
    --nt 5000 \
    --propt 0.01 0.3 0.6 0.9 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega

python3 src/simulation/others/whole_chr/visual_chr22_eqtl_simulation.py \
    --base_path "${OUT_DIR}" \
    --out_dir "${IMG_DIR}" \
    --runname mixture_chr22 \
    --metric all \
    --error_unit "${ERROR_UNIT}" \
    --weighting gene_mean
