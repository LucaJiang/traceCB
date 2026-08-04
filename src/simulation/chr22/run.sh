#!/usr/bin/env bash
set -euo pipefail

# Curated whole-chromosome chr22 mixture-architecture simulation and figures.
#
# Usage:
#   bash src/simulation/chr22/run.sh
#
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
setup_simulation_env

OUT_DIR="${OUT_DIR:-bench/result/chr22_eqtl_mixture}"
IMG_DIR="${IMG_DIR:-${REPO_ROOT}/bench/result/img}"
NREP="${NREP:-1}"
MAX_GENES="${MAX_GENES:-0}"
MAX_SNPS_PER_GENE="${MAX_SNPS_PER_GENE:-0}"
ERROR_UNIT="${ERROR_UNIT:-gene}"
RUN_VISUALS="${RUN_VISUALS:-1}"

"${PYTHON}" src/simulation/chr22/simulate.py \
    --run_prefix baseline \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.9 \
    --n1 100 \
    --n2 400 \
    --nt 5000 \
    --propt 0.9 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega \
    --re2

"${PYTHON}" src/simulation/chr22/simulate.py \
    --run_prefix nt_20000 \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.9 \
    --n1 100 \
    --n2 400 \
    --nt 20000 \
    --propt 0.9 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega \
    --re2

"${PYTHON}" src/simulation/chr22/simulate.py \
    --run_prefix n2_1000 \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.9 \
    --n1 100 \
    --n2 1000 \
    --nt 5000 \
    --propt 0.9 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega \
    --re2

"${PYTHON}" src/simulation/chr22/simulate.py \
    --run_prefix h2sq_0.2 \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.2 \
    --gc 0.9 \
    --n1 100 \
    --n2 400 \
    --nt 5000 \
    --propt 0.9 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega \
    --re2

"${PYTHON}" src/simulation/chr22/simulate.py \
    --run_prefix propt_0.01 \
    --out_dir "${OUT_DIR}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.9 \
    --n1 100 \
    --n2 400 \
    --nt 5000 \
    --propt 0.01 \
    --pcausal 0.005 \
    --architecture_probs 0.25 0.25 0.25 0.25 \
    --nrep "${NREP}" \
    --max_snps_per_gene "${MAX_SNPS_PER_GENE}" \
    --max_genes "${MAX_GENES}" \
    --estimate_omega \
    --re2

if [[ "${RUN_VISUALS}" == "1" ]]; then
    "${PYTHON}" src/simulation/chr22/plot_results.py \
        --base_path "${OUT_DIR}" \
        --out_dir "${IMG_DIR}" \
        --summary_out_dir "${OUT_DIR}" \
        --runname mixture_chr22 \
        --metric all \
        --error_unit "${ERROR_UNIT}" \
        --weighting gene_mean \
        --re2 \
        --run_prefix_order baseline propt_0.01 h2sq_0.2 n2_1000 nt_20000 \
        --power_shared_ylim 0.22 0.42 \
        --power_pop1_specific_ylim 0.22 0.42 \
        --alpha_null_ylim 0 0.58 \
        --alpha_pop2_specific_ylim 0 0.65
fi
