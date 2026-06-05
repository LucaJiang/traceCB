#!/usr/bin/env bash
set -euo pipefail

# Run masked-omega comparisons for pop1 target.
#
# Defaults are production-style and align with src/simulation/run_simulation.sh.
# Default rho grid: 0.01, 0.3, 0.7.
# Uses true omega by default, matching simulation.py when --estimate_omega is absent.
# - type I error: h1sq=0.000000000001
#
# Run from the repository root. The script activates conda env py312.
# Override NREP/NSNP/OUT_DIR from the environment for quick checks, e.g.
# NREP=3 NSNP=500 bash src/simulation/others/run_masked_omega.sh
# ALPHA_GCS="0.01 0.3" bash src/simulation/others/run_masked_omega.sh

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

POP1_GENO="${POP1_GENO:-data/simulation/EAS_n5000_chr22_loci29.npy}"
POP2_GENO="${POP2_GENO:-data/simulation/EUR_n20000_chr22_loci29.npy}"
OUT_DIR="${OUT_DIR:-bench/result/masked_omega_compare}"
NREP="${NREP:-100}"
NSNP="${NSNP:-2000}"
MASKED_OMEGA_SEED="${MASKED_OMEGA_SEED:-20260525}"

ALPHA_GCS="${ALPHA_GCS:-0.01 0.3 0.7}"
ALPHA_RUNNAME_GC001="${ALPHA_RUNNAME_GC001:-alpha_masked_omega_trueomega_gc0.01_h2sq_pcausal_propt}"
ALPHA_RUNNAME_GC03="${ALPHA_RUNNAME_GC03:-alpha_masked_omega_trueomega_gc0.3_h2sq_pcausal_propt}"
ALPHA_RUNNAME_GC07="${ALPHA_RUNNAME_GC07:-alpha_masked_omega_trueomega_gc0.7_h2sq_pcausal_propt}"

for GC in ${ALPHA_GCS}; do
    if [[ "${GC}" == "0.01" ]]; then
        ALPHA_RUNNAME="${ALPHA_RUNNAME_GC001}"
    elif [[ "${GC}" == "0.3" ]]; then
        ALPHA_RUNNAME="${ALPHA_RUNNAME_GC03}"
    elif [[ "${GC}" == "0.7" ]]; then
        ALPHA_RUNNAME="${ALPHA_RUNNAME_GC07}"
    else
        echo "Unexpected alpha GC value: ${GC}" >&2
        exit 1
    fi
    python3 src/simulation/others/simulation_masked_omega.py \
        --pop1_geno "${POP1_GENO}" \
        --pop2_geno "${POP2_GENO}" \
        --runname "${ALPHA_RUNNAME}" \
        --h1sq 0.000000000001 \
        --h2sq 0.1 0.2 \
        --gc "${GC}" \
        --n1 100 \
        --n2 400 \
        --nt 5000 \
        --nsnp "${NSNP}" \
        --propt 0.01 0.2 0.4 0.6 0.8 \
        --pcausal 0.005 0.01 0.02 \
        --out_dir "${OUT_DIR}" \
        --nrep "${NREP}" \
        --seed "${MASKED_OMEGA_SEED}"
done

python3 src/simulation/others/visual_masked_omega.py \
    --base_path "${OUT_DIR}" \
    --runname "${ALPHA_RUNNAME_GC001}" "${ALPHA_RUNNAME_GC03}" "${ALPHA_RUNNAME_GC07}" \
    --save_prefix alpha_masked_omega_trueomega_gc_h2sq_pcausal_propt \
    --metric alpha \
    --x propt \
    --row gc \
    --col pcausal \
    --alpha_ymax 0.48
