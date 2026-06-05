#!/usr/bin/env bash
set -euo pipefail

# Rebuttal simulation grids for showing when traceC and traceCB gain power.
#
# Usage:
#   bash src/simulation/others/run_power_gain.sh
#
# Defaults are true omega, h1sq=0.1, h2sq=0.1, gc=0.7, n1=100, n2=400,
# nt=5000, nsnp=2000, propt=0.4, pcausal=0.005, nrep=100.
# This script only varies the parameters needed for the two rebuttal panels.
#
# Run from the repository root. The script activates conda env py312.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

POP1_GENO="${POP1_GENO:-data/simulation/EAS_n5000_chr22_loci29.npy}"
POP2_GENO="${POP2_GENO:-data/simulation/EUR_n20000_chr22_loci29.npy}"
OUT_DIR="${OUT_DIR:-bench/result_power_gain}"
NREP="${NREP:-100}"
NSNP="${NSNP:-2000}"
POWER_GAIN_SEED="${POWER_GAIN_SEED:-20260605}"

TRACEC_RUNNAME="${TRACEC_RUNNAME:-power_gain_tracec_n2_rho}"
TRACECB_RUNNAME="${TRACECB_RUNNAME:-power_gain_tracecb_propt_nt}"
GAIN_OUTPUT_PREFIX="${GAIN_OUTPUT_PREFIX:-${OUT_DIR}/img/power_gain}"

mkdir -p "${OUT_DIR}/img"

echo "Running traceC gain grid: x=N2, group=rho"
python3 src/simulation/simulation.py \
    --pop1_geno "${POP1_GENO}" \
    --pop2_geno "${POP2_GENO}" \
    --runname "${TRACEC_RUNNAME}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.01 0.3 0.6 0.9 \
    --n1 100 \
    --n2 100 200 400 800 \
    --nt 5000 \
    --nsnp "${NSNP}" \
    --propt 0.4 \
    --pcausal 0.005 \
    --out_dir "${OUT_DIR}" \
    --nrep "${NREP}" \
    --seed "${POWER_GAIN_SEED}"

python3 src/simulation/visual_simulation.py \
    --metric power \
    --runname "${TRACEC_RUNNAME}" \
    --row h2sq \
    --col gc \
    --x n2 \
    --omega true \
    --base_path "${OUT_DIR}" \
    --ymin 0 \
    --ymax 1

echo "Running traceCB gain grid: x=propt, group=Nt"
python3 src/simulation/simulation.py \
    --pop1_geno "${POP1_GENO}" \
    --pop2_geno "${POP2_GENO}" \
    --runname "${TRACECB_RUNNAME}" \
    --h1sq 0.1 \
    --h2sq 0.1 \
    --gc 0.7 \
    --n1 100 \
    --n2 400 \
    --nt 1000 5000 10000 \
    --nsnp "${NSNP}" \
    --propt 0.01 0.3 0.6 0.9 \
    --pcausal 0.005 \
    --out_dir "${OUT_DIR}" \
    --nrep "${NREP}" \
    --seed "${POWER_GAIN_SEED}"

python3 src/simulation/visual_simulation.py \
    --metric power \
    --runname "${TRACECB_RUNNAME}" \
    --row h2sq \
    --col nt \
    --x propt \
    --omega true \
    --base_path "${OUT_DIR}" \
    --ymin 0 \
    --ymax 1

python3 src/simulation/others/visual_power_gain.py \
    --base_path "${OUT_DIR}" \
    --tracec_runname "${TRACEC_RUNNAME}" \
    --tracecb_runname "${TRACECB_RUNNAME}" \
    --output_prefix "${GAIN_OUTPUT_PREFIX}"

echo "Done. Gain figure: ${GAIN_OUTPUT_PREFIX}.pdf"
