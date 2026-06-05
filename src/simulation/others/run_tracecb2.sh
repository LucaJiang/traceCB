#!/usr/bin/env bash
set -euo pipefail

# Curated traceCB^2 simulation grid.
#
# Usage:
#   bash src/simulation/others/run_tracecb2.sh
#
# Run from the repository root. The script activates conda env py312.
# Simulations are run first; figures are generated only after all simulations finish.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

POP1_GENO="${POP1_GENO:-data/simulation/EAS_n5000_chr22_loci29.npy}"
POP2_GENO="${POP2_GENO:-data/simulation/EUR_n20000_chr22_loci29.npy}"
OUT_DIR="${OUT_DIR:-bench/result}"
NREP="${NREP:-100}"
NSNP="${NSNP:-2000}"
TRACECB2_SEED="${TRACECB2_SEED:-20260525}"

## estimated Omega, traceCB^2 power: tissue sample sizes in both populations
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname nt1_nt2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 400 --nt1 500 1000 2000 --nt2 1000 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}" --estimate_omega

## estimated Omega, traceCB^2 power: population 1 GWAS and tissue sample sizes
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n1_nt1_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 200 --n2 400 --nt1 500 1000 2000 --nt2 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}" --estimate_omega

## estimated Omega, traceCB^2 power: population 2 GWAS and population 1 tissue sample sizes
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n2_nt1_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 200 400 --nt1 500 1000 2000 --nt2 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}" --estimate_omega

## estimated Omega, traceCB^2 type I error
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}" --estimate_omega

# ---------------------------

## true Omega, traceCB^2 power: tissue sample sizes in both populations
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname nt1_nt2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 400 --nt1 500 1000 2000 --nt2 1000 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}"

## true Omega, traceCB^2 power: population 1 GWAS and tissue sample sizes
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n1_nt1_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 200 --n2 400 --nt1 500 1000 2000 --nt2 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}"

## true Omega, traceCB^2 power: population 2 GWAS and population 1 tissue sample sizes
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n2_nt1_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 200 400 --nt1 500 1000 2000 --nt2 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}"

## true Omega, traceCB^2 type I error
python3 src/simulation/others/simulation_tracecb2.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${TRACECB2_SEED}"

# ---------------------------
# Generate all figures after all simulation commands finish.

python3 src/simulation/others/visual_tracecb2.py --metric power --runname nt1_nt2_propt --omega false --ymax 0.39 --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric power --runname n1_nt1_propt --omega false --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric power --runname n2_nt1_propt --omega false --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric alpha --runname alpha_h2sq_pcausal_propt --omega false --ymin 0 --ymax 0.45 --base_path "${OUT_DIR}"

python3 src/simulation/others/visual_tracecb2.py --metric power --runname nt1_nt2_propt --omega true --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric power --runname n1_nt1_propt --omega true --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric power --runname n2_nt1_propt --omega true --base_path "${OUT_DIR}"
python3 src/simulation/others/visual_tracecb2.py --metric alpha --runname alpha_h2sq_pcausal_propt --omega true --ymin 0 --ymax 0.45 --base_path "${OUT_DIR}"
