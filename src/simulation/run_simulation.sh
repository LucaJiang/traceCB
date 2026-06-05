#!/usr/bin/env bash
set -euo pipefail

# Curated small-window simulation.py grid for the traceCB paper.
#
# Usage:
#   bash src/simulation/run_simulation.sh
#
# Run from the repository root. The script activates conda env py312.
# Each simulation.py command is followed by its visual_simulation.py command so
# the paper figure is regenerated with the same parameters.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

POP1_GENO="${POP1_GENO:-data/simulation/EAS_n5000_chr22_loci29.npy}"
POP2_GENO="${POP2_GENO:-data/simulation/EUR_n20000_chr22_loci29.npy}"
OUT_DIR="${OUT_DIR:-bench/result}"
NREP="${NREP:-100}"
NSNP="${NSNP:-2000}"
SIMULATION_SEED="${SIMULATION_SEED:-20260525}"

## estimated Omega in fig 2
python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 1000 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric power --runname nt_n2_propt  --ymax 0.38 --ymin 0.16 --base_path "${OUT_DIR}" --omega false

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.42 --ymin 0.12 --base_path "${OUT_DIR}" --omega false

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.32 --ymin 0.05 --base_path "${OUT_DIR}" --omega false

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.48 --base_path "${OUT_DIR}" --omega false

# ---------------------------

## true Omega in supplementary

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 1000 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}"
python3 src/simulation/visual_simulation.py --metric power --runname nt_n2_propt --ymax 0.88 --base_path "${OUT_DIR}" --omega true

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}"
python3 src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.88 --base_path "${OUT_DIR}" --omega true

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}"
python3 src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.88 --base_path "${OUT_DIR}" --omega true

python3 src/simulation/simulation.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${SIMULATION_SEED}"
python3 src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.48 --base_path "${OUT_DIR}" --omega true
