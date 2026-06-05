#!/usr/bin/env bash
set -euo pipefail

# Supplementary robustness grids using simulation_robustness.py.
#
# This script keeps add-on architectures separate from the base
# src/simulation/simulation.py driver:
#   1. inaccurate mean cell-type proportion supplied to GMM tissue
#   2. shared-causal-SNP robustness checks via --causal_overlap
#   3. segmented-null A/B checks via --causal_partition_mode pop2_a_shared_b
#
# Usage:
#   bash src/simulation/others/run_robustness.sh
#
# Run from the repository root. The script activates conda env py312.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py312
export PYTHONUNBUFFERED=1

POP1_GENO="${POP1_GENO:-data/simulation/EAS_n5000_chr22_loci29.npy}"
POP2_GENO="${POP2_GENO:-data/simulation/EUR_n20000_chr22_loci29.npy}"
OUT_DIR="${OUT_DIR:-bench/result}"
NREP="${NREP:-100}"
NSNP="${NSNP:-2000}"
GMM_PROPT_SEED="${GMM_PROPT_SEED:-20260525}"

# true omega with inaccurate mean cell-type proportion supplied to GMM tissue
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode exact --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}"
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode underestimate overestimate --gmm_propt_mode_scale 0.1 0.2 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}"
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode normal --gmm_propt_normal_var 0.1 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}"
python3 src/simulation/visual_simulation.py --metric alpha --runname alpha_pcausal_propt_gmmproptmode --ymax 0.4 --base_path "${OUT_DIR}" --omega true

# estimated omega with inaccurate mean cell-type proportion supplied to GMM tissue
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode exact --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode underestimate overestimate --gmm_propt_mode_scale 0.1 0.2 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_pcausal_propt_gmmproptmode --h1sq 0.000000000001 --h2sq 0.1 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --gmm_propt_mode normal --gmm_propt_normal_var 0.1 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha --runname alpha_pcausal_propt_gmmproptmode --ymax 0.4 --base_path "${OUT_DIR}" --omega false

# estimated omega, correlation controlled by shared causal SNP proportion
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname h2sq_causaloverlap_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0 --causal_overlap 0 0.5 0.9 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric power --runname h2sq_causaloverlap_propt --ymax 0.42 --ymin 0.12 --base_path "${OUT_DIR}"

# estimated omega, causal SNPs constrained to max pairwise abs correlation 0.6
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_h2sq_causaloverlap_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --causal_overlap 0.4 0.8 --causal_max_abs_cor 0.6 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_causaloverlap_propt --ymax 0.4 --base_path "${OUT_DIR}"

python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname power_h2sq_causaloverlap_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0 --causal_overlap 0.4 0.8 --causal_max_abs_cor 0.6 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric power --runname power_h2sq_causaloverlap_propt --ymax 0.42 --ymin 0.0 --base_path "${OUT_DIR}"

# segmented A/B setting: A is pop1-null with pop2-only cis-SNPs; B is shared non-null
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_a_h2sq_nullregionprop_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0 --causal_partition_mode pop2_a_shared_b --null_region_prop 0.8 --causal_max_abs_cor 0.6 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha_a --runname alpha_a_h2sq_nullregionprop_propt --ymax 0.4 --ymin 0.0 --base_path "${OUT_DIR}"

# segmented A/B setting without causal SNP max-correlation constraint
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_a_h2sq_nullregionprop_propt_nocormax --h1sq 0.1 --h2sq 0.1 0.2 --gc 0 --causal_partition_mode pop2_a_shared_b --null_region_prop 0.2 0.5 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha_a --runname alpha_a_h2sq_nullregionprop_propt_nocormax --ymax 0.4 --ymin 0.0 --base_path "${OUT_DIR}"

# segmented A/B setting without causal SNP max-correlation constraint; B effect correlation controlled by gc
python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_a_pA02_h2sq_gc_propt_nocormax --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.3 0.5 0.9 --causal_partition_mode pop2_a_shared_b --null_region_prop 0.2 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha_a --runname alpha_a_pA02_h2sq_gc_propt_nocormax --ymax 0.4 --ymin 0.0 --base_path "${OUT_DIR}"
python3 src/simulation/visual_simulation.py --metric power --runname alpha_a_pA02_h2sq_gc_propt_nocormax --ymax 0.6 --ymin 0.0 --base_path "${OUT_DIR}"

python3 src/simulation/others/simulation_robustness.py --pop1_geno "${POP1_GENO}" --pop2_geno "${POP2_GENO}" --runname alpha_a_pA05_h2sq_gc_propt_nocormax --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.3 0.5 0.9 --causal_partition_mode pop2_a_shared_b --null_region_prop 0.5 --n1 100 --n2 400 --nt 5000 --nsnp "${NSNP}" --propt 0.01 0.3 0.6 --pcausal 0.005 --out_dir "${OUT_DIR}" --nrep "${NREP}" --seed "${GMM_PROPT_SEED}" --estimate_omega
python3 src/simulation/visual_simulation.py --metric alpha_a --runname alpha_a_pA05_h2sq_gc_propt_nocormax --ymax 0.4 --ymin 0.0 --base_path "${OUT_DIR}"
python3 src/simulation/visual_simulation.py --metric power --runname alpha_a_pA05_h2sq_gc_propt_nocormax --ymax 0.6 --ymin 0.0 --base_path "${OUT_DIR}"
