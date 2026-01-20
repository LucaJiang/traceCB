# script to run all simulation experiments
python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 1000 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100 --estimate_omega
python src/simulation/visual_simulation.py --metric power --runname nt_n2_propt  --ymax 0.38 --ymin 0.16

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100 --estimate_omega
python src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.42 --ymin 0.12

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100 --estimate_omega
python src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.32

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100 --estimate_omega
python src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.48


# true Omega runs ----------------

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname nt_n2_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 100 --n2 100 200 400 --nt 1000 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100
python src/simulation/visual_simulation.py --metric power --runname nt_n2_propt --ymax 0.88

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname h2sq_gc_propt --h1sq 0.1 --h2sq 0.1 0.2 --gc 0.01 0.5 0.9 --n1 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 --out_dir bench/result --nrep 100
python src/simulation/visual_simulation.py --metric power --runname h2sq_gc_propt --ymax 0.88

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname n1_pcausal_propt --h1sq 0.1 --h2sq 0.1 --gc 0.7 --n1 50 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100
python src/simulation/visual_simulation.py --metric power --runname n1_pcausal_propt --ymax 0.88

python src/simulation/simulation.py --pop1_geno data/simulation/EAS_n5000_chr22_loci29.npy --pop2_geno data/simulation/EUR_n20000_chr22_loci29.npy --runname alpha_h2sq_pcausal_propt --h1sq 0.000000000001 --h2sq 0.1 0.2 --gc 0 --n1 100 --n2 400 --nt 5000 --nsnp 2000 --propt 0.01 0.2 0.4 0.6 0.8 --pcausal 0.005 0.01 0.02 --out_dir bench/result --nrep 100
python src/simulation/visual_simulation.py --metric alpha --runname alpha_h2sq_pcausal_propt --ymax 0.48
