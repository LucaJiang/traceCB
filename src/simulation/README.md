# Simulation Code Guide

All commands below are intended to be run from the repository root with conda
env `py312`.

```bash
conda activate py312
```

The simulation outputs are written under `bench/result` by default. For smoke
tests or temporary reruns, set `OUT_DIR` or pass `--out_dir` to a temporary
folder and delete it after inspection.

## Main Small-Window Simulations

- Simulation code: `simulation.py`
- Shell entry point: `run_simulation.sh`
- Visualization code: `visual_simulation.py`
- Default output: `bench/result/<runname>/<parameter-setting>/simulation_<rep>.csv`
- Aggregated output: `bench/result/<runname>/result_df*.csv`
- Figure output: `bench/result/img/<runname>*.pdf`

`simulation.py` generates one cis-window of summary statistics for population 1
single-cell, population 2 single-cell, and population 2 bulk tissue. It then
evaluates original population 1 testing, traceC, traceCB, and meta-analysis
variants. The `--estimate_omega` flag estimates omega from summary statistics;
without it, the true simulated omega is used.
When a runname contains both estimated-omega and true-omega settings, pass
`--omega false` or `--omega true` to `visual_simulation.py` so the two modes are
aggregated separately.

Run the curated paper grid and its paired visualizations:

```bash
bash src/simulation/run_simulation.sh
```

For a temporary smoke run:

```bash
python3 src/simulation/simulation.py \
  --runname smoke_small_window \
  --h1sq 0.1 \
  --h2sq 0.1 \
  --gc 0.7 \
  --n1 20 \
  --n2 30 \
  --nt 40 \
  --nsnp 80 \
  --propt 0.2 \
  --pcausal 0.05 \
  --nrep 1 \
  --out_dir bench/result/_smoke_test

python3 src/simulation/visual_simulation.py \
  --base_path bench/result/_smoke_test \
  --runname smoke_small_window \
  --metric power
```

Delete `bench/result/_smoke_test` after confirming the run.

## Supplementary Robustness Simulations

- Simulation code: `others/simulation_robustness.py`
- Shell entry point: `others/run_robustness.sh`
- Visualization code: `visual_simulation.py`

This script keeps add-on experiments out of `simulation.py`, including
inaccurate GMM cell-type proportion, shared-causal-SNP overlap checks, and the
segmented-null `pop2_a_shared_b` architecture.

```bash
bash src/simulation/others/run.sh robustness
```

## Masked-Omega Comparisons

- Simulation code: `others/simulation_masked_omega.py`
- Shell entry point: `others/run_masked_omega.sh`
- Visualization code: `others/visual_masked_omega.py`
- Default output: `bench/result/masked_omega_compare`
- Figure output: `bench/result/masked_omega_compare/img`

This experiment compares original traceC/traceCB against two masked-input
variants: population 1 single-cell plus population 2 bulk, and population 2
single-cell plus population 2 bulk. The default shell script runs both power and
type I error grids and then generates the paired figures.

```bash
bash src/simulation/others/run_masked_omega.sh
```

Use environment variables such as `NREP=3`, `NSNP=500`, or `OUT_DIR=...` for
quick local checks.

## traceCB^2 Simulations

- Simulation code: `others/simulation_tracecb2.py`
- Shell entry point: `others/run_tracecb2.sh`
- Visualization code: `others/visual_tracecb2.py`
- Default output: `bench/result/<runname>`
- Figure output: `bench/result/img/*_tracecb2*.pdf`

This experiment adds a population 1 tissue panel and evaluates the four-source
traceCB^2 method against traceCB and traceC.

```bash
bash src/simulation/others/run_tracecb2.sh
```

## mashr Benchmarks

- Simulation code: `others/simulation_mashr.py`
- Shell entry point: `others/run_mashr.sh`
- Benchmark code: `others/run_mashr_benchmark.R`
- Visualization code: `others/visual_mashr_true_fsr_curve.py`
- Default output: `bench/result_mashr`

```bash
bash src/simulation/others/run_mashr.sh
```

## Power-Gain Panels

- Shell entry point: `others/run_power_gain.sh`
- Visualization code: `others/visual_power_gain.py`
- Default output directory: `bench/result_power_gain`
- Default output: `bench/result_power_gain/img/power_gain.pdf`

```bash
bash src/simulation/others/run_power_gain.sh
```

## Whole-Chromosome chr22 Simulations

- Preparation code: `others/whole_chr/prepare_chr22_annotation.py`
- Input harmonization: `others/whole_chr/prepare_chr22_eqtl_inputs.py`
- Simulation and grid runner: `others/whole_chr/simulation_chr22.py`
- Shell entry point: `others/whole_chr/simulation_chr22.sh`
- Visualization code: `others/whole_chr/visual_chr22_eqtl_simulation.py`
- Default output: `bench/result/chr22_eqtl_mixture`
- Figure output: `bench/result/chr22_eqtl_mixture/img`

These scripts use PLINK genotype inputs and gene-window SNP annotations to run
gene-level chr22 mixture-architecture simulations. See
`bench/result/chr22_eqtl_mixture/img/README.md` for the exact commands used to
generate the existing chr22 figures.

The whole-chromosome visualizer uses gene-level 95% confidence intervals by
default, treating genes as repeat units. Use `--error_unit replicate` only when
`NREP>1` and you want simulation-replicate uncertainty instead.

## Result Hygiene

- Keep long-running production results in named directories under
  `bench/result`.
- Put temporary runs under `bench/result/_smoke_test` or another ignored local
  directory and remove them after validation.
- Do not commit `__pycache__`, `.pyc`, `.DS_Store`, or ad hoc temporary figures.
- Each figure directory should include a `README.md` that records what the
  figures show, which script generated them, and the parameters used.
