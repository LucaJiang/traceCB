# traceCB

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE)
[![Open Tutorial In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lucajiang/traceCB/blob/master/docs/tutorial/run_traceCB_colab.ipynb)

traceCB maps trans-ancestry cell-type-specific eQTL effects by integrating
single-cell and bulk-tissue summary statistics. This repository contains the
Python package, full-data workflows, simulations, and manuscript analyses.

![traceCB workflow](docs/img/traceCB.jpg)

## Repository layout

- `src/traceCB/`: installable traceCB package.
- `scripts/`: full-data preprocessing, LD-score, traceCB, and colocalization workflows.
- `src/simulation/`: main, supplementary, and chromosome 22 simulations.
- `src/enrichment/`: reproducible ancestry-matched S-LDSC and
  pathway-enrichment analyses; see
  [`src/enrichment/README.md`](src/enrichment/README.md).
- `src/figures/`: manuscript figure scripts and shared metadata.
- `src/preprocess/` and `src/coloc/`: workflow implementations used by `scripts/`.
- `tests/`: unit and CLI contract tests.
- `data/toy_example/`: small public inputs for the tutorial.

Large input data and generated results are intentionally excluded from Git.
By default, workflows read from `data/` and write to `results/`; both locations
can be overridden with environment variables documented in `scripts/config.sh`.

## Installation

Python 3.10 or newer is required. For repository-level reproduction, create the
validated Python 3.12 reference environment from the tracked specification:

```bash
git clone https://github.com/lucajiang/traceCB.git
cd traceCB
conda env create -f environment.yml
conda activate py312
```

For a lightweight library-only installation, use `pip install -e .` in any
supported Python environment.

When using a lightweight or custom environment, add the corresponding optional
dependencies for the local notebook and manuscript analyses:

```bash
pip install -e '.[tutorial]'            # local Jupyter tutorial
pip install -e '.[enrichment,figures]'  # enrichment and figure scripts
```

## Quick start

The tutorial notebooks are available at
[`docs/tutorial/run_traceCB.ipynb`](docs/tutorial/run_traceCB.ipynb) and on
[Google Colab](https://colab.research.google.com/github/lucajiang/traceCB/blob/master/docs/tutorial/run_traceCB_colab.ipynb).
They use the tracked files in `data/toy_example/` and demonstrate the model on a
single gene.

## Full-data workflow

External inputs such as population-specific eQTLs, tissue eQTLs, and 1000
Genomes reference panels are not redistributed here. Set their locations in the
environment or edit the repository-relative defaults in `scripts/config.sh`.
The BBJ cell-type eQTL data used for the EAS analysis are available from
[Human Database of Japan: hum0099-v1](https://humandbs.dbcls.jp/en/hum0099-v1).
See the [pipeline guide](https://lucajiang.github.io/traceCB/pipeline/) for all
data sources, expected input formats, and preprocessing commands.

```bash
# Optional examples
export TRACECB_DATA_ROOT=/path/to/input-data
export TRACECB_OUTPUT_ROOT=/path/to/results
export SLDXR_DIR=/path/to/s-ldxr
export PLINK_BIN=/path/to/plink

bash scripts/prepare_inputs.sh
bash scripts/run_ld_scores.sh
bash scripts/run_gmm.sh
bash scripts/run_colocalization.sh  # optional
```

See the [simulation guide](src/simulation/README.md) for manuscript simulation
entry points and the [enrichment guide](src/enrichment/README.md) for the S-LDSC
and pathway-enrichment analyses. Figure-specific dependencies and invocation
patterns are listed in the [figure guide](src/figures/README.md).

## Tests

```bash
pip install -e '.[test]'
conda run -n py312 pytest -q
```

## Citation

If you use traceCB, please cite the paper:

```txt
@article{jiang2026tracecb,
  title={{traceCB}: Trans-ancestry cell-type-specific {eQTLs} mapping by integrating {scRNA-seq} and bulk data},
  author={Jiang, Wenxin and Xiao, Jiashun and Cai, Mingxuan},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.06.20.733502},
  url={https://doi.org/10.64898/2026.06.20.733502},
  publisher={Cold Spring Harbor Laboratory}
}
```

## License

traceCB is distributed under the [GPL-3.0 license](LICENSE).
