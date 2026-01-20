# traceCB: Trans-ancestry cell-type-specific eQTLs mapping by integrating scRNA-seq and bulk data

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains code for the **traceCB** paper, featuring the main algorithm and a complete pipeline for trans-ancestry cell-type-specific eQTL mapping.

![traceCB_workflow](docs/img/traceCB.png)

## Repository Structure

- `src/traceCB`: The main source code for the Python package.
- `src/coloc`: Scripts for colocalization analysis.
- `src/visual`: Visualization scripts for GMM results.
- `shell`: Shell scripts for running the pipeline steps (preprocessing, LDSC, GMM, etc.).
- `docs`: Documentation and tutorials.
- `data`: Folder for storing input/output data (see `docs/pipeline.md` for structure).

## Installation

### Prerequisites
- Python >= 3.8
- `numba`, `pyarrow`, `numpy`, `pandas`, `scipy`, 

### Install from source

Clone the repository and install the package using pip:

```bash
git clone https://github.com/lucajiang/traceCB.git
cd traceCB
pip install -e .
```

Installation in editable mode (`-e`) allows you to import the `traceCB` module in your scripts while keeping the ability to modify the source code if needed.

## Tutorial

A step-by-step tutorial notebook is provided at `docs/tutorial/run_traceCB.ipynb`. This tutorial guides you through running the traceCB algorithm on a single gene example.

It is highly recommended to run this tutorial first to understand the input data format and model outputs.

## Usage Pipeline

For full-scale analysis, we provide a structured shell-script pipeline. Detailed preprocessing steps are described in [Pipeline Documentation](docs/pipeline.md).

### 1. Configuration
Modify `shell/setting.sh` to specify your paths and parameters according to your environment.

### 2. Run Pipeline Steps
The analysis is divided into sequential modules:

```bash
# 1. Merge GWAS summary statistics
source shell/run_merge.sh

# 2. Calculate LD scores (LDSC)
source shell/run_ld.sh

# 3. Run Gaussian Mixture Model (GMM)
source shell/run_gmm.sh

# 4. Colocalization Analysis (Optional)
source shell/run_coloc.sh
```

### 3. Visualization
Scripts for visualization are provided in `src/visual/`.

## Citation

If you use **traceCB** in your research, please cite our paper:

> *Citation pending...*

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For any questions or issues, please contact [wx.jiang@my.cityu.edu.hk](mailto:wx.jiang@my.cityu.edu.hk).
