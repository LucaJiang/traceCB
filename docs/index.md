# Welcome to traceCB

**traceCB** is a Python package for **Tra**ns-ancestry **c**ell-type-specific **e**QTL effects mapping by integrating scRNA-seq and **B**ulk data ("Method of Moments").

![traceCB Workflow](img/traceCB.png)

## Overview

traceCB enables researchers to:
- Integrate single-cell RNA-seq data with bulk GWAS summary statistics.
- Estimate cell-type-specific eQTL effects across different populations.
- Perform trans-ancestry analysis to improve resolution and power.

## Key Features

- **Integration**: Seamlessly combines single-cell precision with bulk data scale.
- **Trans-ancestry**: Leverages diversity to find robust signals.
- **Efficiency**: Optimized using `numba` for high-performance computing.

## Documentation

Explore our documentation to learn how to use traceCB:

- **[Pipeline Workflow](pipeline.md)**: Detailed steps for data preprocessing and analysis modules.
- **[API Reference](api_reference.md)**: Comprehensive documentation for core functions (`GMM`, `LDSC`).
- **[Tutorial](tutorial/run_traceCB.ipynb)**: Step-by-step Jupyter notebook example.

## Installation

### Prerequisites
- Python >= 3.8

### Quick Start

Clone the repository and install the package using pip:

```bash
git clone https://github.com/lucajiang/traceCB.git
cd traceCB
pip install -e .
```

## Citation

If you use **traceCB** in your research, please cite our paper:

<!-- TODO -->

## Support

For any questions or issues, please contact [wx.jiang@my.cityu.edu.hk](mailto:wx.jiang@my.cityu.edu.hk) or open an issue on GitHub.
