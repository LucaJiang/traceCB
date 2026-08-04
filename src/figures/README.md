# Manuscript figure scripts

Run figure scripts from the repository root. Install the Python dependencies
declared for this analysis surface:

```bash
pip install -e '.[figures]'
PYTHONPATH=src python -m figures.egene_counts
```

Using `PYTHONPATH=src` keeps the manuscript scripts importable without placing
the `figures` package in the distributable traceCB wheel.

The R figures use these CRAN packages: `colorspace`, `cowplot`, `data.table`,
`dplyr`, `gggenes`, `ggplot2`, `ggtern`, `ggtext`, `gridExtra`, `jsonlite`,
`locuszoomr`, `patchwork`, and `readr`. The locus-zoom scripts additionally use
the Bioconductor packages `EnsDb.Hsapiens.v75`, `GenomicFeatures`, and
`rtracklayer`. The mashr comparison in `src/simulation/experiments/` requires
the R package `mashr`.

Shared labels and colors are stored in `metadata.json`. Common input and output
locations can be overridden without changing the scripts:

- `TRACECB_STUDY_DIR`: GMM study results.
- `TRACECB_FIGURE_DIR`: generated figure directory.
- `TRACECB_ONEK1K_FILE`: OneK1K replication input.
- `TRACECB_GTEX_LOOKUP`: GTEx variant lookup table.
- `TRACECB_GTEX_GENE_ANNOTATION`: GTEx gene annotation.
- `TRACECB_OASIS_DIR`: OASIS replication directory.

Other analysis-specific paths remain at the top of the corresponding script
and can be changed by the user when preparing the protected or external data.
