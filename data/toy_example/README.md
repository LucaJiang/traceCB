# traceCB toy example

This directory contains the processed one-gene example used by
`docs/tutorial/run_traceCB.ipynb` and its Colab counterpart. The target gene is
`ENSG00000025708`.

- `eas_summary_statistics.csv`, `eur_summary_statistics.csv`, and
  `tissue_summary_statistics.csv` contain aligned summary statistics.
- `eas_ld.csv`, `eur_ld.csv`, and `cross_ld.csv` contain the corresponding
  within- and cross-population LD scores.
- `celltype_proportion.csv` contains the cell-type proportion used by the
  tissue-aware model.

These are processed tutorial inputs, not the full study datasets. Upstream data
sources and access conditions are described in `docs/pipeline.md`. Verify the
bundled inputs before a reproducibility run with:

```bash
sha256sum --check data/toy_example/SHA256SUMS
```
