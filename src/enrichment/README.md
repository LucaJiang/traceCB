# S-LDSC eGene-Interval Enrichment Pipeline

This directory contains the reproducible EUR-reference rerun of the
eGene-interval stratified LD score regression (S-LDSC) analysis. Generated
reference files, annotations, logs, and results are written under
`output/sldsc_gsea_eur_release_matched/` by default and are intentionally not
tracked by Git.

## Scientific definition

The LDSC annotation unit is a SNP, while the biological regions are defined by
EAS/BBJ-derived eGenes:

1. eGenes are called from each study's `GMM/chr*/summary.csv`.
   - `original`: `TAR_SeSNP > 0`
   - `traceC_increment`: `TAR_CeSNP > 0` and not `original`
   - `traceCB_increment`: `TAR_TeSNP > 0` and not `original`
   - Overall annotations are unions across the 10 studies.
2. Each eGene interval is the minimum and maximum tested SNP position in the
   corresponding study's `INFO/chr*.csv` (GRCh37/hg19).
3. The eGene membership and interval definitions are held fixed when changing
   the S-LDSC reference ancestry. Their complete serialized definitions are
   locked by SHA-256 values in
   `config/annotation_definition_checksums.tsv`; a rerun stops if they drift.
4. A 1000 Genomes Project Phase 3 EUR SNP is annotated when its GRCh37 position
   overlaps at least one selected interval.

Each of the 33 custom annotations is tested separately as `EUR baseline-LD
v2.2 + one custom annotation`, using `--overlap-annot`. The analysis therefore
contains 198 models (33 annotations x 6 traits); the three custom annotation
families are not fitted jointly. `Enrichment_FDR` is the Benjamini-Hochberg
adjustment across all 198 enrichment P values.

## Reference provenance

All S-LDSC reference components come from *S-LDSC reference files*, version 4
(Zenodo DOI: [10.5281/zenodo.10515792](https://doi.org/10.5281/zenodo.10515792)).
The preparation script verifies the MD5 values published by Zenodo before
extracting any archive and records SHA-256 values for the exact local files.

| File | Published MD5 |
|---|---|
| `1000G_Phase3_baselineLD_v2.2_ldscores.tgz` | `b261e0caf06a003e7522938e01b3d349` |
| `1000G_Phase3_plinkfiles.tgz` | `a7773ab485827b533cb300c76356d76b` |
| `1000G_Phase3_frq.tgz` | `ac29686ffd5b6378789857a522ebca77` |
| `1000G_Phase3_weights_hm3_no_MHC.tgz` | `a98ac0f089ee285177544a3e6e721ca3` |
| `hm3_no_MHC.list.txt` | `65a34c68833eb4a764d0707b5505b508` |

The compatible regression-SNP list is the exact intersection of the published
HapMap3 non-MHC list, EUR PLINK panel, baseline-LD rows, and EUR regression
weights. Baseline-LD rows are then copied into that exact SNP order; annotation
values, LD-score values, M counts, and the 97 baseline categories are not
changed.

Reference archives are downloaded at runtime and are not redistributed by this
repository. Users remain responsible for the source resources' access and use
terms and for the Pan-UK Biobank data-use requirements.

## Configuration and command

The tracked trait configuration is `config/traits.tsv`. Its summary-statistic
paths are relative to `GWAS_ROOT`. The pipeline requires:

- `STUDY_DIR`: directory containing the 10 `QTD*` study folders
- `GWAS_ROOT`: directory containing the Pan-UK Biobank sumstats subdirectories
- `LDSC_DIR`: checkout containing `ldsc.py`
- an LDSC-compatible environment (the stage scripts activate a Conda
  environment named `ldsc` unless `SKIP_CONDA_ACTIVATE=1`)

From the repository root:

```bash
STUDY_DIR=/path/to/EAS_eQTLGen \
GWAS_ROOT=/path/to/pan_ukb \
LDSC_DIR=/path/to/ldsc \
PYTHON_BIN=python \
RESULT_DIR="$PWD/output/sldsc_gsea_eur_release_matched" \
LDSC_L2_MAX_JOBS=72 \
LDSC_H2_MAX_JOBS=72 \
OVERWRITE=1 \
bash src/enrichment/run_all_sldsc_gsea.sh
```

Use a new, explicitly EUR-labelled result directory for a clean rerun. Set
`OVERWRITE=0` to reuse complete outputs already present in that directory.

## Pipeline stages

- `00_prepare_eur_reference.sh`: downloads, verifies, and stages the EUR
  reference stack.
- `00_build_compatible_hm3_list.py`: builds the exact common HapMap3 non-MHC
  SNP list.
- `01_prepare_annotations.py`: creates the 33 thin EUR-SNP annotations and
  auditable gene/interval manifests.
- `00_validate_eur_stack.py`: checks ancestry/build provenance, full annotation
  definition hashes, identifiers, coordinates, and reference compatibility.
- `02_compute_ldscores.sh`: computes 33 x 22 = 726 custom LD-score files.
- `00_filter_baseline_to_regression_snps.py`: makes the row-restricted
  baseline-LD copy required by the exact common SNP order.
- `00_validate_custom_ldscores.py`: checks all 726 files and every SNP row.
- `03_run_h2.sh`: runs the 198 univariate S-LDSC models.
- `04_aggregate_results.py`: validates and aggregates the result tables and
  applies global BH correction.
- `05_visualize_results.py`: writes the enrichment heatmap. Its fixed 0.9--1.3
  color scale uses colorbar extensions for clipped cells; printed cell values
  retain the estimates.
- `06_write_manuscript_sections.py`: writes a data-driven manuscript draft
  after checking the validation report and FDR values. It never modifies the
  analysis tables.

## Auditable outputs

- `reference/reference_manifest.tsv`: DOI source URLs, published MD5 values,
  and observed SHA-256 values.
- `metadata/annotation_manifest.tsv`: annotation labels, counts, definition
  hashes, ancestry/build, and prefixes.
- `metadata/annotation_egene_sets.tsv.gz` and
  `metadata/annotation_cis_intervals.tsv.gz`: exact biological definitions.
- `metadata/trait_manifest.tsv`: resolved GWAS paths and provenance.
- `metadata/validation_summary.json` and
  `metadata/custom_ldscore_validation.json`: reference and custom-LD audits.
- `metadata/annotation_prefixes.tsv` and `metadata/h2_jobs.tsv`: headered,
  deterministic execution manifests.
- `results/summary/master_sldsc_gsea_results.csv`: all 198 model results.
- `visualization/`: PDF/PNG figure and figure-specific notes.
- `manuscript/sldsc_eur_reference_revision.md`: generated text draft. Its data
  availability section deliberately contains a pre-publication archival-deposit
  reminder rather than claiming that a local filesystem path is public.
- `logs/run_all_sldsc_gsea.*.log`: resolved configuration, traceCB/LDSC Git
  revisions, traceCB dirty-worktree flag, Python version, stage output, and
  timestamps. Publication reruns should use `TRACECB_GIT_DIRTY=0`.

The ORA workflow remains separate in `run_ora_pipeline.sh`; it consumes the
S-LDSC results but is not part of the EUR reference-stack rerun above.
