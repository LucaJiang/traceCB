# S-LDSC eGene-Interval SNP Enrichment Pipeline

This folder contains the reproducible pipeline for the eGene-interval SNP
S-LDSC screen saved under:

`/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea`

## Annotation Definition

The LDSC annotation unit is a SNP, but the selected regions are defined by
eGenes.

1. eGenes are called from each study's `GMM/chr*/summary.csv`:
   - `original`: `TAR_SeSNP > 0`
   - `traceC_increment`: `TAR_CeSNP > 0` and not `original`
   - `traceCB_increment`: `TAR_TeSNP > 0` and not `original`
   - `original_overall`, `traceC_overall`, `traceCB_overall`: union of the
     corresponding eGenes across all 10 studies
2. Each eGene interval is defined as the min/max tested SNP position in the
   corresponding study's `INFO/chr*.csv`.
3. A 1000G EAS reference SNP is annotated as 1 if its BIM position overlaps at
   least one selected eGene interval.

Each annotation is run separately as `baselineLD + one custom annotation`.
The three annotations are not fit as one three-column custom joint model.

## Scripts

- `01_prepare_annotations.py`: creates single-column `.annot.gz` files and
  metadata manifests.
- `02_compute_ldscores.sh`: computes custom LD scores with `ldsc.py --l2`.
- `03_run_h2.sh`: runs partitioned heritability with `ldsc.py --h2`.
- `04_aggregate_results.py`: parses `.results` files into summary CSV files.
- `05_visualize_results.py`: writes heatmaps and overall 95% CI plots.
- `ora_utils.py`: shared GSEApy ORA utilities, gene-group definitions, GMT
  filtering, and publication-style dotplot helpers.
- `06_ora_significant_pathways.py`: ORA for heritability- and
  correlation-significant gene sets.
- `07_ora_two_tier_enrichment.py`: two-tier ORA for original, traceC
  incremental, traceCB incremental, and full original/traceC/traceCB eGenes.
- `run_ora_pipeline.sh`: runs both ORA scripts and overwrites the ORA output
  folders.
- `run_all_sldsc_gsea.sh`: runs all steps in order.

## Command

```bash
RESULT_DIR=/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea \
LDSC_L2_MAX_JOBS=72 \
LDSC_H2_MAX_JOBS=72 \
OVERWRITE=1 \
bash /home/wjiang49/traceCB/src/enrichment/run_all_sldsc_gsea.sh
```

Run the ORA analysis:

```bash
RESULT_ROOT=/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/sldsc_gsea \
GMT_DIR=/home/wjiang49/group/wjiang49/data/gsea_gmt \
WORKERS=24 \
bash /home/wjiang49/traceCB/src/enrichment/run_ora_pipeline.sh
```

`run_ora_pipeline.sh` also prepares the required GMT files before ORA. Existing
GMT files are reused by default. Useful switches:

- `PREPARE_GMT_ONLY=1`: download/build GMT files and stop before ORA.
- `SKIP_GMT_PREP=1`: skip the GMT preparation check.
- `FORCE_GMT=1`: re-download/rebuild GMT files.

## Outputs

- `metadata/annotation_manifest.tsv`: annotation IDs, eGene counts, interval
  counts, annotated SNP counts, and LD score prefixes.
- `metadata/trait_manifest.tsv`: six GWAS traits used in this screen.
- `annotations/ldscores/`: `.annot.gz` and `.l2.ldscore.gz` files.
- `results/raw/`: one `.results` file per trait and custom annotation.
- `results/summary/master_sldsc_gsea_results.csv`: all 198 custom annotation
  results.
- `visualization/`: publication-style PDF/PNG figures.
- `ora/significant_pathways/`: local-GMT ORA for heritability- and
  correlation-significant genes.
- `ora/two_tier_enrichment/`: two-tier local-GMT ORA for incremental and full
  original/traceC/traceCB eGene groups, including tier pathway-count summaries.
