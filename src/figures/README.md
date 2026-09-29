# Manuscript figure scripts

Use the [path configuration guide](../../docs/configuration.md) for a copyable
local configuration and checks that the figures read the intended complete run.

Run figure scripts from the repository root, using a Python environment with
the figure dependencies installed:

```bash
pip install -e '.[figures]'
```

The `figures` extra includes `openpyxl` for the CIMA cell-type annotation Excel
file. The runtime figure overlays one colored point per study and chromosome on
the chromosome mean bars; these points are saved elapsed-time records, not
repeated benchmark runs. Its right bars sum chromosome elapsed times by study.

Figures are saved as PDF without duplicate PNG exports.

The Sankey PDF export uses Plotly and Kaleido. Kaleido 1.x requires Plotly
6.1.1 or newer and a separate Chrome/Chromium installation. The reference
environment pins Plotly 6.5.1 with Kaleido 1.2.0. If a compatible browser is
not already available, install Chrome from the activated figure environment:

```bash
plotly_get_chrome
```

On machines without download access, install Chrome/Chromium separately and
set `BROWSER_PATH` to its executable if automatic detection fails. See the
[Plotly export documentation](https://plotly.com/python/static-image-export/)
for browser setup. HTML export does not verify that PDF export is available.

Before running, **review the shared figure paths in
[`scripts/config.sh`](../../scripts/config.sh) for your own filesystem**. The
input defaults in that file refer to the current machine; the full study data
and GENCODE annotation are not included in this repository.

| Configuration variable | What the user must specify |
| --- | --- |
| `TRACECB_STUDY_DIR` | Study results directory containing `QTD*/GMM/chr*/summary.csv`, for example `/path/to/traceCB/EAS_eQTLGen`. |
| `TRACECB_GTEX_GENE_ANNOTATION` | Path to the uncompressed GENCODE GTF file, for example `/path/to/GTEx/gencode.v26.GRCh38.genes.gtf`. |
| `TRACECB_FIGURE_DIR` | Writable figure output directory; defaults to `${TRACECB_OUTPUT_ROOT}/figures/${TARGET_POPULATION}_${TISSUE_SOURCE}`, normally `results/figures/EAS_eQTLGen` under the repository root. |

Edit the default values after `:-` in `scripts/config.sh`, preserving the
`${VARIABLE:-default}` syntax. Existing environment variables override those
defaults. Set `TARGET_POPULATION` and `TISSUE_SOURCE` to match the study results;
the configured defaults select `EAS_eQTLGen`.

After changing defaults, start a fresh terminal or unset the corresponding
previously exported variables before sourcing again. An IDE's Python/R process
must be started from the configured shell to inherit these settings.

Both `f3cor_density_*.pdf` and the `f3cor_*.pdf` box plots use `COR_X_ORI`, the
original correlation before clipping at ±1. This column is required in the
input summaries; older results containing only clipped `COR_X` must be
regenerated or replaced with results that retain the original correlations.
Finite correlations outside ±1 remain in the plots.

`TARGET_POPULATION` in `scripts/config.sh` also selects the box-plot groups:

- `EAS`: 10 groups, with cut points −0.8, −0.6, −0.4, −0.2, 0, 0.2, 0.4, 0.6, 0.8.
- `AFR`: 6 groups, with cut points −0.8, −0.5, 0, 0.5, 0.8.

Intervals include their right endpoint. The first group contains all values
≤−0.8 and the last contains all values >0.8, including values beyond ±1.
Figure output defaults keep each population/bulk-source configuration in its own
directory. When switching configurations in an existing shell, unset previously
exported `TRACECB_STUDY_DIR` and `TRACECB_FIGURE_DIR` before sourcing again, or
explicitly set them to the desired input and output directories.

Apply the configuration in the current shell, then run the case study:

```bash
source scripts/config.sh
python -m figures.case_study
python -m figures.combine_case_study
```

`source scripts/config.sh` exports the configured paths and population to
Python, sets `PYTHONPATH` and the default headless Matplotlib backend, and
creates the figure output directory. Run it once per terminal session; use
`source` so the settings remain available to subsequent Python commands.
Results are written to `${TRACECB_FIGURE_DIR}/single_study/`, which defaults to
`results/figures/EAS_eQTLGen/single_study/` for EAS + eQTLGen.

For EAS + GTEx, select `TARGET_POPULATION=EAS` and `TISSUE_SOURCE=GTEx` in
`scripts/config.sh`. The study directory then defaults to `EAS_GTEx` and the
figure output directory to `results/figures/EAS_GTEx`. Apply these settings
in a fresh shell so previously exported paths do not override the new defaults.

After the single-study figures finish, `figures.combine_case_study` combines each
study's correlation box plot, eGene UpSet plot, and annotated effective-sample-size
scatter plot. It writes `f3combined_QTD*.pdf` to
`${TRACECB_FIGURE_DIR}/single_study_combined/`, defaulting to
`results/figures/EAS_eQTLGen/single_study_combined/` for EAS + eQTLGen. Run it with
each cohort's configuration to generate separate EAS + eQTLGen, AFR + eQTLGen,
and EAS + GTEx sets.

To run another Python figure script with the same configuration:

```bash
python -m figures.egene_counts
python -m figures.oasis_pathway_egenes
python -m figures.cima_pathway_replication
Rscript src/figures/ternary_bcx.R
```

For `f3egene.pdf`, `TARGET_POPULATION` selects the original cell-type header
layout: EAS adds 10% to the y-axis range and places the colored bar 110 count
units below the new upper limit; AFR adds 24% and uses a 26-unit offset.
These settings apply to both eQTLGen and GTEx. The increment plot retains its
separate spacing based on its y-axis range.

Direct Python commands use the environment configured by
`source scripts/config.sh`. Its `PYTHONPATH` setting keeps the manuscript scripts
importable without placing the `figures` package in the distributable traceCB
wheel.

The R figures use these CRAN packages: `colorspace`, `cowplot`, `data.table`,
`dplyr`, `gggenes`, `ggplot2`, `ggtern`, `ggtext`, `gridExtra`, `jsonlite`,
`locuszoomr`, `patchwork`, and `readr`. The locus-zoom scripts additionally use
the Bioconductor packages `EnsDb.Hsapiens.v75`, `GenomicFeatures`, and
`rtracklayer`. The mashr comparison in `src/simulation/experiments/` requires
the R package `mashr`.

All data paths are configured and exported by `scripts/config.sh`. Python
scripts share `paths.py`; R scripts share `paths.R`. CLI path arguments, where
supported, override these defaults.

`TRACECB_DATA_ROOT` defaults to `/home/wjiang49/group/wjiang49/data` on this
machine. `TRACECB_STUDY_ROOT` defaults to `${TRACECB_DATA_ROOT}/traceCB`, containing
the existing results; set it to `${TRACECB_OUTPUT_ROOT}` to plot newly generated
pipeline results. Figure outputs follow `TRACECB_FIGURE_DIR`.
Pipeline tools are configured through `TRACECB_SOFTWARE_ROOT`, `PLINK_BIN`,
`PLINK2`, and `SLDXR_DIR` in the same file.

| Variable | Default location / purpose |
| --- | --- |
| `TRACECB_OASIS_DIR` | `${TRACECB_DATA_ROOT}/hum0197/eQTL_summary_statistics`; requires the `Mono`, `CD4T`, `CD8T`, `B`, and `NK` files ending in `_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz`. |
| `TRACECB_ONEK1K_FILE` | `${TRACECB_STUDY_ROOT}/onek1k_supp/onek1k_esnp.csv`. |
| `TRACECB_GTEX_LOOKUP` | Full GTEx v8 `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz` under `${TRACECB_DATA_ROOT}/GTEx`. Also used by GTEx preprocessing. |
| `TRACECB_CIMA_DIR` | `${TRACECB_DATA_ROOT}/CIMA`; base for the next two inputs. |
| `TRACECB_CIMA_LEAD_EQTL` | `${TRACECB_CIMA_DIR}/xQTL/CIMA_Lead_cis-xQTL.csv`. |
| `TRACECB_CIMA_CELL_TYPES` | `${TRACECB_CIMA_DIR}/Cell_Atlas/CIMA_Cell_Type_Level_and_Marker.xlsx`. |
| `TRACECB_REPLICATION_EGENES` / `TRACECB_REPLICATION_ESNPS` | `hum0343_eGene.csv` / `hum0343_eSNP.csv` under `${TRACECB_DATA_ROOT}/hum0343`. |
| `TRACECB_CELL_PROPORTIONS` | `${TRACECB_STUDY_ROOT}/cell_type_proportion/ind_celltype_proportion.csv`. |
| `TRACECB_AFR_STUDY_DIR` | `${TRACECB_STUDY_ROOT}/AFR_${TISSUE_SOURCE}`. |
| `TRACECB_AFR_FIGURE_DIR` | `${TRACECB_FIGURE_DIR}/afr`. |
| `TRACECB_ESNP_REPLICATION_DIR` | `${TRACECB_FIGURE_DIR}/esnp_replication`. |
| `TRACECB_TIMING_FILE` | Repository `tmp/timing/summary_timing.csv`. |
| `TRACECB_COLOC_INPUT_DIR` | `${TRACECB_STUDY_ROOT}/coloc`; prepared GWAS/locus inputs. |
| `TRACECB_COLOC_DIR` | `${TRACECB_OUTPUT_ROOT}/${TARGET_POPULATION}_${TISSUE_SOURCE}/coloc`; generated per-study `*_coloc.csv` tables. |
| `TRACECB_COLOC_GENES` | `${TRACECB_COLOC_INPUT_DIR}/bcx/bcx_mon.closest.protein_coding.bed`. |
| `TRACECB_COLOC_REPLICATION` | `${TRACECB_COLOC_DIR}/replication.csv`; separately prepared replication table required by the heatmap. |
| `TRACECB_COLOC_FIGURE_DIR` | `${TRACECB_FIGURE_DIR}/coloc`, shared by Python and R. |
| `TRACECB_LOCUS_GWAS` | `${TRACECB_COLOC_INPUT_DIR}/bcx/bcx_mon_GWAS.csv`. |
| `TRACECB_LOCUS_EQTL_DIR` | `${TRACECB_STUDY_DIR}`; locus scripts require exported per-gene CSVs in their expected layout. |
| `TRACECB_LOCUS_TRACK_DIR` | `${TRACECB_DATA_ROOT}/locuszoom`; external ENCODE bigWig files. |
| `TRACECB_FIGURE_METADATA` | Repository `src/figures/metadata.json`. |

The OASIS plot validates all five input files before reading data or plotting.
A missing file stops execution and points to `source scripts/config.sh`;
the black `x` marker is only used for genes absent from an available summary.
Pathway replication thresholds differ by reference: OASIS uses `+` for
5×10⁻³ > p ≥ 10⁻⁵ and `++` for p < 10⁻⁵, while CIMA uses only `++` for
p < 10⁻⁵, with no 5×10⁻³ tier. In the CIMA plots, `x` means absent from the
available lead cis-eQTL summary for the matching broad cell type.
Configuration does not generate missing colocalization, replication, or locus
inputs; the corresponding scripts report missing files before plotting.
