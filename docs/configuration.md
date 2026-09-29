# Configure paths for your machine

Use this guide before running the full-data workflows or manuscript figures.
The public toy tutorial does not require the external study datasets or this
configuration. Run the commands below in **Bash, from the repository root**.
Replace every `/srv/...` example with a real absolute path on your machine.

`scripts/config.sh` currently contains the authors' data and software roots.
You can edit its defaults or override them with exported environment variables.
Configuration selects files; it does not download data, create missing analysis
inputs, install software, or convert genome builds. See the [pipeline guide](pipeline.md)
for data sources, access restrictions and column formats.

## 1. Choose where inputs and results live

| Setting | What it points to | Example |
| --- | --- | --- |
| `TRACECB_DATA_ROOT` | Root containing external datasets and references. | `/srv/tracecb/input` |
| `TRACECB_SOFTWARE_ROOT` | Parent of external tools; defaults expect `plink`, `plink2`, and `s-ldxr-master/s-ldxr.py` below it. | `/srv/tracecb/software` |
| `TRACECB_OUTPUT_ROOT` | Writable root for pipeline outputs. The population/tissue subdirectory is added automatically. | `/srv/tracecb/runs/run01` |
| `TRACECB_STUDY_ROOT` | Root containing study results to **read**, including population/tissue subdirectories. | `/srv/tracecb/runs/run01` |
| `TRACECB_STUDY_DIR` | One population/tissue directory containing `QTD*` study directories. | `/srv/tracecb/runs/run01/EAS_eQTLGen` |
| `TRACECB_FIGURE_DIR` | Writable output directory for figures in the selected configuration. | `/srv/tracecb/runs/run01/figures/EAS_eQTLGen` |
| `TRACECB_LOG_DIR` | Pipeline logs; default is `TRACECB_OUTPUT_ROOT/logs`. | `/srv/tracecb/runs/run01/logs` |
| `TRACECB_REPO_ROOT` | Repository checkout; normally detected automatically. | `/srv/code/traceCB` |

**Pipeline output and figure input are separate settings.** GMM writes to
`TRACECB_OUTPUT_ROOT/<population>_<tissue>/QTD*/GMM/`. Figures read from
`TRACECB_STUDY_DIR`, whose default is under `TRACECB_DATA_ROOT/traceCB`, not under
the output root. To plot a new run, explicitly set `TRACECB_STUDY_ROOT` to its
output root **before** sourcing the configuration. Do not append `EAS_eQTLGen`
to `TRACECB_OUTPUT_ROOT`; otherwise the wrapper adds that component twice.

Use a new output directory for a new analysis version. Existing files are not
automatically protected from overwriting, and stale outputs can otherwise be
mixed with newly generated results.

## 2. Save a machine-specific configuration

Keep local paths outside the checkout, for example in `../tracecb-local.sh`.
This sample assumes the directory layout in section 3:

```bash
# Contents of ../tracecb-local.sh; replace paths before sourcing.
export TRACECB_DATA_ROOT="/srv/tracecb/input"
export TRACECB_SOFTWARE_ROOT="/srv/tracecb/software"
export TRACECB_OUTPUT_ROOT="/srv/tracecb/runs/run01"
export TRACECB_STUDY_ROOT="$TRACECB_OUTPUT_ROOT"

export TARGET_POPULATION="EAS"   # EAS or AFR
export TISSUE_SOURCE="eQTLGen"   # eQTLGen or GTEx; case matters
export PYTHON_ENV="py312"       # existing Conda environment name
export R_ENV="r4"               # needed for the optional coloc wrapper

# These overrides are optional if the software root has the default layout.
export PLINK_BIN="/srv/tracecb/software/plink"
export PLINK2="/srv/tracecb/software/plink2"
export SLDXR_DIR="/srv/tracecb/software/s-ldxr-master"

export TRACECB_GTEX_GENE_ANNOTATION="/srv/tracecb/input/GTEx/gencode.v26.GRCh38.genes.gtf"
export CELL_TYPE_PROPORTION_FILE="/srv/tracecb/input/GTEx/celltype_proportion.csv"
export MAX_JOBS=8
```

Apply it in a shell that has not already sourced a different traceCB configuration:

```bash
source ../tracecb-local.sh
source scripts/config.sh
```

`source scripts/config.sh` exports figure paths, sets `PYTHONPATH` and the
headless Matplotlib backend, and creates output/log directories. It does not
activate a Conda environment. The pipeline shell wrappers activate `PYTHON_ENV`
or `R_ENV` themselves; direct Python/R plotting commands use the active environment.
Do not run `bash scripts/config.sh`: it cannot configure the calling shell and
the script deliberately rejects that invocation.

Environment values take precedence over defaults written as `${VARIABLE:-default}`.
Use `export` so overrides reach wrappers started with `bash scripts/...`.
Do not just assign `DATA_ROOT` or `OUTPUT_ROOT`; the corresponding user overrides
are `TRACECB_DATA_ROOT` and `TRACECB_OUTPUT_ROOT`.

When changing population, tissue or roots, start from a clean environment and
source the appropriate local settings again. A child `bash` launched from an
already configured shell inherits its exports and is **not** a clean reset.
Re-sourcing only `TARGET_POPULATION=AFR` does not recompute existing exported
`TRACECB_STUDY_DIR`, `TRACECB_FIGURE_DIR`, coloc paths or other derived paths.
If reusing a configured shell, explicitly update or unset every affected override.

## 3. Supply the inputs for the stages you will run

### Harmonization, LD scores and GMM

The defaults expect the following layout. Only the selected target population
and tissue source are needed. `chr1` examples repeat through `chr22`.

```text
input/
├── BBJ_eQTL/by_celltype_chr/Monocytes/chr1.csv
├── popcell/AFB_NS/MONO__NS/eQTL_AFB_chr1_assoc.txt.gz
├── eQTLCatalogue/by_celltype_chr/QTD000021/chr1.csv
├── eQTLGen/chr1.tsv.gz
├── GTEx/
│   ├── GTEx_Whole_Blood_by_chr/chr1.csv
│   ├── celltype_proportion.csv
│   └── gencode.v26.GRCh38.genes.gtf
└── 1000G/
    ├── 1000G_EAS/1000G.EAS.QC.maf.1.bed
    ├── 1000G_EAS/1000G.EAS.QC.maf.1.bim
    ├── 1000G_EAS/1000G.EAS.QC.maf.1.fam
    ├── 1000G_EUR/1000G.EUR.QC.maf.1.{bed,bim,fam}
    └── 1000G_AFR/1000G.AFR.QC.maf.1.{bed,bim,fam}
```

`{bed,bim,fam}` denotes three files, not a literal filename. Add all required
cell types and studies, following `STUDY_IDS` and `CELL_TYPES` in the config.

| Override | Required contents / role |
| --- | --- |
| `BBJ_DIR` | **Formatted** EAS target inputs: `<cell-type>/chr<chr>.csv`, e.g. `Monocytes/chr22.csv`. Defaults to `DATA_ROOT/BBJ_eQTL/by_celltype_chr`. |
| `AFR_DIR` | Formatted AFR target inputs: `<label>__NS/eQTL_AFB_chr<chr>_assoc.txt.gz`. Labels are `MONO`, `B`, `NK`, `T.CD4`, `T.CD8`. |
| `EQTL_CATALOGUE_DIR` | Formatted auxiliary EUR inputs: `<QTD study>/chr<chr>.csv`. |
| `EQTLGEN_DIR` | eQTLGen `chr<chr>.tsv.gz`, tab-delimited and retaining the original headers expected by `load_eQTLGen()`. A whole-genome download is not directly interchangeable with these files. |
| `GTEX_SOURCE_DIR` | Raw GTEx download directory; default parent for GTEx files below. |
| `GTEX_DIR` | Formatted GTEx tissue inputs: `chr<chr>.csv`. |
| `CELL_TYPE_PROPORTION_FILE` | CSV named **`celltype_proportion.csv`**, with `Cell_type,Proportion` columns. Values are **percentages**: `20` means 20%, not `0.2`. Use cell-type names matching the config. The preparation wrapper copies the basename unchanged, and GMM expects this exact basename. |
| `LD_REFERENCE_DIR` | Parent of `1000G_EAS`, `1000G_AFR` and, by default, `1000G_EUR`. |
| `AUX_LD_DIR` | Directory containing the EUR reference BED/BIM/FAM triplets; can override the default `LD_REFERENCE_DIR/1000G_EUR`. |

The loader currently chooses input formats using directory names: the target
path must contain `bbj` or `af`, and the tissue path must contain `gtex` or
`eqtlgen` (case-insensitive). Keep these identifiers when choosing custom paths.

`TARGET_EQTL_DIR`, `TISSUE_DIR`, `AUX_EQTL_DIR`, `TARGET_LD_DIR` and `OUTPUT_DIR`
are derived assignments that the config overwrites. Set their parent settings
from the table instead. In particular, exporting `TARGET_LD_DIR` alone does not
override the selected `LD_REFERENCE_DIR/1000G_<population>` directory. For a
different target layout, arrange the selected subdirectory accordingly or use
the underlying Python tools' explicit path arguments.

Put LD triplets directly in the configured population directory. The loader
expects `1000G.*.QC.maf.<chr>.bim`; there should be one matching reference per
chromosome. A download named `1000G.EUR.QC.<chr>` does not match that pattern.
Confirm ancestry, genome build, QC and allele conventions before adapting a
reference layout; changing a filename does not perform QC or liftover.
`harmonize_inputs.py` matches rsIDs and retains tissue coordinates without liftover.
Coordinate-based coloc/enrichment inputs must use matching builds.

Study IDs, cell types, sample sizes and the chromosome list are Bash arrays or
unconditional assignments in `config.sh`, not environment-variable overrides.
For your own cohorts, edit those coordinated lists and any manuscript figure
metadata that assumes the paper's ten studies. Exporting `CHROMOSOMES=22` will
not restrict the wrappers. For one prepared study/chromosome, use the CLI:

```bash
python -m traceCB.run_gmm \
  --study QTD000021 --cell-type Monocytes --chromosome 22 \
  --data-dir "$TRACECB_OUTPUT_ROOT/${TARGET_POPULATION}_${TISSUE_SOURCE}"
```

This command requires already harmonized inputs, gene LD scores and proportions
under that study directory; it does not prepare them.

### Raw-data preprocessing settings

These settings belong to the preprocessing scripts. Their outputs feed the
formatted-input paths above; setting one side does not always update the other.

| Script | Raw input settings | Output setting and next-stage input |
| --- | --- | --- |
| `preprocess_bbj.sh <cell-type>` | `BBJ_SOURCE_DIR/eQTL_<cell-type>.tar.gz`. | `BBJ_OUTPUT_DIR`; set `BBJ_DIR` to the same directory when overriding. |
| `preprocess_eqtl_catalogue.sh` | `EQTL_CATALOGUE_SOURCE_DIR` containing Catalogue study files. | `EQTL_CATALOGUE_DIR`. |
| `preprocess_gtex.sh` | `GTEX_EQTL_FILE` and `GTEX_LOOKUP_FILE` (both gzip files). `TRACECB_GTEX_LOOKUP` also supplies the lookup default; keep overrides consistent. | `GTEX_OUTPUT_DIR`; set `GTEX_DIR` to the same directory when overriding. |
| `preprocess_1000g.sh` | `REFERENCE_DIR` containing `all_hg38.pgen`, `all_hg38.pvar`, `all_hg38.psam`, or the supported compressed source files. | Population subdirectories under `REFERENCE_DIR`; set `LD_REFERENCE_DIR` to this parent for subsequent stages. Requires `PLINK1` and `PLINK2`. |

The exact default GTEx lookup filename is
`GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz`.
The default association filename is
`GTEx_Analysis_v8_QTLs-GTEx_Analysis_v8_eQTL_all_associations-Whole_Blood.allpairs.txt.gz`.
Override their variables if your downloaded filenames differ.
The repository does not currently provide a dedicated eQTLGen chromosome-splitting
entry point; prepare the specified chromosome files before running harmonization.

### External software

| Setting | Required value |
| --- | --- |
| `PLINK_BIN` | PLINK 1.9 executable path, not its containing directory. A command name on `PATH` also works. |
| `PLINK2` | PLINK 2 executable, used by 1000G preprocessing. |
| `PLINK1` | Optional PLINK 1.9 override specific to 1000G preprocessing; otherwise uses `PLINK_BIN`. |
| `SLDXR_DIR` | Directory containing `s-ldxr.py`; its Python dependencies must be installed in the environment running LD scores. |
| `PYTHON_ENV` / `R_ENV` | Existing Conda environment **names**. The Python `environment.yml` does not create the R environment. |
| `PYTHON_BIN` | Optional Python executable used by the main workflow wrappers after environment activation. Prefer the activated environment's `python`. |
| `MAX_JOBS` | Positive integer for the wrappers that implement this limit: GTEx preprocessing, LD scores, GMM and coloc. Input preparation instead launches one process per configured study. |

See the [pipeline guide](pipeline.md) and the figure README for dependencies.
Record actual PLINK/S-LDXR/R versions and commits with your run; paths alone do
not freeze software versions.

## 4. Check the resolved paths before expensive jobs

After sourcing both configuration files, inspect what the wrappers will use:

```bash
for name in TARGET_POPULATION TISSUE_SOURCE TARGET_EQTL_DIR AUX_EQTL_DIR \
  TISSUE_DIR CELL_TYPE_PROPORTION_FILE TARGET_LD_DIR AUX_LD_DIR OUTPUT_DIR \
  TRACECB_STUDY_DIR TRACECB_FIGURE_DIR PLINK_BIN SLDXR_DIR; do
  printf '%s=%s\n' "$name" "${!name}"
done

# Main prerequisites before input preparation / LD scoring:
(
  set -e
  test -d "$TARGET_EQTL_DIR"
  test -d "$AUX_EQTL_DIR"
  test -d "$TISSUE_DIR"
  test -s "$CELL_TYPE_PROPORTION_FILE"
  command -v "$PLINK_BIN"
  test -s "$SLDXR_DIR/s-ldxr.py"

  # Require exactly one correctly named reference and its complete triplet.
  shopt -s nullglob
  for ref_dir in "$TARGET_LD_DIR" "$AUX_LD_DIR"; do
    for chr in "${CHROMOSOMES[@]}"; do
      bims=("$ref_dir"/1000G.*.QC.maf."$chr".bim)
      if (( ${#bims[@]} != 1 )); then
        printf 'Expected one BIM for chr%s in %s; found %s\n' \
          "$chr" "$ref_dir" "${#bims[@]}" >&2
        exit 1
      fi
      for ext in bed bim fam; do
        test -s "${bims[0]%.bim}.$ext" || exit 1
      done
    done
  done
)
```

These checks cover paths and filenames, not scientific validity or full CSV
schemas. Resolve failures before continuing. A pipeline directory may exist
because the config created it even when no analyses have run.

Then run the needed stages in order:

```bash
bash scripts/prepare_inputs.sh
bash scripts/run_ld_scores.sh
bash scripts/run_gmm.sh
```

Before making full-study figures, check completeness explicitly. The plotting
loader can otherwise read a subset of chromosomes without treating it as an error:

```bash
(
  missing=0
  for study in "${STUDY_IDS[@]}"; do
    for chr in "${CHROMOSOMES[@]}"; do
      file="$TRACECB_STUDY_DIR/$study/GMM/chr$chr/summary.csv"
      if [[ ! -s "$file" ]]; then
        printf 'Missing or empty: %s\n' "$file" >&2
        missing=1
      fi
    done
  done
  exit "$missing"
)
```

This only checks presence and size. Also inspect logs, summary columns/rows and
per-gene Parquet outputs. A header-only summary is not evidence of a successful
full analysis. Label deliberately partial analyses as partial.

## 5. Configure manuscript figures

For case-study plots of a newly generated run, the local example already points
the study input at the output root. Install the figure dependencies and activate
the intended Python environment before running:

```bash
conda activate py312
pip install -e '.[figures]'
source ../tracecb-local.sh
source scripts/config.sh
python -m figures.case_study
python -m figures.combine_case_study
```

For **existing** results, use these overrides instead, in a clean shell:

```bash
export TRACECB_DATA_ROOT="/srv/tracecb/input"
export TRACECB_OUTPUT_ROOT="/srv/tracecb/figure-runs/run01"
export TRACECB_STUDY_ROOT="/srv/tracecb/archive"
export TARGET_POPULATION="EAS"
export TISSUE_SOURCE="eQTLGen"
export TRACECB_GTEX_GENE_ANNOTATION="/srv/tracecb/input/GTEx/gencode.v26.GRCh38.genes.gtf"
source scripts/config.sh
# Reads /srv/tracecb/archive/EAS_eQTLGen/QTD*/GMM/...
# Writes /srv/tracecb/figure-runs/run01/figures/EAS_eQTLGen/...
```

`TRACECB_GTEX_GENE_ANNOTATION` is a **file**, not a directory; the current reader
expects an uncompressed GTF. It is also required for gene annotation when the
tissue source is eQTLGen. `TRACECB_FIGURE_METADATA` defaults to the repository's
`src/figures/metadata.json` and normally should not be changed.

Only configure the following extra inputs for figures that consume them:

| Setting | Required input / use |
| --- | --- |
| `TRACECB_AFR_STUDY_DIR` | AFR cohort root with `QTD*/GMM/chr*/summary.csv`; defaults to `STUDY_ROOT/AFR_<tissue>`. |
| `TRACECB_ONEK1K_FILE` | Processed `onek1k_esnp.csv`, not an arbitrary raw OneK1K download. |
| `TRACECB_OASIS_DIR` | OASIS directory with all five `<cell>_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz` files; `<cell>` is `Mono`, `CD4T`, `CD8T`, `B`, `NK`. |
| `TRACECB_GTEX_LOOKUP` | Full GTEx variant lookup gzip; used for SNP matching in replication plots. |
| `TRACECB_CIMA_DIR` | Parent of the CIMA files below; override the individual files if their layout differs. |
| `TRACECB_CIMA_LEAD_EQTL` | `CIMA_Lead_cis-xQTL.csv`. |
| `TRACECB_CIMA_CELL_TYPES` | `CIMA_Cell_Type_Level_and_Marker.xlsx`; the `openpyxl` Excel reader is included in the `figures` extra. |
| `TRACECB_REPLICATION_EGENES` / `TRACECB_REPLICATION_ESNPS` | Processed `hum0343_eGene.csv` / `hum0343_eSNP.csv`. |
| `TRACECB_CELL_PROPORTIONS` | Individual-level `ind_celltype_proportion.csv` for the proportion figure; distinct from GMM's mean `celltype_proportion.csv`. |
| `TRACECB_TIMING_FILE` | `summary_timing.csv` in the schema consumed by `runtime_comparison.py`; a raw GMM text log is not a substitute. |
| `TRACECB_COLOC_INPUT_DIR` | Prepared GWAS/locus directory containing `bcx/` and `bbj/`. Keep `COLOC_INPUT_DIR` consistent if also set. |
| `TRACECB_COLOC_DIR` | Generated `*_coloc.csv` result directory. Setting it selects **figure input**; the coloc shell wrapper writes to `OUTPUT_DIR/coloc`. |
| `TRACECB_COLOC_GENES` | Closest-gene BED file, e.g. `bcx_mon.closest.protein_coding.bed`. |
| `TRACECB_COLOC_REPLICATION` | Separately generated replication table, default `COLOC_DIR/replication.csv`; the main coloc wrapper does not create it. |
| `TRACECB_LOCUS_GWAS` | Prepared locus GWAS CSV, e.g. `bcx_mon_GWAS.csv`. |
| `TRACECB_LOCUS_EQTL_DIR` | Root of exported per-gene CSVs in the layout used by `manhattan_locuszoom.R`; GMM Parquet files alone are insufficient. |
| `TRACECB_LOCUS_TRACK_DIR` | ENCODE bigWig track directory, with the filenames expected by the locus scripts. |

When `TRACECB_STUDY_ROOT` points to a new run, unchanged external resources such
as `onek1k_supp/`, `cell_type_proportion/` and prepared `coloc/` inputs may still
be stored elsewhere. Override their individual paths; do not assume GMM creates
them in the new study root.

Optional output overrides are `TRACECB_AFR_FIGURE_DIR`,
`TRACECB_ESNP_REPLICATION_DIR` and `TRACECB_COLOC_FIGURE_DIR`. Their defaults are
`afr/`, `esnp_replication/` and `coloc/` under `TRACECB_FIGURE_DIR`. The full script
inventory and R packages are documented in `src/figures/README.md`.

## 6. Simulation and enrichment use additional settings

`scripts/config.sh` does **not** configure every simulation/enrichment path.

For the small-window simulations, use their own overrides explicitly:

```bash
SIM_DATA_DIR=/srv/tracecb/input/simulation \
OUT_DIR=/srv/tracecb/runs/run01/simulation \
IMG_DIR=/srv/tracecb/runs/run01/figures/simulation \
NREP=100 NSNP=2000 OMEGA_MODE=both \
bash src/simulation/run_all.sh all
```

`SIM_DATA_DIR` must contain `EAS_n5000_chr22_loci29.npy` and
`EUR_n20000_chr22_loci29.npy`, or set `POP1_GENO` and `POP2_GENO` to individual files.
The matrices are not distributed with the repository. Reproducing the exact
paper inputs additionally requires their data-generation/selection provenance.
Set `CONDA_ENV` for these wrappers; they do not use `PYTHON_ENV`.
Set `SKIP_CONDA=1` only if the active environment is already suitable.

The `all` target does not include the whole-chromosome suite. Its
`src/simulation/chr22/run.sh` entry point also needs PLINK, SNP-gene and LD inputs;
`SIM_DATA_DIR` does not redirect its hard-coded defaults. Use the underlying
`chr22/simulate.py --help` path arguments for a custom input layout. Explicitly
set `NREP` and `IMG_DIR`: the shared setup currently supplies `NREP=100`, so the
later `NREP=1` fallback in `chr22/run.sh` does not take effect. See
`src/simulation/README.md` for the separate grids and plotting commands.

S-LDSC enrichment uses `STUDY_DIR`, `RESULT_DIR`, `GWAS_ROOT`, `LDSC_DIR` and LD
reference **prefixes**, described in `src/enrichment/README.md`. `LDSC_DIR` is the
LDSC checkout, distinct from `SLDXR_DIR`. Its EUR reference stack uses GRCh37;
setting a path does not make GRCh38 intervals compatible.

For ORA, set both `GMT_DIR` (shell download location) and `TRACECB_GMT_DIR`
(Python reader location) to the same directory. Set `RESULT_ROOT` for the ORA
wrapper and provide its study/GTF inputs. Neither this wrapper nor the S-LDSC
wrapper automatically sources `scripts/config.sh`.

## 7. Record what another researcher needs

Alongside each run, save the code commit and local diff, resolved configuration
without credentials, software/environment versions, data source/access details,
genome builds, input checksums, simulation seeds and exact commands. Include a
manifest connecting each paper figure to its input tables and producing command.
Keep API tokens out of config files committed to Git. Complete path configuration
is necessary, but does not replace missing inputs, scientific validation or
the provenance required for exact paper reproduction.
