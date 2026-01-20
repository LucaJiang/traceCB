# HowTO: Run traceCB pipeline

This page describes how to preprocess data, run traceC and traceCB's gmm model, and visualize the results. 
All logs will be saved according to your `log_path`. Please check the log files carefully.
If you use other eQTL data, please make sure the input data format is the same as described below.

- [xpmm](#xpmm)
  - [Preprocesssing](#preprocesssing)
    - [Format data by CHR](#format-data-by-chr)
    - [eQTLCatalogue Infomation](#eqtlcatalogue-infomation)
    - [Cell type proportion](#cell-type-proportion)
    - [Alignment](#alignment)
    - [Annotate LD](#annotate-ld)
    - [Run s-ldxr](#run-s-ldxr)
  - [GMM](#gmm)
  - [Visualization](#visualization)
  - [COLOC](#coloc)

## Preprocesssing

- eQTLCatalogue: [tabix](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv)
- GTEx: [gcloud](https://console.cloud.google.com/storage/browser/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations;tab=objects?inv=1&invt=Ab037A&prefix=&forceOnObjectsSortingFiltering=true) and [gtex](https://www.gtexportal.org/home/downloads/adult-gtex/qtl)
- BBJ: [bbj](http://jenger.riken.jp/en/result)
- 1000G: [s-ldsc](https://zenodo.org/records/7768714) or [plink2](https://www.cog-genomics.org/plink/1.9/resources#phase1)
- S-LDXR: [github repo](https://github.com/huwenboshi/s-ldxr/tree/master)

### Format data by CHR

To accelerate python loading, we format the data by chromosome.
Preprocess data to `chr@.csv` format

- GTEx Whole Blood: shell/GTEx_Preprocess.sh

Input:

```txt
GTEx_Analysis_v8_QTLs-GTEx_Analysis_v8_eQTL_all_associations-Whole_Blood.allpairs.txt.gz:
gene_id variant_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se
ENSG00000227232.5 chr1_13550_G_A_b38 -16003 19 19 0.0141791 0.734151 0.0587242 0.172837
ENSG00000227232.5 chr1_14671_G_C_b38 -14882 17 17 0.0126866 0.876478 -0.0282343 0.181569

GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table2017.48.22.txt.gz:
variant_id chr variant_pos ref alt num_alt_per_site rs_id_dbSNP151_GRCh38p7 variant_id_b37
chr1_13526_C_T_b38 chr1 13526 C T 1 rs1209314672 1_13526_C_T_b37
chr1_13550_G_A_b38 chr1 13550 G A 1 rs554008981 1_13550_G_A_b37
```

Output:

```txt
chr22.csv
GENE,RSID,CHR,POS,TSS_DISTANCE,A1,A2,MAF,PVAL,BETA,SE
ENSG00000008735,rs117049661,22,49600902,-999783,T,C,0.00671642,0.613759,0.133079,0.263534
ENSG00000008735,rs77177971,22,49600956,-999729,T,C,0.0141791,0.400642,-0.160801,0.191186
```

- BBJ cell type: shell/BBJ_Preprocess.sh

Input example: `eQTL_B_cells.tar.gz` unzipped to `B_cells`

```txt
chr22_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz:
SNP POS REF ALT gene beta t-stat p-value
chr22:16201313:I 16201313 A AG ENSG00000100181.17 0.114427884180853 0.355477022866622 .722966322335031
chr22:16201313:I 16201313 A AG ENSG00000198062.10 -0.0856479858741761 -0.265998487382756 .790777303314436
```

Output example:

```txt
chr22.csv
CHR,RSID,POS,A2,A1,GENE,BETA,Z,PVAL
22,rs1000427,36890105,G,A,ENSG00000100055,-0.152166404682243,-0.858684787452969,0.392527866749798
22,rs1000427,36890105,G,A,ENSG00000100060,-0.039387484946029,-0.221527620769525,0.825124526212667
```

- eQTLCatalog: shell/eQTLCatalog_Preprocess.md

Input example: `QTD000031.all.tsv.gz`

```txt
molecular_trait_id  chromosome  position  ref alt variant ma_samples  maf pvalue  beta  se  ype ac  an  r2  molecular_trait_object_id gene_id median_tpm  rsid
ENSG00000187583 1 14464 A T chr1_14464_A_T  39  0.11976 0.0987921 0.307854  0.185352  NP  40  334 NA  ENSG00000187583 ENSG00000187583 0.512 rs546169444
ENSG00000187608 1 14464 A T chr1_14464_A_T  39  0.11976 0.0987921 0.307854  0.185352  NP  40  334 NA  ENSG00000187608 ENSG00000187608 69.84 rs546169444
```

Output example:

```txt
chr22.csv
CHR,RSID,GENE,POS,A1,A2,BETA,SE,PVAL,Z,N
22,rs5747203,ENSG00000015475,17493644,A,G,-0.0400902,0.148042,0.78691,-0.2708028802637089,167
22,rs2110439,ENSG00000015475,16858862,A,G,0.0271006,0.100051,0.78686,0.27086785739272967,167
```

### eQTLCatalogue Infomation

BLUEPRINT: {QTD000021, QTD000031}
CEDASR index: {QTD000066, QTD000067, QTD000069, QTD000073}
Fair_kasela index: {QTD000081, QTD000371, QTD000372}
Gilchrist_2021 index: {QTD000115}

QTD000021: Monocytes, sample size: 191
QTD000031: CD4+ T cells, sample size: 167

QTD000066: CD8+T cells, sample size: 277
QTD000067: CD4+T cells, sample size: 290
QTD000069: Monocytes, sample size: 286
QTD000073: B cells, sample size: 262

QTD000081: Monocytes, sample size: 420
QTD000371: CD4+ T cells, sample size: 280
QTD000372: CD8+ T cells, sample size: 269

QTD000115: NK_cells, sample size: 247

### Cell type proportion

Note that you can use any other method to obtain cell type proportion.

Use GTEx whole blood tpm file (need individual level data to calculate proportion) and `Cibersortx` online software to calculate cell type proportion.

```csv
Cell_type,Proportion
B_cells,0.8069001744805046
PCs,0.30700693083630154
CD4+T_cells,10.401165801186464
T cells gamma delta,0.4175659034025441
CD8+T_cells,7.140400706490821
NK_cells,9.722658552454147
Monocytes,14.755361092635985
Macrophages,0.8641379801294939
DC,0.34964559978643117
Mast_cells,2.114845815363139
Granulocyte,53.120311443234144
```

### Alignment

Run `src/preprocess/MergeChr.py` by `shell/run_merge.sh` to align all input data files. Result will be saved to `save_path_main`. Cell type proportion file will also be copied to the same directory.

Output example:

```shell
ls <save_path_main>/EAS_GTEx/QTD000021
AUX_Monocytes  INFO  TAR_Monocytes  Tissue

head <save_path_main>/EAS_GTEx/QTD000021/AUX_Monocytes/chr22.csv
RSID,GENE,BETA,SE,Z,PVAL,N
rs4386418,ENSG00000015475,-0.172559,0.107485,-1.6054240126529282,0.110194,191
rs2385713,ENSG00000015475,-0.00438877,0.0981174,-0.0447297828927387,0.964373,191

head <save_path_main>/EAS_GTEx/QTD000021/INFO/chr22.csv
RSID,GENE,A1,A2
rs4386418,ENSG00000015475,A,C
rs2385713,ENSG00000015475,A,G
```

### Annotate LD

Run `src/preprocess/create_ld_annot.py` with `shell/run_ld.sh` to annotate LD information, which will later run `plink` to process 1000G data. This creates files for s-ldxr to use.

Input 1000G data:

```txt
1000G.<pop>.QC.maf.@.bed
1000G.<pop>.QC.maf.@.bim
1000G.<pop>.QC.maf.@.fam
```

Output example:

```txt
<save_path_main>/EAS_GTEx/QTD000021/LDSC/TAR:
1000G.TAR.QC.maf.1.bed
1000G.TAR.QC.maf.1.bim
1000G.TAR.QC.maf.1.fam

<save_path_main>/EAS_GTEx/QTD000021/LDSC/LD_annotation:
1.print_snps.txt
1.annot.gz
```

### Run s-ldxr

Use `shell/run_ld.sh` to run `s-ldxr` for each study. This will calculate the gene level LD scores and save the results in the `save_path_main` directory.

Output example:

```txt
<save_path_main>/EAS_GTEx/QTD000021/LDSC/LDSC_gene/TAR_AUX_std_chr22_pop1.gz
CHR SNP BP  base  ENSG00000015475 ...
22  rs2005631 16869565  9.745773315429688 1.0578601360321045 ...
```

Note if occurs error like:

```text
    self._col = self._internal.col[self._col_indexer]
IndexError: index 16464 is out of bounds for axis 0 with size 16464
```

check if you really aligned all the input correctly.

## GMM

After preprocessing, the folder for a study should look like this:

```shell
├── AUX_Monocytes
│   ├── chr10.csv
│   ├── ...
├── celltype_proportion.csv
├── gene_snp_count
│   ├── chr10.csv
│   ├── ...
│   └── count.csv
├── INFO
│   ├── chr10.csv
│   ├── ...
├── LDSC
│   ├── AUX
│   │   ├── 1000G.AUX.QC.maf.10.bed
│   │   ├── ...
│   ├── LD_annotation
│   │   ├── 10.print_snps.txt
│   │   ├── ...
│   ├── LDSC_gene
│   │   ├── TAR_AUX_std_chr10_pop1.gz
│   │   ├── ...
│   └── TAR
│       ├── 1000G.TAR.QC.maf.10.bed
│       ├── ...
├── TAR_Monocytes
│   ├── chr10.csv
│   ├── ...
└── Tissue
    ├── chr10.csv
    ├── ...
```

Use `shell/run_gmm.sh` to run `src/traceCB/run_gmm.py` for each study. This will calculate the GMM results and save them in the `save_path_main` directory.

Output example:

```shell
<save_path_main>/EAS_GTEx/QTD@/GMM/chr@/ENSG@.parquet: # for each gene
RSID  TAR_SBETA  TAR_CBETA  TAR_TBETA   TAR_SSE   TAR_CSE   TAR_TSE    TAR_SZ    TAR_CZ    TAR_TZ  TAR_SPVAL  TAR_CPVAL     TAR_TPVAL  AUX_SBETA  AUX_CBETA  AUX_TBETA   AUX_SSE   AUX_CSE   AUX_TSE    AUX_SZ    AUX_CZ    AUX_TZ  AUX_SPVAL  AUX_CPVAL     AUX_TPVAL  TISSUE_BETA  TISSUE_SE  TISSUE_Z   TISSUE_PVAL
rs10985869  -0.513371  -0.515194  -0.637880  0.150474  0.143661  0.127948 -3.411698 -3.586188 -4.985484   0.000924   0.000336  6.180674e-07  -0.160632  -0.184845  -0.228108  0.144584  0.053387  0.044848 -1.110994 -3.462335 -5.086225   0.267231   0.000536  3.652615e-07    -0.048964   0.005903   -8.2955  1.080400e-16

<save_path_main>/EAS_GTEx/QTD@/GMM/chr@/summary.csv # for all genes in one chromosome
GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COR_X_ORI,COV_PVAL,COR_X,SIGMAO,RUN_GMM,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
ENSG00000197563,2321,7.496220311530559e-05,2.011502227540292e-05,0.0002531782351655968,2.4300689095640065e-05,1.6474267061036236,2.2415841571456486e-10,0.9899999942102872,6.525303922270945e-06,True,269.0556459941143,933.1492802187638,1548.5439947504565,87,199.0,237.0,1725.9942415217718,1771.6907173545303,2944.1095306644656,536,540.0,718.0,44210.94968964821,1012
```

where 'TAR_\*' means the target population (here is EAS), 'AUX_\*' means the auxiliary population (here is EUR), and '\*_S\*', '\*_C\*', '\*_T\*' means the summary statistics, cross-population enhancement, and cross-population with tissue enhancement results, respectively.

Note that our software use the jit and nogil features of `numba` package to accelerate computation. Feel free to comment out the `numba` decorators (beginning with `@` before functions) if you encounter any issues, which will makes debugging easier but slower.

Also, the gene-level results are saved in `parquet` format to reduce disk usage and speed up loading time. You can use `pandas` to read the files:

```python
import pandas as pd
df = pd.read_parquet("ENSG00000197563.parquet")
```

Or, you can modify the code in `src/traceCB/run_gmm.py` to save as `csv` format if you prefer.

## Visualization

Codes in `src/visual/` are used to visualize the GMM results. The `src/visual/utils.py` controls which study to plot for Python script. And `src/visual/meta_data.json` defines the color, label, and etc. for each study. For visualization of fig 3 or fig 5 in the paper, please choose the corresponding study in `src/visual/utils.py`.

## Colocalization

0. Prerequisites & Dependencies

    *   **LDlinkR API token**: [LDLinkR](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html)
    *   **Bedtools**: The `closestBed` binary is required.
    *   **Reference Files**:
        *   **Cytoband File**: hg19/GRCh37 cytoband definitions.
        *   **Gene Annotation BED**: A BED file containing protein-coding gene coordinates (typically from GENCODE).

Use you own token to access the `ldlinkr` API. You must fill in your token in `src/coloc/run_ldlink.R` before running the scripts.

1. Run `find_leadingSNP.py` to find and annotate leading SNPs

```bash
python src/coloc/find_leadingSNP.py \
  --gwas_sumstats_path <PATH_TO_GWAS_FILE> \
  --save_path <OUTPUT_DIRECTORY> \
  --save_prefix <FILE_PREFIX> \
  [optional arguments]
```

The script used in this paper are commented at the top of `find_leadingSNP.py`.

Output files:

* **`{prefix}_GWAS.csv`**: Standardized GWAS data containing columns `SNP`, `PVALUE`, `CHR`, `POS`, `REF`, `ALT`.
* **`{prefix}_leadingSNP.csv`**: A list of the most significant SNP (leading SNP) for each chromosomal band.
* **`{prefix}_ldlink.csv`**: Output from the R script containing LD information for the identified loci.
* **`{prefix}_GWAS_index_snps.bed`**: A BED file of the index SNPs derived from the LD results.
* **`{prefix}.closest.protein_coding.bed`**: Result from `closestBed`, mapping SNPs to their nearest protein-coding genes.
* **`{prefix}_loci.csv`**: **Final Result**. A summary table merging SNP positions with their closest Ensembl gene IDs and distance.

2. Run COLOC analysis

Use `shell/run_coloc.sh` to run code in `src/coloc/run_coloc.py` for each study. This will calculate the COLOC results and save them in the `save_path_main` directory.
