# How to run traceCB for BBJ, eQTL Catalogue and GTEx data

This page describes how to preprocess data, run traceC and traceCB's gmm model, and visualize the results. If you use other eQTL data, please make sure the input data format is the same as described below.

## Table of Contents

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

## Requirements and setup

1. Install the required python packages:

```shell
pip install numpy pandas scipy pysnptools numba matplotlib seaborn sklearn statsmodels
```

2. Install `plink` and [s-ldxr](https://huwenboshi.github.io/s-ldxr/) software.

3. Fill in the paths in `shell/settings.sh`.

## Preprocesssing

Where to download data and software:

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

- eQTLCatalog: shell/eQTLCatalog_Preprocess.sh

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

<!-- TODO: details -->

Use GTEx whole blood tpm file and [Cibersortx](https://cibersortx.stanford.edu/) online software to calculate cell type proportion.

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

Run `src/preprocess/merge.py` by `shell/run_merge.sh` to align all input data files. Result will be saved to `save_path_main`. Cell type proportion file will also be copied to the same directory.

Output example:

```shell
ls <save_path_main>/EAS_GTEx/QTD000021
AUX_Monocytes  INFO  TAR_Monocytes  Tissue

head <save_path_main>/EAS_GTEx/QTD000021/AUX_Monocytes/chr22.csv
RSID,GENE,BETA,SE,Z,PVAL,N
rs4386418,ENSG00000015475,-0.172559,0.107485,-1.6054240126529282,0.110194,191
rs2385713,ENSG00000015475,-0.00438877,0.0981174,-0.0447297828927387,0.964373,191

head <save_path_main>/EAS_GTEx/QTD000021/INFO/chr22.csv
GENE,RSID,POS,A1,A2
ENSG00000141956,rs72613628,40881223,T,C
ENSG00000141956,rs7279498,40881645,A,G
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

Use `shell/run_gmm.sh` to run `src/run_gmm.py` for each study. This will calculate the GMM results and save them in the `save_path_main` directory.
Output example:

```shell
<save_path_main>/EAS_GTEx/QTD@/GMM/chr@/ENSG@.csv:
RSID,TAR_SBETA,TAR_CBETA,TAR_TBETA,TAR_SSE,TAR_CSE,TAR_TSE,TAR_SZ,TAR_CZ,TAR_TZ,TAR_SPVAL,TAR_CPVAL,TAR_TPVAL,AUX_SBETA,AUX_CBETA,AUX_TBETA,AUX_SSE,AUX_CSE,AUX_TSE,AUX_SZ,AUX_CZ,AUX_TZ,AUX_SPVAL,AUX_CPVAL,AUX_TPVAL,TISSUE_BETA,TISSUE_SE,TISSUE_Z,TISSUE_PVAL
rs10738337,-0.218573418389311,-0.008294419021041131,-0.04038206057243096,0.2431349430015113,0.10591163828453204,0.10244247745371493,-0.898979865629402,-0.07831451911600254,-0.3941925417674147,0.37076020961439,0.9375778683499774,0.6934388916548393,0.0452562,0.03067393532694538,-0.014346591263839595,0.111361,0.11026942057576058,0.10432265251450933,0.4063918247860562,0.27817263541228876,-0.13752134285354997,0.684949,0.7808798399460992,0.890618723856331,-0.0587699,0.0475156,-1.2368548434619369,0.216624
rs10961155,0.0736274810073055,0.09801075539874109,0.08612261969328758,0.164901922904843,0.15878379026850212,0.15775429789552664,0.446492555758688,0.6172592002811225,0.5459288326351818,0.656178806660539,0.537063773310086,0.5851148580952852,0.216722,0.16922020336632115,0.12220997962162766,0.310475,0.21246559055358924,0.20394416913335534,0.6980336581045172,0.7964593368997297,0.5992325259454553,0.486077,0.4257651013263122,0.5490178357069655,-0.0632752,0.107333,-0.5895223277090923,0.555733

<save_path_main>/EAS_GTEx/QTD@/GMM/chr@/summary.csv
GENE,NSNP,H1SQ,H1SQSE,H2SQ,H2SQSE,COV,COV_PVAL,TAR_SNEFF,TAR_CNEFF,TAR_TNEFF,TAR_SeSNP,TAR_CeSNP,TAR_TeSNP,AUX_SNEFF,AUX_CNEFF,AUX_TNEFF,AUX_SeSNP,AUX_CeSNP,AUX_TeSNP,TISSUE_SNEFF,TISSUE_SeSNP
ENSG00000005238,1822.0,2.7857883284876096e-05,1.6657978687744558e-05,1e-12,7.118640817536115e-06,5.2780567716556795e-09,0.9994814551731197,152.86558783402612,,,0.0,,,497416611.8914097,,,0.0,,,638662843.4278672,0.0
ENSG00000011454,1724.0,5.3081856569818305e-05,2.569967695537234e-05,0.00015629291835265752,2.6459813162442765e-05,5.53790512393402e-05,0.022731906153504067,411.55410082730054,1098.801367801141,1211.63562949458,0.0,194.0,197.0,634.9874418978754,730.7892366877431,971.1327181950111,174.0,198.0,207.0,210.78945770297386,54.0
ENSG00000023318,1600.0,1e-12,8.69308557471805e-06,4.829093055210316e-06,6.516001141075617e-06,1e-12,0.9999998243731301,-4382085988.71508,,,0.0,,,813.5234655127368,,,0.0,,,-774.9106238664701,0.0
```

## Visualization

<!-- ```shell
# download visualization results
QTDids=("QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069" "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372")
for j in "${QTDids[@]}"; do
    rsync -avz --progress \
        --include="*.png" --exclude="*" \
        -e ssh \
        wjiang49@burgundy.hpc.cityu.edu.hk:/gpfs1/scratch/wjiang49/xpmm/EAS_GTEx/${j}/results/visualization/ \
        /Users/lucajiang/learn/CityU/xpmm/docs/EAS_GTEx/
done
``` -->

## COLOC

<!-- ```shell
QTDids=("QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069" "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372")
for j in "${QTDids[@]}"; do
    rsync -avz --progress \
        --include="*_coloc.csv" --exclude="*" \
        -e ssh \
        wjiang49@burgundy.hpc.cityu.edu.hk:/gpfs1/scratch/wjiang49/xpmm/EAS_GTEx/${j}/results/coloc/ \
        /Users/lucajiang/learn/CityU/xpmm/coloc/data/
done
``` -->
