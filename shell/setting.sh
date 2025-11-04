# control data path for shell scripts

## Global Settings
use_tissue="GTEx"       # Options: GTEx, eQTLGen
target_population="EAS" # Options: EAS, AFR
chrs=($(seq 1 22))

## File Paths
load_path_main="/home/group1/wjiang49/data"
GTEx_path="${load_path_main}/GTEx/GTEx_Whole_Blood_by_chr"
eQTLGen_path="${load_path_main}/eQTLGen"
celltype_proportion_path="${load_path_main}/GTEx/celltype_proportion.csv"
eQTLCatalog_path="${load_path_main}/eQTLCatalogue/by_celltype_chr"
BBJ_path="${load_path_main}/BBJ_eQTL/by_celltype_chr"
AFR_path="${load_path_main}/popcell/AFB_NS"
LDSC_path="${load_path_main}/1000G"
aux_LDSC_path="${LDSC_path}/1000G_EAS_EUR/EUR"
coloc_data_path="/home/group1/wjiang49/data/xpmm/coloc"
save_path_main="/home/group1/wjiang49/data/traceCB" # add suffix later based on tissue and population used

## Code Paths
main_repo="/home/wjiang49/traceCB"
src_path="${main_repo}/src"
shell_path="${main_repo}/shell"
log_path="${main_repo}/log"

## Environment and Software
python_env="py312"
r_env="r412"
plink_path="/home/group1/wjiang49/software/plink"
sldxr_path="/home/group1/wjiang49/software/s-ldxr-master"

## Dataset Information
QTDids=("QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069" "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372")
Celltypes=("Monocytes" "CD4+T_cells" "CD8+T_cells" "CD4+T_cells" "Monocytes" "B_cells" "Monocytes" "NK_cells" "CD4+T_cells" "CD8+T_cells")
GTEx_Samplesize=670
eQTLGen_Samplesize=30000 # aready in eQTLGen files, this value use for plotting
BBJ_Samplesize=(105 103 103 103 105 104 105 104 103 103)
eQTLCatalogue_Samplesize=(191 167 277 290 286 262 420 247 280 269)
AFR_Samplesize=(80 80 80 80 80 80 80 80 80 80)

aux_population_samplesize=("${eQTLCatalogue_Samplesize[@]}") # eQTL Catalogue as auxiliary population
## Set Tissue and Population based on use_tissue and target_population
aux_population_data_path="$eQTLCatalogue_path"
if [ "$target_population" == "EAS" ]; then
    target_population_data_path="$BBJ_path"
    target_population_samplesize=("${BBJ_Samplesize[@]}")
    tar_LDSC_path="${LDSC_path}/1000G_EAS_EUR/EAS"
    save_path_main="${save_path_main}/EAS_"
elif [ "$target_population" == "AFR" ]; then
    target_population_data_path="$AFR_path"
    target_population_samplesize=("${AFR_Samplesize[@]}")
    tar_LDSC_path="${LDSC_path}/1000G_AFR/AFR"
    save_path_main="${save_path_main}/AFR_"
else
    echo "Invalid target_population value. Please set it to either 'EAS' or 'AFR'."
    exit 1
fi
if [ "$use_tissue" == "GTEx" ]; then
    tissue_data_path="$GTEx_path"
    tissue_samplesize=($(for i in {1..10}; do echo "$GTEx_Samplesize"; done))
    save_path_main="${save_path_main}GTEx"
elif [ "$use_tissue" == "eQTLGen" ]; then
    tissue_data_path="$eQTLGen_path"
    tissue_samplesize=($(for i in {1..10}; do echo "$eQTLGen_Samplesize"; done))
    save_path_main="${save_path_main}eQTLGen"
else
    echo "Invalid use_tissue value. Please set it to either 'GTEx' or 'eQTLGen'."
    exit 1
fi
## Create necessary directories
if [ ! -d "$log_path" ]; then
    mkdir -p "$log_path"
fi
if [ ! -d "$save_path_main" ]; then
    mkdir -p "$save_path_main"
fi
# End of file
