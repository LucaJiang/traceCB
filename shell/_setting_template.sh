# control data path for shell scripts

# Global Settings, set a name for easy identification, you can modify them as needed
use_tissue="GTEx"
target_population="EUR"
chrs=($(seq 1 22)) # chromosomes to be analyzed

## File Paths
load_path_main="" # please set your main data path here
GTEx_path="${load_path_main}/GTEx/GTEx_Whole_Blood_by_chr" # path to tissue eQTL data
celltype_proportion_path="${load_path_main}/GTEx/celltype_proportion.csv" # path to cell type proportion data
eQTLCatalog_path="${load_path_main}/eQTLCatalogue/by_celltype_chr" # path to auxiliary population eQTL data
BBJ_path="${load_path_main}/BBJ_eQTL/by_celltype_chr" # path to target population eQTL data
LDSC_path="${load_path_main}/1000G" # path to LDSC reference data
aux_LDSC_path="${LDSC_path}/1000G_EAS_EUR/EUR" # path to auxiliary population LDSC reference data
tar_LDSC_path="${LDSC_path}/1000G_EAS_EUR/EAS" # path to target population LDSC reference data
coloc_data_path="" # path to coloc data used in COLOC analysis
save_path_main="" # save path for traceCB results

## Code Paths
main_repo="/home/wjiang49/traceCB" # main repository path
src_path="${main_repo}/src"
shell_path="${main_repo}/shell"
log_path="${main_repo}/log"

## Environment and Software
python_env="py312"
r_env="r412"
plink_path="./software/plink" 
sldxr_path="./software/s-ldxr-master"

## Dataset Information
QTDids=("QTD000021" "QTD000031" "QTD000066" "QTD000067" "QTD000069" "QTD000073" "QTD000081" "QTD000115" "QTD000371" "QTD000372")
Celltypes=("Monocytes" "CD4+T_cells" "CD8+T_cells" "CD4+T_cells" "Monocytes" "B_cells" "Monocytes" "NK_cells" "CD4+T_cells" "CD8+T_cells")
GTEx_Samplesize=670
BBJ_Samplesize=(105 103 103 103 105 104 105 104 103 103)
eQTLCatalogue_Samplesize=(191 167 277 290 286 262 420 247 280 269)

aux_population_samplesize=("${eQTLCatalogue_Samplesize[@]}") # eQTL Catalogue as auxiliary population
## Set Tissue and Population based on use_tissue and target_population
aux_population_data_path="$eQTLCatalogue_path"
target_population_data_path="$BBJ_path"
target_population_samplesize=("${BBJ_Samplesize[@]}")
tar_LDSC_path="${LDSC_path}/1000G_EAS_EUR/EAS"
save_path_main="${save_path_main}/EAS_"

tissue_data_path="$GTEx_path"
tissue_samplesize=($(for i in {1..10}; do echo "$GTEx_Samplesize"; done))
save_path_main="${save_path_main}GTEx"

## Create necessary directories
if [ ! -d "$log_path" ]; then
    mkdir -p "$log_path"
fi
if [ ! -d "$save_path_main" ]; then
    mkdir -p "$save_path_main"
fi
# End of file
