#!/bin/bash
# 1000G preprocessing script for EAS, EUR, and AFR populations  

# Set path
PLINK1=/home/group1/wjiang49/software/plink
PLINK2=/home/group1/wjiang49/software/plink2
DATA_ROOT=/home/group1/wjiang49/data/1000G

# Populations to process
POPS=("AFR" "EAS" "EUR")

# ------------------------------------------------------------------
# 1. Global Preparation (Run once in DATA_ROOT)
# ------------------------------------------------------------------
echo "=================================================="
echo "Preparing global files in $DATA_ROOT..."
echo "=================================================="

cd "$DATA_ROOT" || exit 1

# Decompress pgen if needed
if [ -f "all_hg38.pgen.zst" ] && [ ! -f "all_hg38.pgen" ]; then
    echo "Decompressing all_hg38.pgen.zst..."
    $PLINK2 --zst-decompress all_hg38.pgen.zst all_hg38.pgen
fi

# Decompress pvar if needed
if [ -f "all_hg38_rs.pvar.zst" ] && [ ! -f "all_hg38.pvar" ]; then
    echo "Decompressing all_hg38_rs.pvar.zst..."
    $PLINK2 --zst-decompress all_hg38_rs.pvar.zst all_hg38.pvar
fi

# Rename psam file (hg38_corrected.psam -> all_hg38.psam)
if [ -f "hg38_corrected.psam" ] && [ ! -f "all_hg38.psam" ]; then
    echo "Renaming hg38_corrected.psam to all_hg38.psam..."
    mv hg38_corrected.psam all_hg38.psam
fi

# Check if essential files exist
if [ ! -f "all_hg38.pgen" ] || [ ! -f "all_hg38.pvar" ] || [ ! -f "all_hg38.psam" ]; then
    echo "Error: Missing required files (pgen, pvar, or psam) in $DATA_ROOT. Exiting."
    exit 1
fi


# ------------------------------------------------------------------
# Loop through populations
# ------------------------------------------------------------------
for POP in "${POPS[@]}"; do
    echo "=================================================="
    echo "Processing Population: $POP"
    echo "=================================================="
    
    # Define output directory for this population
    OUT_DIR="${DATA_ROOT}/1000G_${POP}"
    mkdir -p "$OUT_DIR"
    
    echo "Output directory: $OUT_DIR"
    
    # ------------------------------------------------------------------
    # 2. Filter samples by SuperPop column 5
    # ------------------------------------------------------------------
    echo "[2/4] Filtering samples for $POP (SuperPop column 5)..."
    
    # Generate sample list in the output directory
    awk -v p="$POP" '!/^#/ && $5 == p {print $1}' "${DATA_ROOT}/all_hg38.psam" > "${OUT_DIR}/${POP}_samples.txt"
    
    SAMPLE_COUNT=$(wc -l < "${OUT_DIR}/${POP}_samples.txt")
    echo "  > Found $SAMPLE_COUNT samples for $POP."

    # ------------------------------------------------------------------
    # 3. Create Dataset and QC
    # ------------------------------------------------------------------
    echo "[3/4] Creating dataset and performing QC..."
    
    # Create raw population dataset in OUT_DIR
    # Input files are in DATA_ROOT
    if [ ! -f "${OUT_DIR}/all_hg38_${POP}.bed" ]; then
        echo "  > Creating all_hg38_${POP}..."
        $PLINK2 \
          --pgen "${DATA_ROOT}/all_hg38.pgen" \
          --psam "${DATA_ROOT}/all_hg38.psam" \
          --pvar "${DATA_ROOT}/all_hg38.pvar" \
          --snps-only just-acgt \
          --max-alleles 2 \
          --chr 1-22 \
          --keep "${OUT_DIR}/${POP}_samples.txt" \
          --make-bed \
          --out "${OUT_DIR}/all_hg38_${POP}"
    fi
     
    # Perform QC
    # Input is the file created in OUT_DIR
    echo "  > Performing QC..."
    $PLINK1 \
      --bfile "${OUT_DIR}/all_hg38_${POP}" \
      --maf 0.05 \
      --geno 0.02 \
      --hwe 1e-6 \
      --make-bed \
      --out "${OUT_DIR}/all_hg38_${POP}.QC"

    # ------------------------------------------------------------------
    # 4. Split by chromosome
    # ------------------------------------------------------------------
    echo "[4/4] Splitting by chromosome..."
    
    # Loop through chromosomes 1-22
    for chr in {1..22}; do
      # echo "  > Processing Chr $chr..."
      $PLINK1 \
        --bfile "${OUT_DIR}/all_hg38_${POP}.QC" \
        --chr ${chr} \
        --make-bed \
        --out "${OUT_DIR}/1000G.${POP}.QC.${chr}"
    done
    
    echo "Success: $POP preprocessing complete, saved in ${OUT_DIR}/1000G.${POP}.QC.*"
    echo ""
done

echo "All populations processed."
