#!/bin/bash
# 1000G preprocessing script for EAS, EUR, and AFR populations  

# Set path
PLINK2=~/ResearchGroup/software/plink2
DATA_ROOT=~/ResearchGroup/data/1000G

# Populations to process
POPS=("AFR" "EAS" "EUR")

for POP in "${POPS[@]}"; do
    echo "=================================================="
    echo "Processing Population: $POP"
    echo "=================================================="
    
    # Define working directory for this population
    WORK_DIR="${DATA_ROOT}/1000G_${POP}"
    
    if [ ! -d "$WORK_DIR" ]; then
        echo "Error: Directory $WORK_DIR does not exist. Skipping $POP."
        continue
    fi
    
    cd "$WORK_DIR" || exit 1
    
    # ------------------------------------------------------------------
    # 1. prepare files
    # ------------------------------------------------------------------
    echo "[1/4] Preparing files..."

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
        echo "Error: Missing required files (pgen, pvar, or psam) in $WORK_DIR. Skipping..."
        continue
    fi

    # ------------------------------------------------------------------
    # 2. Filter samples by SuperPop column 5
    # ------------------------------------------------------------------
    # Convert POP uppercase to lowercase for filenames (AFR -> afr)
    POP_LOWER=$(echo "$POP" | tr '[:upper:]' '[:lower:]')
    
    echo "[2/4] Filtering samples for $POP (SuperPop column 5)..."
    awk -v p="$POP" '$5 == p {print $1, $2}' all_hg38.psam > "${POP_LOWER}_samples.txt"
    
    SAMPLE_COUNT=$(wc -l < "${POP_LOWER}_samples.txt")
    echo "  > Found $SAMPLE_COUNT samples for $POP."

    # ------------------------------------------------------------------
    # 3. Create Dataset and QC
    # ------------------------------------------------------------------
    echo "[3/4] Creating dataset and performing QC..."
    
    # Create raw population dataset
    # --snps-only just-acgt: keep only ACGT SNPs
    # --max-alleles 2: keep only biallelic variants
    if [ ! -f "all_hg38_${POP_LOWER}.bed" ]; then
        echo "  > Creating all_hg38_${POP_LOWER}..."
        $PLINK2 \
          --pgen all_hg38.pgen \
          --psam all_hg38.psam \
          --pvar all_hg38.pvar \
          --snps-only just-acgt \
          --max-alleles 2 \
          --chr 1-22 \
          --keep "${POP_LOWER}_samples.txt" \
          --make-bed \
          --out "all_hg38_${POP_LOWER}"
    fi
     
    # Perform QC
    # MAF > 0.05, GENO < 0.02, HWE > 1e-6
    echo "  > Performing QC..."
    $PLINK2 \
      --bfile "all_hg38_${POP_LOWER}" \
      --maf 0.05 \
      --geno 0.02 \
      --hwe 1e-6 \
      --make-bed \
      --out "all_hg38_${POP_LOWER}.QC"

    # ------------------------------------------------------------------
    # 4. Split by chromosome
    # ------------------------------------------------------------------
    echo "[4/4] Splitting by chromosome..."
    
    # Loop through chromosomes 1-22
    for chr in {1..22}; do
      # echo "  > Processing Chr $chr..."
      $PLINK2 \
        --bfile "all_hg38_${POP_LOWER}.QC" \
        --chr ${chr} \
        --make-bed \
        --out "1000G.${POP}.QC.${chr}"
    done
    
    echo "Success: $POP preprocessing complete, saved in ${WORK_DIR}/1000G.${POP}.QC.*"
    echo ""
done

echo "All populations processed."