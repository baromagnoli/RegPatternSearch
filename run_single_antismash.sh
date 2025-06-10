#!/bin/bash -f
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -J antismash_dsl2
#SBATCH --time=10:00:00


source /temporario2/9877294/anaconda3/etc/profile.d/conda.sh
conda activate

# Diretório de saída
OUTPUT_DIR="/temporario2/9877294/antismash_output_all/${BASENAME}"
mkdir -p "$OUTPUT_DIR"

# Executa antiSMASH
antismash \
    --output-dir "$OUTPUT_DIR" \
    --genefinding-tool prodigal \
    --cb-general \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --pfam2go \
    "$GENOME"
