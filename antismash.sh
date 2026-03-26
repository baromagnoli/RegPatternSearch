#!/bin/bash -f
#SBATCH --partition=SP2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -J antismash_single
#SBATCH --time=10:00:00

source /temporario2/9877294/anaconda3/etc/profile.d/conda.sh
conda activate

GENOME="/temporario2/9877294/extracted_gbff/GCF_045039635.1_genomic.gbff"

BASENAME=$(basename "$GENOME" .gbff)

OUTPUT_DIR="/temporario2/9877294/antismash_output/${BASENAME}"
mkdir -p "$OUTPUT_DIR"
    --output-dir "$OUTPUT_DIR" \ #Define o diretório onde os resultados serão gravados
    --genefinding-tool prodigal \ #Usa o Prodigal como ferramenta de predição de genes
    --cb-general \ #Ativa predição de cluster biossintético geral
    --cb-knownclusters \ #Compara clusters previstos com clusters conhecidos: MIBiG database
    --cb-subclusters \ #Detecta subclusters: motivos recorrentes em clusters
    --asf \ #Ativa a detecção automática de famílias de cluster com base em similaridade
    --pfam2go \ #Adiciona anotações funcionais: Pfam - Gene Ontology
    "$GENOME" 
