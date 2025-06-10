#!/bin/bash
# Diretório onde estão os subdiretórios de genomas
BASE_DIR="/temporario2/9877294/ncbi_dataset/data"

# Para cada subdiretório dentro do diretório base
for dir in "$BASE_DIR"/*; do
    if [ -d "$dir" ]; then
        # Procura um arquivo .gbff dentro do subdiretório
        gbff=$(find "$dir" -maxdepth 1 -type f -name "*.gbff" | head -n 1)
        if [ -n "$gbff" ]; then
            base=$(basename "$gbff" .gbff)
            echo "Submetendo $base"
            sbatch --export=GENOME="$gbff",BASENAME="$base" run_single_antismash.sh
        else
            echo "Nenhum arquivo .gbff encontrado em $dir"
        fi
    fi
done