#!/bin/bash

GBFF_DIR="/temporario2/9877294/ncbi_dataset/data"
OUTPUT_FILE="genomas_encontrados_ncbidataset.txt"
> "$OUTPUT_FILE"

species_list=(
    "Streptomyces albidoflavus"
    "Streptomyces anulatus"
    "Streptomyces microflavus"
    "Streptomyces griseus"
    "Streptomyces rimosus"
    "Kitasatospora papulosa"
    "Streptomyces rochei"
    "Streptomyces mirabilis"
    "Streptomyces scabiei"
)

for species in "${species_list[@]}"; do
    echo "🔍 Procurando: $species"
    
    # Procura todos os arquivos .gbff que contenham a espécie
    matches=$(grep -ril "$species" "$GBFF_DIR" --include="*.gbff")

    found_refseq_file=""

    # Para cada arquivo encontrado, verifica se tem DBSOURCE RefSeq
    while IFS= read -r file; do
        if grep -q "^DBSOURCE\s*RefSeq" "$file"; then
            found_refseq_file="$file"
            break  # Pega o primeiro arquivo refseq encontrado e sai do loop
        fi
    done <<< "$matches"

    if [[ -n "$found_refseq_file" ]]; then
        echo "Encontrado RefSeq: $found_refseq_file"
        echo -e "$species\t$found_refseq_file" >> "$OUTPUT_FILE"
    else
        echo " Nenhum genoma RefSeq encontrado."
        echo -e "$species\tNÃO ENCONTRADO" >> "$OUTPUT_FILE"
    fi

    echo "-----------------------------"
done

echo "Caminhos salvos em: $OUTPUT_FILE"
