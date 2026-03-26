#!/bin/bash

#Caminho base onde estão os subdiretórios com os arquivos .gbff
GBFF_DIR="/temporario2/9877294/ncbi_dataset_refseq/data"

#Arquivo de saída com os caminhos encontrados
OUTPUT_FILE="caminho_genomas_refseq.txt"
> "$OUTPUT_FILE"  #Limpa o arquivo se ele já existir

#Lista de espécies a buscar o caminho para o genoma refseq 
species_list=(
    "Streptomyces albidoflavus"
    "Streptomyces anulatus"
    "Streptomyces griseus"

)

# Para cada espécie, procurar o arquivo correspondente
for species in "${species_list[@]}"; do
    echo "Procurando: $species"
    # Procura arquivos que contenham a espécie e que tenham "refseq" no caminho
    result=$(grep -ril "$species" "$GBFF_DIR" --include="*.gbff" | head -n 1)

    if [[ -n "$result" ]]; then
        echo "Encontrado: $result"
        echo -e "$species\t$result" >> "$OUTPUT_FILE"
    else
        echo "Nenhum genoma encontrado."
        echo -e "$species\tNÃO ENCONTRADO" >> "$OUTPUT_FILE"
    fi

    echo "-----------------------------"
done

echo "📄 Caminhos salvos em: $OUTPUT_FILE"

