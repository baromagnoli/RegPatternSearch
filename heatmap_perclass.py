import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Caminho para o arquivo
file_path = "/home/barbara/documents/RegPatternSearch/output_bigscape_most_freq3/output_files/2025-08-07_09-03-35_c0.3/record_annotations.tsv"

# Ler o TSV
df = pd.read_csv(file_path, sep='\t')

# Mapeamento manual para nomes corretos
name_map = {
    "Streptomyces olivaceus": "S. olivaceus",
    "Streptomyces plicatus": "S. rochei",
    "Streptomyces griseus": "S. griseus",
    "Streptomyces mirabilis": "S. mirabilis",
    "Streptomyces graminilatus": "S. scabiei",
    "Streptomyces microflavus": "S. microflavus",
    "Streptomyces paromomycinus": "S. rimosus",
    "Streptomyces anulatus": "S. anulatus",
    "Streptomyces albidoflavus": "S. albidoflavus",
    "Streptomyces virginiae": "S. virginiae"
}
df['Organism'] = df['Organism'].map(name_map)

# Manter apenas as 10 espécies corretas
species_order = list(name_map.values())
df = df[df['Organism'].isin(species_order)]

# Criar tabela pivotada usando **Class**
pivot_df = df.pivot_table(index='Class', columns='Organism', aggfunc='size', fill_value=0)
pivot_df = pivot_df.loc[:, species_order]  # garante a ordem das espécies

# Plotar heatmap com eixos invertidos e sem anotação
plt.figure(figsize=(12, 8))
sns.heatmap(pivot_df, annot=False, cmap='Blues', cbar_kws={'label': 'Número de BGCs'})

# Melhorar layout
plt.ylabel('Classe de BGC', fontsize=12)
plt.xlabel('Espécie', fontsize=12)
plt.title('Distribuição de Classes de BGC por espécie (Heatmap)', fontsize=14)
plt.tight_layout()

# Salvar figura
plt.savefig('bgcs_por_especie_classe_heatmap_horizontal.png', dpi=300)
plt.show()
