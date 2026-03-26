import json
import re
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns  # import seaborn

arquivo = "/home/barbara/documents/RegPatternSearch/assembly_data_report.jsonl"

def extrair_genero_especie(nome):
    """Extrai gênero e espécie do nome completo."""
    match = re.match(r"^(Streptomyces [a-z]+)", nome)
    if match:
        return match.group(1)
    return None

contador_especies = Counter()

with open(arquivo, "r", encoding="utf-8") as f:
    for linha in f:
        registro = json.loads(linha)
        
        nome = None
        
        # 1) averageNucleotideIdentity.bestAniMatch.organismName
        ani = registro.get("averageNucleotideIdentity", {})
        best_ani = ani.get("bestAniMatch", {})
        nome = best_ani.get("organismName")
        
        # 2) assemblyInfo.biosample.description.organism.organismName (fallback)
        if not nome:
            biosample = registro.get("assemblyInfo", {}).get("biosample", {})
            descricao = biosample.get("description", {})
            organismo = descricao.get("organism", {})
            nome = organismo.get("organismName")
        
        # Limpar e extrair só Streptomyces spp
        if nome:
            nome_limpo = nome.strip().replace('[', '').replace(']', '')
            especie = extrair_genero_especie(nome_limpo)
            if especie:
                contador_especies[especie] += 1

# Mostrar os 10 mais comuns
top10 = contador_especies.most_common(10)
print("\nTop 10 espécies Streptomyces mais frequentes:")
for especie, contagem in top10:
    print(f"{especie}: {contagem}")

total_genomas = sum(contador_especies.values())
print(f"\nTotal de genomas Streptomyces analisados: {total_genomas}")

# Preparar dados para plot
especies = [item[0] for item in top10]
contagens = [item[1] for item in top10]

plt.figure(figsize=(12, 6))

# Gerar paleta de cores com seaborn
cores = sns.color_palette("husl", len(especies))

# Criar gráfico de barras com cores diferentes
barras = plt.bar(especies, contagens, color=cores)

# Título e eixos
plt.title('10 Espécies de Streptomyces mais frequentes')
plt.xlabel(f'Espécie (Total: {total_genomas} genomas Streptomyces)')
plt.ylabel('Contagem')
plt.xticks(rotation=45, ha='right')

# Adicionar valores nas barras
for barra in barras:
    altura = barra.get_height()
    plt.annotate(f'{altura}',
                 xy=(barra.get_x() + barra.get_width() / 2, altura),
                 xytext=(0, 3),  # deslocamento vertical
                 textcoords="offset points",
                 ha='center', va='bottom')

plt.tight_layout()
plt.savefig("top10_streptomyces_species_colored.png")
plt.show()
