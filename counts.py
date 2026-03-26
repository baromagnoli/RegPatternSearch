import json  #Biblioteca para manipular arquivos JSON
import re  #Biblioteca para trabalhar com expressões regulares
from collections import Counter  #Para contar ocorrências de itens
import matplotlib.pyplot as plt  #Biblioteca para plotagem de gráficos
import seaborn as sns  #Biblioteca para gráficos estatísticos mais bonitos
import pandas as pd  #Biblioteca para manipulação de dados (não usada diretamente aqui, mas importada)

#Caminho do arquivo JSONL exportado do NCBI Datasets
jsonl_path = "/home/barbara/documents/RegPatternSearch/assembly_data_report.jsonl"

#Abrir o arquivo JSONL e carregar todas as linhas, convertendo cada linha para objeto Python com json.loads()
with open(jsonl_path, 'r') as f:
    linhas = [json.loads(linha) for linha in f]

#Função para extrair o nome do organismo de forma segura, tentando dois locais diferentes no JSON
def extrair_nome_organismo(linha):
    try:
        return linha["organism"]["organismName"]  #Primeiro tenta extrair aqui
    except KeyError:
        try:
            #Se falhar, tenta extrair em outro caminho dentro do JSON
            return linha["assemblyInfo"]["biosample"]["description"]["organism"]["organismName"]
        except KeyError:
            return None  #Se não encontrar, retorna None

#Lista para armazenar nomes de organismos que começam com "Streptomyces"
strepto_nomes = []
for linha in linhas:
    nome = extrair_nome_organismo(linha)  #Extrai nome do organismo para cada linha
    if nome and nome.startswith("Streptomyces"):  #Filtra somente organismos que começam com "Streptomyces"
        strepto_nomes.append(nome.strip())  #Remove espaços extras e adiciona na lista

total_genomas = len(strepto_nomes)  #Conta total de genomas Streptomyces encontrados

#Função para extrair gênero + espécie usando expressão regular
def extrair_genero_especie(nome):
    #Regex que captura 'Streptomyces sp' ou 'Streptomyces' seguido de uma palavra com letras minúsculas e hífens opcionais
    match = re.match(r"^(Streptomyces (?:sp|[a-z]+(?:-[a-z]+)?))$", nome)
    if match:
        return match.group(1)  #Retorna o nome da espécie
    return None  #Caso não corresponda ao padrão, retorna None

#Contador para acumular frequência de todas as espécies (incluindo "sp")
contador_todas = Counter()
for nome in strepto_nomes:
    especie = extrair_genero_especie(nome)  #Extrai gênero + espécie
    if especie:
        contador_todas[especie] += 1  #Incrementa contador para a espécie correspondente

#Imprime o total de genomas Streptomyces analisados
print(f"Total de genomas Streptomyces analisados: {total_genomas}\n")

#Imprime a contagem completa por espécie, incluindo "sp"
print("Contagem por espécie (incluindo 'sp'):")
for especie, contagem in contador_todas.most_common():
    print(f"{especie}: {contagem}")  # Imprime cada espécie e sua contagem, ordenado do mais comum ao menos comum

#Remove da contagem a espécie "Streptomyces sp" para gerar o top 10 real
contador_filtrado = Counter({k: v for k, v in contador_todas.items() if k != "Streptomyces sp"})
top10 = contador_filtrado.most_common(10)  # Obtém as 10 espécies mais frequentes (excluindo "Streptomyces sp")

#Imprime o top 10 das espécies mais frequentes, excluindo "Streptomyces sp"
print("\nTop 10 espécies Streptomyces mais frequentes (excluindo 'Streptomyces sp'):")
for especie, contagem in top10:
    print(f"{especie}: {contagem}")

#Preparar listas para plotagem: espécies e suas contagens correspondentes
especies = [e for e, _ in top10]
contagens = [c for _, c in top10]

#Criar figura para o gráfico com tamanho personalizado
plt.figure(figsize=(12, 6))

#Definir paleta de cores usando seaborn para número de barras igual ao número de espécies
paleta = sns.color_palette("tab10", n_colors=len(especies))

#Criar gráfico de barras, usando a lista de espécies para o eixo x e contagens para y
#Usa o argumento hue para colorir barras separadamente, e legend=False para não exibir legenda extra
ax = sns.barplot(x=especies, y=contagens, hue=especies, palette=paleta, legend=False)

#Definir título e rótulos dos eixos do gráfico
plt.title(f'As 10 espécies de Streptomyces mais frequentes')
plt.xlabel('Espécie')
plt.ylabel('Contagem')

#Rotacionar os nomes no eixo x para melhor visualização
plt.xticks(rotation=45, ha='right')

#Adicionar valores numéricos no topo de cada barra (a contagem)
for i, v in enumerate(contagens):
    ax.text(i, v + max(contagens)*0.01, str(v), ha='center')  # Pequeno deslocamento para ficar acima da barra

#Ajustar layout para evitar corte de elementos no gráfico
plt.tight_layout()

#Salvar a figura em arquivo PNG
plt.savefig("top10_streptomyces_species.png")

#Mostrar o gráfico na tela
plt.show()
