import os
import csv
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

# Dicionário de GCFs e suas pastas correspondentes
gcf_to_species_folder = {
    "NZ_CP108647.1": "Streptomyces_Albidoflavus_GenRef",
    "NZ_JBFBHO010000011.1": "Streptomyces_Griseus_GenRef",
    "NZ_CM003601.1": "Streptomyces_Anulatus_GenRef",
    "NZ_CP107957.1": "Streptomyces_Mirabilis_GenRef",
    "NZ_CP054926.1": "Streptomyces_Microflavus_GenRef",
    "NZ_BHZD01000001.1": "Streptomyces_Rimosus_GenRef",
    "NZ_BMUJ01000001.1": "Streptomyces_Rochei_GenRef",
    "NZ_LIQQ01000001.1": "Streptomyces_Scabiei_GenRef",
    # ... (adicione os demais conforme necessário)
}

BASE_RESULTS_PATH = "/home/barbara/documents/RegPatternSearch/"

# Mapeamento nome do regulador (lowercase) → classe
mapa_regulador_classe = {
    "whib family transcriptional regulator": "Whib",
    "luxr c-terminal-related transcriptional regulator": "LuxR",
    "padr family transcriptional regulator": "Padr",
    "laci family dna-binding transcriptional regulator": "LacI",
    "helix-turn-helix transcriptional regulator": "HTH",
    "afsr/sarp family transcriptional regulator": "AfsR/Sarp",
    "marr family winged helix-turn-helix transcriptional regulator" : "MarR",

    # Acrescente mais conforme necessário
}

def carregar_bgcs_csv(csv_path):
    df = pd.read_csv(csv_path, sep="\t")
    required_cols = ["GBK", "Record", "Class"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"O arquivo TSV precisa conter a coluna '{col}'.")
    df["GCF"] = df["Record"].apply(lambda x: x.split(".region")[0])
    return df

def construir_caminho_gbk(gcf, bgc):
    gcf_prefix = gcf.split(".region")[0]
    pasta = gcf_to_species_folder.get(gcf_prefix)
    if not pasta:
        print(f"Pasta para {gcf_prefix} não encontrada.")
        return None
    arquivo_gbk = f"{bgc}.gbk"
    caminho = os.path.join(BASE_RESULTS_PATH, pasta, arquivo_gbk)
    if not os.path.exists(caminho):
        print(f"Arquivo .gbk não encontrado: {caminho}")
        return None
    return caminho

def buscar_reguladores_em_gbk(caminho_gbk):
    reguladores = []
    try:
        for record in SeqIO.parse(caminho_gbk, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    qualifiers = feature.qualifiers
                    produto = qualifiers.get("product", [""])[0].lower().strip()
                    gene = qualifiers.get("gene", [""])[0].lower().strip()
                    if "regulator" in produto or "regulator" in gene:
                        reguladores.append(produto)
    except Exception as e:
        print(f"Erro lendo {caminho_gbk}: {e}")
    return reguladores

def processar_bgcs(df_bgcs):
    # Estrutura: { classe_bgc: Counter({classe_regulador: count, ...}), ... }
    contagem_classes = defaultdict(Counter)
    print(f"🔎 Processando {len(df_bgcs)} BGCs...")
    for idx, row in df_bgcs.iterrows():
        bgc = row["GBK"]
        gcf = row["GCF"]
        classe_bgc = row["Class"]
        caminho_gbk = construir_caminho_gbk(gcf, bgc)
        if caminho_gbk is None:
            continue
        reguladores = buscar_reguladores_em_gbk(caminho_gbk)
        if reguladores:
            for reg in reguladores:
                classe_reg = mapa_regulador_classe.get(reg, "Desconhecida")
                contagem_classes[classe_bgc][classe_reg] += 1
    if not contagem_classes:
        print("Nenhum regulador encontrado após processar todos os BGCs.")
    return contagem_classes

def salvar_csv(contagem_classes, output_path):
    # Salvar o total de reguladores por classe de BGC e classe de regulador
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Classe_BGC", "Classe_Regulador", "Quantidade"])
        for classe_bgc, counter_reg in contagem_classes.items():
            for classe_reg, count in counter_reg.items():
                writer.writerow([classe_bgc, classe_reg, count])

def gerar_grafico(contagem_classes):
    if not contagem_classes:
        print("Dados vazios, gráfico não será gerado.")
        return

    # Montar DataFrame: linhas = classe_bgc, colunas = classe_regulador
    df = pd.DataFrame(contagem_classes).T.fillna(0).astype(int)
    df = df.sort_index()

    # Plot com barras empilhadas
    ax = df.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 8),
        colormap="tab20"
    )

    plt.xlabel("Classe do BGC")
    plt.ylabel("Número de Reguladores")
    plt.title("Distribuição de Classes de Reguladores por Classe de BGC")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title="Classe Regulador", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Adiciona valores nas barras empilhadas
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}',
                        (p.get_x() + p.get_width() / 2, p.get_y() + height / 2),
                        ha='center', va='center',
                        fontsize=8, color='white')

    output_grafico = os.path.join(BASE_RESULTS_PATH, "reguladores_por_classe_bgc_bigscape.jpg")
    plt.savefig(output_grafico, dpi=300)
    plt.close()
    print(f"Gráfico empilhado com rótulos salvo em: {output_grafico}")

def main():
    csv_bgcs_path = "/home/barbara/documents/RegPatternSearch/reguladores_top15classes.tsv"
    df_bgcs = carregar_bgcs_csv(csv_bgcs_path)
    contagem_classes = processar_bgcs(df_bgcs)
    output_csv = os.path.join(BASE_RESULTS_PATH, "reguladores_por_classe_bgc_bigscape.csv")
    print(f"Salvando resultados em {output_csv} ...")
    salvar_csv(contagem_classes, output_csv)
    gerar_grafico(contagem_classes)

if __name__ == "__main__":
    main()
