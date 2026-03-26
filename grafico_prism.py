import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# === DADOS COMPLETOS COM 10 ESPÉCIES ===
dados = {
    "Streptomyces mirabilis": {
        "NRPS": 6, "PKS": 6, "BUTYROLACTONE": 3, "COMPLETE_BACTERIOCIN": 2, "NIS_SYNTHASE": 2,
        "CLASS_III_IV_LANTIPEPTIDE": 1, "MELANIN": 1, "ECTOINE": 1, "LASSO_PEPTIDE": 1,
        "ANGUCYCLINE": 1, "PENTANGULAR_POLYPHENOL": 1, "CLASS_II_LANTIPEPTIDE": 1
    },
    "Streptomyces virginiae": {
        "NRPS": 10, "PKS": 7, "CLASS_III_IV_LANTIPEPTIDE": 3, "NIS_SYNTHASE": 2,
        "BUTYROLACTONE": 2, "ECTOINE": 1, "PENTANGULAR_POLYPHENOL": 1, "NULL": 1,
        "PHOSPHONATE": 1, "BENZOISOCHROMANEQUINONE": 1, "THIOPEPTIDE": 1, "TERPENE": 1,
        "MELANIN": 1, "NYH": 1
    },
    "Streptomyces griseus": {
        "NRPS": 10, "PKS": 6, "LASSO_PEPTIDE": 3, "ECTOINE": 3, "CLASS_III_IV_LANTIPEPTIDE": 3,
        "CLASS_II_LANTIPEPTIDE": 2, "NIS_SYNTHASE": 2, "BUTYROLACTONE": 2, "PHENAZINE": 1,
        "NULL": 1, "MELANIN": 1, "COMPLETE_BACTERIOCIN": 1, "THIOPEPTIDE": 1, "TERPENE": 1,
        "CLASS_I_LANTIPEPTIDE": 1
    },
    "Streptomyces rochei": {
        "NRPS": 6, "PKS": 5, "BUTYROLACTONE": 2, "CLASS_III_IV_LANTIPEPTIDE": 2,
        "AUREOLIC_ACID": 1, "ANGUCYCLINE": 1, "ECTOINE": 1, "MELANIN": 1, "NIS_SYNTHASE": 1,
        "PENTANGULAR_POLYPHENOL": 1, "CLASS_I_LANTIPEPTIDE": 1
    },
    "Streptomyces olivaceus": {
        "NRPS": 14, "PKS": 6, "CLASS_I_LANTIPEPTIDE": 2, "ANGUCYCLINE": 1,
        "CLASS_III_IV_LANTIPEPTIDE": 1, "PENTANGULAR_POLYPHENOL": 1, "CLASS_II_LANTIPEPTIDE": 1,
        "LASSO_PEPTIDE": 1, "MELANIN": 1, "ECTOINE": 1
    },
    "Streptomyces microflavus": {
        "NRPS": 11, "PKS": 6, "ECTOINE": 2, "BUTYROLACTONE": 2, "CLASS_III_IV_LANTIPEPTIDE": 2,
        "NIS_SYNTHASE": 2, "CLASS_II_LANTIPEPTIDE": 1, "THIOPEPTIDE": 1, "ANGUCYCLINE": 1,
        "LASSO_PEPTIDE": 1, "NULL": 1, "LINARIDIN": 1, "COMPLETE_BACTERIOCIN": 1,
        "MELANIN": 1, "NYH": 1, "GOUGEROTIN": 1
    },
    "Streptomyces anulatus": {
        "NRPS": 14, "PKS": 8, "BUTYROLACTONE": 2, "ECTOINE": 2, "CLASS_III_IV_LANTIPEPTIDE": 2,
        "NIS_SYNTHASE": 2, "MELANIN": 2, "ASCAMYCIN": 1, "NULL": 1, "TERPENE": 1,
        "CLASS_II_LANTIPEPTIDE": 1, "A_201A": 1, "THIOPEPTIDE": 1, "ANGUCYCLINE": 1,
        "LINARIDIN": 1, "COMPLETE_BACTERIOCIN": 1, "ENEDIYNE_9_MEMBERED": 1
    },
    "Streptomyces rimosus": {
        "NRPS": 23, "PKS": 11, "BUTYROLACTONE": 3, "CLASS_III_IV_LANTIPEPTIDE": 3,
        "NIS_SYNTHASE": 2, "LASSO_PEPTIDE": 2, "PHOSPHONATE": 1, "ECTOINE": 1, "ISONITRILE": 1,
        "CLASS_I_LANTIPEPTIDE": 1, "MELANIN": 1, "NULL": 1, "TETRACYCLINE": 1
    },
    "Streptomyces scabiei": {
        "NRPS": 7, "PKS": 6, "CLASS_III_IV_LANTIPEPTIDE": 3, "BUTYROLACTONE": 3,
        "NIS_SYNTHASE": 2, "MELANIN": 2, "CLASS_I_LANTIPEPTIDE": 1, "PENTANGULAR_POLYPHENOL": 1,
        "BOTTROMYCIN": 1, "ECTOINE": 1, "LINARIDIN": 1
    },
    "Streptomyces albidoflavus": {
        "NRPS": 10, "PKS": 5, "NIS_SYNTHASE": 2, "CLASS_III_IV_LANTIPEPTIDE": 2,
        "THIOPEPTIDE": 1, "ECTOINE": 1, "COMPLETE_BACTERIOCIN": 1, "CLASS_II_LANTIPEPTIDE": 1
    }
}

# === CRIAÇÃO DO DATAFRAME ===
df = pd.DataFrame(dados).fillna(0).astype(int).T
df.index = df.index.str.title()  # Formatar nomes das espécies

# === GERAR 25 CORES HSV ===
num_cores = 25
hues = np.linspace(0, 1, num_cores, endpoint=False)
colors = [mcolors.hsv_to_rgb([h, 0.65, 0.85]) for h in hues]

# Mapear cores para colunas (classes)
classes = df.columns.tolist()
cor_mapeamento = {cls: colors[i % num_cores] for i, cls in enumerate(classes)}

#  GRÁFICO DE BARRAS HORIZONTAIS EMPILHADAS 
fig, ax = plt.subplots(figsize=(12, 8))

df.plot(kind="barh", stacked=True, ax=ax, color=[cor_mapeamento[c] for c in df.columns], width=0.7)

# Títulos e rótulos
ax.set_title("Distribuição BGCs por espécie", fontsize=16)
ax.set_xlabel("Número de BGCs", fontsize=14)
ax.set_ylabel("", fontsize=13)

# Legenda ao lado direito
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Classe de BGC")

# Ajustar layout
plt.tight_layout()

# Salvar em PNG
plt.savefig("bgc_frequencia_10_especies_hsv_prism_final.png", dpi=300, bbox_inches="tight")

plt.show()
