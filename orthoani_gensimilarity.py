import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyorthoani
import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#Função: carregar e concatenar genoma
def load_genome_concat(fasta_file):
    """
    Lê um arquivo FASTA que pode conter múltiplas sequências.
    Todas as sequências são concatenadas em uma única sequência contínua,
    o que é necessário para que o cálculo de ANI seja feito corretamente.
    """
    sequences = SeqIO.parse(fasta_file, "fasta")  #lê as sequências do arquivo
    full_seq = Seq("")  #inicia uma sequência vazia
    for seq_record in sequences:
        full_seq += seq_record.seq  #concatena cada contig
    #retorna um SeqRecord com a sequência total concatenada
    return SeqRecord(full_seq, id=os.path.basename(fasta_file), description="")



#Função principal: calcular matriz e pares ANI
def compute_ani_matrix_and_pairs(genome_folder, matrix_csv="ani_matrix.csv", pairs_csv="ani_pairs.csv"):
    """
    Percorre a pasta fornecida, encontra todos os arquivos de genomas em formato FASTA,
    calcula o ANI entre cada par de genomas usando pyOrthoANI,
    e salva:
      - matriz ANI em formato CSV (em porcentagem)
      - lista de pares comparados com valores ANI
    """
    #extensões aceitas como arquivos de genoma
    valid_exts = (".fa", ".fna", ".fasta", ".FA", ".FNA", ".FASTA")

    #Busca por arquivos FASTA dentro da pasta
    genome_files = []
    for root, dirs, files in os.walk(genome_folder):
        for f in files:
            if f.endswith(valid_exts):
                genome_files.append(os.path.join(root, f))

    #organiza para manter ordem fixa
    genome_files.sort()
    n = len(genome_files)

    if n == 0:
        print("Nenhum arquivo FASTA encontrado.")
        return [], None, None

    #matriz ANI NxN
    ani_matrix = np.zeros((n, n))
    ani_pairs = []  #lista de pares ANI

    #Loop duplo para calcular ANI
    for i in range(n):
        genome_i = load_genome_concat(genome_files[i])  #carrega e concatena genoma i
        for j in range(i, n):
            genome_j = load_genome_concat(genome_files[j])  #carrega e concatena genoma j

            if i == j:
                ani_value = 1.0  #o ANI consigo mesmo é sempre 100%
            else:
                ani_value = pyorthoani.orthoani(genome_i, genome_j)  #cálculo ANI

            #registra na matriz (simétrica)
            ani_matrix[i, j] = ani_value
            ani_matrix[j, i] = ani_value

            #guarda o par se forem genomas diferentes
            if i != j:
                ani_pairs.append((genome_files[i], genome_files[j], ani_value))

            #imprime progresso já convertido para %
            print(f"{os.path.basename(genome_files[i])} vs {os.path.basename(genome_files[j])}: {ani_value*100:.2f}%")

    #Salva a matriz ANI em porcentagem no CSV
    with open(matrix_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([""] + [os.path.basename(g) for g in genome_files])  # cabeçalho
        for i, row in enumerate(ani_matrix):
            writer.writerow([os.path.basename(genome_files[i])] + list(row * 100))

    #Salva os pares ANI (%)
    with open(pairs_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Genome 1", "Genome 2", "ANI (%)"])
        for g1, g2, ani in ani_pairs:
            writer.writerow([os.path.basename(g1), os.path.basename(g2), ani * 100])

    return genome_files, ani_matrix, ani_pairs


#Heatmap
def plot_ani_heatmap(genome_files, ani_matrix, output_file="ani_heatmap_hierarchical.jpg"):
    """
    Gera um heatmap hierárquico usando clustermap,
    que agrupa genomas automaticamente com base na similaridade ANI.
    """
    if ani_matrix is None or len(genome_files) == 0:
        print("Matriz vazia ou lista inválida.")
        return

    labels = [os.path.basename(g) for g in genome_files]  #nomes curtos para exibição

    #clustermap inclui dendrogramas automáticos (hierarchical clustering)
    g = sns.clustermap(
        ani_matrix * 100,     #converte valores para porcentagem
        cmap="RdBu_r",        
        vmin=0, vmax=100,     #escala fixa
        figsize=(12, 10),
        xticklabels=labels,
        yticklabels=labels,
        annot=True, fmt=".1f",   #mostrar valores numéricos
        linewidths=0.5
    )

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    #Título acima
    g.fig.suptitle("Hierarchical OrthoANI Heatmap (%)", y=1.05)

    #Salva figura final
    g.savefig(output_file, dpi=300)
    plt.close()
    print(f"Heatmap hierárquico salvo em: {output_file}")


if __name__ == "__main__":
    #Pasta contendo os genomas baixados do NCBI
    genome_folder = "/home/barbara/documents/RegPatternSearch/ncbi_dataset/data"

    #Calcula ANI
    genomes, ani_mat, ani_pairs = compute_ani_matrix_and_pairs(
        genome_folder,
        matrix_csv="ani_matrix.csv",
        pairs_csv="ani_pairs.csv"
    )

    #Plota o heatmap hierárquico
    plot_ani_heatmap(genomes, ani_mat, output_file="ani_heatmap_hierarchical.jpg")
