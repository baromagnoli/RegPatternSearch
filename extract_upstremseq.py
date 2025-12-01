import os  # módulo para manipulação de arquivos e diretórios
import csv  # módulo para leitura de arquivos CSV/TSV
from Bio import SeqIO  # Biopython, usado para ler arquivos GenBank (.gbk) e manipular sequências


UPSTREAM_LEN = 250  # tamanho da região promotora upstream a ser extraída (250 nucleotídeos)
ANTISMASH_DIR = "/home/barbara/documents/RegPatternSearch/Resultados_AntiSmash_GenRef"  
# diretório contendo arquivos .gbk gerados pelo antiSMASH

RECORD_ANNOTATIONS = "/home/barbara/documents/RegPatternSearch/output_bigscape_most_freq3/output_files/2025-08-07_09-03-35_c0.3/record_annotations.tsv"  
# arquivo TSV do BiG-SCAPE contendo, para cada cluster, o arquivo GBK de origem e a classe química

OUTPUT_FASTA_ALL = "upstream_promoters_all.fasta"  
# arquivo FASTA geral que irá conter todas as sequências upstream extraídas

OUTPUT_DIR_CLASSES = "/home/barbara/documents/RegPatternSearch/promoters_per_class"  
# diretório onde serão salvos os FASTAs separados por classe química

def extract_all_upstream_from_gbk(gbk_path, upstream_len):
    """
    Extrai sequências upstream (promotoras) de todos os CDS em um arquivo .gbk.
    Considera a orientação da fita (reverse complement para genes na fita negativa).
    Retorna lista de tuplas: (locus_tag, sequência).
    """
    try:
        record = SeqIO.read(gbk_path, "genbank")  # lê o arquivo GenBank
    except Exception as e:
        print(f" Erro lendo {gbk_path}: {e}")  # em caso de erro na leitura, exibe mensagem
        return []  # retorna lista vazia

    upstream_seqs = []  # lista para armazenar sequências upstream

    # percorre todas as features do arquivo GenBank
    for feature in record.features:
        if feature.type == "CDS":  # considera apenas genes codificantes
            strand = feature.location.strand  # orientação da fita: +1 ou -1
            start = int(feature.location.start)  # posição inicial do gene
            end = int(feature.location.end)  # posição final do gene

            if strand == 1:  # fita positiva
                up_start = max(0, start - upstream_len)  # garante que não saia do início do cromossomo
                seq = record.seq[up_start:start]  # extrai a região upstream
            elif strand == -1:  # fita negativa
                up_end = min(len(record), end + upstream_len)  # garante que não ultrapasse o tamanho do cromossomo
                seq = record.seq[end:up_end].reverse_complement()  # extrai e obtém complementar reversa
            else:
                continue  # ignora genes sem orientação definida

            # obtém locus_tag do gene, se não existir, usa "unknown"
            locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
            upstream_seqs.append((locus, seq))  # adiciona tupla à lista

    return upstream_seqs  # retorna lista de tuplas (locus_tag, sequência)



def process_records(record_annotations, antismash_dir, upstream_len, output_fasta_all, output_dir_classes):
    """
    Processa o arquivo TSV de anotações (BiG-SCAPE), extrai promotores dos arquivos .gbk correspondentes,
    e salva em dois tipos de saída: geral e separados por classe química.
    """
    os.makedirs(output_dir_classes, exist_ok=True)  # cria diretório para arquivos por classe, se não existir

    class_files = {}  # dicionário para armazenar arquivos FASTA por classe
    count = 0  # contador de promotores extraídos

    # abre o TSV e o FASTA geral
    with open(record_annotations, newline='') as tsvfile, open(output_fasta_all, "w") as out_fasta_all:
        reader = csv.DictReader(tsvfile, delimiter="\t")  # lê TSV como dicionário

        for row in reader:  # percorre cada linha do TSV
            gbk_base = row["GBK"]  # nome base do arquivo GBK
            classe = row.get("Class", "Unknown").replace(" ", "_")  # classe química do cluster
            gbk_filename = f"{gbk_base}.gbk"  # monta o nome completo do arquivo GBK

            # procura o arquivo GBK nas subpastas do diretório antiSMASH
            gbk_path = None
            for root, dirs, files in os.walk(antismash_dir):
                if gbk_filename in files:
                    gbk_path = os.path.join(root, gbk_filename)
                    break

            if gbk_path:  # se o arquivo foi encontrado
                upstream_seqs = extract_all_upstream_from_gbk(gbk_path, upstream_len)  # extrai sequências upstream
                if upstream_seqs:
                    # abre arquivo FASTA por classe, se ainda não estiver aberto
                    if classe not in class_files:
                        path_classe = os.path.join(output_dir_classes, f"promoters_{classe}.fasta")
                        class_files[classe] = open(path_classe, "w")

                    # escreve cada sequência no FASTA geral e por classe
                    for locus, seq in upstream_seqs:
                        header = f">{gbk_base}_{locus}_upstream"  # cabeçalho FASTA
                        out_fasta_all.write(f"{header}\n{seq}\n")  # escreve no geral
                        class_files[classe].write(f"{header}\n{seq}\n")  # escreve no por classe
                        count += 1
                else:
                    print(f" Nenhuma CDS encontrada em {gbk_path}")  # aviso se não há CDS
            else:
                print(f" Arquivo não encontrado: {gbk_filename}")  # aviso se GBK não existe

    # fecha todos os arquivos por classe
    for f in class_files.values():
        f.close()

    print(f"\n Extração concluída: {count} promotores salvos em {output_fasta_all} e em arquivos por classe na pasta {output_dir_classes}")


def main():
    # chama a função de processamento com todos os parâmetros definidos
    process_records(RECORD_ANNOTATIONS, ANTISMASH_DIR, UPSTREAM_LEN, OUTPUT_FASTA_ALL, OUTPUT_DIR_CLASSES)



if __name__ == "__main__":
    main()  # inicia a execução do script
