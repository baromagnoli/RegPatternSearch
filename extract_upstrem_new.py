import os  
import csv 
import re  #expressões regulares (usado para extrair Orig. start)
from Bio import SeqIO 


UPSTREAM_LEN = 1000  
#tamanho da regiao upstream (promotora) que será extraida (1000 pb)

ANTISMASH_DIR = "/home/barbara/documents/RegPatternSearch/Resultados_AntiSmash_GenRef"
#dir onde estão os arquivos .gbk do antiSMASH

GENOMES_DIR = "/home/barbara/documents/RegPatternSearch/ncbi_dataset/data"
#dir com os genomas completos (.fna)

RECORD_ANNOTATIONS = "/home/barbara/documents/RegPatternSearch/output_bigscape_most_freq3/output_files/2025-08-07_09-03-35_c0.3/record_annotations.tsv"
#arquivo TSV contendo os nomes dos arquivos GBK e suas classes

OUTPUT_FASTA_ALL = "upstream_promoters_all.fasta"
#arquivo de saida

OUTPUT_DIR_CLASSES = "/home/barbara/documents/RegPatternSearch/promoters_per_class_test"
#dir onde serão salvos FASTAs separados por classe


def load_all_genomes(genomes_dir):
    """
    Lê todos os arquivos .fna dentro do diretório de genomas
    e armazena as sequências em um dicionário:
    chave = ID do contig
    valor = sequência do contig
    """
    genome_dict = {}  #dicionario vazio para armazenar seqs

    #percorre todas as subpastas do diretório
    for root, _, files in os.walk(genomes_dir):
        for f in files:
            if f.endswith(".fna"):  #considera apenas arquivos FASTA
                path = os.path.join(root, f)  #caminho completo do arquivo

                #percorre cada sequência dentro do FASTA
                for record in SeqIO.parse(path, "fasta"):
                    genome_dict[record.id.split(".")[0]] = record.seq  
                    #armazena sequência com ID como chave

    print(f"{len(genome_dict)} contigs carregados.")  
    #informa quantos contigs foram carregados

    return genome_dict  #retorna o dicionário com todos os genomas      


def get_region_start_from_comment(record):
    """
    Extrai o valor de 'Orig. start' do campo COMMENT do .gbk
    Esse valor indica onde a região começa no genoma completo
    """
    comments = record.annotations.get("comment", "")  
    #pega o campo COMMENT do arquivo GenBank

    #procura padrão "Orig. start: numero"
    match = re.search(r"Orig\. start\s*::\s*(\d+)", comments)

    if match:
        return int(match.group(1)) - 1  
        #converte para inteiro e ajusta para índice 0-based

    return 0  
    #se não encontrar, assume 0 (fallback)


#extract

def extract_upstream(gbk_path, genome_dict, upstream_len):
    """
    Para um arquivo .gbk:
    - encontra o contig correto
    - corrige coordenadas usando Orig.start
    - extrai upstream de cada CDS
    """

    try:
        record = SeqIO.read(gbk_path, "genbank")  
        #le o arquivo .gbk
    except Exception as e:
        print(f"Erro lendo {gbk_path}: {e}")  
        return []  #retorna lista vazia em caso de erro

    #pega ID do contig
    if "accessions" in record.annotations:
        contig_id = record.annotations["accessions"][0]
    else:
        contig_id = record.id

    #verifica se o contig existe no genoma carregado
    if contig_id not in genome_dict:
        print(f"Contig não encontrado: {contig_id}")
        return []

    genome_seq = genome_dict[contig_id]  
    #sequência do genoma correspondente do .gbk

    #obtém deslocamento real da região
    region_start = get_region_start_from_comment(record)

    upstream_seqs = []  #lista para armazenar sequências extraídas

    #percorre todas as features do arquivo
    for feature in record.features:

        if feature.type != "CDS":  
            continue  #ignora  o que não for CDS

        strand = feature.location.strand  #orientação (+1 ou -1)
        start = int(feature.location.start)  #inicio no GBK
        end = int(feature.location.end)  # fim no GBK

        #converte para as coordenadas "reais" no genoma
        real_start = start + region_start
        real_end = end + region_start

        if strand == 1:
            #fita positiva: upstream está antes do start codon
            up_start = max(0, real_start - upstream_len)
            seq = genome_seq[up_start:real_start]

        elif strand == -1:
            #fita negativa: upstream está depois do gene
            up_end = min(len(genome_seq), real_end + upstream_len)
            seq = genome_seq[real_end:up_end].reverse_complement()

        else:
            continue  #ignora se indefinida

        #pega o identificador do gene
        locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]

        #adiciona na lista
        upstream_seqs.append((locus, seq))

    return upstream_seqs  #retorna a lista de tuplas


#para o debug

def debug_coordinate_conversion(gbk_path, genome_dict, upstream_len):
    """
    Função para testar se as coordenadas estão corretas.
    Imprime informações de um CDS.
    """

    print("\n=== TESTE DE COORDENADAS ===")

    record = SeqIO.read(gbk_path, "genbank")

    #identifica a contig
    if "accessions" in record.annotations:
        contig_id = record.annotations["accessions"][0]
    else:
        contig_id = record.id

    print(f"Contig ID: {contig_id}")

    if contig_id not in genome_dict:
        print("Contig não encontrado no genoma!")
        return

    genome_seq = genome_dict[contig_id]

    #pega o Orig.start
    comments = record.annotations.get("comment", "")
    match = re.search(r"Orig\. start\s*::\s*(\d+)", comments)

    if not match:
        print("Orig.start não encontrado!")
        return

    region_start = int(match.group(1)) - 1

    print(f"Orig.start (0-based): {region_start}")

    #pega apenas o primeiro CDS para teste
    for feature in record.features:
        if feature.type == "CDS":

            strand = feature.location.strand
            start = int(feature.location.start)
            end = int(feature.location.end)

            print(f"\nCDS encontrado:")
            print(f"  Strand: {strand}")
            print(f"  Coordenadas no GBK: {start}..{end}")

            #calcula as coordenadas "reais"
            real_start = start + region_start
            real_end = end + region_start

            print(f"  Coordenadas no genoma: {real_start}..{real_end}")

            if strand == 1:
                up_start = max(0, real_start - upstream_len)
                seq = genome_seq[up_start:real_start]
                print(f"  Upstream: {up_start}..{real_start}")

            elif strand == -1:
                up_end = min(len(genome_seq), real_end + upstream_len)
                seq = genome_seq[real_end:up_end].reverse_complement()
                print(f"  Upstream (rev): {real_end}..{up_end}")

            print(f"  Tamanho upstream: {len(seq)} bp")
            print(f"  Seq (50 nt): {seq[:50]}")

            break  #testa apenas um CDS

    print("=== FIM DO TESTE ===\n")


def process_records(genome_dict):
    """
    Percorre o TSV do BiG-SCAPE, encontra os GBKs,
    extrai upstream e salva em arquivos FASTA
    """

    os.makedirs(OUTPUT_DIR_CLASSES, exist_ok=True)  
    #cria o diretório de saída se não existir

    class_files = {}  #armazena arquivos FASTA por classe
    count = 0  #contador de sequências

    #abre TSV e FASTA geral
    with open(RECORD_ANNOTATIONS) as tsvfile, open(OUTPUT_FASTA_ALL, "w") as out_all:

        reader = csv.DictReader(tsvfile, delimiter="\t")

        for row in reader:

            gbk_base = row["GBK"]  
            classe = row.get("Class", "Unknown").replace(" ", "_")
            gbk_filename = f"{gbk_base}.gbk"

            gbk_path = None

            #procura o arquivo GBK
            for root, _, files in os.walk(ANTISMASH_DIR):
                if gbk_filename in files:
                    gbk_path = os.path.join(root, gbk_filename)
                    break

            if not gbk_path:
                print(f"GBK não encontrado: {gbk_filename}")
                continue

            #extrai as sequencias upstream
            upstream_seqs = extract_upstream(gbk_path, genome_dict, UPSTREAM_LEN)

            if not upstream_seqs:
                continue

            #cria os arquivos por classes se necessário
            if classe not in class_files:
                path = os.path.join(OUTPUT_DIR_CLASSES, f"promoters_{classe}.fasta")
                class_files[classe] = open(path, "w")

            #escreve as sequências
            for locus, seq in upstream_seqs:
                header = f">{gbk_base}_{locus}_upstream"

                out_all.write(f"{header}\n{seq}\n")
                class_files[classe].write(f"{header}\n{seq}\n")

                count += 1

    #fecha os arquivos
    for f in class_files.values():
        f.close()

    print(f"\n✔ {count} promotores extraídos com sucesso!")




if __name__ == "__main__":

    print("Carregando genomas...")
    genome_dict = load_all_genomes(GENOMES_DIR)

    #teste
    #test_gbk = "/home/barbara/documents/RegPatternSearch/Resultados_AntiSmash_GenRef/Streptomyces_Albidoflavus_GenRef/NZ_CP108647.1.region001.gbk"
    #debug_coordinate_conversion(test_gbk, genome_dict, UPSTREAM_LEN)

    #execução principal
    process_records(genome_dict)