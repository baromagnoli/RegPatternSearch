import os
import csv
from Bio import SeqIO

#Caminhos de entrada
gbk_base_path = "/home/barbara/documents/RegPatternSearch/Resultados_AntiSmash_GenRef"
top15_tsv_path = "/home/barbara/documents/RegPatternSearch/record_annotations_all_ref_top15_10.tsv"
output_file = "reguladores_top15classes_10.tsv"

#Palavras-chave para identificar reguladores
keywords = ["regulator", "transcription factor"]

#Obter BGCs das top15 classes
top15_records = set()
with open(top15_tsv_path, "r") as tsv_in:
    reader = csv.DictReader(tsv_in, delimiter="\t")
    for row in reader:
        record = row["GBK"]
        top15_records.add(record.strip())

#Percorrer os subdiretórios procurando arquivos .gbk
results = []
for root, dirs, files in os.walk(gbk_base_path):
    for file in files:
        if file.endswith(".gbk") or file.endswith(".gbff"):
            full_path = os.path.join(root, file)
            record_id = file.replace(".gbk", "").replace(".gbff", "")
            if record_id not in top15_records:
                continue

            for record in SeqIO.parse(full_path, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        product = feature.qualifiers.get("product", [""])[0].lower()
                        if any(k in product for k in keywords):
                            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                            protein_id = feature.qualifiers.get("protein_id", [""])[0]
                            start = int(feature.location.start) + 1
                            end = int(feature.location.end)
                            strand = "+" if feature.location.strand == 1 else "-"

                            results.append([
                                full_path, product, locus_tag, protein_id, start, end, strand
                            ])

#Escreve saída
with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Arquivo", "Produto", "Locus_tag", "Protein_ID", "Start", "End", "Strand"])
    writer.writerows(results)

print(f"{len(results)} reguladores encontrados. Resultado salvo em: {output_file}")
