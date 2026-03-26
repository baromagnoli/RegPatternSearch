from Bio import SeqIO
import os

gbff_root_dir = "/home/barbara/documents/RegPatternSearch/ncbi_dataset/data/"
upstream_len = 300
output_fasta = "/home/barbara/documents/RegPatternSearch/background_regions.fasta"

with open(output_fasta, "w") as out_f:
    for root, dirs, files in os.walk(gbff_root_dir):
        for filename in files:
            if filename.endswith(".gbff") or filename.endswith(".gbk") or filename.endswith(".gbff.gz"):
                gbff_path = os.path.join(root, filename)
                print(f"Processando {gbff_path} ...")
                try:
                    for record in SeqIO.parse(gbff_path, "genbank"):
                        seq_len = len(record.seq)
                        for feature in record.features:
                            if feature.type == "CDS":
                                gene_id = feature.qualifiers.get("gene", ["unknown"])[0]
                                strand = feature.location.strand
                                if strand == 1:
                                    start = int(feature.location.start)
                                    upstream_start = max(0, start - upstream_len)
                                    upstream_end = start
                                    seq_upstream = record.seq[upstream_start:upstream_end]
                                elif strand == -1:
                                    end = int(feature.location.end)
                                    upstream_start = end
                                    upstream_end = min(seq_len, end + upstream_len)
                                    seq_upstream = record.seq[upstream_start:upstream_end].reverse_complement()
                                else:
                                    continue
                                
                                if len(seq_upstream) == upstream_len:
                                    header = f">{filename}_{gene_id}_{strand}_{upstream_start}_{upstream_end}"
                                    out_f.write(f"{header}\n{seq_upstream}\n")
                except Exception as e:
                    print(f"Erro ao processar {gbff_path}: {e}")

print(f"Regiões upstream extraídas em: {output_fasta}")
