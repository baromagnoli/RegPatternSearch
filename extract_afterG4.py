from Bio import SeqIO

# Caminhos dos arquivos
coord_file = "/home/barbara/documents/RegPatternSearch/G4PromFinder/Promoter coordinates.txt"
fasta_file = "/home/barbara/documents/RegPatternSearch/ncbi_dataset/data/GCF_000091305.1/GCF_000091305.1_ASM9130v1_genomic.fna"
output_file = "promoter_sequences_afterG4PROM.fasta"

# Ler a sequência do genoma (assumindo apenas 1 contig)
record = SeqIO.read(fasta_file, "fasta")
genome_seq = record.seq

# Abrir arquivo de saída
with open(output_file, "w") as out:

    with open(coord_file, "r") as f:
        next(f)  # pular cabeçalho

        for line in f:
            if line.strip() == "":
                continue

            parts = line.split()

            region = parts[0]
            strand = parts[1]
            start = int(parts[2])
            end = int(parts[3])

            # Ajuste para Python (0-based)
            seq = genome_seq[start-1:end]

            # Se quiser considerar fita negativa:
            if strand.lower() in ["negativo", "minus", "-"]:
                seq = seq.reverse_complement()

            # Escrever no FASTA
            out.write(f">{region}_{start}_{end}_{strand}\n")
            out.write(str(seq) + "\n")

print("Extração concluída! Arquivo salvo como:", output_file)