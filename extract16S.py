import os
import glob
from Bio import SeqIO
from datetime import datetime

def extract_16S_rRNA(genbank_file, output_fasta, create_dir=False, log_file=None):
    """
    Extracts 16S rRNA sequences from a GenBank file and saves them in FASTA format.
    Only includes features with 'product=16S ribosomal RNA', avoiding duplicates.
    """

    #Check if input GenBank file exists
    if not os.path.exists(genbank_file):
        msg = f"[ERROR] File not found: {genbank_file}"
        print(msg)
        if log_file:
            with open(log_file, "a") as log:
                log.write(msg + "\n")
        return 0

    #Create output directory if needed
    if create_dir:
        output_dir = os.path.dirname(output_fasta)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"[INFO] Directory created: {output_dir}")

    count = 0
    extracted_locations = set()  #Track extracted positions to avoid duplicates

    #Open GenBank and output FASTA files
    with open(genbank_file, "r") as gb, open(output_fasta, "w") as fasta_out:
        for record in SeqIO.parse(gb, "genbank"):
            for feature in record.features:
                if feature.type == "rRNA" and "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].strip().lower()
                    if product == "16s ribosomal rna":
                        location = str(feature.location)
                        if location in extracted_locations:
                            continue  #Skip duplicate
                        extracted_locations.add(location)

                        #Write extracted sequence to FASTA
                        seq = feature.extract(record.seq)
                        header = f">{record.id}_16S_rRNA_{location}"
                        fasta_out.write(f"{header}\n{seq}\n")
                        count += 1

    #Handle case where no sequences were found
    if count > 0:
        msg = f"[OK] {count} 16S sequence(s) extracted from: {genbank_file}"
    else:
        msg = f"[WARNING] No 16S sequence found in: {genbank_file}"
        if os.path.exists(output_fasta):
            os.remove(output_fasta)

    print(msg)
    if log_file:
        with open(log_file, "a") as log:
            log.write(msg + "\n")

    return count


def process_genomes_in_subfolders(root_dir, output_folder="16S_extracted_GENOMASPRINCIPAIS", log_path="log.txt"):
    """
    Searches subdirectories for GenBank files, extracts 16S sequences, and logs the results.
    Skips files if output already exists.
    """

    total_files = 0
    skipped_files = 0
    total_16S = 0

    #Create or reset the log file
    with open(log_path, "w") as log:
        log.write("=== 16S rRNA EXTRACTION LOG ===\n")
        log.write(f"Start: {datetime.now()}\n\n")

    #Find all GenBank files recursively
    genbank_files = glob.glob(os.path.join(root_dir, '**', '*.gb*'), recursive=True)

    for filepath in genbank_files:
        total_files += 1

        #Generate unique output name using parent folder + filename
        parent = os.path.basename(os.path.dirname(filepath))
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        output_name = f"{parent}_{base_name}"

        output_fasta = os.path.join(output_folder, f"{output_name}_16S.fasta")

        #Skip if output already exists
        if os.path.exists(output_fasta):
            msg = f"[SKIPPED] Already exists: {output_fasta}"
            print(msg)
            with open(log_path, "a") as log:
                log.write(msg + "\n")
            skipped_files += 1
            continue

        #Extract 16S sequences
        num_16S = extract_16S_rRNA(
            filepath,
            output_fasta,
            create_dir=True,
            log_file=log_path
        )

        total_16S += num_16S

    #Final summary log
    with open(log_path, "a") as log:
        log.write("\n=== FINAL SUMMARY ===\n")
        log.write(f"Total files found: {total_files}\n")
        log.write(f"Files skipped (already existed): {skipped_files}\n")
        log.write(f"Total 16S sequences extracted: {total_16S}\n")
        log.write(f"End: {datetime.now()}\n")
        log.write("==============================\n")


def main():
    #Input root directory with GenBank files
    root_dir = '/home/barbara/documents/RegPatternSearch/ncbi_dataset/data'
    #Output directory for extracted FASTA files
    output_dir = "16S_resultados"
    #Log file path
    log_file = "log.txt"

    #Start processing
    process_genomes_in_subfolders(root_dir, output_dir, log_file)


if __name__ == "__main__":
    main()
