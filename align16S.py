import os
import subprocess
from Bio import SeqIO

def concatenate_fasta(fasta_dir, output_fasta):

    #List all FASTA files in the directory
    fasta_files = [os.path.join(fasta_dir, filename) for filename in os.listdir(fasta_dir) if filename.endswith('.fasta')]

    #Open the output file to write the concatenated sequences
    with open(output_fasta, "w") as output_file:
        for fasta_file in fasta_files:
            #Open each individual FASTA file
            with open(fasta_file, "r") as input_file:
                #Write the contents of each FASTA file to the output file
                output_file.write(input_file.read())

    #Confirmation message after the concatenated file is created
    print(f"Concatenated file generated: {output_fasta}")

def align_sequences_with_mafft(input_fasta, output_fasta):
    """
    Function to align sequences using MAFFT.
    - Calls MAFFT via command line.
    - Produces an output file with the aligned sequences.
    """
    #Command to run MAFFT with the --auto option, which automatically selects parameters
    mafft_command = f"mafft --auto {input_fasta} > {output_fasta}"
    
    #Execute the MAFFT command using subprocess
    subprocess.run(mafft_command, shell=True, check=True)
    
    #Confirmation message after the aligned file is created
    print(f"Aligned sequences generated: {output_fasta}")

def main():
    """
    - Concatenates FASTA files.
    - Aligns the sequences using MAFFT.
    """
    #Path to the directory containing the FASTA files
    fasta_dir = "/home/barbara/documents/RegPatternSearch/16S_extracted_principaisgenomas"

    #Path to the output file that will store all concatenated sequences
    concatenated_fasta = "/home/barbara/documents/RegPatternSearch/16S_alinhadas_principaisgenomas/all_sequences10.fasta"
    
    #Path to the output file that will store the aligned sequences
    aligned_fasta = "/home/barbara/documents/RegPatternSearch/16S_alinhadas_principaisgenomas/aligned_sequences10.fasta"

    #Step 1: Concatenate sequences from all FASTA files in the directory
    concatenate_fasta(fasta_dir, concatenated_fasta)

    #Step 2: Align the concatenated sequences using MAFFT
    align_sequences_with_mafft(concatenated_fasta, aligned_fasta)

#Ensure that main() is executed only when the script is run directly
if __name__ == "__main__":
    main()
