# Usage: python translate_CDS.py -i <input_dir> [-o <output_dir>]


import re
import os
import argparse
from subprocess import call
from Bio import SeqIO

# Function to pad sequences to same length
def pad_sequences(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    max_length = max(len(record.seq) for record in sequences)
    
    for record in sequences:
        if len(record.seq) < max_length:
            record.seq = record.seq + "-" * (max_length - len(record.seq))
        
    # Overwrite the original FASTA file
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

# Function to translate and modify the sequences
def process_fasta(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over all CDS.fasta files in input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith("CDS.fasta"):
            input_fasta = os.path.join(input_dir, file_name)
            output_fasta = os.path.join(output_dir, f"{file_name[:-9]}_AA.fasta")

            # Translate sequences using seqkit
            call(["seqkit", "translate", "-x", "-w", "0", input_fasta, "-o", output_fasta])

            with open(output_fasta, "r") as fasta_file:
                lines = fasta_file.readlines()

            with open(output_fasta, "w") as fasta_file:
                for line in lines:
                    if not line.startswith(">"):
                        # Replace/remove everything after the first stop codon "*"
                        line = re.sub(r"\*.*", "*", line)
                    fasta_file.write(line)
            
            # Pad sequences to ensure uniform length
            pad_sequences(output_fasta)
            print(f"Padded sequences saved to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate CDS nucleotide alignments to amino acid and pad sequences to uniform length.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing CDS nucleotide alignments (FASTA format).")
    parser.add_argument("-o", "--output_dir", default="protein_aln", help="Output directory for amino acid alignment files (FASTA format).")
    
    args = parser.parse_args()

    process_fasta(args.input_dir, args.output_dir)
