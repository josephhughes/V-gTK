# Usage: python trim_CDS.py --config_file <config_file> -i <input_dir> [-o <output_dir>]


import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_ranges(ranges, alignment_length):
    trimmed_ranges = []
    for r in ranges:
        start_end = r.split(":")
        start = int(start_end[0]) if start_end[0] else 1
        end = int(start_end[1]) if len(start_end) > 1 and start_end[1] else alignment_length
        
        trimmed_ranges.append((start - 1, end))
    return trimmed_ranges

def trim_alignment(input_fasta, output_fasta, ranges):
    trimmed_seqs = []
    
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        trimmed_seq_parts = []
        alignment_length = len(seq_record.seq)

        trimmed_ranges = parse_ranges(ranges, alignment_length)

        for start, end in trimmed_ranges:
            trimmed_seq_parts.append(str(seq_record.seq[start:end]))

        trimmed_seq = Seq("".join(trimmed_seq_parts))
        
        trimmed_record = SeqRecord(trimmed_seq, id=seq_record.id, description=seq_record.description)
        trimmed_seqs.append(trimmed_record)

    SeqIO.write(trimmed_seqs, output_fasta, "fasta")

def process_config_file(config_file):
    with open(config_file, mode='r') as file:
        csv_reader = csv.DictReader(file, delimiter='\t')
        for row in csv_reader:
            input_fasta = row['input_fasta']
            output_fasta = row['output_fasta']
            ranges = row['ranges'].split(",")  # Split ranges by comma
            
            trim_alignment(input_fasta, output_fasta, ranges)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim sequences in a FASTA file based on specified ranges in the config file.")
    parser.add_argument("--config_file", required=True, help="TSV config file with input and output file names and CDS positions.")

    args = parser.parse_args()

    process_config_file(args.config_file)
