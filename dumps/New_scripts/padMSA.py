#python padMSA.py -r nxt_test/ -i nextalign_cluster/ -o padMSA/

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Insert gaps from the reference alignment into the subalignment sequences.
def insert_gaps(reference_aligned, subalignment_seqs):
    ref_with_gaps_list = list(reference_aligned)
    updated_sequences = []
    
    for seq_record in subalignment_seqs:
        sequence = list(str(seq_record.seq))
        gapped_sequence = []
        seq_idx = 0
        
        # Insert gaps from the reference sequence
        for char in ref_with_gaps_list:
            if char == '-':
                gapped_sequence.append('-')
            else:
                if seq_idx < len(sequence):
                    gapped_sequence.append(sequence[seq_idx])
                    seq_idx += 1

        gapped_seq_str = ''.join(gapped_sequence)
        seq_record.seq = Seq(gapped_seq_str)
        updated_sequences.append(seq_record)
    
    return updated_sequences

# Process a single reference alignment file and its subalignments. 
def process_master_alignment(reference_alignment_file, input_dir, output_dir, keep_intermediate_files=False):
    print("testing")
    master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
    merged_sequences = []  # List to store sequences for the merged output

    for ref_id, ref_record in master_alignment.items():
        ref_aligned = ref_record.seq
        subalignment_file = os.path.join(input_dir, f"{ref_id}/{ref_id}.aligned.fasta")

        if os.path.exists(subalignment_file):
            print(f"Processing subalignment for {ref_id} using {subalignment_file}")
            subalignment_seqs = list(SeqIO.parse(subalignment_file, "fasta"))
            # Insert the gaps from the reference into the subalignment sequences
            updated_seqs = insert_gaps(ref_aligned, subalignment_seqs)

            os.makedirs(output_dir, exist_ok=True)
            
            output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
            with open(output_file, "w") as output_handle:
                SeqIO.write(updated_seqs, output_handle, "fasta")
            print(f"Saved updated alignment to {output_file}")
            
            # Add sequences to merged list
            merged_sequences.extend(updated_seqs)
        else:
            pass
            #print(f"Subalignment file {subalignment_file} not found. Skipping {ref_id}.")

    print ("Length of merged seq", len(merged_sequences))
    if merged_sequences:
        # Concatenate all subalignment sequences into one file
        merged_output_file = os.path.join(output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
        with open(merged_output_file, "w") as merged_output_handle:
            SeqIO.write(merged_sequences, merged_output_handle, "fasta")
        print(f"Saved merged alignment to {merged_output_file}")
    else:
        print(f"No sequences found for {reference_alignment_file}, skipping merged file creation.")

    # Remove intermediate _aligned_padded.fasta files (optional)
    if not keep_intermediate_files:
        for ref_id in master_alignment:
            padded_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
            if os.path.exists(padded_file):
                os.remove(padded_file)
                print(f"Deleted intermediate file {padded_file}")

def process_all_fasta_in_directory(reference_dir, input_dir, output_dir, keep_intermediate_files):
    for file_name in os.listdir(reference_dir):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):
            reference_alignment_file = os.path.join(reference_dir, file_name)
            print(f"Processing master alignment file: {reference_alignment_file}")
            process_master_alignment(reference_alignment_file, input_dir, output_dir, keep_intermediate_files)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
    parser.add_argument("-r", "--reference_alignment", required=True, help="Directory containing master alignment files (FASTA format).")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing subalignment files (Nextalign output).")
    parser.add_argument("-o", "--output_dir", default="tmp/Pad_alignment", help="Directory to save padded subalignments and merged files. Default: Pad_alignment.")
    parser.add_argument("--keep_intermediate_files", action="store_true", help="Flag to keep intermediate files (padded subalignment). Default: remove")

    args = parser.parse_args()

    process_all_fasta_in_directory(args.reference_alignment, args.input_dir, args.output_dir, args.keep_intermediate_files)
