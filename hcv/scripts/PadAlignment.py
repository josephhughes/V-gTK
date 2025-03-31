# Usage: python pad_subMSA.py -r <reference_alignment_dir> -i <input_dir> [-o <output_dir>] [--keep_intermediate_files]

import os
import argparse
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq

class PadAlignment:
	def __init__(self, reference_alignment, input_dir, base_dir, output_dir, keep_intermediate_files):
		self.reference_alignment = reference_alignment
		self.input_dir = input_dir
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.keep_intermediate_files = keep_intermediate_files
 
	# Insert gaps from the reference alignment into the subalignment sequences.
	def insert_gaps(self, reference_aligned, subalignment_seqs):
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
	def process_master_alignment(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False):
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		merged_sequences = []  # List to store sequences for the merged output

		for ref_id, ref_record in master_alignment.items():
			ref_aligned = ref_record.seq
			ref_id = ref_id.split('.')[0]  # Remove .1 or any other suffix after the dot
			subalignment_file = os.path.join(input_dir, f"{ref_id}/{ref_id}.aligned.fasta")

			if os.path.exists(subalignment_file):
				print(f"Processing subalignment for {ref_id} using {subalignment_file}")
				subalignment_seqs = list(SeqIO.parse(subalignment_file, "fasta"))
				# Insert the gaps from the reference into the subalignment sequences
				updated_seqs = self.insert_gaps(ref_aligned, subalignment_seqs)

				os.makedirs(output_dir, exist_ok=True) 
				output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")

				with open(output_file, "w") as output_handle:
					SeqIO.write(updated_seqs, output_handle, "fasta")
					print(f"Saved updated alignment to {output_file}")
            
					# Add sequences to merged list
					merged_sequences.extend(updated_seqs)
			else:
				print(f"Subalignment file {subalignment_file} not found. Skipping {ref_id}.")

		if merged_sequences:
			# Concatenate all subalignment sequences into one file
			merged_output_file = os.path.join(base_dir, output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
			os.makedirs(join(base_dir, output_dir), exist_ok=True)

			with open(merged_output_file, "w") as merged_output_handle:
				SeqIO.write(merged_sequences, merged_output_handle, "fasta")
				print(f"Saved merged alignment to {merged_output_file}")
		else:
			print(f"No sequences found for {reference_alignment_file}, skipping merged file creation.")

		# Remove intermediate _aligned_padded.fasta files (optional)
		if not keep_intermediate_files:
			for ref_id in master_alignment:
				padded_file = os.path.join(output_dir, f"{ref_id.split('.')[0]}_aligned_padded.fasta")
				if os.path.exists(padded_file):
					os.remove(padded_file)
					print(f"Deleted intermediate file {padded_file}")

	def process_all_fasta_in_directory(self):
		# Process all FASTA files in the specified reference directory.
		for file_name in os.listdir(self.reference_alignment):
			if file_name.endswith(".fasta"):
				reference_alignment_file = os.path.join(self.reference_alignment, file_name)
				print(f"Processing master alignment file: {reference_alignment_file}")
				self.process_master_alignment(reference_alignment_file, self.input_dir, self.base_dir, self.output_dir, self.keep_intermediate_files)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
	parser.add_argument("-r", "--reference_alignment", help="Directory containing master alignment files (FASTA format).", default="tmp/Nextalign/reference_aln/NC_001542")
	parser.add_argument("-i", "--input_dir", help="Directory containing subalignment files (Nextalign output).", default="tmp/Nextalign/query_aln")
	parser.add_argument("-d", "--base_dir", help="Base directory.", default="tmp")
	parser.add_argument("-o", "--output_dir", help="Directory to save padded subalignments and merged files.", default="Pad-alignment")
	parser.add_argument("--keep_intermediate_files", action="store_true", help="Keep intermediate files (padded subalignment). Default: disabled (files will be removed).")

	args = parser.parse_args()

	processor = PadAlignment(args.reference_alignment, args.input_dir, args.base_dir, args.output_dir, args.keep_intermediate_files)
	processor.process_all_fasta_in_directory()
