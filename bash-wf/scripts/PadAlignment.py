import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join

class PadAlignmentSequences:
	def __init__(self, reference_sequence, input_dir, output_dir, keep_intermediate_files, output_file):
		self.reference_sequence = reference_sequence
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.keep_intermediate_files = keep_intermediate_files
		self.output_file = output_file

	@staticmethod
	def insert_gaps(reference_aligned, subalignment_seqs):
		ref_with_gaps_list = list(reference_aligned)
		updated_sequences = []
		seq_id = [] 
		for seq_record in subalignment_seqs:
			sequence = list(str(seq_record.seq))
			gapped_sequence = []
			seq_idx = 0
            
			for char in ref_with_gaps_list:
				if char == '-':
					gapped_sequence.append('-')
				else:
					if seq_idx < len(sequence):
						gapped_sequence.append(sequence[seq_idx])
						seq_idx += 1

			gapped_seq_str = ''.join(gapped_sequence)
			seq_record.seq = Seq(gapped_seq_str)
			if seq_record.id not in seq_id:
				updated_sequences.append(seq_record)
				seq_id.append(seq_record.id)
		return updated_sequences

	def process_master_alignment(self, reference_alignment_file):
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		merged_sequences = []

		for ref_id, ref_record in master_alignment.items():
			ref_aligned = ref_record.seq
			subalignment_file = os.path.join(self.input_dir, f"{ref_id}/{ref_id}.aligned.fasta")

			if os.path.exists(subalignment_file):
				print(f"Processing subalignment for {ref_id} using {subalignment_file}")
				subalignment_seqs = list(SeqIO.parse(subalignment_file, "fasta"))
				updated_seqs = self.insert_gaps(ref_aligned, subalignment_seqs)

				os.makedirs(self.output_dir, exist_ok=True)
                
				output_file = os.path.join(self.output_dir, f"{ref_id}_aligned_padded.fasta")
				with open(output_file, "w") as output_handle:
					SeqIO.write(updated_seqs, output_handle, "fasta")
				print(f"Saved updated alignment to {output_file}")
                
				merged_sequences.extend(updated_seqs)

		if merged_sequences:
			merged_output_file = os.path.join(self.output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
			merged_output_file = join(self.output_dir, self.output_file)
			with open(merged_output_file, "w") as merged_output_handle:
				SeqIO.write(merged_sequences, merged_output_handle, "fasta")
			print(f"Saved merged alignment to {merged_output_file}")
		else:
			print(f"No sequences found for {reference_alignment_file}, skipping merged file creation.")

		if self.keep_intermediate_files=="N":
			for ref_id in master_alignment:
				padded_file = os.path.join(self.output_dir, f"{ref_id}_aligned_padded.fasta")
				if os.path.exists(padded_file):
					os.remove(padded_file)
					print(f"Deleted intermediate file {padded_file}")
		else:
			print(f"Stored all the intermediate files at location {self.output_dir}")

	def process(self):
		reference_alignment_file = os.path.join(self.reference_sequence)
		print(f"Processing master alignment file: {reference_alignment_file}")
		self.process_master_alignment(reference_alignment_file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
	parser.add_argument("-r", "--reference_sequence", default="tmp/Sequences/ref_seq.fa", help="FASTA file containing reference sequences")
	parser.add_argument("-i", "--input_dir", default="tmp/Nextalign/query_aln", help="Directory containing Nextalign output alignments for each reference sequences.")
	parser.add_argument("-o", "--tmp_dir", default="tmp/Pad-Alignment", help="Directory to save padded subalignments and merged files.")
	parser.add_argument("-f", "--output_file", default="paded-query-alignment.fa", help="Output file name to store paded alignment file")
	parser.add_argument("--keep_intermediate_files", default="N", help="Flag to keep intermediate files (padded subalignment). Default: remove")
	args = parser.parse_args()
	pad_aln = PadAlignmentSequences(args.reference_sequence, args.input_dir, args.tmp_dir, args.keep_intermediate_files, args.output_file)
	pad_aln.process()

