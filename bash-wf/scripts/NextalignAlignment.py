import os
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser

class NextalignAlignment:
	def __init__(self, query_dir, ref_dir, tmp_dir, ref_fa_file, master_seq_dir, master_ref):
		self.query_dir = query_dir
		self.ref_dir = ref_dir
		self.ref_fa_file = ref_fa_file
		self.master_seq_dir = master_seq_dir
		self.master_ref = master_ref
		self.tmp_dir = tmp_dir
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"

	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]


	def nextalign_master(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(ref_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
			]

		command_str = " ".join(command)
		print(f"Executing command: {command_str}")

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def nextalign(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(query_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]
        
		command_str = " ".join(command)
		print(f"Executing command: {command_str}")
        
		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def process(self):
		query_aln_output_dir = join(self.tmp_dir, "query_aln")
		ref_aln_output_dir = join(self.tmp_dir, "reference_aln")

		#align query against reference sequence alignment
		for each_query_file in os.listdir(self.query_dir):
			ref_file = each_query_file
			self.nextalign(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file), 
					query_aln_output_dir
			)
		
		#align master and reference sequence alignment
		for each_master_ref in os.listdir(self.master_seq_dir):
			self.nextalign_master(
					self.ref_fa_file,
          join(self.master_seq_dir, each_master_ref),
					ref_aln_output_dir 
		)
		
if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-q', '--query_dir', help='Query file directory.', default="tmp/Blast/grouped_fasta")
	parser.add_argument('-r', '--ref_dir', help='Reference fasta directory', default="tmp/Blast/ref_seqs")
	parser.add_argument('-f', '--ref_fa_file', help='Reference fasta file combined, this file is will used to perform Nextalign against the master reference sequence', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-ms', '--master_seq_dir', help='Master sequence directory', default="tmp/Blast/master_seq")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp/Nextalign")
	parser.add_argument('-m', '--master_ref', help='Master reference accession. Generally, the Ref Seq accession. In case of Rabies it is NC_001542', required=True) 
	args = parser.parse_args()

	processor = NextalignAlignment(args.query_dir, args.ref_dir, args.tmp_dir, args.ref_fa_file, args.master_seq_dir, args.master_ref)
	processor.process()
