import os
from os.path import join
from argparse import ArgumentParser

class NextalignAlignment:
	def __init__(self, query_dir, ref_dir, tmp_dir):
		self.query_dir = query_dir
		self.ref_dir = ref_dir
		self.tmp_dir = tmp_dir
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"

	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def nextalign(self, query_acc_path, ref_acc_path):
		accession = self.path_to_basename(query_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(self.tmp_dir, f'{accession}'),
			'--output-basename', f'{accession}',
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
		for each_query_file in os.listdir(self.query_dir):
			ref_file = each_query_file
			self.nextalign(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file)
			)

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-q', '--query_dir', help='Query file directory.', default="tmp/Blast/grouped_fasta")
	parser.add_argument('-r', '--ref_dir', help='Reference fasta directory', default="tmp/Blast/ref_seqs")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp/Nextalign")
	args = parser.parse_args()

	processor = NextalignAlignment(args.query_dir, args.ref_dir, args.tmp_dir)
	processor.process()
