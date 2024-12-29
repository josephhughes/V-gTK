import os
import csv
import subprocess
from os.path import join
from argparse import ArgumentParser

class BlastProcessor:
	def __init__(self, query_fa, ref_fa, tmp_dir, output_file):
		self.query_fa = query_fa
		self.ref_fa = ref_fa
		self.tmp_dir = tmp_dir
		self.output_file = output_file
		self.dbtype = "nucl"
		os.makedirs(self.tmp_dir, exist_ok=True)

	def read_ref_list(self):
		ref_list = []
		with open(self.ref_fa) as file:
			for each_ref in file:
				ref_list.append(each_ref.strip())
		return ref_list

	@staticmethod
	def check_blast_exists(command):
		try:
			subprocess.run([command, '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			print(f"{command} is working.")
			return True
		except subprocess.CalledProcessError:
			print(f"{command} is not installed correctly.")
			return False
		except FileNotFoundError:
			print(f"{command} is not found on the system.")
			return False

	@staticmethod
	def run_makeblastdb(ref_fa, dbtype, output_dir):
		os.makedirs(output_dir, exist_ok=True)
		base_name = os.path.basename(ref_fa)
		command = [
			'makeblastdb',
			'-in', ref_fa,
			'-out', join(output_dir, base_name),
			'-title', "alignment",
			'-dbtype', dbtype
			]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {ref_fa}.")
		except subprocess.CalledProcessError as e:
			print(f"Error running makeblastdb: {e}")

	@staticmethod
	def run_blastn(query_fa, db_path, output_file):
		command = [
			'blastn',
			'-query', query_fa,
			'-db', db_path,
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
			'-out', output_file,
			'-outfmt', "6 qacc sacc pident sstrand"
		]
		try:
			print(" ".join(command))
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {output_file}.")
		except subprocess.CalledProcessError as e:
			print(f"Error running blastn: {e}")

	def process(self):
		#query_seq_path, ref_seq_path = self.create_query_and_ref_seq()
		if self.check_blast_exists("blastn"):
			db_dir = join(self.tmp_dir, 'DB')
			#self.run_makeblastdb(self.ref_fa, 'nucl', db_dir)
			db_path = join(db_dir, os.path.basename(self.ref_fa))
			self.run_blastn(self.query_fa, db_path, join(self.tmp_dir, self.output_file))

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the BLAST alignment of query sequences against the given reference')
	parser.add_argument('-q', '--query_fa', help='Query FASTA sequence file', default="tmp/Sequences/query_seq.fa")
	parser.add_argument('-r', '--ref_fa', help='Reference FASTA sequence file', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp/Blast")
	parser.add_argument('-o', '--output_file', help='Name of the alignment file', default='blast_alignment.tsv')
	args = parser.parse_args()

	processor = BlastProcessor(args.query_fa, args.ref_fa, args.tmp_dir, args.output_file)
	processor.process()

