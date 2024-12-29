import os
import csv
import subprocess
from os.path import join
from argparse import ArgumentParser
import read_file

class BlastAlignment:
	def __init__(self, query_fasta, db_fasta, tmp_dir, output_file, is_segmented_virus, segment_file=None):
		self.query_fasta = query_fasta
		self.db_fasta = db_fasta
		self.tmp_dir = tmp_dir
		self.output_file = output_file
		self.is_segmented_virus = is_segmented_virus
		self.segment_file = segment_file
		self.db_file_name = os.path.basename(db_fasta)

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

	def run_makeblastdb(self):
		os.makedirs(join(self.tmp_dir, 'DB'), exist_ok=True)
		command = [
			'makeblastdb',
			'-in', self.db_fasta,
			'-out', join(self.tmp_dir, 'DB', self.db_file_name),
			'-title', "alignment",
			'-dbtype', 'nucl'
		]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {self.db_fasta}")
		except subprocess.CalledProcessError as e:
			print(f"Error running makeblastdb: {e}")

	def run_blastn(self):
		command = [
			'blastn',
			'-query', self.query_fasta,
			'-db', join(self.tmp_dir, "DB", self.db_file_name),
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
			'-out', join(self.tmp_dir, self.output_file),
			'-outfmt', "6 qacc sacc pident sstrand"
			]
		try:
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {self.output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running blastn: {e}")

	def write_tophits(self):
		input_file = join(self.tmp_dir, "query_tophits.tsv")
		query_tophit_uniq = join(self.tmp_dir, "query_tophits_uniq.tsv")
		grouped_fasta = join(self.tmp_dir, "grouped_fasta")
		ref_seq_dir = join(self.tmp_dir, "ref_seqs")

		os.makedirs(grouped_fasta, exist_ok=True)
		os.makedirs(ref_seq_dir, exist_ok=True)

		records = {}
		values = {}
	
		with open(input_file, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				col1, col2, col3, col4 = row[0], row[1], float(row[2]), row[3]
				if col1 in values:
					existing_value = values[col1]
					if col3 > existing_value:
						records[col1] = [col1, col2, col3, col4]
						values[col1] = col3
				else:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

		with open(query_tophit_uniq, 'w', newline='') as file:
			writer = csv.writer(file, delimiter='\t')
			for record in records.values():
				writer.writerow(record)
	
		'''	
		grouped_dict = {}
		for each_line in open(query_tophit_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			grouped_dict.setdefault(ref_acc, []).append(query_acc)
			seq_dicts = {rows[0].strip(): rows[1].strip() for rows in read_file.fasta(self.query_fasta)}
		'''

		grouped_dict = {}
		for each_line in open(query_tophit_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if ref_acc not in grouped_dict:
				grouped_dict[ref_acc] = [query_acc]
			else:
				grouped_dict[ref_acc].append(query_acc)

		seq_dicts = {}
		query_seqs = read_file.fasta(self.query_fasta)
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + '\n')
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + '\n')
		
		#writing reference sequences into individual fasta files
		ref_seqs = read_file.fasta(self.db_fasta)
		for rows in ref_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()
			
		for ref_accs in grouped_dict.keys():
			seqs = seq_dicts[ref_accs]
			with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
				write_file.write(">" + ref_accs + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')
			
		
	def process_segment_virus(self, input_file, uniq_hit_output, segment_file, annotated_output):
		uniq_hits = {}
		segments = {}
		segment_assigned = {}

		os.makedirs(join(self.tmp_dir, "segment_sorted"), exist_ok=True)
		os.makedirs(join(self.tmp_dir, "segment_merged_fasta"), exist_ok=True)

		segment_sorted = join(self.tmp_dir, "segment_sorted")
		segment_merged = join(self.tmp_dir, "segment_merged_fasta")

		with open(input_file, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
				if col1 in uniq_hits:
					if float(col3) > float(uniq_hits[col1][2]):
						uniq_hits[col1] = [col1, col2, col3, col4]
					else:
						uniq_hits[col1] = [col1, col2, col3, col4]

		with open(uniq_hit_output, 'w') as write_uniq_hits:
			for k, v in uniq_hits.items():
				write_uniq_hits.write('\t'.join(v) + '\n')
				for line in open(segment_file):
					accession, segment = line.strip().split('\t')
					segments[accession] = segment

		with open(annotated_output, 'w') as write_segment:
			with open(uniq_hit_output, newline='') as file:
				reader = csv.reader(file, delimiter='\t')
				for row in reader:
					col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
					if col2 in segments:
						val = [col1, col2, col3, col4, segments[col2]]
						segment_assigned[col1] = val
						write_segment.write('\t'.join(val) + '\n')
					else:
						print("Could not find the segment for :" + col2)

		fasta_seqs = read_file.fasta(self.query_fasta)
		for each_seq in fasta_seqs:
			header = each_seq[0].strip()
			sequence = each_seq[1].strip()

			if header in segment_assigned:
				accession, reference, score, strand, segment = segment_assigned[header]
				if strand == "plus":
					with open(join(segment_sorted, f"seg_{segment}_plus.fa"), 'a') as write_file:
						write_file.write(f">{header}\n{sequence}\n")
				else:
						with open(join(segment_sorted, f"seg_{segment}_minus.fa"), 'a') as write_file:
							write_file.write(f">{header}\n{sequence}\n")

		for each_seg in os.listdir(segment_sorted):
			name, seg_num, strand = each_seg.split('.')[0].split('_')
			print(f"Processing {name}-{strand}")
			if strand == "plus":
				command = ["seqkit", "grep", "-n", "-f", join(segment_sorted, each_seg), ">>", join(segment_merged, each_seg)]
			else:
				command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(segment_sorted, each_seg), ">>", join(segment_merged, each_seg)]


	def process(self):
		if self.is_segmented_virus == 'Y' and not self.segment_file:
			raise ValueError("Segment file is required for segmented viruses.")

		if self.is_segmented_virus == 'Y':
			self.run_makeblastdb()
			self.run_blastn()
			self.process_segment_virus(
				join(self.tmp_dir, "query_tophits.tsv"),
				join(self.tmp_dir, "query_uniq_tophits.tsv"),
				self.segment_file,
				join(self.tmp_dir, "query_uniq_tophit_annotated.tsv")
				)
		else:
			self.check_blast_exists("blastn")
			self.run_makeblastdb()
			self.run_blastn()
			self.write_tophits()

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the BLAST alignment of query sequences against the given reference sequences')
	parser.add_argument('-q', '--query_fa', help='query fasta file', default="tmp/Sequences/query_seq.fa")
	parser.add_argument('-r', '--ref_fa', help='Blast DB fasta file. (Note: program will consider this file as db file and create the blast index for the given file)', default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-t', '--tmp_dir', help='temp directory to process the data', default="tmp/Blast")
	parser.add_argument('-o', '--output_file', help='output file', default='query_tophits.tsv')
	parser.add_argument('-s', '--is_segmented_virus', help='Type Y for segmented virus else N', default='N')
	parser.add_argument('-f', '--segment_file', help='File containing information about the segments')
	args = parser.parse_args()

	processor = BlastAlignment(
		args.query_fa,
		args.ref_fa,
		args.tmp_dir,
		args.output_file,
		args.is_segmented_virus,
		args.segment_file
		)
	processor.process()

