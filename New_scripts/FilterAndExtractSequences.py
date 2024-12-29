import os
import csv
from os.path import join
from argparse import ArgumentParser

genbank_divisions = ['VRL', 'PAT', 'SYN', 'ENV']

class FilterAndExtractSequences:
	def __init__(self, genbank_matrix, genbank_matrix_filtered, ref_file, tmp_dir, total_length, real_length, prop_ambigious, gb_division, seq_type):
		self.genbank_matrix = genbank_matrix
		self.genbank_matrix_filtered = genbank_matrix_filtered
		self.ref_file = ref_file
		self.tmp_dir = tmp_dir
		self.total_length = total_length
		self.real_length = real_length
		self.prop_ambigious = prop_ambigious
		self.gb_division = gb_division
		self.seq_type = seq_type
		os.makedirs(self.tmp_dir, exist_ok=True)

	def read_ref_list(self):
		ref_list = []
		with open(self.ref_file) as f:
			for each_ref in f:
				ref_list.append(each_ref.strip())
		return ref_list

	def check_gb_division(self):
		if self.gb_division is not None:
			for each_val in self.gb_division:
				if each_val not in genbank_divisions:
					return False
			return True
		return True

	def filter_columns(self):
		if not self.check_gb_division():
			print("Error: Invalid GenBank division specified.")
			return

		ref_list = self.read_ref_list()

		query_seq_file = join(self.tmp_dir, 'query_seq.fa')
		ref_seq_file = join(self.tmp_dir, 'ref_seq.fa')

		with open(query_seq_file, 'w') as write_query_seq, open(ref_seq_file, 'w') as write_ref_seq:
			with open(self.genbank_matrix) as file, open(join(self.genbank_matrix_filtered, 'gB_matrix.tsv'), 'w', newline='') as filtered_matrix:
				csv_reader = csv.DictReader(file, delimiter='\t')
				fieldnames = csv_reader.fieldnames
				writer = csv.DictWriter(filtered_matrix, fieldnames=fieldnames, delimiter='\t')
				writer.writeheader()

				for row in csv_reader:
					gb_division = row['division']
					seq_len = int(row['length'])
					seq_len_without_n = seq_len - int(row['n'])
					accession = row['primary_accession']
					sequence = row['sequence']
					if (	
								(self.gb_division is None or gb_division not in self.gb_division) and
								seq_len >= self.total_length and 
								seq_len_without_n >= self.real_length
								):
						filtered_row = {key: value for key, value in row.items() if key != 'sequence'}
						writer.writerow(filtered_row)
						if accession in ref_list:
							write_ref_seq.write(f">{accession}\n{sequence}\n")
						else:
							write_query_seq.write(f">{accession}\n{sequence}\n")

	def process(self):
		self.filter_columns()

if __name__ == "__main__":
	parser = ArgumentParser(description='Filter sequences and prepare sequences for BLAST alignment')
	parser.add_argument('-g', '--genbank_matrix', help='GenBank matrix file', required=True)
	parser.add_argument('-f', '--genbank_matrix_filtered', help='Filtered GenBank matrix file storage directory', default='tmp/GenBank-matrix')
	parser.add_argument('-r', '--ref_file', help='Text file containing list of reference sequence accessions', required=True)
	parser.add_argument('-t', '--tmp_dir', help='Temporary directory to process the data', default="tmp/Sequences")
	parser.add_argument('-l', '--total_length', help='Total length of the sequence should be more or equal to the length provided', default=1, type=int)
	parser.add_argument('-n', '--real_length', help='Length of sequences without N', default=1, type=int)
	parser.add_argument('-a', '--prop_ambigious_data', help='Proportion of ambiguous data to be excluded.', nargs='+', default=None, type=str)
	parser.add_argument(
		'-d', '--genbank_division',
		help=(
				'GenBank division to be excluded. General divisions are:\n'
				'    "VRL" = Viral sequences\n'
				'    "ENV" = Environmental sequences\n'
				'    "PAT" = Patented sequences\n'
				'    "SYN" = Synthetic and chimeric sequences\n'
			),
				nargs='+',
				default=None,
				type=str
    )
	parser.add_argument('-s', '--seq_type', help='Sequence type', default=None, type=str)
	args = parser.parse_args()

	processor = FilterAndExtractSequences(
		genbank_matrix=args.genbank_matrix,
		genbank_matrix_filtered=args.genbank_matrix_filtered,
		ref_file=args.ref_file,
		tmp_dir=args.tmp_dir,
		total_length=args.total_length,
		real_length=args.real_length,
		prop_ambigious=args.prop_ambigious_data,
		gb_division=args.genbank_division,
		seq_type=args.seq_type
	)
	processor.process()
