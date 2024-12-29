import os
import re
import csv
import time
import read_file
import subprocess
from os.path import join
from argparse import ArgumentParser
from collections import defaultdict


class SequenceProcessor:
	def __init__(self, genbank_matrix, tmp_dir, blast_hits, pad_aln):
		self.genbank_matrix = genbank_matrix
		self.tmp_dir = tmp_dir
		self.blast_hits = blast_hits
		self.pad_aln = pad_aln
		self.output_dir = tmp_dir
		os.makedirs(self.output_dir, exist_ok=True)

	def load_gb_matrix(self):
		cds_keys = ['start', 'end', 'gene', 'protein_id', 'product']
		meta_file_path = join(self.output_dir, "meta_data.tsv")
		cds_info_file_path = join(self.output_dir, "features.tsv")

		with open(meta_file_path, 'w') as meta_file, open(cds_info_file_path, 'w') as cds_info_file:
			cds_header = ["ref_seq_name", "feature_name", "ref_start", "ref_end", "product", "protein_id"]
			cds_info_file.write("\t".join(cds_header) + '\n')

			with open(self.genbank_matrix) as file:
				csv_reader = csv.DictReader(file, delimiter='\t')
				headers = [field for field in csv_reader.fieldnames if field != 'sequence']
				meta_file.write('\t'.join(headers) + '\n')

				for row in csv_reader:
					pop_sequence = row.pop('sequence', None)
					line = '\t'.join(row[field] for field in headers)
					meta_file.write(line + '\n')

					cds_info = self.parse_cds_info(row['cds_info'])
					gene_info = self.parse_gene_locations(row['genes'])

					for each_gene_info_block in gene_info:
						cds_info_file.write(row['primary_accession'] + '\t' +
						each_gene_info_block['feature'] + '\t' +
						each_gene_info_block['start'] + '\t' +
						each_gene_info_block['end'] + '\t' +
						each_gene_info_block['product'] + '\t' +
						each_gene_info_block['protein_id'] + '\n')

					for each_cds_info_block in cds_info:
						cds_info_file.write(row['primary_accession'] + '\t' +
						each_cds_info_block['feature'] + '\t' +
						each_cds_info_block['start'] + '\t' +
						each_cds_info_block['end'] + '\t' +
						each_cds_info_block['product'] + '\t' +
						each_cds_info_block['protein_id'] + '\n')

	@staticmethod
	def parse_cds_info(cds_info_string):
		cds_pattern = re.compile(r"cds_location:\s*<?(\d+)\.\.>?(\d+);")
		gene_pattern = re.compile(r"gene:\s*([^;]+);")
		protein_id_pattern = re.compile(r"protein_id:\s*([\w\.]+);")
		product_pattern = re.compile(r"product:\s*([^;]+);")

		cds_matches = cds_pattern.findall(cds_info_string)
		gene_matches = gene_pattern.findall(cds_info_string)	
		protein_id_matches = protein_id_pattern.findall(cds_info_string)
		product_matches = product_pattern.findall(cds_info_string)

		cds_data = []
		for i in range(len(cds_matches)):
			start, end = cds_matches[i]
			gene = gene_matches[i] if i < len(gene_matches) else ''
			protein_id = protein_id_matches[i] if i < len(protein_id_matches) else ''
			product = product_matches[i] if i < len(product_matches) else ''
			cds_data.append({
				'start': start,
				'end': end,
				'feature': 'CDS',
				'gene': gene,
				'protein_id': protein_id,
				'product': product
				})

		return cds_data

	@staticmethod
	def parse_gene_locations(gene_locations):
		pattern = re.compile(r"([A-Z])\((\d+)\.\.>?(\d+)\);?")
		matches = pattern.findall(gene_locations)
		blank = ''
		parsed_data = []
		for match in matches:
			gene, start, end = match
			parsed_data.append({'feature': 'gene',
				'start': start,
				'end': end,
				'product': gene,
				'protein_id': blank
			})
		return parsed_data

	@staticmethod
	def load_blast_hits(blast_hit_file):
		acc_dict = {}
		for i in open(blast_hit_file):
			query, ref, score, strand = i.strip().split('\t')
			acc_dict[query] = ref
		return acc_dict

	def created_alignment_table(self, blast_dict):
		header = ["sequence_id", "alignment_name", "alignment"]
		write_file = open(join(self.output_dir, "sequence.tsv"), 'w')
		write_file.write("\t".join(header) + "\n")
		rds = read_file.fasta(self.pad_aln)
		for rows in rds:
			if rows[0] in blast_dict:
				write_file.write(rows[0].strip() + '\t' +
					blast_dict[rows[0]].strip() + '\t' + rows[1] + '\n')

		write_file.close()

	def process(self):
		blast_dictionary = self.load_blast_hits(self.blast_hits)
		self.load_gb_matrix()
		self.created_alignment_table(blast_dictionary)

if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-g', '--genbank_matrix', help='Genbank matrix table', default="tmp/GenBank-matrix/gB_matrix.tsv")
	parser.add_argument('-t', '--tmp_dir', help='tmp directory to store all the db-ready tsv files', default="tmp/Tables/")
	parser.add_argument('-b', '--blast_hits', help='BLASTN unique hits', default="tmp/Blast/query_tophits_uniq.tsv")
	parser.add_argument('-p', '--pad_aln', help='Padded alignment file', default="tmp/Pad-Alignment/paded-alignment.fa")
	args = parser.parse_args()

	processor = SequenceProcessor(args.genbank_matrix, args.tmp_dir, args.blast_hits, args.pad_aln)
	processor.process()

