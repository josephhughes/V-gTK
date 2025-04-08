import os
import re
import csv
import time
import read_file
import subprocess
import urllib.error
from Bio import Entrez
from os.path import join
from itertools import islice
from argparse import ArgumentParser
from collections import defaultdict
from FeatureCalculator  import FeatureCordCalculator 
#from AlignmentCordCalc import AlignmentCordCalculator
from tenacity import retry, stop_after_attempt, wait_exponential

class SequenceProcessor:
	def __init__(self, genbank_matrix, tmp_dir, blast_hits, query_aln, host_taxa_file, reference_aln,email, reference_feature_op, query_feature_op, gaps_to_ignore):
		self.genbank_matrix = genbank_matrix
		self.tmp_dir = tmp_dir
		self.blast_hits = blast_hits
		self.query_aln = query_aln
		self.output_dir = tmp_dir
		self.host_taxa_file = host_taxa_file
		self.reference_aln = reference_aln
		self.email = email
		self.reference_feature_op = reference_feature_op
		self.query_feature_op = query_feature_op
		self.gaps_to_ignore = gaps_to_ignore
		os.makedirs(self.output_dir, exist_ok=True)
		self.feature_cords = FeatureCordCalculator()
		#self.feature_cords = AlignmentCordCalculator

	@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=10), retry=(lambda e: isinstance(e, urllib.error.HTTPError)))
	def fetch_taxonomy_details(self, tax_id):
		Entrez.email = self.email
		try:
			handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
			records = Entrez.read(handle)
			time.sleep(1)
			handle.close()

			if records:
				tax_record = records[0]
				taxonomy_info = {
					"Scientific Name": tax_record.get("ScientificName", "N/A"),
					"Taxonomy ID": tax_record.get("TaxId", "N/A"),
					"Rank": tax_record.get("Rank", "N/A"),
					"Lineage": tax_record.get("Lineage", "N/A"),
					"Other Names": tax_record.get("OtherNames", {}).get("Synonym", []),
					}
				return taxonomy_info
			else:
				return "No taxonomy details found for the given ID."

		except urllib.error.HTTPError as e:
			print(f"HTTPError: {e}. Retrying...")
			raise

	def host_taxa_file_check(self):
		host_tax_id_list = []
		host_file = join(self.tmp_dir, self.host_taxa_file)
		if os.path.exists(host_file):
			with open(host_file) as f:
				for each_line in f:
					tax_id, scientific_name, rank, lineage = each_line.split('\t')
					host_tax_id_list.append(tax_id)
		
		return host_tax_id_list

	def host_table(self):
		existing_host_taxa_id = self.host_taxa_file_check()
		taxa_file = join(self.output_dir, self.host_taxa_file)
		write_file = open(taxa_file, 'a')
		header = ['scientific_name', 'taxonomy_id', 'rank', 'lineage']
		write_file.write('\t'.join(header))
		write_file.write('\n')
		host_taxa_list = []

		with open(self.genbank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				if row['host_taxa_id'] not in host_taxa_list:
					if 'NA' not in row['host_taxa_id']:
						if row['host_taxa_id'] not in existing_host_taxa_id:
							host_taxa_list.append(row['host_taxa_id'])
	
		print(f"\nTotal {len(host_taxa_list)} taxa id informations to fetch\n")	
		for idx, each_host_taxaid in enumerate(host_taxa_list, start=1):
			print(f" - Fetching taxa details for {each_host_taxaid} which is {idx} of {len(host_taxa_list)}")
			host_taxa = self.fetch_taxonomy_details(each_host_taxaid)
			tax_id = host_taxa['Taxonomy ID'] if host_taxa['Taxonomy ID'] is not None else 'NA'
			scientific_name = host_taxa['Scientific Name'] if host_taxa['Scientific Name'] is not None else 'NA'
			rank = host_taxa['Rank']	if host_taxa['Rank'] is not None else 'NA'
			lineage = host_taxa['Lineage'] if host_taxa['Lineage'] is not None else 'NA'
			write_file.write(tax_id + '\t' + scientific_name + '\t' + rank + '\t' + lineage + '\n')
		write_file.close() 
	
	def load_gb_matrix(self):
		cds_keys = ['start', 'end', 'gene', 'protein_id', 'product']
		meta_file_path = join(self.output_dir, "meta_data.tsv")
		cds_info_file_path = join(self.output_dir, "gB_features.tsv")

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

	def create_gff_dict(gff_file):
		gff_dict = {}
		with open(gff_file) as f:
			for each_line in f:
				if not each_line.startswith('#'):
					split_line = each_line.strip().split('\t')
					if len(split_line) >= 5:
						feature = split_line[2]
						start = split_line[3]
						end = split_line[4]
						if feature not in gff_dict:
							gff_dict[feature] = []
							gff_dict[feature].append({'start': start, 'end': end})
						else:
							gff_dict[feature].append({'start': start, 'end': end})
		return gff_dict

	def query_feature_table(self):
		write_file = open(join(self.tmp_dir, self.query_feature_op), 'w')
		header = ['query_acc', 'start', 'end']
		write_file.write("\t".join(header))
		write_file.write("\n")
		sequence_object = read_file.fasta(self.query_aln)

		for row in sequence_object:
			content = self.feature_cords.calculate_alignment_coords(row[0], row[1], self.gaps_to_ignore)
			for each_cords in content['aligned']:
				start, end = each_cords[0], each_cords[1]
				data = [row[0], str(start), str(end)]	
				write_file.write('\t'.join(data))
				write_file.write("\n")
		write_file.close()

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
		rds = read_file.fasta(self.query_aln)
		for rows in rds:
			print(rows[0])
			if rows[0] in blast_dict:
				write_file.write(rows[0].strip() + '\t' +
					blast_dict[rows[0]].strip() + '\t' + rows[1] + '\n')

		write_file.close()

	def process(self):
		blast_dictionary = self.load_blast_hits(self.blast_hits)
		self.load_gb_matrix()
		self.created_alignment_table(blast_dictionary)
		self.host_table()
		self.query_feature_table()
# create the nextalign master alignment, add the function to create master alignment in Alignment.py tmp/Nextalign/

if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-g', '--genbank_matrix', help='Genbank matrix table', default="tmp/Validate-matrix/gB_matrix_validated.tsv")
	parser.add_argument('-t', '--tmp_dir', help='tmp directory to store all the db-ready tsv files', default="tmp/Tables/")
	parser.add_argument('-f', '--host_taxa', help='Host taxa file name, default="host_taxa.tsv"', default="host_taxa.tsv")
	parser.add_argument('-b', '--blast_hits', help='BLASTN unique hits', default="tmp/Blast/query_uniq_tophits.tsv")
	parser.add_argument('-p', '--query_aln', help='Padded alignment file', default="tmp/Pad-Alignment/paded-query-alignment.fa")
	parser.add_argument('-m', '--reference_aln', help='Master alingment file', default="tmp/Pad-Alignment/paded-reference-alignment.fa")
	parser.add_argument('-e', '--email', help='Email id', default='your-email@example.com')
	parser.add_argument('-rf', '--reference_feature_op', help='Output reference feature table file name', default='reference_features.tsv')
	parser.add_argument('-qf', '--query_feature_op', help='Output query feature table file name', default='query_features.tsv') 
	parser.add_argument('-gp', '--gaps_to_ignore', help='Lenght of gaps to ignore', default=30)
	args = parser.parse_args()

	try:
		processor = SequenceProcessor(args.genbank_matrix, args.tmp_dir, args.blast_hits, args.query_aln, args.host_taxa, args.reference_aln, args.email, args.reference_feature_op, args.query_feature_op, args.gaps_to_ignore)
		processor.process()
	except urllib.error.HTTPError as e:
		print(f"HTTPError: {e}. Retrying...")
