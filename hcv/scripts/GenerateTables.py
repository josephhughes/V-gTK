import os
import re
import csv
import time
import shutil
import read_file
import subprocess
import urllib.error
from Bio import Entrez
from Bio.Seq import Seq
from os.path import join
from itertools import islice
from argparse import ArgumentParser
from collections import defaultdict

class GenerateTables:
	def __init__(self, genbank_matrix, tmp_dir, blast_hits, paded_aln, host_taxa_file, nextalign_dir, email):
		self.genbank_matrix = genbank_matrix
		self.tmp_dir = tmp_dir
		self.blast_hits = blast_hits
		self.paded_aln = paded_aln
		self.output_dir = tmp_dir
		self.host_taxa_file = host_taxa_file
		self.nextalign_dir = nextalign_dir
		self.email = email
		os.makedirs(self.output_dir, exist_ok=True)

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
			with open(host_file, newline='') as f:
				reader = csv.DictReader(f, delimiter='\t', fieldnames=["taxonomy_id"])
				for row in reader:
					host_tax_id_list.append(row["taxonomy_id"])

		return host_tax_id_list

	def host_table(self):
		existing_host_taxa_id = self.host_taxa_file_check()
		taxa_file = join(self.output_dir, self.host_taxa_file)
		file_exists = os.path.isfile(taxa_file)
		write_file = open(taxa_file, 'a')
		if not file_exists or os.path.getsize(taxa_file) == 0:
				header = ['taxonomy_id', 'scientific_name', 'rank', 'lineage']
				write_file.write('\t'.join(header) + '\n')

		host_taxa_list = []

		with open(self.genbank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				if row['host_taxa_id'] not in host_taxa_list:
					if 'NA' not in row['host_taxa_id']:
						if row['host_taxa_id'] not in existing_host_taxa_id:
							host_taxa_list.append(row['host_taxa_id'])

		host_taxa_list = [item for item in host_taxa_list if item]
		host_taxa_list = [x for x in host_taxa_list if x != 'nan']
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

	@staticmethod
	def load_blast_hits(blast_hit_file):
		acc_dict = {}
		for i in open(blast_hit_file):
			query, ref, score, strand = i.strip().split('\t')
			acc_dict[query] = ref
		return acc_dict

	def created_alignment_table(self, blast_dict):
		accessions = {}
		header = ["sequence_id", "alignment_name", "alignment"]
		write_file = open(join(self.output_dir, "sequence_alignment.tsv"), 'w')
		write_file.write("\t".join(header) + "\n")
		rds = read_file.fasta(self.paded_aln)
		for rows in rds:
			if rows[0] not in accessions:
				accessions[rows[0]] = 1
				if rows[0] in blast_dict:
					write_file.write(rows[0].strip() + '\t' +
						blast_dict[rows[0]].strip() + '\t' + rows[1] + '\n')

		write_file.close()

	def create_sequence_table(self):
		shutil.copy(self.sequences, self.tmp_dir)

	def create_insertion_table(self):
		write_file = open(join(self.tmp_dir, "insertions.tsv"), 'w')
		header = ["accession", "reference", "insertion"]
		write_file.write("\t".join(header) + "\n")
		for aln_dir in os.listdir(self.nextalign_dir):
			for each_aln_dir in os.listdir(join(self.nextalign_dir, aln_dir)):
				with open(join(self.nextalign_dir, aln_dir, each_aln_dir, each_aln_dir + ".insertions.csv")) as f:
					for each_line in islice(f, 1, None):
						accession, insertion, aa_insertion = each_line.strip().split(",")
						if len(insertion) > 0:
							data = [accession, each_aln_dir, insertion]
							write_file.write("\t".join(data) + "\n")

		write_file.close()
					
	def process(self):
		blast_dictionary = self.load_blast_hits(self.blast_hits)
		#self.load_gb_matrix()
		self.created_alignment_table(blast_dictionary)
		self.host_table()
		self.create_insertion_table()

if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-g', '--genbank_matrix', help='Genbank matrix table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-t', '--tmp_dir', help='tmp directory to store all the db-ready tsv files', default="tmp/Tables/")
	parser.add_argument('-f', '--host_taxa', help='Host taxa file name, default="host_taxa.tsv"', default="host_taxa.tsv")
	parser.add_argument('-b', '--blast_hits', help='BLASTN unique hits', default="tmp/Blast/query_uniq_tophits.tsv")
	parser.add_argument('-p', '--paded_aln', help='Paded alignment file', default="tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta")
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign aligned directory', default="tmp/Nextalign/")
	parser.add_argument('-e', '--email', help='Email id', default='your-email@example.com')
	args = parser.parse_args()
	
	processor = GenerateTables(args.genbank_matrix, args.tmp_dir, args.blast_hits, args.paded_aln, args.host_taxa, args.nextalign_dir, args.email)
	processor.process()
