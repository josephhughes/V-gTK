import io
import os
import requests
import tarfile
import pandas as pd
from os.path import join
from itertools import islice
from datetime import datetime
from argparse import ArgumentParser

class ValidateMatrix:
	def __init__(self, url, taxa_path, tmp_dir, gb_matrix, country):
		self.url = url
		self.taxa_path = taxa_path
		self.tmp_dir = tmp_dir
		self.gb_matrix = gb_matrix
		self.country = country
		self.dump_file = "taxadump.tar.gz"
		self.count = 0
		self.count_set = 100
		os.makedirs(self.tmp_dir, exist_ok=True)
		os.makedirs(self.taxa_path, exist_ok=True)

	def download(self):
		print(f"Downloading {self.url}/{self.dump_file}")
		response = requests.get(self.url)
		response.raise_for_status()
		with open(join(self.taxa_path, self.dump_file), 'wb') as file:
			file.write(response.content)
			print("Download complete.")

	def read_tar(self):
		print(f"Reading file {self.taxa_path}/{self.dump_file}")
		output_file = join(self.taxa_path, "taxa-name-dump.txt")
		with open(output_file, 'w') as write_taxa, tarfile.open(join(self.taxa_path, self.dump_file), 'r') as tar:
			names_dmp = tar.extractfile('names.dmp')
			if names_dmp:
				for line in names_dmp:
					split_line = line.decode('utf-8').strip().split('|')
					#write_taxa.write(split_line[0].strip() + '\t' + split_line[1].strip() + '\t' + split_line[2].strip() + "\n")
					write_taxa.write("\t".join(split_line) + '\n')
		print("Extraction complete.")

	def taxa_name_dump_to_dict(self):
		taxa_dump = {}
		print(f"Reading {self.taxa_path}/taxa-name-dump.txt")
		file_name = join(self.taxa_path, "taxa-name-dump.txt")
		with open(file_name) as f:
			for each_line in f:
				split_line = each_line.strip().split('\t')
				if 'genbank common name' in split_line or 'common name' in split_line or 'scientific name' in split_line:
					taxa_id = split_line[0]
					name = split_line[3]
					if name not in taxa_dump:
						taxa_dump[name] = [taxa_id]
					else:
						taxa_dump[name].append(taxa_id)
		return taxa_dump

	@staticmethod
	def validate_date(date):
		if pd.isna(date):
			return date
		for fmt in ('%d-%b-%Y', '%b-%Y', '%Y'):
			try:
				return "Yes" #int(datetime.strptime(date, fmt).year)
			except ValueError:
				pass
		return None

	@staticmethod
	def country_to_dict(infile):
		country_dict = {}
		with open(infile) as f:
			for i in islice(f, 1, None):
				split_line = i.strip().split(',')
				country = split_line[1]
				country_dict[country] = i.strip()
		return country_dict

	def validate_country(self, country, country_dict):
		if pd.isna(country):
			return "NA"
		country = country.split(':')[0]
		if country in country_dict:
			split_line = country_dict[country].split('\t')
			alpha3_code = split_line[0]
			region = split_line[3]
			sub_region = split_line[3]
			return "Yes"
		else:
			return "NA"

	def validate_host(self, host_name, taxa_dict):
		if pd.isna(host_name):
			return "NA"
		if host_name in taxa_dict.keys():
			return "Yes"
		else:
			pass #print(f"{host_name}")

	def read_meta_sheet(self, country_dict, taxa_dict):
		print(f"Taxa dictionary has {len(taxa_dict.keys())} taxa names") 
		print(f"Reading meta data {self.gb_matrix}")
		df = pd.read_csv(self.gb_matrix, sep='\t')
		df['collection_date_validated'] = df['collection_date'].apply(lambda x: self.validate_date(x))
		df['country_validated'] = df['country'].apply(lambda x: self.validate_country(x, country_dict))
		df['host_validated'] = df['host'].apply(lambda x: self.validate_host(x, taxa_dict))
		filtered_df = df[df[['collection_date_validated', 'host_validated', 'country_validated']].isnull().any(axis=1)]
		filtered_df = filtered_df[['primary_accession', 'collection_date_validated', 'host_validated', 'country_validated']]

		filtered_df = filtered_df.fillna("NA")
		failed_df = filtered_df.merge(df[['primary_accession', 'host', 'country', 'collection_date']],
                   on='primary_accession', how='left')
		validated_df = df
		validated_df.to_csv(join(self.tmp_dir, 'gB_matrix_validated.tsv'), sep='\t', index=False)
		failed_df.to_csv(join(self.tmp_dir, 'gB_matrix_failed_validation.tsv'), sep='\t', index=False)
		print("######---Validation summary---######")
		print(f"Total accession: {len(df)}")
		diff = len(df) - len(filtered_df)
		print(f"Overall {diff}  accessions were validated")
		print(f"Missing information: {len(filtered_df)}")
		print(f"Validated and failed files saved to {self.tmp_dir}")
		#print("Reading meta data complete.")

	def process(self):
		country_dict = self.country_to_dict(self.country)
		if os.path.exists(join(self.taxa_path, self.dump_file)):
			self.read_tar()
		else:
			self.download()
			self.read_tar()
		taxa_dict = self.taxa_name_dump_to_dict()
		#for k, v in taxa_dict.items():
		#	print(k, v)
		self.read_meta_sheet(country_dict, taxa_dict)
		#print("Validation complete")

if __name__ == "__main__":
	parser = ArgumentParser(description='Validated the gB_matrix file based on the country, date and host columns. The validation process involves in making sure the country is matched with m_49 list, date is in standard format, and host name is valid as per the NCBI taxa names')
	parser.add_argument('-u', '--url', help='URL to download taxa file', default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
	parser.add_argument('-t', '--taxa_path', help='Path to save taxa dump names, generally "tmp" preferred ', default='tmp/taxa')
	parser.add_argument('-d', '--tmp_dir', help='Path to save validated gB_matrix, generally "tmp/Validate-matrix" preferred ', default='tmp/Validate-matrix')
	parser.add_argument('-g', '--gb_matrix', help='Genbank matrix sheet generated by genbank_to_tsv.py', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-c', '--country', help='m49 country file', default='assets/m49_country.csv')
	args = parser.parse_args()
	
	validator = ValidateMatrix(args.url, args.taxa_path, args.tmp_dir, args.gb_matrix, args.country)
	validator.process()

