#Validates the TSV file information including country, host and date
#Validate date, country, host and add species scientific name and taxa id to the validated sheet
import io
import os
import requests
import tarfile
import pandas as pd
from os.path import join
from itertools import islice
from datetime import datetime
from argparse import ArgumentParser

dump_file = "taxdump.tar.gz"
# Function to download the file and read names.dmp
def download(url, output_path):
	print(f"Downloading {url}/{dump_file}")
	os.makedirs(output_path, exist_ok=True)
	response = requests.get(url)
	response.raise_for_status()
	with open(join(output_path, dump_file), 'wb') as file:
		file.write(response.content)
	#tar_data = io.BytesIO(response.content)

def read_tar(file_path, taxa_path):
	print(f"reading file {file_path}/{dump_file}")
	write_taxa = open(join(taxa_path, "taxa-name-dump.txt"), 'w')
	with tarfile.open(join(file_path, dump_file), 'r') as tar:
		names_dmp = tar.extractfile('names.dmp')
		if names_dmp:
			for line in names_dmp:
				split_line = line.decode('utf-8').split('|')
				write_taxa.write(split_line[0].strip() + '\t' + split_line[1].strip() + '\t' + split_line[2].strip() + "\n")
		else:
			print("names.dmp not found in the archive.")
	print(f"taxa name dump written to {taxa_path}/taxa-name-dump.txt")
	write_taxa.close()

def taxa_name_dump_to_dict(file_path):
	taxa_dump	= {}
	print("reading_file taxa-dump")
	file_name = join(file_path, "taxa-name-dump.txt")
	with open(file_name) as f:
		for each_line in f:
			split_line = each_line.strip().split('\t')
			species_name = split_line[1] + "|" + split_line[2] if len(split_line) > 2 else split_line[1]
			if split_line[0] not in taxa_dump:
				taxa_dump[split_line[0]] = [species_name]
			else: 
				taxa_dump[split_line[0]].append(species_name)

	return taxa_dump

def validate_date(date):
	if pd.isna(date):
		return date
	try:
		for fmt in ('%d-%b-%Y', '%b-%Y', '%Y'):
			try:
				return int(datetime.strptime(date, fmt).year)
			except ValueError:
				pass
	except Exception as e:
		return np.nan

def country_to_dict(infile):
	country_dict = {}
	with open(infile) as f:
		for i in islice(f, 1, None):
			split_line = i.strip().split('\t')
			country = split_line[1]
			country_dict[country] = i.strip()
	return country_dict
			
#taxa_dicts = taxa_name_dump_to_dict("tmp")
def validate_host(host):
	host = host.strip()
	if pd.isna(host):
		return host
	try:
		if taxa_dicts[host]:
				return "yes"
		elif host in taxa_dicts.keys():
			return "yes"
		else:
			return "no"
	except Exception as e:
		pass 

def validate_country(country, country_dict):
	if pd.isna(country):
		return ["NA", "NA", "NA"]
	else:
		country = country.split(':')[0]
		if country in country_dict:
			#alpha3_code, country, region, sub_region, intermediate_sub_region, ldc, lldc, sids, developed_developing_countries 
			split_line = country_dict[country].split('\t')
			alpha3_code = split_line[0]
			region = split_line[3]
			sub_region = split_line[3]
			return [alpha3_code, region, sub_region]
		else:
			return ["NA", "NA", "NA"]

def read_meta_sheet(infile, country_dict):
	df = pd.read_csv(infile, sep='\t')
	df['Collection Date Validated'] = df['Collection Date'].apply(validate_date)
	#df['Host Validated'] = df['Host'].apply(validate_host) 
	#df['Collection Date Formatted'] = [int(year) if not pd.isna(year) else year for year in df['Collection Date Formatted'].tolist()]
	#print (df['Host Validated'].tolist())
	country_details = df['Country'].apply(lambda x: validate_country(x, country_dict))
	df[['alpha3_code', 'region', 'sub_region']] = pd.DataFrame(country_details.tolist(), index=df.index)
	filtered_df = df[df[['Collection Date Validated', 'alpha3_code', 'region', 'sub_region']].isnull().any(axis=1)]
	filtered_df = filtered_df[['Locus', 'Primary Accession', 'Collection Date Validated', 'alpha3_code', 'region', 'sub_region']]
	filtered_df = filtered_df.fillna("NA")
	validated_df = df
	validated_df.to_csv('tmp/gB_matrix_validated.tsv', sep='\t', index=False)
	filtered_df.to_csv('missing_data.tsv', sep='\t', index=False)
		
def process(args):
	country_dict = country_to_dict('m49_countries.txt')
	read_meta_sheet(args.gb_matrix, country_dict)
	if os.path.exists(join(args.path, dump_file)):
		pass
		#read_tar(args.path, args.taxa_path)
		#	taxa_name_dump_to_dict(args.taxa_path)
	else:
		pass
		#download(args.url, args.path)
		#read_tar(args.path, args.taxa_path)
		#taxa_name_dump_to_dict(args.taxa_path)

if __name__ == "__main__":
	parser = ArgumentParser(description='Download and update GenBank XML files for a given species')
	parser.add_argument('-u', '--url', help='URL to download taxa file', default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
	parser.add_argument('-p', '--path', help='Path to download taxa file, generally "tmp" prefered ', default='tmp')
	parser.add_argument('-t', '--taxa_path', help='Path to save taxa dump names, generally "tmp" prefered ', default='tmp')
	parser.add_argument('-g', '--gb_matrix', help='Genbank matrix sheet generated by genbank_to_tsv.py', default="tmp/gB_matrix.tsv")
	#parser.add_argument('-e', '--email', help='email id', default='your_email@example.com')
	args = parser.parse_args()
	process(args)
