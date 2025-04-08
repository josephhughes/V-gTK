import os
import csv
import requests
from time import sleep
from os.path import join
from argparse import ArgumentParser

class GenBankFetcher:
	def __init__(self, taxid, base_url, email, output_dir, batch_size, sleep_time, base_dir, update_file):
		self.taxid = taxid
		self.base_url = base_url
		self.email = email
		self.output_dir = output_dir
		self.batch_size = batch_size
		self.sleep_time = sleep_time
		self.base_dir = base_dir
		self.update_file = update_file

	def get_record_count(self):
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmode=json&email={self.email}"
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return int(data['esearchresult']['count'])

	def fetch_ids(self):
		retmax = self.get_record_count()
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmax={retmax}&usehistory=y&email={self.email}&retmode=json"
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return data["esearchresult"]["idlist"]

	def fetch_accs(self):
		retmax = self.get_record_count()
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmax={retmax}&idtype=acc&usehistory=y&email={self.email}&retmode=json"
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return data["esearchresult"]["idlist"]

	def fetch_genbank_data(self, ids):
		batch_count = self.batch_size
		for start in range(0, len(ids), self.batch_size):
			batch_ids = ids[start:start + self.batch_size]
			ids_str = ",".join(batch_ids)
			fetch_url = f"{self.base_url}efetch.fcgi?db=nucleotide&id={ids_str}&retmode=xml&email={self.email}"
			response = requests.get(fetch_url)
			response.raise_for_status()
			self.save_data(response.text, batch_count)
			batch_count += self.batch_size
			sleep(self.sleep_time)

	def save_data(self, data, batch_size):
		os.makedirs(self.output_dir, exist_ok=True)
		os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

		base_filename = join(self.output_dir, self.base_dir, f"batch-{batch_size}")
		filename = f"{base_filename}.xml"
		counter = 1
		while os.path.exists(filename):
			filename = f"{base_filename}_{counter}.xml"
			counter += 1

		with open(filename, "w") as file:
				file.write(data)
		print(f"Data written to: {filename}")

	def download(self):
		ids = self.fetch_ids()
		print(f"Found {len(ids)} IDs")
		self.fetch_genbank_data(ids)

	def update(self, update_file):
		with open(update_file, 'r') as file:
			reader = csv.DictReader(file, delimiter='\t')
			if 'accession_version' not in reader.fieldnames:
					raise ValueError("Expecting a column called 'accession_version'")
			accession_versions = [row['accession_version'] for row in reader]
			print(f"Found {len(accession_versions)} accession_versions")
			ids = self.fetch_accs()
			missing_ids = [id for id in ids if id not in accession_versions]
			print(f"Found {len(missing_ids)} missing IDs")
			self.fetch_genbank_data(missing_ids)

if __name__ == "__main__":
	parser = ArgumentParser(description='Download and update GenBank XML files for a given species')
	parser.add_argument('-t', '--taxid', help='TaxID example: 11292', required=True)
	parser.add_argument('-o', '--tmp_dir', help='Output directory where all the XML files are stored', default='tmp')
	parser.add_argument('-b', '--batch_size', help='Max number of XML files to pull and merge in a single file', default=100, type=int)
	parser.add_argument('--update', help='Run the script in update mode to download only new sequences, expecting a TSV with a column called Accession Version')
	parser.add_argument('-u', '--base_url', help='Base URL to download the XML files', default='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/')
	parser.add_argument('-e', '--email', help='Email ID', default='your_email@example.com')
	parser.add_argument('-s', '--sleep_time', help='Delay after each set of information fetch', default=2, type=int)
	parser.add_argument('-d', '--base_dir', help='Directory where all the XML files are stored', default='GenBank-XML')
	args = parser.parse_args()

	fetcher = GenBankFetcher(
		taxid=args.taxid,
		base_url=args.base_url,
		email=args.email,
		output_dir=args.tmp_dir,
		batch_size=args.batch_size,
		sleep_time=args.sleep_time,
		base_dir = args.base_dir,
		update_file = args.update
	)

	if args.update:
		fetcher.update(args.update)
	else:
		fetcher.download()

