#conda install python
#conda install requests
# check for max number of records and add that in the max_rets

import os
import requests
from time import sleep
from argparse import ArgumentParser
from os.path import join

# Fetch accession ID's
def fetch_ids(search_term, base_url, email):
	retmax = get_record_count(search_term, email)
	search_url = base_url + "esearch.fcgi?db=nucleotide&term=txid" + search_term + "[Organism:exp]&retmax=" + str(retmax) + "&usehistory=y&email=" + email + "&retmode=json"
	print (search_url)
	response = requests.get(search_url)
	response.raise_for_status()
	data = response.json()
	ids = data["esearchresult"]["idlist"]
	return ids

def get_record_count(search_term, email):
	base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
	search_url = f"{base_url}esearch.fcgi?db=nucleotide&term=txid{search_term}[Organism:exp]&retmode=json&email={email}"
	
	response = requests.get(search_url)
	data = response.json()
	
	count = int(data['esearchresult']['count'])
	return count

# Fetch GenBank XML data in batches
def fetch_genbank_data(ids, batch_size, base_url, email, output_dir):
	batch_count = batch_size
	for start in range(0, len(ids), batch_size):
		batch_ids = ids[start:start + batch_size]
		ids_str = ",".join(batch_ids)
		fetch_url = base_url + "efetch.fcgi?db=nucleotide&id=" + ids_str + "&retmode=xml&email=" + email
		response = requests.get(fetch_url)
		response.raise_for_status()
		save_data(response.text, batch_count, output_dir)
		batch_count = batch_count + batch_size
		sleep(2) 

# Save XML to file as batch
def save_data(data, batch_size, output_dir):
	os.makedirs(output_dir, exist_ok=True)
	filename = join(output_dir, "batch-" + str(batch_size) + ".xml")

	with open(filename, "w") as file:
		file.write(data)
	print(f"Data written to: "+ filename)

# Download from scratch
def download(args):
	ids = fetch_ids(args.search_term, args.base_url, args.email)
	print("Found " + str(len(ids)) + " IDs")
	fetch_genbank_data(ids, args.batch_size, args.base_url,args.email, args.output_dir)

if __name__ == "__main__":
	parser = ArgumentParser(description='Download and update GenBank XML files for a given species')
	parser.add_argument('-s', '--search_term', help='Search term, generally a taxaid example: 11292', required=True)
	parser.add_argument('-o', '--output_dir', help='Output directory where all the XML files are stored', default='GenBank')
	#parser.add_argument('-r', '--max_ret', help='Max number of records to be retrieved. Default is 100000', default=100000, type=int)
	parser.add_argument('-b', '--batch_size', help='Max number of XML files to pull and merge in to a single file', default=100, type=int)
	parser.add_argument('--update', action='store_true', help='Run the script in update mode to download only new sequences')
	parser.add_argument('-u', '--base_url', help='base url to downlaod the XML files', default='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/')
	parser.add_argument('-e', '--email', help='email id', default='your_email@example.com')
	args = parser.parse_args()
	if args.update:
		update(args)
	else:
		download(args)
