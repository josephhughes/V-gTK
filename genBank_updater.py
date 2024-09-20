#conda install python
#conda install requests
# check for max number of records and add that in the max_rets (completed)



import os
import requests
import csv
from time import sleep
from argparse import ArgumentParser
from os.path import join

# Fetch accession ID's
def fetch_ids(search_term, base_url, email):
	retmax = get_record_count(search_term, email)
	search_url = base_url + "esearch.fcgi?db=nucleotide&term=txid" + search_term + "[Organism:exp]&retmax=" + str(retmax) + "&usehistory=y&email=" + email + "&retmode=json"
	#print (search_url)
	response = requests.get(search_url)
	response.raise_for_status()
	data = response.json()
	ids = data["esearchresult"]["idlist"]
	return ids

def fetch_accs(search_term, base_url, email):
	retmax = get_record_count(search_term, email)
	search_url = base_url + "esearch.fcgi?db=nucleotide&term=txid" + search_term + "[Organism:exp]&retmax=" + str(retmax) + "&idtype=acc&usehistory=y&email=" + email + "&retmode=json"
	#print (search_url)
	response = requests.get(search_url)
	response.raise_for_status()
	data = response.json()
	accs = data["esearchresult"]["idlist"]
	return accs


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
		
	base_filename = join(output_dir, "batch-" + str(batch_size))
	filename = base_filename + ".xml"
	counter = 1
	while os.path.exists(filename):
		filename = f"{base_filename}_{counter}.xml"
		counter += 1

	with open(filename, "w") as file:
		file.write(data)
	print(f"Data written to: "+ filename)

# Download from scratch
def download(args):
	ids = fetch_ids(args.taxid, args.base_url, args.email)
	print("Found " + str(len(ids)) + " IDs")
	fetch_genbank_data(ids, args.batch_size, args.base_url,args.email, args.output_dir)

def update(args):
	# check if the output directory exists and throw an error if it does not
	if os.path.exists(args.output_dir):
		raise ValueError("Output directory exists, it is recomended to create a new dated directory for the updated files")
	with open(args.update, 'r') as file:
		reader = csv.DictReader(file, delimiter='\t')
		if 'Accession Version' not in reader.fieldnames:
			raise ValueError("Expecting a column called 'Accession Version'")
		accession_versions = [row['Accession Version'] for row in reader]
		print("Found " + str(len(accession_versions)) + " Accession Versions")
		ids = fetch_accs(args.taxid, args.base_url, args.email)
		missing_ids = [id for id in ids if id not in accession_versions]
		print("Found " + str(len(missing_ids)) + " missing IDs")
		fetch_genbank_data(missing_ids, args.batch_size, args.base_url, args.email, args.output_dir)



if __name__ == "__main__":
	parser = ArgumentParser(description='Download and update GenBank XML files for a given species')
	parser.add_argument('-t', '--taxid', help='taxid example: 11292', required=True)
	parser.add_argument('-o', '--output_dir', help='Output directory where all the XML files are stored', default='GenBank')
	#parser.add_argument('-r', '--max_ret', help='Max number of records to be retrieved. Default is 100000', default=100000, type=int)
	parser.add_argument('-b', '--batch_size', help='Max number of XML files to pull and merge in to a single file', default=100, type=int)
	parser.add_argument('--update', help='Run the script in update mode to download only new sequences, expecting a tsv with a column called Accession Version',default="gB_matrix.tsv")
	parser.add_argument('-u', '--base_url', help='base url to downlaod the XML files', default='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/')
	parser.add_argument('-e', '--email', help='email id', default='your_email@example.com')
	args = parser.parse_args()
	if args.update:
		update(args)
	else:
		download(args)
