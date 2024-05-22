import os
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from os.path import join as join

Entrez.email = "your_email@example.com"

#Fetch GenBank id's
def fetch_rabies_genbank_ids(search_term, max_ret):
	search_handle = Entrez.esearch(
		db="nucleotide", 
		term=search_term, 
		retmax=max_ret,
		retmode="xml")
	search_results = Entrez.read(search_handle)
	search_handle.close()
	return search_results['IdList']

# Extract accession version from GenBank XML data
def extract_accession_version(xml_data):
	root = ET.fromstring(xml_data)
	accession_version = root.findtext(".//GBSeq_accession-version")
	return accession_version

# Download GenBank XML for given IDs and save with accession version filename
def download_genbank_xml(ids, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for genbank_id in ids:
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            id=genbank_id,
            rettype="gb",
            retmode="xml"
        )
        xml_data = fetch_handle.read()
        fetch_handle.close()

        accession_version = extract_accession_version(xml_data)
        if accession_version:
            filename = join(output_dir,accession_version + ".xml")
        else:
            filename = join(output_dir, genbank_id + ".xml")

        with open(filename, "w") as file:
            file.write(xml_data.decode('utf-8'))

def run(args):
	print("Fetching GenBank IDs for: " + args.search_term)
	genbank_ids = fetch_rabies_genbank_ids(args.search_term, args.max_ret)
	id_count = len(genbank_ids)
	print("Found " + str(id_count) + " records. Downloading GenBank XML files...")
	download_genbank_xml(genbank_ids, args.output_dir)
	print("Download completed.")
	
# Main script execution
if __name__ == "__main__":
	parser = ArgumentParser(description='Download GenBank XML file for a given species')
	parser.add_argument('-s', '--search_term', help='Search term, example: rabies virus[Organism]', default='rabies virus[Organism]')
	parser.add_argument('-o', '--output_dir', help='Output directory where all the XAL files are stored', default='tmp')
	parser.add_argument('-r', '--max_ret', help='Max number of records to be retrived. Default is 100000', default=100000)
	args = parser.parse_args()
	run(args)

