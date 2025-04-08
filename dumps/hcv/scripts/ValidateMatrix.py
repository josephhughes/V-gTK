import io
import os
import requests
import tarfile
import hashlib
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
        self.md5_url = f"{self.url}.md5"
        os.makedirs(self.tmp_dir, exist_ok=True)
        os.makedirs(self.taxa_path, exist_ok=True)

    def get_remote_md5(self):
        print(f"Fetching MD5 checksum from {self.md5_url}")
        response = requests.get(self.md5_url)
        response.raise_for_status()
        return response.text.split()[0]

    def get_local_md5(self, file_path):
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def download(self):
        print(f"Downloading {self.url}")
        response = requests.get(self.url, stream=True)
        response.raise_for_status()
        with open(join(self.taxa_path, self.dump_file), 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        print("Download complete.")

    def verify_and_download(self):
        local_file = join(self.taxa_path, self.dump_file)
        if os.path.exists(local_file):
            remote_md5 = self.get_remote_md5()
            local_md5 = self.get_local_md5(local_file)
            if local_md5 == remote_md5:
                print("Checksum matches. Using existing taxadump file.")
                return
            else:
                print("Checksum mismatch. Downloading new taxadump file.")
        else:
            print("Local taxadump file not found. Downloading...")
        self.download()

    def read_tar(self):
        print(f"Extracting {self.dump_file}")
        output_file = join(self.taxa_path, "taxa-name-dump.txt")
        with open(output_file, 'w') as write_taxa, tarfile.open(join(self.taxa_path, self.dump_file), 'r') as tar:
            names_dmp = tar.extractfile('names.dmp')
            if names_dmp:
                for line in names_dmp:
                    split_line = line.decode('utf-8').strip().split('|')
                    write_taxa.write("\t".join(split_line) + '\n')
        print("Extraction complete.")

    def taxa_name_dump_to_dict(self):
        taxa_dump = {}
        file_name = join(self.taxa_path, "taxa-name-dump.txt")
        with open(file_name) as f:
            for each_line in f:
                split_line = each_line.strip().split('\t')
                if 'genbank common name' in split_line or 'common name' in split_line or 'scientific name' in split_line:
                    taxa_id = split_line[0]
                    name = split_line[3]
                    taxa_dump.setdefault(name, []).append(taxa_id)
        return taxa_dump

    def validate_date(self, date):
        if pd.isna(date):
            return date
        for fmt in ('%d-%b-%Y', '%b-%Y', '%Y'):
            try:
                return "Yes"
            except ValueError:
                pass
        return None

    def country_to_dict(self, infile):
        country_dict = {}
        with open(infile) as f:
            for i in islice(f, 1, None):
                split_line = i.strip().split(',')
                country_dict[split_line[1]] = i.strip()
        return country_dict

    def validate_country(self, country, country_dict):
        if pd.isna(country):
            return "NA"
        country = country.split(':')[0]
        return "Yes" if country in country_dict else "NA"

    def validate_host(self, host_name, taxa_dict):
        return "Yes" if not pd.isna(host_name) and host_name in taxa_dict else "NA"

    def read_meta_sheet(self, country_dict, taxa_dict):
        print(f"Reading meta data {self.gb_matrix}")
        df = pd.read_csv(self.gb_matrix, sep='\t')
        gb_matrix_df = df.copy()
        df['collection_date_validated'] = df['collection_date'].apply(self.validate_date)
        df['country_validated'] = df['country'].apply(lambda x: self.validate_country(x, country_dict))
        df['host_validated'] = df['host'].apply(lambda x: self.validate_host(x, taxa_dict))
        #df['host_taxa_id'] = df['host'].apply(lambda x: taxa_dict.get(x, ['NA'])[0])
        df['host_taxa_id'] = df['host'].apply(lambda x: str(taxa_dict.get(x, ['NA'])[0]))

        failed_df = df[df[['collection_date_validated', 'host_validated', 'country_validated']].isnull().any(axis=1)]
        failed_df = failed_df[['primary_accession', 'collection_date_validated', 'host_validated', 'country_validated']]
        failed_df = failed_df.fillna("NA")
        failed_df = failed_df.merge(df[['primary_accession', 'host', 'country', 'collection_date']], on='primary_accession', how='left')
        
        df.to_csv(join(self.tmp_dir, 'gB_matrix_validated.tsv'), sep='\t', index=False)
        failed_df.to_csv(join(self.tmp_dir, 'gB_matrix_failed_validation.tsv'), sep='\t', index=False)
				
        merge_df = gb_matrix_df.merge(df[['gi_number', 'collection_date_validated', 'country_validated', 'host_validated', 'host_taxa_id']],on='gi_number',how='right')
        merge_df.to_csv(self.gb_matrix, sep="\t", index=False)
        print("######---Validation summary---######")
        print(f"Total accession: {len(df)}")
        print(f"Validated: {len(df) - len(failed_df)}")
        print(f"Missing information: {len(failed_df)}")
        print(f"Results saved to {self.tmp_dir}")

    def process(self):
        country_dict = self.country_to_dict(self.country)
        self.verify_and_download()
        self.read_tar()
        taxa_dict = self.taxa_name_dump_to_dict()
        self.read_meta_sheet(country_dict, taxa_dict)

if __name__ == "__main__":
    parser = ArgumentParser(description='Validate the gB_matrix file based on country, date, and host columns.')
    parser.add_argument('-u', '--url', help='URL to download taxa file', default="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
    parser.add_argument('-t', '--taxa_path', help='Path to save taxa dump names', default='tmp/taxa')
    parser.add_argument('-d', '--tmp_dir', help='Path to save validated gB_matrix', default='tmp/Validate-matrix')
    parser.add_argument('-g', '--gb_matrix', help='Genbank matrix file', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
    parser.add_argument('-c', '--country', help='m49 country file', default='assets/m49_country.csv')
    args = parser.parse_args()
    
    validator = ValidateMatrix(args.url, args.taxa_path, args.tmp_dir, args.gb_matrix, args.country)
    validator.process()

