import os
import argparse
import subprocess
from os.path import join

class NCBI_GFF_Downloader:
	def __init__(self, accession_ids, tmp_dir):
		self.accession_ids = accession_ids
		self.tmp_dir = tmp_dir

	def download_gff(self):
		os.makedirs(join(self.tmp_dir), exist_ok=True)

		for each_accession in self.accession_ids.split(','):
			try:
				command = ["efetch", "-db", "nuccore", "-id", each_accession, "-format", "gff3", ">", join(self.tmp_dir, each_accession + ".gff3")]
				with open(join(self.tmp_dir, each_accession + ".gff3"), "w") as output:
					subprocess.run(" ".join(command), shell=True, check=True, stdout=output, stderr=subprocess.PIPE)
					print(f"Successfully downloaded GFF3 file: {each_accession}")
			except subprocess.CalledProcessError as e:
				print(f"Error downloading GFF3 file: {e.stderr.decode().strip()}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Download GFF3 file from NCBI using efetch")
	parser.add_argument("-id", "--accession_ids", required=True, help="NCBI accession ID")
	parser.add_argument("-o", "--tmp_dir", help="Directory where the GFF files are saved", default="tmp/Gff")
    	
	args = parser.parse_args()
	downloader = NCBI_GFF_Downloader(args.accession_ids, args.tmp_dir)
	downloader.download_gff()
