import os
import csv
import subprocess
from os.path import join
from argparse import ArgumentParser

def read_ref_list(ref_list_file):
	ref_list = []
	for each_ref in open(ref_list_file):
		ref_list.append(each_ref.strip())

	return ref_list

def create_query_and_ref_seq(index_list, genBank_matrix, ref_list_file, tmp_dir):
	ref_list = read_ref_list(ref_list_file)
	os.makedirs(tmp_dir, exist_ok=True)

	write_query_seq = open(join(tmp_dir, 'query_seq.fa'), 'w')
	write_ref_seq = open(join(tmp_dir, 'ref_seq.fa'), 'w')

	with open(genBank_matrix) as file:
		csv_reader = csv.DictReader(file, delimiter='\t')
		for row in csv_reader:
			accession = row["Primary Accession"]
			sequence = row["Sequence"]
			if accession in ref_list:
				write_ref_seq.write(">" + accession + '\n' + sequence + '\n')
			else:
				write_query_seq.write(">" + accession + '\n' + sequence + '\n')
		
		print ("database and query file generated successfully")

	write_query_seq.close()
	write_ref_seq.close()

def check_blast_exists(command):
    try:
        subprocess.run([command, '-version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"{command} is working.")
        return True
    except subprocess.CalledProcessError:
        print(f"{command} is not installed correctly.")
        return False
    except FileNotFoundError:
        print(f"{command} is not found on the system.")
        return False

def run_makeblastdb(ref_file, dbtype, output_dir):
	os.makedirs(output_dir, exist_ok=True)
	base_name = os.path.basename(ref_file)
	command = [
			'makeblastdb',
			'-in', ref_file,
			'-out', output_dir + "/" + base_name,
			'-title', "alignment",
			'-dbtype', 'nucl'
		]
	try:
		subprocess.run(command, check=True)
		print(f"makeblastdb ran successfully on {ref_file}.")
	except subprocess.CalledProcessError as e:
		print(f"Error running makeblastdb: {e}")

def run_blastn(query_file, db_path, output_file):
	command = [
			'blastn', 
			'-query', query_file,
			'-db', db_path,
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
	 		'-out', 'query_tophit.tsv',
			'-outfmt', "6 qacc sacc pident sstrand"
		]
	try:
		subprocess.run(command, check=True)
		print(f"blastn ran successfully. Results saved in {output_file}.")
	except subprocess.CalledProcessError as e:
		print(f"Error running blastn: {e}")

def process(args):
	tmp_dir = args.tmp_dir
	ref_file = args.ref_file
	out_file = args.output_dir

	ref_acc_list = read_ref_list(args.ref_file)
	create_query_and_ref_seq(ref_acc_list, args.genBank_matrix, args.ref_file, tmp_dir)
	check_blast_exists("blastn")
	run_makeblastdb(join(tmp_dir, 'ref_seq.fa'), 'nucl', join(tmp_dir,'db'))
	run_blastn(join(tmp_dir, 'query_seq.fa'), join(tmp_dir, "db", "ref_seq.fa"), out_file)

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the BLAST alignment of query sequences against the given reference')
	parser.add_argument('-g', '--genBank_matrix', help='genBank matrix file', required=True)
	parser.add_argument('-r', '--ref_file', help='Text file containing list of reference sequence accessions', required=True)
	parser.add_argument('-t', '--tmp_dir', help='temp directory to process the data', default="tmp/blast")
	parser.add_argument('-o', '--output_dir', help='output file', default='blast_alignment.tsv')
	args = parser.parse_args()
	process(args)
