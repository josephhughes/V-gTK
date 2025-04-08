import os
import csv
import subprocess
from os.path import join
from argparse import ArgumentParser
import read_file

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
		print(f"makeblastdb ran successfully on {ref_file}")
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
		print(f"blastn ran successfully. Results saved in {output_file}")
	except subprocess.CalledProcessError as e:
		print(f"Error running blastn: {e}")

def write_tophits(input_file, query_tophit_uniq, tmp_dir, query_fasta):

	merged_fasta = join(tmp_dir, "merged_fasta")
	os.makedirs(merged_fasta, exist_ok=True)
	output_file = "query_tophit_unique.tsv"

	records = {}
	values = {}
	uniq_values = {}

	with open(input_file, newline='') as file:
		reader = csv.reader(file, delimiter='\t')
		for row in reader:
			col1, col2, col3, col4 = row[0], row[1], float(row[2]), row[3]

			if col1 in values:
				existing_value = values[col1]
				if col3 > existing_value:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

			else:
				records[col1] = [col1, col2, col3, col4]
				values[col1] = col3

	with open(query_tophit_uniq, 'w', newline='') as file:
		writer = csv.writer(file, delimiter='\t')
		for record in records.values():
			writer.writerow(record)
	
	#writing merge fasta files
	write_fasta = open(join(merged_fasta, "seqs.fa"), 'w')
	fasta_seqs = read_file.fasta(query_fasta)
	for each_seq in fasta_seqs:
		header = each_seq[0].strip()
		sequence = each_seq[1].strip()
		if header in records:
			write_fasta.write(">" + header + "\n" + sequence + "\n")
	write_fasta.close()
		
# create a function to process segmented virus
# rename the query_seq.fa to seg_output for next script processing
def process_segment_virus(input_file, uniq_hit_output, segment_file, annotated_output, query_fasta, tmp_dir):	
	uniq_hits = {}
	segments = {}
	segment_assigned = {}

	os.makedirs(join(tmp_dir,"segment_sorted"), exist_ok=True)
	os.makedirs(join(tmp_dir,"segment_merged_fasta"), exist_ok=True)

	segment_sorted = join(tmp_dir, "segment_sorted")
	segment_merged = join(tmp_dir, "segment_merged_fasta")

	with open(input_file, newline='') as file:
		reader = csv.reader(file, delimiter='\t')
		for row in reader:
			col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
			if col1 in uniq_hits:
				if float(col3) > float(uniq_hist[col1][2]):
					uniq_hits[col1] = [col1, col2, col3, col4]
			else:
					uniq_hits[col1] = [col1, col2, col3, col4]

	#writing uniq hits
	write_uniq_hits = open(uniq_hit_output, 'w')
	for k, v in uniq_hits.items():
		write_uniq_hits.write('\t'.join(v) + '\n')
	write_uniq_hits.close()

	for line in open(segment_file):
		accession, segment = line.strip().split('\t')
		segments[accession] = segment

	#writng annotated file
	write_segment = open(annotated_output, 'w')
	with open(uniq_hit_output, newline='') as file:
		reader = csv.reader(file, delimiter='\t')
		for row in reader:
			col1, col2, col3, col4 = row[0], row[1], row[2], row[3]
			if col2 in segments:
				val = [col1, col2, col3, col4, segments[col2]]
				segment_assigned[col1] = val
				write_segment.write('\t'.join(val) + '\n')
			else:
				print("Could not find the segment for :" + col2)
	write_segment.close()

	#writing segment sorted fasta files
	fasta_seqs = read_file.fasta(query_fasta)
	for each_seq in fasta_seqs:
		header = each_seq[0].strip()
		sequence = each_seq[1].strip()

		if header in segment_assigned:
			accession, reference, score, strand, segment = segment_assigned[header]
			if strand == "plus":
				write_file = open(join(segment_sorted, "seg_" + str(segment) + "_" + strand + ".fa"), 'a')
				write_file.write(">" + each_seq[0].strip() + '\n' + sequence + '\n')
			else:
				write_file = open(join(segment_sorted, "seg_" + str(segment) + "_" + strand + ".fa"), 'a')
				write_file.write(">" + each_seq[0].strip() + '\n' + sequence + '\n')			
	
	#combine segments
	for each_seg in os.listdir(segment_sorted):
		name, seg_num, strand = each_seg.split('.')[0].split('_')	
		print ("Processing {name}-{strand}")
		if strand == "plus":
			command = ["seqkit", "grep", "-n", "-f", join(segment_sorted, each_seg), ">>", join(segment_merged, each_seg)]
		else:
			command = ['seqkit', 'seq', '-r', '-p', '-v', '-t', 'dna', join(segment_sorted, each_seg), '>>', join(segment_merged, each_seg)]
			
def execute_seqkit(command):
	try:
		subprocess.run(command, shell=True, check=True)
		print(f"Command executed successfully for seg{i}")
	except subprocess.CalledProcessError as e:
		print(f"Error occurred while executing the command for seg{i}: {e}")
			
def process(args):
	tmp_dir = args.tmp_dir
	db_file = args.db_fasta
	query_file = args.query_fasta
	out_file = args.output_file
	seg_file = args.segment_file

	if args.is_segmented_virus == 'Y' and not args.segment_file:
		parser.error("Argument -f is required when -s is 'Y' for segmented virus.")
		exit(0)
	else:
		if args.is_segmented_virus == 'Y':
			run_makeblastdb(join(tmp_dir, 'ref_seq.fa'), 'nucl', join(tmp_dir,'db'))
			run_blastn(join(tmp_dir, 'query_seq.fa'), join(tmp_dir, "db", "ref_seq.fa"), out_file)
			process_segment_virus(join(args.tmp_dir, "query_tophit.tsv"), join(args.tmp_dir, "query_uniq_tophit.tsv"), seg_file, join(tmp_dir, "query_uniq_tophit_annotated.tsv"), join(tmp_dir, 'query_seq.fa'), tmp_dir)
		else:
			check_blast_exists("blastn")
			run_makeblastdb(join(tmp_dir, 'ref_seq.fa'), 'nucl', join(tmp_dir,'db'))
			run_blastn(join(tmp_dir, 'query_seq.fa'), join(tmp_dir, "db", "ref_seq.fa"), out_file)
			write_tophits(join("tmp/blast", "query_tophit.tsv"), join(tmp_dir, "query_tophits_uniq.txt"), tmp_dir, query_file)

if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the BLAST alignment of query sequences against the given reference')
	parser.add_argument('-q', '--query_fasta', help='query fasta file', required=True, default="tmp/blast/query_seq.fa")
	parser.add_argument('-r', '--db_fasta', help='Blast DB fasta file. (Note: program will consider this file as db file and create the blast index for the given file)', required=True, default="tmp/blast/ref_seq.fa")
	parser.add_argument('-t', '--tmp_dir', help='temp directory to process the data', default="tmp/blast")
	parser.add_argument('-o', '--output_file', help='output file', default='tmp/blast/query_tophits.tsv')
	parser.add_argument('-s', '--is_segmented_virus', help='Type Y for segmented virus else N', default='N')
	parser.add_argument('-f', '--segment_file', help='File containing information about the segments')
	args = parser.parse_args()
	process(args)
