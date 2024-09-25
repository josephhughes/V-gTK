import os
import csv
from time import sleep
from os.path import join
import xml.etree.ElementTree as ET
from argparse import ArgumentParser

genbank_divisions = ['VRL','PAT','SYN','ENV']

def read_ref_list(ref_list_file):
	ref_list = []
	for each_ref in open(ref_list_file):
		ref_list.append(each_ref.strip())

	return ref_list

def check_gb_division(gb_div_arg_val):
	return_val = []
	if gb_div_arg_val != None:
		for each_val in gb_div_arg_val:
			if each_val in genbank_divisions:
				return_val.append(1)
			else: 
				return_val.append(0)
		return 0 if 0 in return_val else 1
	return 1
	
def filter_columns(genBank_matrix, ref_list_file, tmp_dir, total_length, real_length, prop_ambigious, gb_div, seq_type):
	#add sequence type and prop_ambegious_function values here
	if check_gb_division(gb_div) != 0:
		ref_list = read_ref_list(ref_list_file)
		os.makedirs(tmp_dir, exist_ok=True)

		write_query_seq = open(join(tmp_dir, 'query_seq.fa'), 'w')
		write_ref_seq = open(join(tmp_dir, 'ref_seq.fa'), 'w')

		with open(genBank_matrix) as file:
			csv_reader = csv.DictReader(file, delimiter='\t')
			for row in csv_reader:
				gb_division = row['Division']
				seq_len = row['Length']
				seq_len_without_n = int(row['Length']) - int(row['N'])
				accession = row["Primary Accession"]
				sequence = row["Sequence"]
				if (gb_div is None or gb_division not in gb_div) and int(seq_len) >= total_length and int(seq_len_without_n) >= real_length:
					if accession in ref_list:
						write_ref_seq.write(">" + accession + '\n' + sequence + '\n')
					else:
						write_query_seq.write(">" + accession + '\n' + sequence + '\n')
	
		write_query_seq.close()
		write_ref_seq.close()
	else:
		pass
		# display error message

def process(args):
	gb_division = args.genbank_division
	check_gb_div = check_gb_division(gb_division)
	filter_columns(args.genBank_matrix, args.ref_file, args.tmp_dir,args.total_length, args.real_length, args.prop_ambigious_data, args.genbank_division, args.seq_type)
		
if __name__ == "__main__":
	parser = ArgumentParser(description='Filter sequences and prepare seuences for blast alignment')
	parser.add_argument('-g', '--genBank_matrix', help='genBank matrix file', required=True)
	parser.add_argument('-r', '--ref_file', help='Text file containing list of reference sequence accessions', required=True)
	parser.add_argument('-t', '--tmp_dir', help='temp directory to process the data', default="tmp/blast")
	parser.add_argument('-l', '--total_length', help='Total length of the sequence should be more or equal to the length privided', default=1, type=int)
	parser.add_argument('-n', '--real_length', help='Length of sequences without N', default=1, type=int)

	parser.add_argument('-a', '--prop_ambigious_data', 
		help='Proportion of ambiguous data to be excluded.', 
		nargs='+', 
		default=None, 
		type=str)

	parser.add_argument(
    '-d', '--genbank_division', 
    help=(
        'Genbank division to be excluded. General divisions are:\n'
        '    "VRL" = Viral sequences\n'
        '    "ENV" = Environmental sequences\n'
        '    "PAT" = Patented sequences\n'
        '    "SYN" = Synthetic and chimeric sequences\n'
    ),
    nargs='+',  # Accept multiple values
    default=None,
    type=str
	)
	parser.add_argument('-s', '--seq_type', help='', default=None, type=str)
	args = parser.parse_args()
	process(args)
