import os
import read_file
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser

def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

def mafft(input_seq, output):
		accession = path_to_basename(input_seq)
		command = [
			'mafft',
			f'{input_seq}', '>',
			join(f'{output}', f'{accession}' + '_aln.fa')
		]
        
		command_str = " ".join(command)
		print(f"Executing command: {command_str}")
    
		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

def run_mafft(input_seq_dir, output_dir):
    for each_input_seq in os.listdir(input_seq_dir):
        mafft(join(input_seq_dir, each_input_seq), output_dir)


input_dir = "tmp/MafftSeqs"
output = "tmp/MafftAlignment"
os.makedirs(output, exist_ok=True)
run_mafft(input_dir, output)
