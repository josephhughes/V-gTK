import os
import argparse
import read_file
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join

'''
script only handles single master sequence, multiple master seqeunce feature is not implemented
'''
class PadAlignmentSequences:
	def __init__(self, input_dir, reference_aln_table, output_dir, output_file, master_seq_dir):
		self.input_dir = input_dir
		self.reference_aln_table = reference_aln_table
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_seq_dir = master_seq_dir

	def load_ref_aln_table(self):
		ref_cords = {}
		for i in open(self.reference_aln_table):
			acc, aln_start, aln_end = i.strip().split('\t')
			if acc not in ref_cords:
				ref_cords[acc] = [(aln_start, aln_end)]
			else:
				ref_cords[acc].append((aln_start, aln_end))
		return ref_cords 

	def read_master_seq(self):
		for each_fasta in os.listdir(self.master_seq_dir):
			for record in SeqIO.parse(join(self.master_seq_dir, each_fasta), "fasta"):
				return len(record.seq)
	
	def pad_align(self):
		master_seq_len = self.read_master_seq()
		query_acc_lst = []
		ref_coords = self.load_ref_aln_table()
		write_file = open(join(self.output_dir, self.output_file), 'w')
		for each_query_alignment in os.listdir(self.input_dir):
			ref_acc = each_query_alignment
			fasta_seqs_obj = read_file.fasta(join(self.input_dir, each_query_alignment, ref_acc + ".aligned.fasta"))
			for row in fasta_seqs_obj:
				header = row[0]
				seq = row[1]
				start, end = ref_coords[ref_acc][0][0], ref_coords[ref_acc][0][1]
				seq_len = len(seq)
				if header not in query_acc_lst:
					query_acc_lst.append(header)
					if len(ref_coords[ref_acc]) <= 1:
						ref_cord_diff = abs(int(ref_coords[ref_acc][0][0]) - int(ref_coords[ref_acc][0][1]))
						if seq_len != ref_cord_diff:
							start, end = ref_coords[ref_acc][0][0], seq_len
						#11809 11801 71 11880 71 11801
						#if header == "KR337512": print(ref_cord_diff, seq_len, ref_coords[ref_acc][0][0], ref_coords[ref_acc][0][1], start, end) ;l
						prefix_char = "-" * (int(start)-1)
						suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)#"-" * (abs(int(end)-11932))
						write_file.write('>' + header + "\n")#+ "|" + start + "|" + str(end) + "\n")
						write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
					else:
						ref_cord_diff = abs(int(ref_coords[ref_acc][0][0]) - int(ref_coords[ref_acc][-1][1]))
						if seq_len != ref_cord_diff + 1:
							start, end = ref_coords[ref_acc][0][0], seq_len
						prefix_char = "-" * (int(start)-1)
						suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)
						write_file.write('>' + header + "\n")#"|" + start + "|" + str(end) + "\n")
						write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
		write_file.close()

	def process(self):
		os.makedirs(self.output_dir, exist_ok=True)
		self.pad_align()			

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
	parser.add_argument("-i", "--input_dir", default="tmp/Nextalign/query_aln", help="Directory containing Nextalign output alignments for each reference sequences.")
	parser.add_argument("-r", "--reference_aln_table", default="tmp/Tables/reference_features.tsv", help="reference alignment table contains calculated alignment coordinates")
	parser.add_argument("-o", "--tmp_dir", default="tmp/Pad-Alignment", help="Directory to save padded subalignments and merged files.")
	parser.add_argument("-f", "--output_file", default="paded-query-alignment.fa", help="Output file name to store paded alignment file")
	parser.add_argument("-m", "--master_seq_dir", help="Master sequence directory, it can be one or more master sequence directory", default="tmp/Blast/master_seq/")
	args = parser.parse_args()
	pad_aln = PadAlignmentSequences(args.input_dir, args.reference_aln_table, args.tmp_dir, args.output_file, args.master_seq_dir)
	pad_aln.process()	
