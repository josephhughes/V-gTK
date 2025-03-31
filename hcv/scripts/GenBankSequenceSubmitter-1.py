#python scripts/GenBankSequenceSubmitter.py -q test_data/fasta-seq/ -m test_data/metadata.tsv -n test_data/template.sbt
#python scripts/GenBankSequenceSubmitter-1.py -q test_data/example-2/fasta-seq/ -m test_data/example-2/metadata.tsv -n test_data/example-2/template.sbt -gff tmp/Gff/NC_001542.gff3 -db gdb.db

import os
import re
import csv
import uuid
import glob
import shutil
import sqlite3
import requests
import read_file
import subprocess
from Bio import SeqIO
from time import sleep
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
from PadAlignment import PadAlignment
from FeatureCalculator  import FeatureCordCalculator 

class GenBankSequenceSubmitter:
	def __init__(self, sequence_dir, tmp_dir, output_dir, metadata, ncbi_submission_template, gff_file, db, gaps_to_ignore):
		self.sequence_dir = sequence_dir
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.metadata = metadata
		self.ncbi_submission_template = ncbi_submission_template
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"
		self.gff_file = gff_file
		self.db = db
		self.gaps_to_ignore = gaps_to_ignore
		self.feature_cords = FeatureCordCalculator()
		self.analysis_dir = "analysis_dir" #str(uuid.uuid4())

	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def table2asn(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		command = [
			'table2asn',
			'-indir', join(self.tmp_dir, self.output_dir, "tmp"),
			'-t', self.ncbi_template
		]
		try:
			subprocess.run(command, check=True)
			print(f"Table2asn ran successfully on {self.squence_dir}")
		except subprocess.CalledProcessError as e:
			print(f"Error running Table2asn: {e}")

	def run_makeblastdb(self, db_fasta):
		db_file_name = db_fasta
		command = [
			'makeblastdb',
			'-in', db_fasta,
			'-out', join(db_file_name),
			'-title', "alignment",
			'-dbtype', 'nucl'
		]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {db_file_name}")
		except subprocess.CalledProcessError as e:
			print(f"Error running makeblastdb: {e}")

	def blastn(self, query_path, db_path):
		db_file_name = os.path.basename(db_path)
		output_file = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits.tsv")
		command = [
			'blastn',
			'-query', query_path,
			'-db', join(self.tmp_dir, "analysis", self.analysis_dir, "DB", db_file_name),
			'-task', 'blastn',
			'-max_target_seqs', '1',
			'-max_hsps', '1',
			'-out', output_file,
			'-outfmt', "6 qacc sacc pident sstrand"
		]
		try:
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running blastn: {e}")

	def metadata_header(self):
		headers = ["SeqID", "organism", "genotype", "host", "isolate", "isolation_source", "collection-date", "note", "country", "note"]
		#headers = ["FastaID", "sample_id", "isolation_source", "country" ,"place_sampled" ,"host_species" ,"collection_year","notes", "collection_date"]
		return headers

	def load_metadata(self):
		headers = self.metadata_header() 
		df = {}
		with open(self.metadata, mode="r", encoding="utf-8") as file:
			reader = csv.DictReader(file, delimiter="\t")

			missing_headers = [col for col in headers if col not in reader.fieldnames]
			if missing_headers:
				raise ValueError(f"Missing columns in header: {', '.join(missing_headers)}")

			for row in reader:
				seq_id = row["SeqID"]
				if seq_id in df:
					raise ValueError(f"Duplicate SeqID found: {SeqID}")

				df[seq_id] = {key: row[key] for key in headers if key in row}

		return df

	def load_fasta(self):
		fasta_dict = {}

		for filename in os.listdir(self.sequence_dir):
			if filename.endswith(".fasta") or filename.endswith(".fa"):
				filepath = os.path.join(self.sequence_dir, filename)
				with open(filepath, "r") as file:
					for record in SeqIO.parse(file, "fasta"):
						if record.id in fasta_dict:
							raise ValueError(f"Duplicate accession found: {accession}. Ensure unique accessions.")
						fasta_dict[record.id] = str(record.seq)

		return fasta_dict

	def data_prep(self):
		os.makedirs(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp"), exist_ok=True)
		header = self.metadata_header()
		meta_data = self.load_metadata()
		sequence_obj = self.load_fasta()
		
		missing_in_fasta = [seq_id for seq_id in sequence_obj.keys() if seq_id not in meta_data.keys()]
		if len(missing_in_fasta) < 1:
			write_query_seq = open(join(self.tmp_dir, "analysis", self.analysis_dir, "query.fa"), 'w')
			#creating meta-data file
			for each_acc, data in meta_data.items():
				write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", each_acc + ".src"), "w")
				write_file.write("\t".join(header))
				write_file.write("\n")
				write_file.write("\t".join(data.values()))
				write_file.close()

			#creating fasta file
			for acc, seq in sequence_obj.items():
				write_query_seq.write(">" + acc + "\n" + seq + "\n")
				write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", acc + ".fsa"), "w")
				write_file.write(">" + acc)
				write_file.write("\n")
				write_file.write(seq)
				write_file.write("\n")
				write_file.close()
			write_query_seq.close()
		else:
			raise ValueError(f"Following sequences accssions are not found on meta data {missing_in_fasta}")		

	def extract_ref_seq(self):
		db_path = join(self.tmp_dir, "analysis", self.analysis_dir, "DB")
		os.makedirs(db_path, exist_ok=True)
		write_file = open(join(db_path, "db.fa"), 'w')
		conn = sqlite3.connect(self.db)
		cursor = conn.cursor()

		cursor.execute("SELECT primary_accession, accession_type FROM meta_data where accession_type='reference' OR accession_type = 'master'")
		ref_accs = cursor.fetchall()
		ref_accs_list = [(item[0], item[1]) for item in ref_accs]
		for each_acc, accession_type in ref_accs_list:
			cursor.execute("SELECT sequence FROM sequences WHERE header = ?", (each_acc,))
			if each_acc == "NC_001542": accession_type="master"
			result = cursor.fetchone()
			if result:
				sequence = result[0]
				acc_type = "|" + accession_type
				write_file.write(">" + each_acc + acc_type)
				write_file.write("\n")
				write_file.write(sequence)
				write_file.write("\n")
		write_file.close()

	def blast_analysis(self):

		values = {}
		records = {}

		query_path = join(self.tmp_dir, "analysis", self.analysis_dir, "query.fa")
		ref_path = join(self.tmp_dir, "analysis", self.analysis_dir, "ref.fa")
		db_path = join(self.tmp_dir, "analysis", self.analysis_dir, "DB", "db.fa")
		query_tophits = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits.tsv")
		query_tophits_uniq = join(self.tmp_dir, "analysis", self.analysis_dir, "query_tophits_uniq.tsv")
		sorted_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "sorted_fasta")
		merged_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "merged_fasta")
		sorted_all = join(self.tmp_dir, "analysis", self.analysis_dir, "sorted_all")
		grouped_fasta = join(self.tmp_dir, "analysis", self.analysis_dir, "grouped_fasta")
		ref_seq_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "reference_sequences")
		master_seq_dir = join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_sequences')

		self.run_makeblastdb(db_path)
		self.blastn(query_path, db_path)

		with open(query_tophits, newline='') as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				#query, ref, identity, strand,
				col1, col2, col3, col4 = row[0], row[1].split("|")[0], float(row[2]), row[3]
				if col1 in values:
					existing_value = values[col1]
					if col3 > existing_value:
						records[col1] = [col1, col2, col3, col4]
						values[col1] = col3
				else:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

		with open(query_tophits_uniq, 'w', newline='') as file:
			writer = csv.writer(file, delimiter='\t')
			for record in records.values():
				writer.writerow(record)

		seq_dicts = {}
		query_seqs = read_file.fasta(query_path)
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		#Seperate plus and minus strand sequences
		for each_line in open(query_tophits_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if strand == "plus":
				with open(join(sorted_fasta, "plus.fa"), 'a') as file_plus:
					file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
			else:
				with open(join(sorted_fasta, "minus.fa"), 'a') as file_minus:
					file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

		for each_file in os.listdir(sorted_fasta):
			if "minus" in each_file:
				command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
			else:
				command = ["cp", join(sorted_fasta, each_file), join(merged_fasta, each_file)]

			try:
				print(' '.join(command))
				os.system(' '.join(command))
				print(f"seqkit ran successfully for {each_file}")
			except subprocess.CalledProcessError as e:
				print(f"Error running seqkit: {e}")

		file_list = []
		for each_file in os.listdir(merged_fasta):
			prefix = merged_fasta + "/"
			file_list.append(prefix + each_file)
		
		command = ["cat", " ".join(file_list), ">", join(sorted_all, "query_seq.fa")]
		print(' '.join(command))
		try:
			os.system(' '.join(command))
			print(f"concatenation sucessful")
		except subprocess.CalledProcessError as e:
			print(f"Error in concatenation: {e}")

		grouped_dict = {}
		for each_line in open(query_tophits_uniq, 'r'):
			query_acc, ref_acc, identity, strand = each_line.strip().split('\t')
			if ref_acc not in grouped_dict:
				grouped_dict[ref_acc] = [query_acc]
			else:
				grouped_dict[ref_acc].append(query_acc)

		seq_dicts = {}
		master_seq_dict = {}
		query_seqs = read_file.fasta(join(sorted_all, "query_seq.fa"))
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			with open(join(grouped_fasta, each_ref_acc + '.fasta'), 'a') as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + '\n')
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + '\n')
	
		ref_seqs = read_file.fasta(db_path)
		for rows in ref_seqs:
			accs, accs_type = rows[0].split("|")
			seq_dicts[accs] = rows[1].strip()
			if accs_type == "master":
				master_seq_dict[accs] = rows[1].strip()

		for ref_acc_info in grouped_dict.keys():
			if "|" in ref_acc_info:
				split_info = ref_acc_info.split("|")
				ref_accs = split_info[0]
				acc_type = split_info[1]
				seqs = seq_dicts[ref_accs]
			else:
				ref_accs = ref_acc_info
				acc_type = ""
				seqs = seq_dicts[ref_acc_info]

			with open(join(ref_seq_dir, ref_accs + '.fasta'), 'w') as write_file:
				write_file.write(">" + ref_accs + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')

		#writing master sequences
		for master_acc, seqs in master_seq_dict.items():
			with open(join(master_seq_dir, master_acc + '.fasta'), 'w') as write_file:
				write_file.write(">" + master_acc + '\n')
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + '\n')

		#merge reference sequence
		write_file = open(ref_path, 'w')
		for each_ref_fa in os.listdir(ref_seq_dir):
			df = read_file.fasta(join(ref_seq_dir, each_ref_fa))
			for rows in df:
				write_file.write(">" + rows[0].strip() + "\n" + rows[1].strip() + "\n")
		write_file.close()
			
			
	def nextalign(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(query_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]
        
		command_str = " ".join(command)
		print(f"Executing command: {command_str}")
        
		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def nextalign_analysis(self):
		query_dir = join(self.tmp_dir, "analysis", self.analysis_dir,  "grouped_fasta")
		ref_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "reference_sequences") 
		query_aln_output_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "query_aln")
		reference_aln_output_dir = join(self.tmp_dir, "analysis", self.analysis_dir,  "ref_aln")
		master_seq_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "master_sequences")
		ref_path = join(self.tmp_dir, "analysis", self.analysis_dir, "ref.fa")
		
		#query alignment
		for each_query_file in os.listdir(query_dir):
			ref_file = each_query_file
			self.nextalign(
				join(query_dir, each_query_file),
				join(ref_dir, ref_file),
				query_aln_output_dir
			)

		#reference alignment
		for each_master_file in os.listdir(master_seq_dir):
			self.nextalign(
			join(ref_path),
			join(master_seq_dir, each_master_file),
			reference_aln_output_dir
			)

	def read_master_seq(self):
		for each_fasta in os.listdir(self.master_seq_dir):
			for record in SeqIO.parse(join(self.master_seq_dir, each_fasta), "fasta"):
				return len(record.seq)
	
	def load_ref_aln_table(self):
		ref_cords = {}

		db_path = join(self.tmp_dir, "analysis", self.analysis_dir, "DB")
		os.makedirs(db_path, exist_ok=True)
		write_file = open(join(db_path, "db.fa"), 'w')
		conn = sqlite3.connect(self.db)
		cursor = conn.cursor()

		cursor.execute("SELECT * FROM reference_features")
		reference_features = cursor.fetchall()

		for each_item in reference_features:
			acc, aln_start, aln_end, cds_start, cds_end, product = each_item
			if acc not in ref_cords:
				ref_cords[acc] = [(aln_start, aln_end)]
			else:
				ref_cords[acc].append((aln_start, aln_end))
		return ref_cords


	'''
	def pad_alignment(self):
		nextalign_op_path = join(self.tmp_dir, "analysis", self.analysis_dir, "query_aln")
		master_seq_len = 11932
		query_acc_lst = []
		ref_coords = self.load_ref_aln_table()
		write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, "pad_alignment", "paded_query.fa"), 'w')
		for each_query_alignment in os.listdir(nextalign_op_path):
			ref_acc = each_query_alignment
			fasta_seqs_obj = read_file.fasta(join(nextalign_op_path, each_query_alignment, ref_acc + ".aligned.fasta"))
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
						prefix_char = "-" * (int(start)-1)
						suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)
						write_file.write('>' + header + "\n")
						write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
					else:
						ref_cord_diff = abs(int(ref_coords[ref_acc][0][0]) - int(ref_coords[ref_acc][-1][1]))
						if seq_len != ref_cord_diff + 1:
							start, end = ref_coords[ref_acc][0][0], seq_len
						prefix_char = "-" * (int(start)-1)
						suffix_char = "-" * abs(len(prefix_char) + seq_len - master_seq_len)
						write_file.write('>' + header + "\n")
						write_file.write(prefix_char + seq.strip() + suffix_char + "\n")
		write_file.close()
	'''
	def pad_alignment(self):
		reference_alignment_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "ref_aln")
		input_dir = join(self.tmp_dir, "analysis", self.analysis_dir, "query_aln")
		base_dir = self.tmp_dir
		output_dir = self.analysis_dir
		keep_intermediate = False
		print(reference_alignment_dir, input_dir, base_dir, output_dir);l
		processor = PadAlignment(reference_alignment_dir, input_dir, base_dir, output_dir, keep_intermediate)
		processor.process_all_fasta_in_directory()

	def create_gff_dict(self):
		gff_dict = {}
		with open(self.gff_file) as f:
			for each_line in f:
				if not each_line.startswith('#'):
					split_line = each_line.strip().split('\t')
					if len(split_line) >= 5:
						feature = split_line[2]
						start = split_line[3]
						end = split_line[4]
						description = split_line[8]
						match = re.search(r'product=([^;]+)', description)
						product = match.group(1) if match else None

						if feature not in gff_dict:
							gff_dict[feature] = []
							gff_dict[feature].append({'start': start, 'end': end, 'product': product})
						else:
							gff_dict[feature].append({'start': start, 'end': end, 'product':product})
							
		return gff_dict

	def find_cds_for_coordinates(self, gff_dict, query_start, query_end):
		matching_cds = []
		for feature, cds_list in gff_dict.items():
			if feature=='CDS':
				for cds in cds_list:
					cds_start = int(cds['start'])
					cds_end = int(cds['end'])
					if query_end >= cds_start and query_start <= cds_end:
						overlap_start = max(query_start, cds_start)
						overlap_end = min(query_end, cds_end)
						matching_cds.append({
							'start': overlap_start,
							'end': overlap_end,
							'product': cds['product']
					})

		return matching_cds

	def table2asn_coordinates(self, gff_dict, start, end):
		for cds_entry in gff_dict['CDS']:
			cds_start = int(cds_entry['start'])
			cds_end = int(cds_entry['end'])
			product = cds_entry['product']
		
			if cds_start <= start and end <= cds_end:
				if cds_start == start and cds_end == end:
					return [str(start), str(end), product]
				elif cds_start == start:
					return [start, ">" + str(end), product]
				elif cds_end == end:
					return ["<" + str(start), str(end), product]
				else:
					print(start, end, cds_start, cds_end, type(start))
					return ["<" + str(start), ">" + str(end), product]

	def extract_gene_name(self, protein_name):
		words = protein_name.split()
		return words[-1] if words[-1].isalpha() else words[-2]


	def translate_six_frames(self, seq):
		frames = [seq[i:].translate(to_stop=False) for i in range(3)]
		rev_seq = seq.reverse_complement()
		frames += [rev_seq[i:].translate(to_stop=False) for i in range(3)]
		return frames

	def find_orfs(self, protein_seq):
		orfs = []
		seq_str = str(protein_seq)
		start = 0
		while start < len(seq_str):
			start = seq_str.find("M", start)
			if start == -1:
				break
			end = seq_str.find("*", start)
			if end != -1:
				orfs.append(seq_str[start:end+1])
			start += 1
		return orfs

	def select_best_orf(self, orfs):
		return max(orfs, key=len) if orfs else None

	def recalculate_cds_cords(self, dna_sequence: Seq, start: int, end: int):
		best_orf = None
		best_length = 0
		dna_seq = Seq(dna_sequence.replace("-", ""))
		for frame in self.translate_six_frames(Seq(dna_seq)):
			orfs = self.find_orfs(frame)
			longest_orf = self.select_best_orf(orfs)
			if longest_orf and len(longest_orf) > best_length:
				best_orf = longest_orf
				best_length = len(longest_orf)

		if best_orf:
			recalculated_end = start + (best_length * 3) - 1
			return start, recalculated_end
		else:
			return start, end

	def query_feature_table(self):
		gff_dict = self.create_gff_dict()
		sequence_object = read_file.fasta(join(self.tmp_dir, "analysis", self.analysis_dir, "pad_alignment", "paded_query.fa"))

		for row in sequence_object:
			content = self.feature_cords.calculate_alignment_coords(row[0], row[1], self.gaps_to_ignore)
			for each_cords in content['aligned']:
				#write_file = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", row[0] + ".tablel"), 'w')
				write_tbl = open(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp", row[0] + ".tbl"), 'w')
				start, end = each_cords[0], each_cords[1]
				gff_info = self.find_cds_for_coordinates(gff_dict, start, end)
				write_tbl.write(">Feature " + row[0] + "\n")
				for match in gff_info:
					cds_start = match['start']
					cds_end = match['end']
					
					recalc_cords = self.recalculate_cds_cords(row[1][match['start']-1:match['end']], match['start'], match['end'])
					table2asn_cords = self.table2asn_coordinates(gff_dict, recalc_cords[0], recalc_cords[1])
					gene = self.extract_gene_name(table2asn_cords[2])
					data = [row[0], str(start), str(end), str(recalc_cords[0]), str(recalc_cords[1]), match['product']]					
					data_tbl = [str(table2asn_cords[0]), str(table2asn_cords[1])]

					write_tbl.write("\t".join(data_tbl) + "\t"  + "gene" + "\n")
					write_tbl.write("\t" + "\t" + "\t" + "gene" + "\t" + gene + "\n")
					write_tbl.write("\t".join(data_tbl) + "\t"  + "CDS" + "\n")
					write_tbl.write("\t" + "\t" + "\t" + "product" + "\t" + table2asn_cords[2] + "\n")
					write_tbl.write("\t" + "\t" + "\t" + "gene" + "\t" + gene + "\n")
					#write_file.write('\t'.join(data))
					#write_file.write("\n")
		#write_file.close()
		write_tbl.close()
	def table2asn(self, input_dir, ncbi_template):
		command = [
			'table2asn' ,
			'-indir', input_dir,
			'-t', ncbi_template
		]

		command_str = " ".join(command)
		print(f"Executing command: {command_str}")

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{input_dir} completed successfully.")
		else:
			print(f"{input_dir} failed with return code {return_code}")

	def run_table2asn(self):
		self.table2asn(join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp"), self.ncbi_submission_template)

	def move_sqn_to_output_dir(self):
		output_dir = join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "SQN_files")
		tmp_dir = join(self.tmp_dir, "analysis", self.analysis_dir, self.output_dir, "tmp")
		os.makedirs(output_dir, exist_ok=True)
		sqn_files = glob.glob(os.path.join(tmp_dir, '*.sqn'))
		if sqn_files:
			for file_path in sqn_files:
				shutil.copy(file_path, output_dir)
			print(f"GenBank SQN files copied sucessfully to location : {output_dir}")
		else:
				print("No *.SQN file exists")

	def process(self):
		os.makedirs(join(self.tmp_dir, 'analysis'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'merged_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'sorted_all'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'grouped_fasta'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'reference_sequences'), exist_ok=True)
		query_aln_output_dir = join(self.tmp_dir, self.analysis_dir, "query_aln")
		os.makedirs(join(self.tmp_dir, 'analysis', self.analysis_dir, 'master_sequences'), exist_ok=True)
		os.makedirs(join(self.tmp_dir, "analysis", self.analysis_dir, "pad_alignment"), exist_ok=True)
		
		#self.data_prep()
		#self.extract_ref_seq()
		#self.blast_analysis()
		#self.nextalign_analysis()
		self.pad_alignment()
		#self.query_feature_table()
		#self.run_table2asn()
		#self.move_sqn_to_output_dir()
		'''
		1. Read all the sequences from the directory (Done)
		2. Read the meta data (Done)
			1. Check if the sequence header is matching to the meta data provided (Done)
		3. Prepare the meta data file suitable to table2asn (Done but needs minor change in header) 
		4. Add them to tmp dir (Done)
			1. Create column seq_type in gB_matrix (Done)
			2. Seperate sequences from gB_matrix (Done)
			3. Method to extract only reference sequeces from DB using sqlite (Done)
		5. Blast the in-house sequences against the reference sequence (Done)


		6. Nextalign the reference and in-house sequences (Done)
		7. Based on the Blast result fetch only the required reference sequences and master reference (On-going)
		7. Nextalign for reference sequences against master reference (On-going)
		7. Pad alignment
		8. Add insertions
		7. Generate feature coordinates for query sequences(Done)
		8. perform table2asn on the tmp directory (Done)
		'''
	
if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')
	parser.add_argument('-q', '--sequence_dir', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data', default="tmp")
	parser.add_argument('-o', '--output_dir', help='Output directory where processed data and results are stored', default='Table2asn')
	parser.add_argument('-m', '--metadata', help='A tab limited meta data file name', required=True)
	parser.add_argument('-n', '--ncbi_submission_template', help='NCBI submission template file generated from "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/"', required=True)
	parser.add_argument('-gff', '--gff_file', help='Master reference GFF3 file', required=True)
	parser.add_argument('-db', '--vgtk-db', help='VGTK Database', required=True)
	parser.add_argument('-gp', '--gaps_to_ignore', help='Lenght of gaps to ignore', default=30)
	
	args = parser.parse_args()

	processor = GenBankSequenceSubmitter(args.sequence_dir, args.tmp_dir, args.output_dir, args.metadata, args.ncbi_submission_template, args.gff_file, args.vgtk_db, args.gaps_to_ignore)
	processor.process()
	
