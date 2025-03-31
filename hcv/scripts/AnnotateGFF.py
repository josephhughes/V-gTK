import re
import os
from os.path import join

# python annotate_gff.py --reference_list ref_list.txt --tmp_dir /tmp --ref_gff_dir /path/to/gffs --output_dir output --master_gff /path/to/master.gff3 --diff_window 15

'''
NCBI-downloaded GFF files often have different gene/CDS names. This script compares the gene and CDS names to those in the master reference GFF and rewrites them to match the master reference
'''

class AnnotateGFF:
	def __init__(self, reference_list, tmp_dir, ref_gff_dir, output_dir, annotate_gtf, master_gff, diff_window):
		self.reference_list = reference_list
		self.tmp_dir = tmp_dir
		self.ref_gff_dir = ref_gff_dir
		self.output_dir = output_dir
		self.annotate_gtf = annotate_gtf
		self.master_gff = master_gff
		self.diff_window = diff_window
		self.feature_type = "CDS"

	def load_master_gff(self):
		gff_dict = {}
		with open(self.master_gff) as f:
			for each_line in f:
				if not each_line.startswith('#'):
					seqid,source,feature,start,end,score,strand,phase,attributes = each_line.strip().split('\t')
					match = re.search(r'product=([^;]+)', attributes)
					product = match.group(1) if match else None
					if feature not in gff_dict:
						gff_dict[feature] = []
						gff_dict[feature].append({'start': int(start), 'end': int(end), 'product': product})
					else:
						gff_dict[feature].append({'start': int(start), 'end': int(end), 'product':product})
		return gff_dict

	def replace_product_line(self, line, replace_word):
		match = re.search(r'product=([^;]+)', line)
		if match:
			product_desc = match.group(1).strip()
			last_word = product_desc.split()[-1]
			new_line = re.sub(r'product=[^;]+', f'product={replace_word}', line)
			return new_line
		return line

	def annotate_gff(self):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)
		master_gff = self.load_master_gff()
		for each_gff in os.listdir(self.ref_gff_dir):
			print(f"Processing: {each_gff}")
			write_file = open(join(self.tmp_dir, self.output_dir, each_gff), 'w')
			with open(join(self.ref_gff_dir, each_gff)) as f:
				for each_line in f:
					if not each_line.startswith('#'):
						seqid,source,feature,start,end,score,strand,phase,attributes = each_line.strip().split('\t')
						cds_info_partial = [seqid,source,feature,start,end,score,strand,phase]
						if feature == self.feature_type:
							for each_cds in master_gff['CDS']:
								master_cds_start = each_cds['start']
								master_cds_end = each_cds['end']
								if end >= master_cds_start and start <= master_cds_end:
									replace_string = self.replace_product_line(attributes, each_cds['product'])
									write_file.write("\t".join(cds_info_partial) + "\t" + replace_string + "\n")
									break
						else:
							write_file.write(each_line.rstrip() + "\n")
					else:
						write_file.write(each_line.rstrip() + "\n")
			write_file.close()

	def process(self):
		self.annotate_gff()

#1. Load the reference_accession list (Done)
#2. Load the master_gff (Done)
#3. Load each of the reference gff if annotate_gff is disabled (Done)
#4. Load each gff and annotate (Done)
#5. Test it (on-going)
