import os
from time import sleep
from argparse import ArgumentParser
from os.path import join
import xml.etree.ElementTree as ET
import pandas as pd

# add source example NCBI or GISAID, or user define, add temp sequences 
# temp sequences which should available for temp purpose 
class GenBankParser:
	def __init__(self, input_dir, output_dir, ref_list, exclusion_list):
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.ref_list = ref_list
		self.exclusion_list = exclusion_list
		os.makedirs(self.output_dir, exist_ok=True)

	def count_ATGCN(self, sequence):
		nucl_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
		sequence = sequence.upper()
		for each_nucl in sequence:
			if each_nucl in nucl_dict:
				nucl_dict[each_nucl] += 1
		return nucl_dict['A'], nucl_dict['T'], nucl_dict['G'], nucl_dict['C'], nucl_dict['N']

	def load_ref_list(self, acc_list_file):
		if acc_list_file is None:
			return []
		ref_dict = {}
		try:
			with open(acc_list_file) as file:
				for line in file:
					accession, accession_type = line.strip().split("\t")
					if accession not in ref_dict:
						ref_dict[accession] = accession_type#.append(trimmed_line)
		except FileNotFoundError:
			print(f"Warning: File {acc_list_file} not found. Skipping.")
		return ref_dict

	'''
	def load_ref_list(self, acc_list_file):
		ref_list = []
		for each_line in open(acc_list_file):
			trimmed_line = each_line.strip()
			if trimmed_line not in ref_list:
				ref_list.append(trimmed_line)
	
		return ref_list
	'''

	def xml_to_tsv(self, xml_file):
		tree = ET.parse(xml_file)
		root = tree.getroot()
		data = []
		
		ref_seq_dict = self.load_ref_list(self.ref_list)
		exclusion_acc_list = self.load_ref_list(self.exclusion_list)
		
		
		for gbseq in root.findall('GBSeq'):
			content = {}
			content['locus'] = gbseq.find('GBSeq_locus').text
			content['length'] = gbseq.find('GBSeq_length').text
			content['strandedness'] = gbseq.find('GBSeq_strandedness').text if gbseq.find('GBSeq_strandedness') is not None else None
			content['molecule_type'] = gbseq.find('GBSeq_moltype').text
			content['topology'] = gbseq.find('GBSeq_topology').text
			content['division'] = gbseq.find('GBSeq_division').text
			content['update_date'] = gbseq.find('GBSeq_update-date').text
			content['create_date'] = gbseq.find('GBSeq_create-date').text
			content['definition'] = gbseq.find('GBSeq_definition').text
			content['primary_accession'] = gbseq.find('GBSeq_primary-accession').text
			content['accession_version'] = gbseq.find('GBSeq_accession-version').text
			content['gi_number'] = content['primary_accession']
			content['source'] = gbseq.find('GBSeq_source').text
			content['organism'] = gbseq.find('GBSeq_organism').text
			content['taxonomy'] = gbseq.find('GBSeq_taxonomy').text
			
			if content['primary_accession'] in ref_seq_dict:
				content['accession_type'] = ref_seq_dict[content['primary_accession']]
			else:
				content['accession_type'] = 'query'
	
			if content['primary_accession'] in exclusion_acc_list:
				content['accession_type'] = 'excluded'
				content['exclusion_criteria'] = 'excluded by the user'
				content['exclusion_status'] = '1'
			else:
				content['exclusion_criteria'] = ''
				content['exclusion_status'] = '0'

			mol_type = ''
			strain = ''
			isolate = ''
			isolation_source = ''
			db_xref = ''
			country = ''
			host = ''
			collection_date = ''
			segment = ''
			serotype = ''
			genes = []
			cds = []

			for gb_feature in gbseq.findall('GBSeq_feature-table/GBFeature'):
				if gb_feature.find('GBFeature_key').text == 'source':
					for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
						name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
						value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
						if name == 'mol_type':
							mol_type = value
						if name == 'strain':
							strain = value
						elif name == 'isolate':
							isolate = value
						elif name == 'isolation_source':
							isolation_source = value
						elif name == 'db_xref':
							db_xref = value
						elif name == 'country' or name == "geo_loc_name":
							country = value
						elif name == 'host':
							host = value
						elif name == 'collection_date':
							collection_date = value
						elif name == 'segment':
							segment = value
						elif name == 'serotype':
							serotype = value
				elif gb_feature.find('GBFeature_key').text == 'gene':
						gene_info = {}
						gene_info['gene_location'] = gb_feature.find('GBFeature_location').text
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
							value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
							if name == 'gene':
								gene_info['gene_name'] = value
								genes.append(gene_info)
				elif gb_feature.find('GBFeature_key').text == 'CDS':
						cds_info = {}
						cds_info['cds_location'] = gb_feature.find('GBFeature_location').text if gb_feature.find('GBFeature_location') is not None else None
						for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
							name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
							value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
							cds_info[name] = value
						cds.append(cds_info)

			pubmed_ids = []
			for reference in gbseq.findall('GBSeq_references/GBReference'):
				pubmed_tag = reference.find("GBReference_pubmed")
				if pubmed_tag is not None:
					pubmed_id = pubmed_tag.text
					if pubmed_id:
						pubmed_ids.append(pubmed_id)

			content['pubmed_id'] = '; '.join(pubmed_ids)
			content['mol_type'] = mol_type
			content['strain'] = strain
			content['isolate'] = isolate
			content['isolation_source'] = isolation_source
			content['db_xref'] = db_xref
			
			if ":" in country:
				tmp_country = country.split(":")
				content['country'] = tmp_country[0]
				content['geo_loc'] = tmp_country[1] if len(tmp_country) > 0 else ""
			else:
				content['country'] = country
				content['geo_loc'] = ""
		
			content['host'] = host
			content['collection_date'] = collection_date
			content['segment'] = segment
			content['serotype'] = serotype
			sequence = gbseq.find('GBSeq_sequence')
			content['sequence'] = sequence.text if sequence is not None else ''
			content['genes'] = '; '.join([f"{gene['gene_name']}({gene['gene_location']})" for gene in genes])
			content['cds_info'] = '; '.join([f"{k}: {v}" for each_cds in cds for k, v in each_cds.items()])
			nucl_count = self.count_ATGCN(content['sequence'])
			content['a'] = nucl_count[0]
			content['t'] = nucl_count[1]
			content['g'] = nucl_count[2]
			content['c'] = nucl_count[3]
			content['n'] = nucl_count[4]
			atgc_list = [int(content['a']), int(content['t']), int(content['g']), int(content['c'])]
			content['real_length'] = sum(atgc_list)

			references = []
			for reference in gbseq.findall('GBSeq_references/GBReference'):
				ref = {}
				content['reference_number'] = reference.find('GBReference_reference').text
				content['position'] = reference.find('GBReference_position').text
				authors = [author.text for author in reference.findall('GBReference_authors/GBAuthor')]
				content['authors'] = ', '.join(authors)
				content['title'] = reference.find('GBReference_title').text if reference.find('GBReference_title') is not None else None
				content['journal'] = reference.find('GBReference_journal').text
				data.append(content)
	
		return data

	def process(self):
		merged_data = []
		for each_xml in os.listdir(self.input_dir):
			print("Parsing: " + each_xml)
			data = self.xml_to_tsv(join(self.input_dir, each_xml))
			merged_data.extend(data)

		df = pd.DataFrame(merged_data)
		df = df.drop_duplicates(subset='locus', keep="last")
		with open(join(self.output_dir, "sequences.fa"), "w") as fasta_file:
			for _, row in df.iterrows():
				fasta_file.write(f">{row['primary_accession']}\n{row['sequence']}\n")

		df.to_csv(join(self.output_dir, "gB_matrix_raw.tsv"), sep="\t", index=False)
		
if __name__ == "__main__":
	parser = ArgumentParser(description='Extract GenBank XML files to a TSV table')
	parser.add_argument('-d', '--input_dir', help='Input directory', default='tmp/GenBank-XML')
	parser.add_argument('-o', '--output_dir', help='Output directory', default='tmp/GenBank-matrix')
	parser.add_argument('-r', '--ref_list', help='Set of reference accessions', default='generic/ref_list.txt')
	parser.add_argument('-e', '--exclusion_list', help='Set of sequence accssions to be excluded')
	args = parser.parse_args()

	parser = GenBankParser(args.input_dir, args.output_dir, args.ref_list, args.exclusion_list)
	parser.process()

