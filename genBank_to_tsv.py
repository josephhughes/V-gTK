import os
from time import sleep
from argparse import ArgumentParser
from os.path import join
import xml.etree.ElementTree as ET
import pandas as pd

# Parse the XML and extract data
def xml_to_tsv(xml_file, output_dir):
	tree = ET.parse(xml_file)
	root = tree.getroot()
	data = []
	for gbseq in root.findall('GBSeq'):
		content = {}
		content['Locus'] = gbseq.find('GBSeq_locus').text
		content['Length'] = gbseq.find('GBSeq_length').text
		content['Strandedness'] = gbseq.find('GBSeq_strandedness').text
		content['Molecule Type'] = gbseq.find('GBSeq_moltype').text
		content['Topology'] = gbseq.find('GBSeq_topology').text
		content['Division'] = gbseq.find('GBSeq_division').text
		content['Update Date'] = gbseq.find('GBSeq_update-date').text
		content['Create Date'] = gbseq.find('GBSeq_create-date').text
		content['Definition'] = gbseq.find('GBSeq_definition').text
		content['Primary Accession'] = gbseq.find('GBSeq_primary-accession').text
		content['Accession Version'] = gbseq.find('GBSeq_accession-version').text
		content['Source'] = gbseq.find('GBSeq_source').text
		content['Organism'] = gbseq.find('GBSeq_organism').text
		content['Taxonomy'] = gbseq.find('GBSeq_taxonomy').text

		mol_type = ''
		isolate = ''
		isolation_source = ''
		db_xref = ''
		country = ''
		host = ''
		collection_date = ''
		segment = ''
		genes = []
		cds = []

		for gb_feature in gbseq.findall('GBSeq_feature-table/GBFeature'):
			if gb_feature.find('GBFeature_key').text == 'source':
				for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
					name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
					value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
					if name == 'mol_type':
						mol_type = value
					elif name == 'isolate':
						isolate = value
					elif name == 'isolation_source':
						isolation_source = value
					elif name == 'db_xref':
						db_xref = value
					elif name == 'country':
						country = value
					elif name == 'host':
						host = value
					elif name == 'collection_date':
						collection_date = value
					elif name == 'segment':
						segment = value
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
				cds_info['cds_location'] = gb_feature.find('GBFeature_location').text if qualifier.find('GBFeature_location') is not None else None
				for qualifier in gb_feature.findall('GBFeature_quals/GBQualifier'):
					name = qualifier.find('GBQualifier_name').text if qualifier.find('GBQualifier_name') is not None else None
					value = qualifier.find('GBQualifier_value').text if qualifier.find('GBQualifier_value') is not None else None
					cds_info[name] = value
				cds.append(cds_info)
						
        
		content['Mol Type'] = mol_type
		content['Isolate'] = isolate
		content['Isolation Source'] = isolation_source
		content['DB Xref'] = db_xref
		content['Country'] = country
		content['Host'] = host
		content['Collection_date'] = collection_date
		content['segment'] = segment
		sequence = gbseq.find('GBSeq_sequence')
		content['Sequence'] = sequence.text if sequence is not None else ''
		content['Genes'] = '; '.join([f"{gene['gene_name']}({gene['gene_location']})" for gene in genes])
		content['CDS Info'] = '; '.join([f"{k}: {v}" for each_cds in cds for k, v in each_cds.items()])

		references = []
		for reference in gbseq.findall('GBSeq_references/GBReference'):
			ref = {}
			content['Reference Number'] = reference.find('GBReference_reference').text
			content['Position'] = reference.find('GBReference_position').text
			authors = [author.text for author in reference.findall('GBReference_authors/GBAuthor')]
			content['Authors'] = ', '.join(authors)
			content['Title'] = reference.find('GBReference_title').text if qualifier.find('GBReference_title') is not None else None
			content['Journal'] = reference.find('GBReference_journal').text
			data.append(content)
	return data

def process(args):
	input_dir = args.input_dir
	output_dir = args.output_dir
	os.makedirs(output_dir, exist_ok=True)
	merged_data = []
	for each_xml in os.listdir(input_dir):
		print("parsing : " + each_xml)
		data = xml_to_tsv(join(input_dir, each_xml), output_dir)
		merged_data.extend(data)

	df = pd.DataFrame(merged_data)
	df = df.drop_duplicates(subset='Locus', keep="last")
	df.to_csv(join(output_dir, "gB_matrix.tsv"), sep="\t", index=False)

if __name__ == "__main__":
	parser = ArgumentParser(description='Extract GenBank XML files to a TSV table')
	parser.add_argument('-d', '--input_dir', help='Input directory', required=True)
	parser.add_argument('-o', '--output_dir', help='Output directory', required=True)
	args = parser.parse_args()
	process(args)

