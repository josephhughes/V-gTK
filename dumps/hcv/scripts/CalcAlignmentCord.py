import os
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from GffToDictionary import GffDictionary

class CalculateAlignmentCoordinates:

	def __init__(self, paded_alignment, master_gff, tmp_dir, output_dir, output_file, master_accession):
		self.paded_alignment = paded_alignment
		self.master_gff = master_gff
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_accession = master_accession

	def find_cds_for_coordinates(self, gff_dict, query_start, query_end):
		query_start = int(query_start)
		query_end = int(query_end)

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

	def count_alignment_coordinates(self):
		#fasta_file, master_gff_file, master_id=None, output_file="features.tsv"
		master_gff_dict = GffDictionary(self.master_gff)
		master_gff_dict = master_gff_dict.gff_dict
		
		records = list(SeqIO.parse(self.paded_alignment, "fasta"))

		if self.master_accession:
			master_seq = next((r for r in records if r.id == self.master_accession), None)
			if master_seq is None:
				raise ValueError(f"Master sequence with ID '{self.master_accession}' not found.")
		else:
			master_seq = records[0]
    
		master_alignment = str(master_seq.seq)
    
		master_coords = []
		master_res_count = 0
		for align_pos, base in enumerate(master_alignment, start=1):
			if base != "-":
				master_res_count += 1
				master_coords.append((align_pos, master_res_count))
    
		base_output_dir = join(self.tmp_dir, self.output_dir)
		os.makedirs(base_output_dir, exist_ok=True)
		output_file = join(base_output_dir, self.output_file)
		header = ["accession", "master_ref_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product"]
		with open(output_file, "w") as out_f:

			out_f.write("\t".join(header))
			out_f.write("\n")

			for record in records:
				if record.id == master_seq.id:
					continue 
            
				seq = str(record.seq)
				start_coord = None
				end_coord = None
            
				for (align_pos, master_coord) in master_coords:
					if seq[align_pos - 1] != "-":
						if start_coord is None:
							start_coord = master_coord
						end_coord = master_coord
            
				if start_coord is None or end_coord is None:
					out_f.write(f"{record.id} NA NA\n")
				else:
					cds_cords = self.find_cds_for_coordinates(master_gff_dict, f"{start_coord}", f"{end_coord}")
					for cds in cds_cords:
						data = [record.id, self.master_accession, str(start_coord), str(end_coord), str(cds['start']), str(cds['end']), cds['product']]
						#out_f.write(f"{record.id} {start_coord} {end_coord} {cds['start']} {cds['end']} {cds['product']} \n")
						out_f.write("\t".join(data))
						out_f.write("\n") 
		print(f"Done! Coordinates written to {output_file}")

if __name__ == "__main__":
	parser = ArgumentParser(description='Calculates the genome and cds coordinates for a given sequences')
	parser.add_argument('-i', '--paded_alignment', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	parser.add_argument('-b', '--tmp_dir', help='Base directory', default="tmp")
	parser.add_argument('-d', '--output_dir', help='Output directory where processed data and results are stored', default='Tables')
	parser.add_argument('-o', '--output_file', help='Output file name', default='features.tsv')
	parser.add_argument('-m', '--master_accession', help='Master accession', required=True)
	parser.add_argument('-g', '--master_gff', help='Master GFF3 file', required=True)
	args = parser.parse_args()

	processor = CalculateAlignmentCoordinates(args.paded_alignment, args.master_gff, args.tmp_dir, args.output_dir, args.output_file, args.master_accession)
	processor.count_alignment_coordinates()

'''
# Example usage:
fasta_file = "test_pad/NC_001542.aligned_merged_MSA.fasta"
master_accession = "NC_001542"  # replace with actual master accession or keep None to take the first sequence
count_alignment_coordinates_with_master(fasta_file, master_id=master_accession)
'''
