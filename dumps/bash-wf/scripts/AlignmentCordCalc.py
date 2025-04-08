from Bio import SeqIO

class AlignmentCordCalculator:
	def __init__(self):
		pass

	def extract_alignment_coordinates(self, aligned_fasta, reference_id):
		alignment_data = {}
		sequences = list(SeqIO.parse(aligned_fasta, "fasta"))
		ref_seq = None
    
		for seq in sequences:
			if seq.id == reference_id:
				ref_seq = str(seq.seq)
				break

		if not ref_seq:
			raise ValueError("Reference sequence not found in the FASTA file!")

		for record in sequences:
			query_seq = str(record.seq)
        
			ref_index = 0
			start = None
			end = None

			for i, (ref_base, query_base) in enumerate(zip(ref_seq, query_seq)):
				if ref_base != "-":
					ref_index += 1
            
				if query_base != "-": 
					if start is None:
						start = ref_index
					end = ref_index 

			alignment_data[record.id] = (start, end)

		return alignment_data

