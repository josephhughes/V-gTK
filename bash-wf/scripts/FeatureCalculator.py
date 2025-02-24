
class FeatureCordCalculator:
	def __init__(self):
		pass

	def gap_index(self, sequence):
		gap_indices = []
		current_gap = []

		for index, char in enumerate(sequence, start=1):
			if char == '-':
				if not current_gap:
					current_gap.append(index)
			else:
				if current_gap:
					current_gap.append(index - 1)
					gap_indices.append(current_gap)
					current_gap = []

		if current_gap:
			current_gap.append(len(sequence))
			gap_indices.append(current_gap)

		return gap_indices

	def calculate_alignment_coords(self, header, seq, threshold):
		gaps = self.gap_index(seq)
		seq_len = len(seq)
		if len(gaps) == 0:
			return {'aligned': [[1, seq_len]], 'unaligned': []}
	
		overall_start = 1
		overall_end = gaps[-1][1]

		large_gaps = [gap for gap in gaps if (gap[1] - gap[0] + 1) >= threshold]
		aligned_segments = []
		current_start = overall_start

		if len(large_gaps) != 0:
			for gap in large_gaps:
				gap_start, gap_end = gap
				if current_start < gap_start:
					aligned_segments.append([current_start, gap_start - 1])
				current_start = gap_end + 1
		else:
			return {'aligned': [[1, seq_len]], 'unaligned': []}

		return {"aligned": aligned_segments, "unaligned": large_gaps}
