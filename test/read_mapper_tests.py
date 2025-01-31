import unittest

class TestReadMapper(unittest.TestCase):

	def find_exact_matches(self, substring, reference_genome):
		suffix_array = sorted(range(len(reference_genome)), key=lambda x: reference_genome[x:])
		left = self.find_lower_bound(substring, reference_genome)
		right = self.find_upper_bound(substring, reference_genome)
		return suffix_array[left:right] if left < right else []
    

	def find_lower_bound(self, substring, reference_genome):
		suffix_array = sorted(range(len(reference_genome)), key=lambda x: reference_genome[x:])
		left, right = 0, len(suffix_array)
		while left < right:
			mid = (left + right) // 2
			suffix = reference_genome[suffix_array[mid] : suffix_array[mid] + len(substring)]
			if suffix < substring:
				left = mid + 1
			else:
				right = mid
		return left
    

	def find_upper_bound(self, substring, reference_genome):
		suffix_array = sorted(range(len(reference_genome)), key=lambda x: reference_genome[x:])
		left, right = 0, len(suffix_array)
		while left < right:
			mid = (left + right) // 2
			suffix = reference_genome[
				suffix_array[mid] : suffix_array[mid] + len(substring)
			]
			if suffix > substring:
				right = mid
			else:
				left = mid + 1
		return left

        
	def divide_read_into_seeds(self, read, num_seeds):
		read_length = len(read)
		seed_length = read_length // num_seeds
		seeds = []
		for i in range(num_seeds):
			start = i * seed_length
			if i == num_seeds - 1:
				end = read_length	# if read doesn't divide evenly, remaining characters all go in last seed
			else:
				end = (i + 1) * seed_length
			seeds.append(read[start:end])
		return seeds


	def align_read_seeding(self, reference_genome, read, num_seeds, max_missing_seeds):

		if len(read) == 0:
			return None
		
		read_length = len(read)
		ref_length = len(reference_genome)

		if read_length > ref_length:
			return None

		
		best_alignment = None
		min_distance = float("inf")
		suffix_array = sorted(range(len(reference_genome)), key=lambda x: reference_genome[x:])

		reverse_complement_read = str(self.compute_reverse_complement(read)).upper()
		for strand, current_read in [("+", read), ("-", reverse_complement_read)]:
			seeds = self.divide_read_into_seeds(current_read, num_seeds)
			seed_lengths = [len(seed) for seed in seeds]
			cumulative_lengths = [0] + [sum(seed_lengths[: i + 1]) for i in range(len(seed_lengths))]

			seed_matches = []
			for seed in seeds:
				matches = self.find_exact_matches(seed, reference_genome)
				seed_matches.append(matches)

			alignment_result = self.align_seeds(seeds, seed_matches, cumulative_lengths, max_missing_seeds, reference_genome)

			if alignment_result:
				alignment_start, alignment_end, distance = alignment_result
				if distance < min_distance:
					min_distance = distance
					best_alignment = (alignment_start, alignment_end, strand)

		if best_alignment is not None:
			return best_alignment
		else:
			return None
		

	def align_seeds(self, seeds, seed_matches, cumulative_lengths, max_missing_seeds, reference_genome):
		complete_read = "".join(seeds)
		alignment_starts = {}

		for seed_idx, matches in enumerate(seed_matches):
			offset = cumulative_lengths[seed_idx]
			for match_pos in matches:
				align_start = match_pos - offset
				alignment_starts.setdefault(align_start, set()).add(seed_idx)

		alignment_distances = {}
		for align_start, matched_seeds in alignment_starts.items():
			if len(matched_seeds) >= len(seeds) - max_missing_seeds:
				align_end = align_start + cumulative_lengths[-1]
				if not (0 <= align_start <= len(reference_genome) and align_end <= len(reference_genome)):
					continue
				distance = self.compute_hamming_distance(reference_genome[align_start:align_end], complete_read)
				alignment_distances[(align_start, align_end)] = distance

		if alignment_distances:
			best_alignment = min(alignment_distances, key=alignment_distances.get)
			return (best_alignment[0], best_alignment[1], alignment_distances[best_alignment])
		return None
	

	def compute_reverse_complement(self, seq):
		complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
		complement_sequence = "".join([complement[base] for base in seq])
		return complement_sequence[::-1]
	
	# Could consider using quality scores to weight the distance calculation
	def compute_hamming_distance(self, ref_seq, read_seq):
		if len(ref_seq) != len(read_seq):
			raise ValueError("Sequences must be of equal length.")
		return sum(base1 != base2 for base1, base2 in zip(ref_seq, read_seq))


	
	# TEST CASES FOR FIND_EXACT_MATCHES() FUNCTION -- returns indices from largest to smallest

	def test_find_exact_matches_empty(self):
		self.assertEqual(self.find_exact_matches("", ""), [])
	
	def test_find_exact_matches_one_match(self):
		self.assertEqual(self.find_exact_matches("CT", "ACTG"), [1])

	def test_find_exact_matches_many_match(self):
		self.assertEqual(self.find_exact_matches("ACT", "ACTGACTGACT"), [8, 4, 0])

	def test_find_exact_matches_full_match(self):
		self.assertEqual(self.find_exact_matches("ACTGACTG", "ACTGACTG"), [0])

	def test_find_exact_matches_no_match(self):
		self.assertEqual(self.find_exact_matches("G", "AAAAA"), [])

	def test_find_exact_matches_substring_too_long(self):
		self.assertEqual(self.find_exact_matches("AA", "A"), [])
	

	# TEST CASES FOR FIND_LOWER_BOUND() FUNCTION
	
	def test_find_lower_bound_empty(self):
		self.assertEqual(self.find_lower_bound("", "ACTG"), 0)
            
	def test_find_lower_bound_one_match(self):
		self.assertEqual(self.find_lower_bound("CT", "ACTG"), 1)
            
	def test_find_lower_bound_many_match(self):
		self.assertEqual(self.find_lower_bound("CT", "GACTGACT"), 2)
            
	def test_find_lower_bound_no_match(self): # if no match, LB is 0
		self.assertEqual(self.find_lower_bound("A", "GGGGG"), 0)

	def test_find_lower_bound_substring_too_long(self): # if substring longer than genome, LB is length of genome
		self.assertEqual(self.find_lower_bound("AAAAA", "A"), 1)

	
	# TEST CASES FOR FIND_UPPER_BOUND() FUNCTION
	
	def test_find_upper_bound_empty(self):
		self.assertEqual(self.find_upper_bound("", "ACTG"), 4)
            
	def test_find_upper_bound_one_match(self):
		self.assertEqual(self.find_upper_bound("CT", "ACTG"), 2)
            
	def test_find_upper_bound_many_match(self):
		self.assertEqual(self.find_upper_bound("CT", "GACTGACT"), 4)
            
	def test_find_upper_bound_no_match(self): # if no match, UB is 0
		self.assertEqual(self.find_upper_bound("A", "GGGGG"), 0)

	def test_find_upper_bound_substring_too_long(self): # if substring longer than genome, UB is length of genome
		self.assertEqual(self.find_upper_bound("AAAAA", "A"), 1)


	# TEST CASES FOR DIVIDE_READ_INTO_SEEDS() FUNCTION

	def test_divide_read_empty(self):
		self.assertEqual(self.divide_read_into_seeds("", 1), [""])

	def test_divide_read_one_seed(self):
		self.assertEqual(self.divide_read_into_seeds("ACTG", 1), ["ACTG"])

	def test_divide_read_n_seeds(self):
		self.assertEqual(self.divide_read_into_seeds("ACTG", 4), ["A", "C", "T", "G"])

	def test_divide_read_even(self):
		self.assertEqual(self.divide_read_into_seeds("ACGTACGT", 4), ["AC", "GT", "AC", "GT"])

	def test_divide_read_uneven(self):
		self.assertEqual(self.divide_read_into_seeds("ACGTACGT", 3), ['AC', 'GT', 'ACGT'])


	# TEST CASES FOR ALIGN_READ_SEEDING() FUNCTION
	def test_align_seed_reading_empty(self):
		self.assertEqual(self.align_read_seeding("", "", 1, 0), None)

	def test_align_seed_reading_forward_strand(self):
		self.assertEqual(self.align_read_seeding("CCCACTGTTT", "ACTG", 2, 0), (3, 7, "+"))

	def test_align_seed_reading_reverse_strand(self):
		self.assertEqual(self.align_read_seeding("CCCACTGTTT", "CAGT", 3, 0), (3, 7, "-"))

	def test_align_seed_reading_allow_missing(self):
		self.assertEqual(self.align_read_seeding("GGGCATGGG", "CGT", 2, 1), (3, 6, "+"))

	def test_align_seed_reading_all_missing(self):
		self.assertEqual(self.align_read_seeding("AAAAA", "GGGGG", 1, 5), None)	# no valid alignment

	def test_align_seed_reading_read_too_long(self):
		self.assertEqual(self.align_read_seeding("ACTG", "ACTGACTG", 2, 0), None) # read longer than ref -> no alignment

	# Test case for empty read (should return None)
	def test_align_read_seeding_empty_read(self):
		self.assertEqual(self.align_read_seeding("ACTG", "", 1, 0), None)

	# Test case where read is longer than genome (should return None)
	def test_align_read_seeding_read_longer_than_genome(self):
		self.assertEqual(self.align_read_seeding("ACT", "ACTGACTG", 1, 0), None)

	def test_align_read_seeding_full_genome_match(self):
		self.assertEqual(self.align_read_seeding("ACTG", "ACTG", 1, 0), (0, 4, "+"))
	
	# TEST CASES FOR ALIGN_SEEDS() FUNCTION

	def test_align_seeds_empty(self):
		self.assertEqual(self.align_seeds([], [], [0], 0, "ACTG"), None)

	def test_align_seeds_one_match(self):
		self.assertEqual(self.align_seeds(["ACT"], [[0]], [0, 3], 0, "ACTG"), (0, 3, 0))

	def test_align_seeds_many_match(self): # multiple matches, choose the one with smaller hamming distance
		self.assertEqual(self.align_seeds(["ACT", "GAC"], [[0, 4], [3, 7]], [0, 3, 6], 0, "ACTGACGACTGA"), (0, 6, 0))

	def test_align_seeds_all_match(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC", "TGA"], [[0], [3], [6]], [0, 3, 6, 9], 0, "ACTGACTGA"), (0, 9, 0))

	def test_align_seeds_one_missing(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC", "TGA"], [[0], [], [6]], [0, 3, 6, 9], 1, "ACTGACTGACTG"), (0, 9, 0))

	def test_align_seeds_many_missing(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC", "TGA"], [[0], [], []], [0, 3, 6, 9], 2, "ACTGACTGACTG"), (0, 9, 0))

	def test_align_seeds_too_many_missing(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC", "TGA"], [[0], [], []], [0, 3, 6, 9], 1, "ACTGACTGACTG"), None)

	def test_align_seeds_out_of_bounds(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC"], [[10], [13]], [0, 3, 6], 1, "ACTGACTG"), None)

	def test_align_seeds_out_of_bounds(self):
		self.assertEqual(self.align_seeds(["ACT", "GAC"], [[10], [13]], [0, 3, 6], 1, "ACTGACTG"), None)


	# TEST CASES FOR COMPUTE_REVERSE_COMPLEMENT() FUNCTION

	def test_reverse_complement_empty_string(self):
		self.assertEqual(self.compute_reverse_complement(""), "")

	def test_reverse_complement_single_char(self):
		self.assertEqual(self.compute_reverse_complement("A"), "T")

	def test_reverse_complement_simple_genome(self):
		self.assertEqual(self.compute_reverse_complement("ATGC"), "GCAT")

	def test_reverse_complement_complex_genome(self):
		self.assertEqual(self.compute_reverse_complement("AAAACCCGGT"), "ACCGGGTTTT")

	def test_reverse_complement_very_complex_genome(self):
		self.assertEqual(self.compute_reverse_complement("CTACAGGATCGGTGAGTTCGATCCCGTCTAGACAACGATCAAAGCTGTGATGGGAATCTGTCCTCCCTTCCACGTACTTAGCTAGTTAAAAGATAATTGCCCACTTGATCGGTCGCCACTGGTGCGCGTGGTTCATGGACTGAAGATACCGATACCGAGGCGGGAACTTGTAGAGCCCCTTATCGGCCATGTGTGGATTAGGAAACGACCTTGGCCCGGCTCATCGGCTGGCAGGAGTGAGAATCACGGGCGCTTACCTTACTTCAATAAGCCGTTAAAACTATATTAGTCGGCTGGATATAGTCACTGCGAACATGATCCAGCTGGAAGTGATGTATGTTATGCTGGCTAAGCCCGGGGCAGAGCGATTCATAGCTTGACCATAAGTAGCTGCGCTATACCCCATACGTGTACGTACAGAACACCTTATACACCACGGAGAGTCCTCGGATATTAACATGGCGTTGCAACACGAGACGGATTAATGGAAGCTTACCTTTCCTCCTTCAATAGCTGGGTAAATCACGTGCTATCCAAGGTGATATTTGCAATTAAACCCTTAACTATAATCGTTAAACATTGGAGTCTCTGTCCTTCAACGGTCCAACCCTTCTTTGCCGGCGCGGCGCTGAGACCCTGATCGACTCGGGCATTACTTAAATAGATGAATCAAGATCGTCGTCCATAGCCCAGGCCCGTGTCCGCCTCCGTAGTGCTGTTCAGGCGGTACCAAAGCTCACGCGCTATCGGGGCACGAACTTTACCAGTCGTTGAATCCCCATTTCGTGCACTATTTGTCTTTACACCTCTCCCCTCATGAGGGCAACATTATGAGCTTGGAGATGCGGGGTCGCGTAAACCCAAGGAACCGAATTCTTTCACTGTAATTGACAGGCTGTACGGCCTGGAAAGACTAAATATTCGAGTAGTAGCATTAATTGAACAAATATTGAATCCGCTCCAGATAATTTCTGGCTGGCATCAAGTGACCTGTATGTCGTACGCGGAAACTTGGATGTCATCTCCAATCAAGGGCCAGACTGGCAGTTCGTCGCGACAAGGACCGCTTCCCGTAGGTAAGTGTCGTCCCCAGCAGTTCGACACCCTTTTGGGCAGCACCGTGCGGGTGAGACTCTGACTCGGAGCGCCCGGAGCGTACGGAGACACGTGACACTTTGGGCACCACCGGGGCGCGTAATGTAAGGACCCCGGGCTTTTACGCTATGTTTGCGGCCGAGGCACTGACGGTAGCTAAGAGTGGAAGTTTGACCTAGCTATGTTTTTATTCGGATAACAATGTAATGACCCTGACGCGTAGTCATCGACTCGTTTCTGTAACGCAGGAACGTACCTCGGTGCTGCCCGGTGACTTCCGACGAGACGTATCGGCATGATCCCGCAATAAGTCCGAGGTTACGGTCGGATAGCACTCGGAGCTATAGCCCACAGATTCGGAACGCGTATAATCTCGGTCTGCGCTGGTAATCTTCTCACAATGATTTCCTCTGGATCGATGTCGTAACCCTTTGCCGGACCTTCCCCCACCCGTGCAGCCTAAGATACGCACTTCCAGCGGCGCTCCTCAATGAATCAGTCTCACTATATCACAGATTATGAACTGCCAAACCTGGCGAGTTATACTCTTAAGTAGATTACTATCAGCGACGGTCGCACCTTGGAACAAGATAGATAGTTAGATGATTTTCCTCATCGTAAGTAAGAATGCGTTGCCTCTTGCCAGACTGGGTTTGCACTGAATCAGCGGCACGGCTGTGACTGACTCCGCATCATTGATCCCGCGTGTGATTCAGCAAAGCATTTTGGCATCCACAACGGCTTCCAAAATCTGGGCAGAATATACCGCCCGCTGCTAAGTGCCGACGTCACGCAAACCTTTTATCTGATCACAGCACGGCGTTGTGCCAACACCCAACCCTAAGATTCTATGCTGCCTGAAGGTGTGACTTCCATCACCTTTAGGAGCGCTAGGAATGATTGGAGGGATTGTGTGATCATCATTTGTGTTTCTTCATTGTCATAGTACTACGACTCCCATACGGCCCCTGAGCATCCTATATTGTGCATATCGGCTGAGGTCTAAGCTGATACCATCATCGTAACCAGACCAGTAAAGTAGTGATAAACAGCAGCGTCAATGTCGAGCCCAAATGCTTAACCCGACCAATAAATGTAGCCCCCAATAGCATAGTCCGTTTCTATAACTTACCACTATCGTCGGCGTGGATCGAACGCCTGCCTAATCGTGTAGCCGGCCCAGGTATGACCCCGTAATCGATCTGACCAAAACGGGTTCATATAGAATAGGGCTGCTTAGGGTTTGTTCGCGGACTGCTCCTAACACCCGTGGTGATCCTCGTGCGGACCTATCCTCCCGATGCATTATTGGTCTCTAAGAGTCCCGACGGTTGCCGATTCGTACGCTTGCACTACTCTCTGGTTTATATAGCTTACTATCTCTCATTGGCCAGCCATAAGCCGGTTAATCAAGGAACCCACCAACATTACCCAAGTTTACCCGTGGACGGCTCGCTCGGCAAGCTACGCATTATTGATGTCAGATCAGGTTACCACTTCTACACAATTCGCGCGTCTATTTAAATGTTTCCAGCCCCCACATCCATACGCCGTTCTTGCGTTACTCTCAAATTATCATGTGTGGGTTCGAAGCTGGTTCACCACAGGCAAGTTTTGGATTAAGCGCTGGTCCTTTCGGGAAGGACTCAGCGACACTCAGAGTAAGGGCCATGACGCCCGTTCGAGGCGATTCGGTCTCTGCCGGTCTATTTGTGCGATTACCCCTAGTCGTCTCAAATCCTGCCTGCGGTTGTACCCACAACTTTGAACAGCCCACGATTAAATCTACAGAGGCTCTTCAGCCAACCCCACATGCTCCGTGGACCCGCAAGCGCGTCTCTCCTTGCTCCGGCGCTAGCGCTATGTCGTTCGTTCGTCGTTTGACATCTTGATCCATCAGCTCCACGGACGGCTTAATATGTGAAGCGCAGGTATAAACTTAAGCCCAGCCATGAAATCAGGGTGGACGGCTTTGAGAGTCTTTGAGTACGCTGGCATCCGATTGGCTGCAGTCTGCGAACCGCAGCCACCTGACTCGAGACAGTGAAGCTGGACTTTCTCGGGAATGGGCGCGCAAAAACCATGGTTGGCATGACCTAAATTGCTGTTACGTACCGAGTCGCGAGGGATAGACCCTTTAGGAAAGTGAGGAGATTAACAGTGCAGGGGCACATAAGACCGGCCGTAGTCCCGGTTGTAAGTGGATGGTCTCGACCATGCTCGTCGATCACTGCTTAGGTCCTGATGAAATTGCCTATTAAACTACTCCGACCGAATTGGTGCACCTGTATCGTCTTGTTTGGTACTGAGTTCACTAACTCATTAGACTCTCGTACGGCCCGTCGGCTGTCTTGGCTGGCTAGACTCTACTCCCGTAGAGAATAACCGAAGACCGGGGGACCATGCCAGGCAGGCCACAAATCAAGACTAAGTCGCGCGCCTAATCCATGTTCCGTGGGCGTTCCTGGCCGGCATACCCACCTATTCAAAGATCGATGAATCTCCGCTTGGGGGCTGATCGTGGAACAAAGAGAGTAGCGCGCCATACACGGTAAGCCTAACTCGTGCGCAGGTATCTCAAACTTCAATATTAGTGTACAATGTGGTGACGGCGCGACCGTGGAACTTGACTCCAGGACACTCGCGCACCCCTCGCAAATCGGGCGAAGCTAAATCGGCGGAATCTAAATGCGTACACCGTATTCACACCGAGCTTAGATTCCTTCCAACTACTACTCCGTTACCGATTGCATATTACATTCTGGCCCACTTACAGTAAATTTAAATAGGGGTACGACGCAGGCAGTACCGGCGTACGTACCTGTCCGTACAGGCTGGACGGACATAGTTAATGAGAGGTGCCCTAGCGGGGTTTACACTGGAACAGTCCTTAATAAACCATCAAATCCTTGAGGTTCGGGTCTTATTTAACCTATAGGACCTTCCAAGATCGAAAGCGGTGGAAGCAAAACTCTAAGACGAAAAGTTCTTAGAAAAGGGCGTGTTTGGAATCGCCGTGCGAGGTGTCGGATCCACGCGTGGCGCTAATTGACGAAATACCGTAGGCATGATACCCGTCCTAATAAGAAGGCTAGGAGCAGTTTCCAGTGGCTCGTGGATGGCGAATTATGAAGCAACGTCAAAGGGGCAAAGTCGGAGACGGAGAACACTGCCCAATGATATCACCACACGATTGTTATCGTTTGTCTCACACCCTAGGGATCATATCGGGTTGGCGTCCAGCACCCTGCGGCTGCGACTCAAGGCCGTGAAGGAGACTTGTATCTGGAACTCGTCGAATGTAAGCGTCTTCTGGTAAGCGGGGGATTCCTCAACCCAAAGAAGGGTCTCAAGACCATGGTGGATCAGTACCGCGTACGGAACTAAGAAGATGGAGAGCCTACTCTAAAATTGTTTGATCTGGTCTCCTAGCATGTTGACCCACGCACCGTCCTTAGTTGCCGCAGGCGAGAGAAAATGGTAGGCATTGGTAACCGCAAGCAGTTAGCTGCCTCATCGTCCGAATCACCGTTTGGCCACTTCCAGTCCCTAGGCTATCCCCTATTTTGGTAAATTCTAGGCGGCCAAGTAATTTGTCTTTAGTCTTGACACCGTTCATCAGTCTTGTTAAAATCAGCGAGCTCAACCCGAAAACACGATAATATCAGGACACAATTTGAGTGCGATATTTCCCTTTGGAAAGGAGAACTGAACCGGTGCTGAACTACATAGAGGGGTGTATGCGTGCATAAGGTAAAATGTGAAAGTCATCTATGAGACCGCTATGAGTTGGGCAGTCCACAAAAGGATTTGTGATGAGGGCTCTACCTCGTTCTATCCCCCCGTTTGTCTGGAAAATTGCGGAGTACTATCTGTGAAATGCTGGTCTCACAGAGGTAAGTACGCTAACATCGATACCGGTCGAATAAGACGATCCAAAAAGGTAAATCTAGACCAGGAATAATGGTAAGGGAGACTTCGATCGCCCCAAACTCAAACCGGTACTCGGCTCCCTTGAGGTCGTGCGCGTATACGTATTAAGAGAGGCCTCGATGGGGTCGATCAGGCCGGGCCACCGCGACGTAATGGTGTTGTCTAGTCAGGTCGCCGTCATCGCAATATGCATCTATGTCACTCAGTCGCCTTCTAGCTGGCGACACAAATATCTAAGTGCGGAGGCGGCGAGCGCAGACCCTCCTAAAGGTTCACGAAATACCATGAGTATTACCATTCTACTACGGGACTCAAGCTTAGGCCTCTAACTGGTGGAAGCGAGTTAGGTCGTTCATCATTCCGAGGAAATGGATCAGTTTGGGAACCCATTTATCAGCCTCCAATCGGTAAGTTTGGCAATTTACAAATAGATCACCACTGCGTACTGCGGGCGGCTTAGGAGCTTCCCTCTGGTTTTAGAGATAGTTCTTTTACCTTACGTAATCGTAGAACCCCATCAACTCGACGTTCGCGGGTGCCAGTACAAAAGTGCACTCCATGGTTGAGTCTACGCACCTTGTATTGGCTGACACTGGGTTGATCATAATCCCTTGCTCATAATAGGTCAGCTCGTGGGCCCCAGACCCATGGTGAGAAGTCAACGCCGAATCCTACGATCTGCTTGTCCAAGGTAGATGCGCGAAGCCAATAGATAAACAGCATTACGCTATAGCACTAATGTTAGCCCTTCTCCTCGGACCGACGCATTGGACCGGAATCGTTTCGGATCCCGCCACTCGCCTGAGAAAGAGTGAGACCCGGCTACGTTAGTGTAGGGCTTTTGCGGCCCGTCGGTTAACATCATTTTTTAGACCAGGCATAGAGAAAAGACCAGTGAACTAACTGTGTGCTAGACACGCAGAACGCGATACGTTGCGCACTCCCGTCAGTGAGTGCCATGTCCTTATATGTAGGTGAGTAAACAGGGGCGTCTCTCGTAGAGCTTTGCATGGACTCAACCGTTCTAAACTAGCGAGCGGATCTATTGTTGTGTGCTAAGCTGCAGTTCTCTGGCAACCCTGAGCTTCGGTAAGCACCATAGACCTCCCGGCCGGTACCCTAACTAGGGTAACACGCATGGCTGGCACGAGACAGGTTCTTTCTAGCTCCCGTTGATGACAATATGGTTGCGCTAGGCCGGGATAATCCGCAGCGTCATGACCGCCGAGTCATGACATGATGACTGTGAATGTGCTAACCCGGCCCGTACCCAGAACTGCCTGGAGGACTGAAACCGTATCAGTGAATGAGCCGGCTCAGTCGACAGATTACAAGAGATGGATCGTAGATAGTGTAGGCTCCAAAAGCAGGGGCGTAACAGCCTTACGTACGCCATAGGGCAAAGCACCGCGGACTCTGAACGTTGAGCATCGCTGACGCACGCGCCTTGTACGCGGAAAAGTCGGGATACAACCAGCACACTGTGAGTGGTTGTTGTGTATAGCTTCTGCCTGCTGCCAACCTTAACCCTGCGGCTGTATTCCCTTCCGATTATACCAGCTTGCGCTGTACGTGACCAAGTGCAGTTTTGGGCGTAAATTCAACGTAAGAACTGAGAGATGTGAGCTATTCAGTCGGTGTCGTAAATTGAGTTTGGGAAGGTCCCGATGTATTACTCCGGTGAGGTGGCATGATTAGTAGCGTCTAGCACTAACTGTACAACGTTGAACTGTGTAGATGTATCCTCACCAGTGATCAAAGACGGTTGGGGGGTGAGAAGCAGATCTACACTGACCTGTTATTTCTTACTTGAGTATTGTCCTGTACACTCCCTCGGCTTAGACTGCCGGGCATCGACCTGCCATAGAGCGTCATCTTAAGGCGACGTTTTGCTAGGATGTATACCGGCATCGTCTAGGTCTTGCGCTATTAACCTTAAAGGCGGATTTCTAAATAAGACTTTTCTCTCAACGTATGGCGGCATCCGAGGACATAGACATTTTTTGCTTCCGACGGCTTGAAAGTGTTGTGGAAATGCCCGCCACTCGACCAGCTAAGCTGCCGACTTATTCAGGTTCGCGTATATTGAAGTATTATGCATTGAAGTTATCTTAACCAAGGCTGATGGCCAGCCTCTCTACAGAAAGACCAAGCAATTAAAATCCCAGCTAGATGCTCTCACCCAAATATGCTTCAACATTCAGCTTGGGCTCCGGATAATTCGGATTAATGCCAGTTTCCGACGATTTGCACTATGCTTATGCTGCTGAAGTTCCTAATCGGGAAGAATTGAAAAGTCTGTAATTGCTTTCCTAGTATGCGTTAGATGTGCGGCACACACCAGGCACTGCCTTGCTCCGGACTTAAGACCGAACCGTGCCCCTTGTATACGAATTATTGTACGGCTGACCCCATTTATCAACAACGTGTAACCAATAGATACCGCTTTTGGTTTGAACAATAGGAGGTTCGTACGACCAACAGTCGAACAACGTAGAGCCTCCGCCTATCGGTCGGCGACTCTGTCGGCTCCTAGGAGAAGGAGTCCCCAAACGCCCACGCAAGCTCTTCGTGGGGTGTCTCAGACCGGGTAGCCGACTGCACCTGCATGGGGTGACTAACAGTGTGGGATTCCACAGTCTCACCGCTAGTATATCGCTGCTGCGGACGGGAATTGCGAAGGCTCAGGCCGACAACTTCACCCCCCCGGGTTCCTGGGCATGGATTGCTGGCCAAACTAACACGCATCTGCATGTACAGACCTTTCTTTACTTATGTGATCACGTGTTCCTAGGCATAAACGTTAGCCCACGTTACTAACTGAGGTGCACCGATCTCTCAGTTTCCCCTTCATGCCACAACTTAAAGTATGTCCACGTGAAAAGCTCTTCTTGGTGCTTTGGGCGACACACCTGACGGCTTTAAGTAGTAACAGTTCTAGTCCGCTAGCAGACCCAAATAAAGGAGGCGCGGGAACTAGCTCTTGGGCTGATTCTCTTGTAAGAACTGCAGGTAGAACAAAACGAATGCTGTGCGCTAGCTGAGAAGAATGGCCGACTTTGATTTGACAGTAGGCCGAGGTATCTTGGCGTACTATAAAACAGCGATATACAGTACCCTGGAAAACCGTGGACACGTGTTCTCCGACAACTCTTCAGAGACTGACGTCCGTAGAGCCAATGGAAACAGTGACGAGTGCCTGTTAGTAATGTGGTTGCGCGACTTGGGTTATATTGCTATCCTTACTTTCTGTATAGCATAGAGAGAAGGATGGTATCCCGTAATCCCCTCCGATTCCAAGGCTGAGAGAATACCCGAGACCACACGGGCTCGCGTCCTTGGAATATTACGTACATCATAACGCACTATTGTTCAGCTTCGCCGGAGATGGCTTGGCCCCCTGTAATACAAGACGAGTTCAAAACCCTTTGTTGTGTTCTCCCCAGGGGGTAGAGACTGATTCGTACTGAGAGTCGCCTACTGTCACCAGTTCCGGTTACTGGAACCACCACTGGTGCACCACGATCTCAGTTATCGCCTTACATGGTGGTAGAAAGTCCTCAAACACTCGCAAAAAATGCGGTCGGAGCGAGACACCTCCTACTCGGCGGAAGTCCGAGGGCCGCCGTGCTAACCGGCGTCGGCTCCCCTTTTCCCCTAGGTTATGATAAAGTCGTCCCATGACGGTGGTGCACCCACCTGCTGGCCATGGTGTCAGTGCGTTATTACCTTATGTCGCTATTGAGTGGTGCTCGTGCGTAAGTCTTGGGATAAGGGCACGGCCGGCTGCCCTCGGCGGGTCCGGAATGGTTTAGACTGAGGCAGTTGGCTGTAGTGATTAGCGCATGCTACGTTCTTTGCTGAAGTCTGGTGAACCCATACGAAGATGATGGTGTAAAGATATTTGGCTGTGCGGGGCACTGGAGCCTAGATCATGACGTTGGTGTGACAGACCGTAATTGGGCTCTCATACGATTTGTTTCTCTCGACCCGGAAGCCAGGTTCTACGGGTTGGCAACAGGGCATCAGCGACGTGTCAGAGGTAAGGTCGCTCATTCTAGGTCGTATTGGAGGATAATAAGGCAGGGTGGCCATGCCGCACAGGCATACATACTTGCACGTGGGTTTTAACCGGACGGGTCTACCCTATCTTCGGGTAGCGTCGCTTTAGCGATACAACGTGGTTAAGTTGGCCCACGTAGCGCCAGGAGCCTTCGAGCATCGTCACCGTAGACTAGTTGCAGATCAAGGACTTTTTCCGATTCCACAAGCACCATGAAGTGATATGGAATGAGTCCTATTGCC"), "GGCAATAGGACTCATTCCATATCACTTCATGGTGCTTGTGGAATCGGAAAAAGTCCTTGATCTGCAACTAGTCTACGGTGACGATGCTCGAAGGCTCCTGGCGCTACGTGGGCCAACTTAACCACGTTGTATCGCTAAAGCGACGCTACCCGAAGATAGGGTAGACCCGTCCGGTTAAAACCCACGTGCAAGTATGTATGCCTGTGCGGCATGGCCACCCTGCCTTATTATCCTCCAATACGACCTAGAATGAGCGACCTTACCTCTGACACGTCGCTGATGCCCTGTTGCCAACCCGTAGAACCTGGCTTCCGGGTCGAGAGAAACAAATCGTATGAGAGCCCAATTACGGTCTGTCACACCAACGTCATGATCTAGGCTCCAGTGCCCCGCACAGCCAAATATCTTTACACCATCATCTTCGTATGGGTTCACCAGACTTCAGCAAAGAACGTAGCATGCGCTAATCACTACAGCCAACTGCCTCAGTCTAAACCATTCCGGACCCGCCGAGGGCAGCCGGCCGTGCCCTTATCCCAAGACTTACGCACGAGCACCACTCAATAGCGACATAAGGTAATAACGCACTGACACCATGGCCAGCAGGTGGGTGCACCACCGTCATGGGACGACTTTATCATAACCTAGGGGAAAAGGGGAGCCGACGCCGGTTAGCACGGCGGCCCTCGGACTTCCGCCGAGTAGGAGGTGTCTCGCTCCGACCGCATTTTTTGCGAGTGTTTGAGGACTTTCTACCACCATGTAAGGCGATAACTGAGATCGTGGTGCACCAGTGGTGGTTCCAGTAACCGGAACTGGTGACAGTAGGCGACTCTCAGTACGAATCAGTCTCTACCCCCTGGGGAGAACACAACAAAGGGTTTTGAACTCGTCTTGTATTACAGGGGGCCAAGCCATCTCCGGCGAAGCTGAACAATAGTGCGTTATGATGTACGTAATATTCCAAGGACGCGAGCCCGTGTGGTCTCGGGTATTCTCTCAGCCTTGGAATCGGAGGGGATTACGGGATACCATCCTTCTCTCTATGCTATACAGAAAGTAAGGATAGCAATATAACCCAAGTCGCGCAACCACATTACTAACAGGCACTCGTCACTGTTTCCATTGGCTCTACGGACGTCAGTCTCTGAAGAGTTGTCGGAGAACACGTGTCCACGGTTTTCCAGGGTACTGTATATCGCTGTTTTATAGTACGCCAAGATACCTCGGCCTACTGTCAAATCAAAGTCGGCCATTCTTCTCAGCTAGCGCACAGCATTCGTTTTGTTCTACCTGCAGTTCTTACAAGAGAATCAGCCCAAGAGCTAGTTCCCGCGCCTCCTTTATTTGGGTCTGCTAGCGGACTAGAACTGTTACTACTTAAAGCCGTCAGGTGTGTCGCCCAAAGCACCAAGAAGAGCTTTTCACGTGGACATACTTTAAGTTGTGGCATGAAGGGGAAACTGAGAGATCGGTGCACCTCAGTTAGTAACGTGGGCTAACGTTTATGCCTAGGAACACGTGATCACATAAGTAAAGAAAGGTCTGTACATGCAGATGCGTGTTAGTTTGGCCAGCAATCCATGCCCAGGAACCCGGGGGGGTGAAGTTGTCGGCCTGAGCCTTCGCAATTCCCGTCCGCAGCAGCGATATACTAGCGGTGAGACTGTGGAATCCCACACTGTTAGTCACCCCATGCAGGTGCAGTCGGCTACCCGGTCTGAGACACCCCACGAAGAGCTTGCGTGGGCGTTTGGGGACTCCTTCTCCTAGGAGCCGACAGAGTCGCCGACCGATAGGCGGAGGCTCTACGTTGTTCGACTGTTGGTCGTACGAACCTCCTATTGTTCAAACCAAAAGCGGTATCTATTGGTTACACGTTGTTGATAAATGGGGTCAGCCGTACAATAATTCGTATACAAGGGGCACGGTTCGGTCTTAAGTCCGGAGCAAGGCAGTGCCTGGTGTGTGCCGCACATCTAACGCATACTAGGAAAGCAATTACAGACTTTTCAATTCTTCCCGATTAGGAACTTCAGCAGCATAAGCATAGTGCAAATCGTCGGAAACTGGCATTAATCCGAATTATCCGGAGCCCAAGCTGAATGTTGAAGCATATTTGGGTGAGAGCATCTAGCTGGGATTTTAATTGCTTGGTCTTTCTGTAGAGAGGCTGGCCATCAGCCTTGGTTAAGATAACTTCAATGCATAATACTTCAATATACGCGAACCTGAATAAGTCGGCAGCTTAGCTGGTCGAGTGGCGGGCATTTCCACAACACTTTCAAGCCGTCGGAAGCAAAAAATGTCTATGTCCTCGGATGCCGCCATACGTTGAGAGAAAAGTCTTATTTAGAAATCCGCCTTTAAGGTTAATAGCGCAAGACCTAGACGATGCCGGTATACATCCTAGCAAAACGTCGCCTTAAGATGACGCTCTATGGCAGGTCGATGCCCGGCAGTCTAAGCCGAGGGAGTGTACAGGACAATACTCAAGTAAGAAATAACAGGTCAGTGTAGATCTGCTTCTCACCCCCCAACCGTCTTTGATCACTGGTGAGGATACATCTACACAGTTCAACGTTGTACAGTTAGTGCTAGACGCTACTAATCATGCCACCTCACCGGAGTAATACATCGGGACCTTCCCAAACTCAATTTACGACACCGACTGAATAGCTCACATCTCTCAGTTCTTACGTTGAATTTACGCCCAAAACTGCACTTGGTCACGTACAGCGCAAGCTGGTATAATCGGAAGGGAATACAGCCGCAGGGTTAAGGTTGGCAGCAGGCAGAAGCTATACACAACAACCACTCACAGTGTGCTGGTTGTATCCCGACTTTTCCGCGTACAAGGCGCGTGCGTCAGCGATGCTCAACGTTCAGAGTCCGCGGTGCTTTGCCCTATGGCGTACGTAAGGCTGTTACGCCCCTGCTTTTGGAGCCTACACTATCTACGATCCATCTCTTGTAATCTGTCGACTGAGCCGGCTCATTCACTGATACGGTTTCAGTCCTCCAGGCAGTTCTGGGTACGGGCCGGGTTAGCACATTCACAGTCATCATGTCATGACTCGGCGGTCATGACGCTGCGGATTATCCCGGCCTAGCGCAACCATATTGTCATCAACGGGAGCTAGAAAGAACCTGTCTCGTGCCAGCCATGCGTGTTACCCTAGTTAGGGTACCGGCCGGGAGGTCTATGGTGCTTACCGAAGCTCAGGGTTGCCAGAGAACTGCAGCTTAGCACACAACAATAGATCCGCTCGCTAGTTTAGAACGGTTGAGTCCATGCAAAGCTCTACGAGAGACGCCCCTGTTTACTCACCTACATATAAGGACATGGCACTCACTGACGGGAGTGCGCAACGTATCGCGTTCTGCGTGTCTAGCACACAGTTAGTTCACTGGTCTTTTCTCTATGCCTGGTCTAAAAAATGATGTTAACCGACGGGCCGCAAAAGCCCTACACTAACGTAGCCGGGTCTCACTCTTTCTCAGGCGAGTGGCGGGATCCGAAACGATTCCGGTCCAATGCGTCGGTCCGAGGAGAAGGGCTAACATTAGTGCTATAGCGTAATGCTGTTTATCTATTGGCTTCGCGCATCTACCTTGGACAAGCAGATCGTAGGATTCGGCGTTGACTTCTCACCATGGGTCTGGGGCCCACGAGCTGACCTATTATGAGCAAGGGATTATGATCAACCCAGTGTCAGCCAATACAAGGTGCGTAGACTCAACCATGGAGTGCACTTTTGTACTGGCACCCGCGAACGTCGAGTTGATGGGGTTCTACGATTACGTAAGGTAAAAGAACTATCTCTAAAACCAGAGGGAAGCTCCTAAGCCGCCCGCAGTACGCAGTGGTGATCTATTTGTAAATTGCCAAACTTACCGATTGGAGGCTGATAAATGGGTTCCCAAACTGATCCATTTCCTCGGAATGATGAACGACCTAACTCGCTTCCACCAGTTAGAGGCCTAAGCTTGAGTCCCGTAGTAGAATGGTAATACTCATGGTATTTCGTGAACCTTTAGGAGGGTCTGCGCTCGCCGCCTCCGCACTTAGATATTTGTGTCGCCAGCTAGAAGGCGACTGAGTGACATAGATGCATATTGCGATGACGGCGACCTGACTAGACAACACCATTACGTCGCGGTGGCCCGGCCTGATCGACCCCATCGAGGCCTCTCTTAATACGTATACGCGCACGACCTCAAGGGAGCCGAGTACCGGTTTGAGTTTGGGGCGATCGAAGTCTCCCTTACCATTATTCCTGGTCTAGATTTACCTTTTTGGATCGTCTTATTCGACCGGTATCGATGTTAGCGTACTTACCTCTGTGAGACCAGCATTTCACAGATAGTACTCCGCAATTTTCCAGACAAACGGGGGGATAGAACGAGGTAGAGCCCTCATCACAAATCCTTTTGTGGACTGCCCAACTCATAGCGGTCTCATAGATGACTTTCACATTTTACCTTATGCACGCATACACCCCTCTATGTAGTTCAGCACCGGTTCAGTTCTCCTTTCCAAAGGGAAATATCGCACTCAAATTGTGTCCTGATATTATCGTGTTTTCGGGTTGAGCTCGCTGATTTTAACAAGACTGATGAACGGTGTCAAGACTAAAGACAAATTACTTGGCCGCCTAGAATTTACCAAAATAGGGGATAGCCTAGGGACTGGAAGTGGCCAAACGGTGATTCGGACGATGAGGCAGCTAACTGCTTGCGGTTACCAATGCCTACCATTTTCTCTCGCCTGCGGCAACTAAGGACGGTGCGTGGGTCAACATGCTAGGAGACCAGATCAAACAATTTTAGAGTAGGCTCTCCATCTTCTTAGTTCCGTACGCGGTACTGATCCACCATGGTCTTGAGACCCTTCTTTGGGTTGAGGAATCCCCCGCTTACCAGAAGACGCTTACATTCGACGAGTTCCAGATACAAGTCTCCTTCACGGCCTTGAGTCGCAGCCGCAGGGTGCTGGACGCCAACCCGATATGATCCCTAGGGTGTGAGACAAACGATAACAATCGTGTGGTGATATCATTGGGCAGTGTTCTCCGTCTCCGACTTTGCCCCTTTGACGTTGCTTCATAATTCGCCATCCACGAGCCACTGGAAACTGCTCCTAGCCTTCTTATTAGGACGGGTATCATGCCTACGGTATTTCGTCAATTAGCGCCACGCGTGGATCCGACACCTCGCACGGCGATTCCAAACACGCCCTTTTCTAAGAACTTTTCGTCTTAGAGTTTTGCTTCCACCGCTTTCGATCTTGGAAGGTCCTATAGGTTAAATAAGACCCGAACCTCAAGGATTTGATGGTTTATTAAGGACTGTTCCAGTGTAAACCCCGCTAGGGCACCTCTCATTAACTATGTCCGTCCAGCCTGTACGGACAGGTACGTACGCCGGTACTGCCTGCGTCGTACCCCTATTTAAATTTACTGTAAGTGGGCCAGAATGTAATATGCAATCGGTAACGGAGTAGTAGTTGGAAGGAATCTAAGCTCGGTGTGAATACGGTGTACGCATTTAGATTCCGCCGATTTAGCTTCGCCCGATTTGCGAGGGGTGCGCGAGTGTCCTGGAGTCAAGTTCCACGGTCGCGCCGTCACCACATTGTACACTAATATTGAAGTTTGAGATACCTGCGCACGAGTTAGGCTTACCGTGTATGGCGCGCTACTCTCTTTGTTCCACGATCAGCCCCCAAGCGGAGATTCATCGATCTTTGAATAGGTGGGTATGCCGGCCAGGAACGCCCACGGAACATGGATTAGGCGCGCGACTTAGTCTTGATTTGTGGCCTGCCTGGCATGGTCCCCCGGTCTTCGGTTATTCTCTACGGGAGTAGAGTCTAGCCAGCCAAGACAGCCGACGGGCCGTACGAGAGTCTAATGAGTTAGTGAACTCAGTACCAAACAAGACGATACAGGTGCACCAATTCGGTCGGAGTAGTTTAATAGGCAATTTCATCAGGACCTAAGCAGTGATCGACGAGCATGGTCGAGACCATCCACTTACAACCGGGACTACGGCCGGTCTTATGTGCCCCTGCACTGTTAATCTCCTCACTTTCCTAAAGGGTCTATCCCTCGCGACTCGGTACGTAACAGCAATTTAGGTCATGCCAACCATGGTTTTTGCGCGCCCATTCCCGAGAAAGTCCAGCTTCACTGTCTCGAGTCAGGTGGCTGCGGTTCGCAGACTGCAGCCAATCGGATGCCAGCGTACTCAAAGACTCTCAAAGCCGTCCACCCTGATTTCATGGCTGGGCTTAAGTTTATACCTGCGCTTCACATATTAAGCCGTCCGTGGAGCTGATGGATCAAGATGTCAAACGACGAACGAACGACATAGCGCTAGCGCCGGAGCAAGGAGAGACGCGCTTGCGGGTCCACGGAGCATGTGGGGTTGGCTGAAGAGCCTCTGTAGATTTAATCGTGGGCTGTTCAAAGTTGTGGGTACAACCGCAGGCAGGATTTGAGACGACTAGGGGTAATCGCACAAATAGACCGGCAGAGACCGAATCGCCTCGAACGGGCGTCATGGCCCTTACTCTGAGTGTCGCTGAGTCCTTCCCGAAAGGACCAGCGCTTAATCCAAAACTTGCCTGTGGTGAACCAGCTTCGAACCCACACATGATAATTTGAGAGTAACGCAAGAACGGCGTATGGATGTGGGGGCTGGAAACATTTAAATAGACGCGCGAATTGTGTAGAAGTGGTAACCTGATCTGACATCAATAATGCGTAGCTTGCCGAGCGAGCCGTCCACGGGTAAACTTGGGTAATGTTGGTGGGTTCCTTGATTAACCGGCTTATGGCTGGCCAATGAGAGATAGTAAGCTATATAAACCAGAGAGTAGTGCAAGCGTACGAATCGGCAACCGTCGGGACTCTTAGAGACCAATAATGCATCGGGAGGATAGGTCCGCACGAGGATCACCACGGGTGTTAGGAGCAGTCCGCGAACAAACCCTAAGCAGCCCTATTCTATATGAACCCGTTTTGGTCAGATCGATTACGGGGTCATACCTGGGCCGGCTACACGATTAGGCAGGCGTTCGATCCACGCCGACGATAGTGGTAAGTTATAGAAACGGACTATGCTATTGGGGGCTACATTTATTGGTCGGGTTAAGCATTTGGGCTCGACATTGACGCTGCTGTTTATCACTACTTTACTGGTCTGGTTACGATGATGGTATCAGCTTAGACCTCAGCCGATATGCACAATATAGGATGCTCAGGGGCCGTATGGGAGTCGTAGTACTATGACAATGAAGAAACACAAATGATGATCACACAATCCCTCCAATCATTCCTAGCGCTCCTAAAGGTGATGGAAGTCACACCTTCAGGCAGCATAGAATCTTAGGGTTGGGTGTTGGCACAACGCCGTGCTGTGATCAGATAAAAGGTTTGCGTGACGTCGGCACTTAGCAGCGGGCGGTATATTCTGCCCAGATTTTGGAAGCCGTTGTGGATGCCAAAATGCTTTGCTGAATCACACGCGGGATCAATGATGCGGAGTCAGTCACAGCCGTGCCGCTGATTCAGTGCAAACCCAGTCTGGCAAGAGGCAACGCATTCTTACTTACGATGAGGAAAATCATCTAACTATCTATCTTGTTCCAAGGTGCGACCGTCGCTGATAGTAATCTACTTAAGAGTATAACTCGCCAGGTTTGGCAGTTCATAATCTGTGATATAGTGAGACTGATTCATTGAGGAGCGCCGCTGGAAGTGCGTATCTTAGGCTGCACGGGTGGGGGAAGGTCCGGCAAAGGGTTACGACATCGATCCAGAGGAAATCATTGTGAGAAGATTACCAGCGCAGACCGAGATTATACGCGTTCCGAATCTGTGGGCTATAGCTCCGAGTGCTATCCGACCGTAACCTCGGACTTATTGCGGGATCATGCCGATACGTCTCGTCGGAAGTCACCGGGCAGCACCGAGGTACGTTCCTGCGTTACAGAAACGAGTCGATGACTACGCGTCAGGGTCATTACATTGTTATCCGAATAAAAACATAGCTAGGTCAAACTTCCACTCTTAGCTACCGTCAGTGCCTCGGCCGCAAACATAGCGTAAAAGCCCGGGGTCCTTACATTACGCGCCCCGGTGGTGCCCAAAGTGTCACGTGTCTCCGTACGCTCCGGGCGCTCCGAGTCAGAGTCTCACCCGCACGGTGCTGCCCAAAAGGGTGTCGAACTGCTGGGGACGACACTTACCTACGGGAAGCGGTCCTTGTCGCGACGAACTGCCAGTCTGGCCCTTGATTGGAGATGACATCCAAGTTTCCGCGTACGACATACAGGTCACTTGATGCCAGCCAGAAATTATCTGGAGCGGATTCAATATTTGTTCAATTAATGCTACTACTCGAATATTTAGTCTTTCCAGGCCGTACAGCCTGTCAATTACAGTGAAAGAATTCGGTTCCTTGGGTTTACGCGACCCCGCATCTCCAAGCTCATAATGTTGCCCTCATGAGGGGAGAGGTGTAAAGACAAATAGTGCACGAAATGGGGATTCAACGACTGGTAAAGTTCGTGCCCCGATAGCGCGTGAGCTTTGGTACCGCCTGAACAGCACTACGGAGGCGGACACGGGCCTGGGCTATGGACGACGATCTTGATTCATCTATTTAAGTAATGCCCGAGTCGATCAGGGTCTCAGCGCCGCGCCGGCAAAGAAGGGTTGGACCGTTGAAGGACAGAGACTCCAATGTTTAACGATTATAGTTAAGGGTTTAATTGCAAATATCACCTTGGATAGCACGTGATTTACCCAGCTATTGAAGGAGGAAAGGTAAGCTTCCATTAATCCGTCTCGTGTTGCAACGCCATGTTAATATCCGAGGACTCTCCGTGGTGTATAAGGTGTTCTGTACGTACACGTATGGGGTATAGCGCAGCTACTTATGGTCAAGCTATGAATCGCTCTGCCCCGGGCTTAGCCAGCATAACATACATCACTTCCAGCTGGATCATGTTCGCAGTGACTATATCCAGCCGACTAATATAGTTTTAACGGCTTATTGAAGTAAGGTAAGCGCCCGTGATTCTCACTCCTGCCAGCCGATGAGCCGGGCCAAGGTCGTTTCCTAATCCACACATGGCCGATAAGGGGCTCTACAAGTTCCCGCCTCGGTATCGGTATCTTCAGTCCATGAACCACGCGCACCAGTGGCGACCGATCAAGTGGGCAATTATCTTTTAACTAGCTAAGTACGTGGAAGGGAGGACAGATTCCCATCACAGCTTTGATCGTTGTCTAGACGGGATCGAACTCACCGATCCTGTAG")


	# TEST CASES FOR COMPUTE_HAMMING_DISTANCE() FUNCTION

	def test_reverse_complement_empty_string(self):
		self.assertEqual(self.compute_hamming_distance("", ""), 0)

	def test_reverse_complement_match(self):
		self.assertEqual(self.compute_hamming_distance("ATGC", "ATGC"), 0)

	def test_reverse_complement_mismatch(self):
		self.assertEqual(self.compute_hamming_distance("ATGC", "AGTC"), 2)

	def test_reverse_complement_complex_genome(self):
		self.assertEqual(self.compute_hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC"), 3)

	def test_reverse_complement_very_complex_genome(self):
		self.assertEqual(self.compute_hamming_distance("GCCCTATGCTTGACCTCAAATGAGTCTACCGTCAAAGTTCCATTAGTGACATTACTGGTGACATCGGTTTCGGGCAGCGACAGTTATCTGTAGACACTTTCCGAATATAGGTAACGTGCGGGTGTTGGAATAGCAAGGTAACCAGATTGATGGCTCTCGTGCTGGTGGATAGCAAAGCCTGGGGTCTTTCACAAACTGTTCCATTTCGCTACACGGAGTCCGTGAGATATTTGGGTGCTATTACACGCTACCTATTGGTGAGCTCCTTCGTGGCAGTATGTCCATTCACAATCTGGCAAAAAGAGGCGAGAAATTTGGGTTCTCCAGAAGATGCTCTTCACAGTTGCGGCAAAGGGAGTTGACGTACAGCCGTTGCATTACTGGTGACGACCCGGATGCGCTCTACCTCCCCCTTTCTAATCTGTGTATAAGACTATAGGCAGATCTATCCCGTCTCCTGGCGCTTGGTCAAACAGAAACTTAGAAGACGTTTGGGATAGCGTTCACTTAAACACAGTACACACCGTCTCGCTGTACACAGGTCCTGACCTCTGGACGCAATTGGAGTCGGTAAGTGTGTCGTGTGGCGATACTTAAGTAACTTTATGGCGACCCGAACATGACGAAGCGCCCGTCGCTCTTGTGATTTTACGTTGGAGCCGGCGTGCCGCACACCGGATTTGGGCAAGGAAGATGTGAAGCGTGAATACTATTCGTACGGCGCGTGGTGACACCGTCTTTACTGCCGAGTAATCAGACGTATCTTTTGGGCGCTGTGAAGGAGGGCCGGCTGAGCACGACACGCTTTGCAGTACATTGCCTAGCGCATACTAGTAGAGACGTAGCGACCACCCGATGACGAAAGAGAAAAACTGTAAGGATATTCTTCTCGTTAAGTTCAATACCTGAAAGAACCGCTGAAGCCTTGGTTGTCTCCGGCCCCTTCCGAACGGTAAATTCGGACACCCGTATCTAGACGCATGTCTCCGACGCATCCTGCTTAGTGTCGTACGAGTGATG", "GGTTCCAGTTATAAAAGCATGATTGCTGAAGACAGAAAAAGACCTAAACCACGCATTCGTGGCAGTTGTTGCCAGAATTCGTTGCTCACAGGGCAAAAGAACCAATGTTAATACCCCCTTCATACTTGGCATGAGAAGCACGCTACCTAACAGCTTGGTTAGGACTCACTGAATGGCGTAACATGGAATACATTCAGATAAATACTGACATAACATGCTGATCATTTCCAGCGTACAGTATATCGCATTAACACAAACCTGGCTTACTACAGTTGAGGCAAGGGACATATCATAACTTTGGCTGCTCTTCCGATCGCGGTAGTTGCCTACCTCAGAATGACACATTGGGTGCTAGATCGGTAGTTTACACGGTGAACACTAGAGCGCCCTGGTATACCGTCCCCTTGACGAATCATGACAAGCGACGCGTCGATTAAGTAAGACCTGTCGGACGTCTCACGTGCAAATGAATTATGTTTAGTTACGGGGGTGCGTGGGTCTCCCAGTCAATGAAGTGGACTCCGGCTGGGGAGCCACTTCGCAGACGCTGAGCCCAACTTCGAAAAGTATTTCACGGACAATAGCTTGGGGGGCCATACAATTCCGGGTCTATTTGTCGCAGTTTAGGTATCATAAGTTGTTACGCGACCACTATACTATTATTCGGCGGTTTTATAGAAGGTGGGCCACAAGTGTGCCGGGCCCCTCTCAGTTAAGGGCTTAGAGCATCTTTAGACAGACCCTAAAGCGGAGACGCTTGATAACGACAGCGCCGGACACAATGAGCAAAGGACTCGGCCTATCAGGAGCACGGTCACCGATACTAGTGAAACACGGATTTTGCGACCTGTCAGGCTTAGTAATTGCAGGTCTGTTTGCTCGGCACTTATTCCGCCAGCTTAGGTACTTAAGTGGGACCTTTTATGTCATCTAATGTAGGTAGAGCAATACGAGGTCTCTGAATTCCAGATATTTAGATTGTGGGCAGTGCTGTGCGCTTGGTTTTAACTAACCGACTTT"), 785)
	


if __name__ == '__main__':
    unittest.main()    