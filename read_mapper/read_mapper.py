import resource
from multiprocessing import Pool

class ReadMapper:
    def __init__(self, reference_genome, suffix_array):
        """
        Initialize ReadMapper with reference genome and suffix array.

        Inputs:
            - reference_genome: Reference genome sequence.
            - suffix_array: Pre-built suffix array for the reference genome.
        """

        self.reference_genome = reference_genome.upper()
        self.suffix_array = suffix_array

    def find_exact_matches(self, substring):
        """
        Find all exact matches of a substring using binary search on the suffix array.

        Inputs:
            - substring: The query string to find matches for.

        Returns:
            - A list of positions where the substring matches in the reference genome.
        """

        # Find lower and upper bounds of the matches
        left = self.find_lower_bound(substring)
        right = self.find_upper_bound(substring)
        # Return matches if found, otherwise return empty list
        return self.suffix_array[left:right] if left < right else []

    def find_lower_bound(self, substring):
        """
        Perform binary search to find the lower bound of matches for a substring.

        Inputs:
            - substring: The query string.

        Returns:
            - The lower bound index in the suffix array.
        """

        # Set search bounds
        left, right = 0, len(self.suffix_array)
        # Perform binary search
        while left < right:
            mid = (left + right) // 2
            # Extract suffix for comparison
            suffix = self.reference_genome[
                self.suffix_array[mid] : self.suffix_array[mid] + len(substring)
            ]
            # Update search bounds based on comparison
            if suffix < substring:
                left = mid + 1
            else:
                right = mid
        # Return lower bound index
        return left

    def find_upper_bound(self, substring):
        """
        Perform binary search to find the upper bound of matches for a substring.

        Inputs:
            - substring: The query string.

        Returns:
            - The upper bound index in the suffix array.
        """

        # Set search bounds
        left, right = 0, len(self.suffix_array)
        # Perform binary search
        while left < right:
            mid = (left + right) // 2
            # Extract suffix for comparison
            suffix = self.reference_genome[
                self.suffix_array[mid] : self.suffix_array[mid] + len(substring)
            ]
            # Update search bounds based on comparison
            if suffix > substring:
                right = mid
            else:
                left = mid + 1
        # Return upper bound index
        return left

    def divide_read_into_seeds(self, read, num_seeds):
        """
        Divide a read into smaller seed substrings.

        Inputs:
            - read: The read sequence to be divided.
            - num_seeds: Number of seeds to divide the read into.

        Returns:
            - A list of seed substrings.
        """

        read_length = len(read)
        seed_length = read_length // num_seeds
        seeds = []
        # Divide the read into seeds of equal length
        for i in range(num_seeds):
            # Calculate start and end indices for each seed
            start = i * seed_length
            if i == num_seeds - 1:
                end = read_length
            else:
                end = (i + 1) * seed_length
            # Append seed substring to the list
            seeds.append(read[start:end])
        return seeds

    def align_read_seeding(
        self,
        read,
        num_seeds,
        max_missing_seeds,
    ):
        """
        Align a read by dividing it into seeds and aligning the seeds.

        Inputs:
            - read: The read sequence to align.
            - num_seeds: Number of seeds to divide the read into.
            - max_missing_seeds: Maximum number of seeds allowed to be unaligned.

        Returns:
            - The best alignment (start, end, strand) or None if no valid alignment.
        """

        best_alignment = None
        min_distance = float("inf")

        # Compute reverse complement of the read
        reverse_complement_read = str(self.compute_reverse_complement(read)).upper()
        # Try both forward and reverse complement of the read
        for strand, current_read in [("+", read), ("-", reverse_complement_read)]:
            # Divide the read into seeds
            seeds = self.divide_read_into_seeds(current_read, num_seeds)
            # Calculate lengths of each seed
            seed_lengths = [len(seed) for seed in seeds]
            # Calculate cumulative lengths for offset calculation
            cumulative_lengths = [0] + [
                sum(seed_lengths[: i + 1]) for i in range(len(seed_lengths))
            ]

            seed_matches = []
            # Find exact matches for each seed by searching the suffix array
            for seed in seeds:
                matches = self.find_exact_matches(seed)
                seed_matches.append(matches)

            # Align seeds by computing the number of matching seeds for each potential alignment
            alignment_result = self.align_seeds(
                seeds, seed_matches, cumulative_lengths, max_missing_seeds
            )

            # Update best alignment if a valid alignment is found (missing seeds <= max_missing_seeds)
            if alignment_result:
                # Extract alignment information
                alignment_start, alignment_end, distance = alignment_result
                # Update best alignment if distance is smaller
                if distance < min_distance:
                    min_distance = distance
                    best_alignment = (alignment_start, alignment_end, strand)

        # Return best alignment or None
        if best_alignment is not None:
            return best_alignment
        else:
            return None

    def align_seeds(
        self,
        seeds,
        seed_matches,
        cumulative_lengths,
        max_missing_seeds,
    ):
        """
        Align seeds to the reference genome by computing how many seeds match for each potential alignment.

        Inputs:
            - seeds: List of seeds from the read.
            - seed_matches: List of match positions for each seed.
            - cumulative_lengths: Cumulative lengths of seeds.
            - max_missing_seeds: Maximum allowed missing seeds.

        Returns:
            - Best alignment (start, end, distance) or None if no valid alignment.
        """

        # Combine seeds to form complete read
        complete_read = "".join(seeds)

        # Maps alignment start to set of seed indices
        alignment_starts = {}

        for seed_idx, matches in enumerate(seed_matches):
            # Calculate offset for alignment start
            offset = cumulative_lengths[seed_idx]
            # For each match position, add the corresponding seed index to its set
            for match_pos in matches:
                align_start = match_pos - offset
                alignment_starts.setdefault(align_start, set()).add(seed_idx)

        # Maps alignment (start, end) to Hamming distance
        alignment_distances = {}

        # For each potential alignment start
        for align_start, matched_seeds in alignment_starts.items():
            # Check if the number of matched seeds for the current alignment exceeds the threshold
            if len(matched_seeds) >= len(seeds) - max_missing_seeds:
                # Calculate alignment end
                align_end = align_start + cumulative_lengths[-1]
                # Check if alignment is within bounds
                if not (
                    0 <= align_start <= len(self.reference_genome)
                    and align_end <= len(self.reference_genome)
                ):
                    continue
                # Compute Hamming distance
                distance = self.compute_hamming_distance(
                    self.reference_genome[align_start:align_end], complete_read
                )
                # Store alignment distance
                alignment_distances[(align_start, align_end)] = distance

        # If there is at least one valid alignment
        if alignment_distances:
            # Find alignment with minimum distance
            best_alignment = min(alignment_distances, key=alignment_distances.get)
            # Return alignment information
            return (
                best_alignment[0],
                best_alignment[1],
                alignment_distances[best_alignment],
            )
        # Return None if no valid alignment found
        return None

    def compute_hamming_distance(self, ref_seq, read_seq):
        """
        Compute the Hamming distance between two sequences.

        Inputs:
            - ref_seq: The reference sequence.
            - read_seq: The read sequence.

        Returns:
            - The Hamming distance between the two sequences.
        """

        # Check if sequences are of equal length
        if len(ref_seq) != len(read_seq):
            raise ValueError("Sequences must be of equal length.")

        # Calculates number of mismatches between the two sequences
        return sum(base1 != base2 for base1, base2 in zip(ref_seq, read_seq))

    def compute_reverse_complement(self, seq):
        """
        Compute the reverse complement of a DNA sequence.

        Inputs:
            - seq: The DNA sequence.

        Returns:
            - The reverse complement of the sequence.
        """

        # Complement base pairs
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        # Convert sequence to complement
        complement_sequence = "".join([complement[base] for base in seq])
        # Reverse the complement sequence
        return complement_sequence[::-1]

    def get_memory_usage(self):
        """
        Get memory usage of the current process.

        Returns:
            - Memory usage in GB.
        """

        # Get memory usage in bytes and convert to GB
        mem_usage_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        mem_usage_gb = mem_usage_bytes / (1024 * 1024 * 1024)
        return mem_usage_gb

    def map_reads(self, reads, read_qualities, num_seeds=6, max_missing_seeds=5):
        """
        Map a set of reads to the reference genome using multiprocessing
        """
        max_memory = 0

        read_items = list(reads.items())

        with Pool(processes=6) as pool:
            results = pool.map(self.map_single_read_wrapper, read_items)

        alignments = {}
        for read_id, alignment, mem_usage in results:
            alignments[read_id] = alignment
            if mem_usage > max_memory:
                max_memory = mem_usage

        return alignments, max_memory

    def map_single_read_wrapper(self, item):
        return self.map_single_read(item)

    def map_single_read(self, item):
        read_id, read = item
        alignment = self.align_read_seeding(read, num_seeds=6, max_missing_seeds=5)
        mem_usage = self.get_memory_usage()
        return (read_id, alignment, mem_usage)
    
    def compute_cigar_string(self, read, start_pos, is_reverse=False):
        """
        Compute the CIGAR string for the alignment of a read to the reference sequence.

        Parameters:
            read (str): The read sequence.
            start_pos (int): Start position of the alignment on the reference.
            is_reverse (bool): Whether the read is reverse-complemented.

        Returns:
            str: The CIGAR string representing the alignment.
        """
        reference = self.reference_genome

        # Reverse complement the read if necessary
        if is_reverse:
            read = self.reverse_complement(read)

        reference_segment = reference[start_pos:start_pos + len(read)]

        # Validate alignment
        # if not self.is_valid_alignment(reference_segment, read):
        #     return None  # Invalid alignment
        if len(reference_segment) != len(read):
            print(f"Reference segment length mismatch: Ref len={len(reference_segment)}, Read len={len(read)}")
            return None

        cigar = []
        operation = None
        operation_length = 0

        for i in range(len(read)):
            ref_base = reference[start_pos + i]
            read_base = read[i]

            # Determine operation: '=' for match, 'X' for mismatch
            current_operation = "=" if ref_base.upper() == read_base.upper() else "X"
            # print(f"Position {i + start_pos}: Ref={ref_base}, Read={read_base}, Op={current_operation}")

            # Extend or finalize the current operation
            if current_operation == operation:
                operation_length += 1
            else:
                if operation is not None:
                    # print(f"Appending to CIGAR: {operation_length}{operation}")
                    cigar.append(f"{operation_length}{operation}")
                operation = current_operation
                operation_length = 1

        # Finalize the last operation
        if operation is not None:
            cigar.append(f"{operation_length}{operation}")

        return "".join(cigar)

    def is_valid_alignment(self, reference_segment, read_sequence):
        mismatches = sum(1 for ref_base, read_base in zip(reference_segment, read_sequence) if ref_base.upper() != read_base.upper())
        return mismatches / len(read_sequence) < 0.1  # Example: Allow up to 70% mismatches
    
    def reverse_complement(seq):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complement[base] for base in reversed(seq))
    
