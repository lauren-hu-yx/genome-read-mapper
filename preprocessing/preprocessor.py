import fastq as fq
from Bio import SeqIO


class Preprocessor:
    def __init__(self, short_reads_ref_file, short_reads_1_file, short_reads_2_file):
        """
        Initialize Preprocessor with reference genome and short reads.

        Inputs:
            - short_reads_ref_file: Path to the reference genome file in FASTA format.
            - short_reads_1_file: Path to the first short reads file in FASTQ format.
            - short_reads_2_file: Path to the second short reads file in FASTQ format.
        """

        # read reference genome from file
        ref_record = next(SeqIO.parse(short_reads_ref_file, "fasta"), None)
        if not ref_record:
            raise Exception("Reference file is empty or improperly formatted.")

        # store reference genome and build suffix array
        self.reference_genome = str(ref_record.seq)
        self.reference_genome_name = ref_record.id
        self.suffix_array = self.build_suffix_array(self.reference_genome)

        # read short reads from files
        self.short_reads_1 = {
            read.getHead()[1:]: (read.getSeq(), read.getQual())
            for read in fq.read(short_reads_1_file)
        }

        self.short_reads_2 = {
            read.getHead()[1:]: (read.getSeq(), read.getQual())
            for read in fq.read(short_reads_2_file)
        }

    def build_suffix_array(self, s):
        """
        Builds a suffix array for the reference genome.

        Inputs:
            - reference_genome: The reference genome sequence as a string.

        Returns:
            - A sorted list of starting indices of all suffixes in the reference genome.
        """
        n = len(s)
        k = 1
        rank = [ord(c) for c in s]
        tmp = [0] * n
        sa = list(range(n))

        while k < n:
            sa.sort(key=lambda x: (rank[x], rank[x + k] if x + k < n else -1))

            tmp[sa[0]] = 0
            for i in range(1, n):
                prev = sa[i - 1]
                curr = sa[i]
                tmp[curr] = tmp[prev] + (
                    rank[prev] != rank[curr]
                    or (rank[prev + k] if prev + k < n else -1)
                    != (rank[curr + k] if curr + k < n else -1)
                )
            rank, tmp = tmp, rank
            k <<= 1

            if rank[sa[-1]] == n - 1:
                break

        return sa

    def get_reads(self):
        """
        Extract and return read sequences and quality scores for the reads.

        Returns:
            - reads_1: A dictionary mapping read IDs to sequences from first set of short reads.
            - quals_1: A dictionary mapping read IDs to quality scores from first set of short reads.
            - reads_2: A dictionary mapping read IDs to sequences from second set of short reads.
            - quals_2: A dictionary mapping read IDs to quality scores from second set of short reads.
        """

        # extract read sequences and quality scores for paired-end reads
        reads_1 = {
            read_id: seq_qual[0] for read_id, seq_qual in self.short_reads_1.items()
        }
        quals_1 = {
            read_id: list(seq_qual[1])
            for read_id, seq_qual in self.short_reads_1.items()
        }
        reads_2 = {
            read_id: seq_qual[0] for read_id, seq_qual in self.short_reads_2.items()
        }
        quals_2 = {
            read_id: list(seq_qual[1])
            for read_id, seq_qual in self.short_reads_2.items()
        }
        return reads_1, quals_1, reads_2, quals_2
