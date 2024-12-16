import math

def write_sam_file(sam_filename, read_1_alignments, read_2_alignments, read_mapper, reads_1, reads_2, quals_1, quals_2, reference_name):
    """
    Write alignments to a SAM file.

    Parameters:
        sam_filename (str): The name of the SAM file to write.
        read_1_alignments (dict): Alignments for reads from the first fastq file.
        read_2_alignments (dict): Alignments for reads from the second fastq file.
    """
 
    with open(sam_filename, "w") as samfile:
        samfile.write("@HD\tVN:1.0\tSO:unsorted\n")

        for read_id, alignment in read_1_alignments.items():

            if alignment != None:
                # Create a new AlignedSegment object
                start_1, end_1, strand_1 = alignment[0], alignment[1], alignment[2]
                read_seq = reads_1[read_id]
                
                qual = "".join(quals_1[read_id])
                #Compute FLAG
                mate_alignment = read_2_alignments.get(read_id)
                mate_unmapped = mate_alignment is None or mate_alignment[0] is None
                mate_reverse = mate_alignment[2] == "-" if mate_alignment else False
                is_paired = True  # Assuming paired-end data
                is_mapped = alignment is not None
                is_reverse = strand_1 == "-"

                is_first = True  # First read in pair
                is_second = False  # Not the second read in pair
                proper_pair = is_paired and is_mapped and not mate_unmapped and not is_reverse and not mate_reverse
                flag1 = compute_flag(
                     is_paired=True, is_mapped=True, is_reverse=is_reverse, is_first=True, is_second=False, mate_unmapped=mate_unmapped, mate_reverse=mate_reverse, proper_pair=True)
                qual_1_score = quals_1[read_id]
                quality_scores = [ord(char) - 33 for char in qual_1_score] # Convert quality string to Phred scores
                alignment_score = sum(quality for ref_base, read_base, quality in zip(read_mapper.reference_genome[start_1:end_1], read_seq, quality_scores) if ref_base == read_base)
                max_alignment_score = sum(quality_scores)
                mapq = compute_mapq(alignment_score, max_alignment_score)
                #Compute CIGAR
                cigar_string = read_mapper.compute_cigar_string(read_seq, start_1)
                #Compute SEQ
                SEQ = read_seq
                if not cigar_string:
                    cigar_string = "*"

                samfile.write(
                    f"{read_id}\t{flag1}\t{reference_name}\t{start_1 + 1}\t{mapq}\t{cigar_string}\t{"="}\t{99}\t{99}\t{SEQ}\t{qual}\n"
                )
            else:
                continue

        for read_id, alignment in read_2_alignments.items():
            if alignment != None:
                # Create a new AlignedSegment object
                start_2, end_2, strand_2 = alignment[0], alignment[1], alignment[2]
                read_seq = reads_2[read_id]
                qual = "".join(quals_2[read_id])

                # Compute FLAG
                mate_alignment = read_1_alignments.get(read_id)
                mate_unmapped = mate_alignment is None or mate_alignment[0] is None
                mate_reverse = mate_alignment[2] == "-" if mate_alignment else False
                is_paired = True  # Assuming paired-end data
                is_mapped = alignment is not None
                is_reverse = strand_2 == "-" 
                is_first = False  # Not the first read in pair
                is_second = True  # Second read in pair
                proper_pair = is_paired and is_mapped and not mate_unmapped and not is_reverse and not mate_reverse

                flag2 = compute_flag(
                     is_paired=True, is_mapped=True, is_reverse=is_reverse, is_first=False, is_second=True, mate_unmapped=mate_unmapped, mate_reverse=mate_reverse, proper_pair=True) 
                # Compute MAPQ
                qual_2_score = quals_2[read_id]
                quality_scores = [ord(char) - 33 for char in qual_2_score] # Convert quality string to Phred scores
                alignment_score = sum(quality for ref_base, read_base, quality in zip(read_mapper.reference_genome[start_2:end_2], read_seq, quality_scores) if ref_base == read_base)
                max_alignment_score = sum(quality_scores)
                mapq = compute_mapq(alignment_score, max_alignment_score)
                #Compute CIGAR
                cigar_string = read_mapper.compute_cigar_string(read_seq, start_2)
                #Compute SEQ
                SEQ = read_seq
                if not cigar_string:
                    # Handle poor alignments
                    cigar_string = "*"

                samfile.write(
                    f"{read_id}\t{flag2}\t{reference_name}\t{start_2 + 1}\t{mapq}\t{cigar_string}\t{"="}\t{99}\t{99}\t{SEQ}\t{qual}\n"
                )
            else:
                continue

"""
Compute MAPQ field 
Parameters:
alignment_score: integer
max_alignment_score: integer

Return:
Mapping quality score
"""
def compute_mapq(alignment_score, max_alignment_score):

	if alignment_score is None or alignment_score <= 0 or max_alignment_score <= 0:
		return 255  # MAPQ not available
	prob_wrong = 1 - (alignment_score / max_alignment_score)
	if prob_wrong <= 0:
		return 99  # Assign MAPQ to 99 if perfectly aligned
	return max(0, min(99, round(-10 * math.log10(prob_wrong))))		# highest MAPQ score without perfect alignment is 33 (lower MAPQ = more possible alignment possibilities)



def compute_flag(is_paired, is_mapped, is_reverse, is_first, is_second, mate_unmapped=False, mate_reverse=False, proper_pair=False):
    """
    Compute the FLAG value for a SAM entry.

    Parameters:
        is_paired (bool): Whether the read is part of a pair.
        is_mapped (bool): Whether the read is mapped.
        is_reverse (bool): Whether the read is mapped to the reverse strand.
        is_first (bool): Whether this is the first read in a pair.
        is_second (bool): Whether this is the second read in a pair.

    Returns:
        int: The computed FLAG value.
    """
    flag = 0
    if is_paired:
        flag |= 0x1  # Template has multiple segments in sequencing
    if proper_pair and is_paired:
        flag |= 0x2  # Each segment properly aligned
    if not is_mapped:
        flag |= 0x4  # Segment unmapped
    if mate_unmapped:
        flag |= 0x8  # Next segment in the template unmapped
    if is_reverse:
        flag |= 0x10  # SEQ is reverse complemented
    if mate_reverse:
        flag |= 0x20  # SEQ of the next segment is reverse complemented
    if is_first:
        flag |= 0x40  # This is the first segment in the template
    if is_second:
        flag |= 0x80  # This is the last segment in the template
    return flag
