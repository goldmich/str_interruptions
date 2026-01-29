import pysam
import argparse
import gzip
import pdb
import pandas as pd
from collections import Counter

def get_bam(sample, manifest):
	"""
	Retrieves the path to the sample's bam file.
	:param sample: Sample to find bam file for.
	:returns: AlignmentFile object for the sample's bam.
	"""
	row = manifest[manifest['sample_id'] == sample]
	bam_path = row.bam.values[0]
	bam = pysam.AlignmentFile(bam_path, "rc")
	return bam

def get_family_bams(manifest, family_id):
	"""
	Reads in all the alignment files (bams or crams) for a family.
	:param manifest: Pandas DataFrame with ssc_id, sample_id, and path to bam or cram file for each sample.
	:param family_id: Integer of the unique ID for each family.
	:return: lLst of (name of member, AlignmentFile for bam) for mom, and dad.
	"""

	mom = ("mom", get_bam("%s.mo" % family_id, manifest))
	dad = ("dad", get_bam("%s.fa" % family_id, manifest))

	return [mom, dad]

def get_header(file):
	"""
	Read the header of the VCF, and add two new lines for the vcf info added in this script.
	:param file: Input VCF file.
	:return: List of each line in the header.
	"""
	header = []
	for line in file:
		if line.startswith("#"):
			if len(header) > 0:
				if header[-1].startswith("##FORMAT") and line.startswith("##INFO"): # adding another format line
					header.append(
						'##FORMAT=<ID=ALLELE_SUPPORT,Number=A,Type=String,Description="Alternate allele read counts">')
				elif header[-1].startswith("##INFO") and line.startswith("##FILTER"): # adding another info line
					header.append(
						'##INFO=<ID=UNCALLABLE_DNMS,Number=A,Type=Integer,Description="Lengths of uncallable dnm alleles">')
			header.append(line.rstrip())
		else:
			if not header:
				raise ValueError(f"{file} has malformed header")
			return header, line

def get_site_info(vcf_row, dad_idx, mom_idx):
	"""
	Extract relevant data about site from VCF line.
	:param vcf_row: String of site data from VCF.
	:param dad_idx: The index of the dad's genotype information.
	:param mom_idx: The index of the mom's genotype information.
	:return: Chromosome, position, motif length (period), motif bases, reference motif count, dad gentoype, mom genotype.
	"""
	info = vcf_row[7].split(';')
	period = int([x for x in info if 'PERIOD=' in x][0].split('=')[1])
	ref = int([x for x in info if 'REF=' in x][0].split('=')[1])
	repeat_unit = [x for x in info if 'RU=' in x][0].split('=')[1].upper()

	gt_idx = vcf_row[8].split(":").index("REPCN")

	dad = vcf_row[dad_idx].split(':')[gt_idx]
	mom = vcf_row[mom_idx].split(':')[gt_idx]

	return vcf_row[0], int(vcf_row[1]), period, repeat_unit, ref, dad, mom

def get_longest_genotype(dad, mom, motif_len):
	"""
	Find the length of longest genotype.
	:param dad: Dad's genotype.
	:param mom: Mom's genotype.
	:param motif_len: Integer length of the repeat unit.
	:return: Integer of the length of longest genotype (STR motif count multiplied by motif length).
	"""
	if '.' in dad:
		return -1
	elif '.' in mom:
		return -1

	genotypes = [int(x) for x in dad.split(',')] + [int(x) for x in mom.split(',')]
	return max(genotypes) * motif_len

def get_soft_clipping(read):
	"""
	Gets the number of soft clipped bases on either end of the read.
	:param read: AlignedSegment of sequencing read.
	:return: Integer lengths of leading and tailing soft clipped bases.
	"""
	cigar = read.cigartuples

	if cigar[0][0] == 4:
		leading_clip = cigar[0][1]
	else:
		leading_clip = 0

	if cigar[-1][0] == 4:
		tailing_clip = cigar[-1][1]
	else:
		tailing_clip = 0

	return leading_clip, tailing_clip

def check_cigar_start(read, str_start, longest_gt, padding):
	"""
	Checks to see if their are bases that do not match the reference before the STR start.
	:param read: AlignedSegment of sequencing read.
	:param str_start: Integer naive prediction of STR starting position in read.
	:param longest_gt: Integer length of longest observed STR genotype.
	:param padding: Integer length of padding to add to either side of STR sequence.
	:return: Integer STR start and end positions adjusted to accommodate for inserted/deleted bases earlier in the read.
	"""
	start_pos = 0
	adjusted_str_start = str_start

	for operation, length in read.cigartuples:
		# Once we reach the start, we no longer need to adjust its location
		if adjusted_str_start - 1 < start_pos + length :
			left_boundary = adjusted_str_start - padding
			right_boundary = adjusted_str_start + longest_gt + padding
			# If there's an insertion or deletion that spans the starting coordinate, include those bases in the sequence
			if operation in (1,2,3):
				left_boundary -= length
				right_boundary += length
			break
		# If there is an insertion (1) or soft clip (4), the position of the STR will be later in the read
		if operation in (1, 4):
			adjusted_str_start += length
		# If there is a deletion (2) or reference skip (3), the position of the STR will be earlier in the read
		elif operation in (2, 3):
			adjusted_str_start -= length

		start_pos += length

	return left_boundary, right_boundary

def get_read_sequence(read, str_coordinate, longest_gt, padding):
	"""
	Checks that reads pass mapq and span the STR, and pulls out the STR sequence with padding.
	:param read: AlignedSegment of sequencing read.
	:param str_coordinate: Integer for STR starting coordinate in reference.
	:param longest_gt: Number of bases in longest observed STR genotype.
	:param padding: Number of bases to add on either side of the STR.
	:return: String for STR sequence, with length of longest_gt + 2*padding, or None if read does not pass filters.
	"""
	# Filter on some of the same basic features that GangSTR uses
	if read.mapping_quality != 60 or read.is_secondary or read.is_supplementary:
		return

	# Find the start and end coordinates of the STR sequence with padding, based on the read alignment and CIGAR string
	naive_start = str_coordinate - read.reference_start - 1
	left, right = check_cigar_start(read, naive_start, longest_gt, padding)
	left_clip, right_clip = get_soft_clipping(read)

	# Check that the read spans the STR sequence and padding, excluding soft-clipped bases
	if left < left_clip or right > read.query_length - right_clip:
		return

	return read.query_sequence[left:right]


def check_interruption(sequence, start_pos, padding, prev_sequence, next_sequence, str_length):
	"""
	Check to see if the sequence following the STR matches the expected reference sequence.
	:param sequence: Read sequence.
	:param start_pos: Position of the base that interrupts the STR. Should be STR end + 1 for uninterrupted STR.
	:param padding: Padding around motif sequence.
	:param prev_sequence: The expected reference sequence preceding the STR.
	:param next_sequence: The expected reference sequence following the STR.
	:return: True if the following sequence matches expectation, False if it does not.
	"""
	pdb.set_trace()
	start, end = True, True
	
	start_sequence = sequence[start_pos-padding-str_length + 1: start_pos-str_length]
	end_sequence = sequence[start_pos: start_pos + padding - 1]
	if len(start_sequence) != padding - 1:
		start = False
	if len(end_sequence) != padding - 1:
		end = False

	if start or end:
		start_differences = sum(c1 != c2 for c1, c2 in zip(start_sequence, prev_sequence[-(padding -1):]))
		end_differences = sum(c1 != c2 for c1, c2 in zip(end_sequence, next_sequence[:padding-1]))

		# Allow for up to 1 difference between expected and observed to account for base calling error around STR.
		# If there is an error in the first motif copy, the first motif will be included in the preceding sequence, so you should never get a match.
		if start_differences > 1:
			start = False
		# If it is an interrupted repeat, the base with the error will always be the first base, so you should never get a match with the subsequent sequence, regardless of similarity.
		if end_differences > 1:
			end = False

	return start, end

def find_motif_count(sequence, motif, longest_gt, padding, prev_seq, next_seq):
	"""
	Counts the number of contiguous repeat motifs in a read.
	:param sequence: String of the sequence from a read.
	:param motif: String of the STR repeat motif.
	:param padding: Integer distance from the STR reference start coordinate to the start of the sequence.
	:param next_seq: The 5 bp of sequence following the end of the STR in the reference.
	:return: Integer count of contiguous repeats of STR motif.
	"""
	start = end = 0
	count = 0

	# extra_padding = len(sequence) - (longest_gt + 2 * padding)
	extra_padding = 1

	while start <= len(sequence) - len(motif) :
		end += 1
		if sequence[end - 1] == motif[end - start - 1]:
			if end - start == len(motif):
				# Found an instance of the motif, count it.
				count += 1
				start = end
		else:
			# This is intended to address any alignment issues, so if an STR begins at a coordinate before the reference, those motifs will still be counted.
			if end - 1 <= padding + extra_padding:
				# Motif in the sequence before STR coordinates begin, and it is not contiguous with the STR.
				start = end = start + 1
				count = 0
			elif count == 0:
				# Keep going because you have not found the STR start yet.
				start = end = start + 1
			elif count != 0:
				# If there is a base error in the STR sequence, we don't want to count the length of that STR.
				start_ok, end_ok = check_interruption(sequence, start, padding, prev_seq, next_seq, count*len(motif))
				if start_ok and end_ok:
					return count
				break
	return -1


def check_reads(bam_info, chrom, pos, longest_gt, motif, prev_seq, next_seq):
	"""
	Iterate through all reads aligned to STR coordinate and find the number of repeated motifs in each read.
	:param bam_info: Tuple of sample name and AlignmentFile object.
	:param chrom: Chromosome name.
	:param pos: STR start coordinate.
	:param longest_gt: Length of the longest observed STR genotype.
	:param motif: Repeated sequence in STR.
	:param next_seq: The sequence prior to the start of the STR in the reference.
	:param next_seq: The sequence following the end of the STR in the reference.
	:return: List of number of motifs observed in each read.
	"""
	sample_name, bam = bam_info
	motif_counts = []

	# padding = 3 * len(motif) + 1
	padding = 5 + 1

	try:
		for aligned_segment in bam.fetch(chrom, pos - 1, pos):
			str_sequence = get_read_sequence(aligned_segment, pos, longest_gt, padding)
			if str_sequence is not None:
				motif_ct = find_motif_count(str_sequence, motif, longest_gt, padding, prev_seq, next_seq)
				if motif_ct != -1 :
					# print(aligned_segment.query_name, " allele count = ", motif_ct)
					motif_counts.append(motif_ct)
				# else:
					# print(aligned_segment.query_name, " interrupted read ")
		return motif_counts

	# if the CRAM is corrupted, return NA
	except OSError:
		return ['NA']

def check_gt_support(dad, mom, dad_alleles, mom_alleles):
	"""
	Check to make sure the called parental alleles are present in read data.
	:param dad: List of dad's two called genotypes.
	:param mom: List of mom's two called genotypes.
	:param dad_alleles: Dictionary of all of dad's observed alleles and the number of reads supporting them.
	:param mom_alleles: Dictionary of all of mom's observed alleles and the number of reads supporting them.
	:return: True if all called alleles have been observed, false if at least one allele not observed.
	"""
	dad_gt1, dad_gt2 = dad
	mom_gt1, mom_gt2 = mom

	if dad_gt1 not in dad_alleles or dad_gt2 not in dad_alleles:
		return False

	if mom_gt1 not in mom_alleles or mom_gt2 not in mom_alleles:
		return False

	return True


def check_dnm_possibility(dad, mom, dad_alleles, mom_alleles):
	"""
	Get all callable DNMs in a 3 motif radius.
	:param dad: List of dad's genotypes.
	:param mom: List of mom's genotypes.
	:param dad_alleles: Dictionary of all of dad's observed alleles and the number of reads supporting them.
	:param mom_alleles: Dictionary of all of mom's observed alleles and the number of reads supporting them.
	:return: Tuple of all callable DNM lengths.
	"""
	observed = set()
	for parent in (dad_alleles, mom_alleles):
		observed.update(gt for gt, count in parent.items() if count > 1)

	possible = set()
	for parent in (dad, mom):
		for gt in parent:
			for mutations in (-3, -2, -1, 1, 2, 3):
				dnm = gt + mutations

				if dnm not in observed:
					possible.add(dnm)

	return tuple(possible)

def get_uncallable_dnms(dad_alleles, mom_alleles):
	"""
	List all alleles observed in at least two reads in one parent.
	:param dad_alleles: Dictionary of all of dad's observed alleles and the number of reads supporting them.
	:param mom_alleles: Dictionary of all of mom's observed alleles and the number of reads supporting them.
	:return: Tuple of all observed DNM lengths.
	"""
	observed = set()
	for parent in (dad_alleles, mom_alleles):
		observed.update(gt for gt, count in parent.items() if count >= 1)

	return tuple(observed)


def get_flanking_sequences(str_chrom, str_pos, motif_len, ref_len, reference):
	"""
	Get the reference base pairs before/after the STR.
	:param str_chrom: Chromosome of the STR.
	:param str_pos: Starting coordinate of the STR.
	:param motif_len: Length of the STR motif.
	:param ref_len: Number of STR copies in the reference genome.
	:param reference: PySam object of reference genome.
	:return: Strings of the 10 base pairs before and 10 base pairs after the STR.
	"""
	# Subtract 1 from the positions to correctly index bam with pysam
	start = str_pos - 1
	end = str_pos + (motif_len * ref_len) - 1
	prev_seq = reference.fetch(str_chrom, start - 10, start)
	next_seq = reference.fetch(str_chrom, end, end + 10)
	return next_seq.upper(), prev_seq.upper()


def get_alleles(allele_counts):
	"""
	Make a string of the observed alleles and how many reads support them.
	:param allele_counts: Dictionary of observed alleles and the number of reads supporting them.
	:return: String of alleles/read count, formatted to append to VCF column.
	"""
	allele_list = []
	for allele in allele_counts:
		allele_list.append(str(allele) + '=' + str(allele_counts[allele]))

	return ','.join(allele_list)

def check_variant(line, dad_idx, mom_idx):
	row = line.rstrip().split('\t')
	if row[6] != 'PASS':
		return None

	chrom, pos, motif_length, motif, ref_length, dad_gt, mom_gt = get_site_info(row, dad_idx, mom_idx)

	next_seq, prev_seq = get_flanking_sequences(chrom, pos, motif_length, ref_length, ref)

	longest_gt = get_longest_genotype(dad_gt, mom_gt, motif_length)

	if longest_gt == -1:
		return None

	dad_reads = Counter(check_reads(dad_bam, chrom, pos, longest_gt, motif, prev_seq, next_seq))
	mom_reads = Counter(check_reads(mom_bam, chrom, pos, longest_gt, motif, prev_seq, next_seq))

	dad_gts = tuple(int(x) for x in dad_gt.split(','))
	mom_gts = tuple(int(x) for x in mom_gt.split(','))

	if check_gt_support(dad_gts, mom_gts, dad_reads, mom_reads):
		# possibilities = check_dnm_possibility(dad_gts, mom_gts, dad_reads, mom_reads)
		impossibilities = get_uncallable_dnms(dad_reads, mom_reads)
		row[7] = row[7] + ";UNCALLABLE_DNMS=" + ",".join(str(x) for x in impossibilities)
		row[8] = row[8] + ":ALLELE_SUPPORT"
		row[dad_idx] = row[dad_idx] + ":" + get_alleles(dad_reads)
		row[mom_idx] = row[mom_idx] + ":" + get_alleles(mom_reads)
		return row

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-v", "--vcf", required=True, help="VCF of STR sites for a family")
	ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	ap.add_argument("-f", "--family", required=False, help="Family ID")
	ap.add_argument("-o", "--output", required=False, help="Output file")
	ap.add_argument("-r", "--reference", required=False, help="Reference genome fasta")

	args = ap.parse_args()
	bams = pd.read_csv(args.bams, sep='\t')
	ref = pysam.FastaFile(args.reference)

	dad_id = bams.loc[bams['sample_id'] == "{}.fa".format(args.family)]['ssc_id'].values[0]
	mom_id = bams.loc[bams['sample_id'] == "{}.mo".format(args.family)]['ssc_id'].values[0]

	mom_bam, dad_bam = get_family_bams(bams, args.family)

	with gzip.open(args.vcf, 'rt') as vcf, open(args.output, 'w') as outfile:
		vcf_header, first_line = get_header(vcf)
		print('\n'.join(vcf_header), file = outfile)
		index_header = vcf_header[-1].split('\t')
		dad_idx = index_header.index(dad_id)
		mom_idx = index_header.index(mom_id)

		out_row = check_variant(first_line, dad_idx, mom_idx)
		if out_row:
			print('\t'.join(out_row), file=outfile)


		for line in vcf:
			out_row = check_variant(line, dad_idx, mom_idx)
			if out_row:
				print('\t'.join(out_row), file=outfile)
