#!/usr/bin/env python

import pdb
import gzip
import random
# import re

gff_file = snakemake.input[0]
outfile = snakemake.output[0]
# autosomes = ['chr' + str(x) for x in range(1, 23)]
autosomes = [str(x) for x in range(1, 23)]
cds_dict = {}

# instead of processing each CDS separately, I now need to keep track of all segments per gene and calculate all reading frames once i'm done with the gene
with gzip.open(gff_file, 'r') as open_gff:
	for l in open_gff.readlines():
		l_split = l.decode('ascii').rstrip('\n').split('\t')
		if l_split[0].startswith("#") or not l_split[0].startswith('NC'):
			continue
		curr_chrom = str(int(l_split[0].split('_')[1].split('.')[0]))
		if curr_chrom not in autosomes:
			continue
		if l_split[2] == 'CDS':
			if curr_chrom not in cds_dict:
				cds_dict[curr_chrom] = []
			curr_info = [x.lstrip(' ') for x in l_split[8].split('"; ')]
			curr_info_dict = dict([tuple(x.split(' "')) for x in curr_info if x])
			try:
				curr_transcript = curr_info_dict['transcript_id'].rstrip('"')
				curr_gene = curr_info_dict['gene_id'].rstrip('"')
			except KeyError:
				pdb.set_trace()
			if not curr_transcript.startswith('NM_'):
				continue
			try:
				curr_tag = curr_info_dict['tag'].rstrip('"')
			except KeyError:
				continue # not a canonical transcript
			# for now i am including all RefSeq and MANE Select and Plus Clinical tags
			if 'Select' not in curr_tag and 'Plus Clinical' not in curr_tag:
				pdb.set_trace()
			if not curr_transcript:
				pdb.set_trace()
			cds_dict[curr_chrom].append(['chr' + curr_chrom, l_split[3], l_split[4], curr_gene, curr_transcript, l_split[6], l_split[7]])

# the segments seem inclusive (i.e. include both start and end base positions)
# this code block seems no longer necessary as column 7 (python) contains the reading frame of the 5' end of the codon

# for chrom in cds_dict:
# 	for curr_transcript in cds_dict[chrom]:
# 		curr_segs = [x[1:3] for x in cds_dict[chrom][curr_transcript]]
# 		is_negative_strand = (cds_dict[chrom][curr_transcript][0][4] == '-')
# 		curr_segs = sorted(curr_segs, key = lambda x: x[0], reverse = is_negative_strand)
# 		curr_frames = []
# 		counter = 0
# 		for seg in curr_segs:
# 			curr_frames.append(counter % 3)
# 			counter += (seg[1] - seg[0]) + 1 # accounting for inclusive segments
# 		if is_negative_strand:
# 			curr_frames.reverse()
# 		cds_dict[chrom][curr_transcript] = [x + [y] for x, y in zip(cds_dict[chrom][curr_transcript], curr_frames)]
# 		pdb.set_trace()
		

with open(outfile, 'w') as open_outfile:
	for chrom in autosomes:
		open_outfile.write('\n'.join(['\t'.join(x) for x in sorted(cds_dict[chrom], key = lambda x: int(x[1]))]))
		open_outfile.write('\n')



# 			[x.split(' "')[1].replace('"', '') for x in curr_info if x.split(' "')[0].startswith('transcript')][0]
# 			curr_gene = [x.split(' "')[1].replace('"', '') for x in curr_info if x.split(' "')[0].startswith('gene')][0]
