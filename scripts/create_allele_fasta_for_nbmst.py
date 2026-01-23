#!/usr/bin/env python

import common
import pdb
import pysam

allele_seq_file = snakemake.input['allele_perfect_seq']
ref_complete = pysam.FastaFile(snakemake.input['ref_genome'])

curr_chr = allele_seq_file.split('/')[-1].split('.')[-2]
ref = ref_complete[curr_chr]

output_list = []
with open(allele_seq_file, 'r') as open_allele_file:
	header = open_allele_file.readline()
	curr_start = -1
	curr_max_af = 0.
	curr_af = 0.
	curr_max_af_line = []
	for l in open_allele_file.readlines():
		l_split = l.rstrip('\n').split('\t')
		if curr_start < 0:
			curr_start = int(l_split[1])
		if int(l_split[1]) != curr_start:
			# process
			# chr|pos|allele_seq|obs_or_perfect
			curr_desc_obs = '>%s|%s|%s|obs' % (curr_chr, curr_max_af_line[1], curr_max_af_line[3])
			curr_desc_perf = '>%s|%s|%s|perf' % (curr_chr, curr_max_af_line[1], curr_max_af_line[3])
			curr_end = int(curr_max_af_line[2])
			upstream_seq = ref[curr_start - 151: curr_start - 1]
			downstream_seq = ref[curr_end - 1: curr_end + 149]
			output_list += [curr_desc_obs, 
							upstream_seq + curr_max_af_line[3] + downstream_seq,
							curr_desc_perf,
							upstream_seq + curr_max_af_line[4] + downstream_seq]
			curr_max_af = 0
			curr_start = int(l_split[1])
		curr_af = float(l_split[5])
		if curr_af > curr_max_af:
			curr_max_af_line = l_split
			curr_max_af = curr_af
	# process once more
	curr_desc_obs = '%s|%s|%s|obs' % (curr_chr, curr_max_af_line[1], curr_max_af_line[3])
	curr_desc_perf = '%s|%s|%s|perf' % (curr_chr, curr_max_af_line[1], curr_max_af_line[3])
	curr_end = int(curr_max_af_line[2])
	upstream_seq = ref[curr_start - 151: curr_start - 1]
	downstream_seq = ref[curr_end - 1: curr_end + 149]
	output_list += [curr_desc_obs, 
					upstream_seq + curr_max_af_line[3] + downstream_seq,
					curr_desc_perf,
					upstream_seq + curr_max_af_line[4] + downstream_seq]
		
with open(snakemake.output[0], 'w') as open_outfile:
	open_outfile.write('\n'.join(output_list))
