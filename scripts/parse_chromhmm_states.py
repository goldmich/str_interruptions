#!/usr/bin/env python

import pdb

locus_characteristics = snakemake.input[0]
output_file = snakemake.output[0]

def parse_chromhmm_state(state_string):
	if state_string == '.':
		return 'NA'
	state_set = set([int(x.lstrip('E')) for x in state_string.split(',')])
	if state_set.intersection(set(range(5, 10))):
		return 'tx'
	elif state_set.intersection(set([1])):
		return 'tss'
	elif state_set.intersection(set(range(2, 5))):
		return 'prom'
	elif state_set.intersection(set(range(13, 16))): # 10-12 may not be relevant
		return 'enh_str'
	elif state_set.intersection(set([10, 11, 12, 16, 17, 18])):
		return 'enh_wk'
	elif state_set.intersection(set([19, 20, 22, 23])):
		return 'other'
	else:
		return 'het'

with open(locus_characteristics, 'r') as open_locus_characteristics:
	with open(output_file, 'w') as open_outfile:
		outlines = ''
		counter = 0
		for l in open_locus_characteristics.readlines():
			counter += 1
			l_split = l.rstrip('\n').split('\t')
			outlines += '\t'.join(l_split + [parse_chromhmm_state(l_split[-1])]) + '\n' # for now i will keep the chromhmm state number
			if counter == 1000:
				open_outfile.write(outlines)
				outlines = ''
				counter = 0
		open_outfile.write(outlines)

