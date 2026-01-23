#!/usr/bin/env python

import pdb
import random

loeuf_file = snakemake.input[0]
ad_gene_list = snakemake.input[1]
locus_characteristics = snakemake.input[2]

loeuf_dict = {}

by_transcript = False

with open(loeuf_file, 'r') as open_loeuf:
	for l in open_loeuf.readlines():
		l_split = l.rstrip('\n').split('\t')
		curr_key = ''
		if by_transcript:
			curr_key = l_split[4]
		else:
			curr_key = l_split[3]
		curr_loeuf = float(l_split[5])
		loeuf_dict[curr_key] = curr_loeuf

ad_genes = []
with open(ad_gene_list, 'r') as open_ad_genes:
	ad_genes = open_ad_genes.read().split('\n')
	ad_genes.pop(-1) # empty line


outfile = snakemake.output[0]
match_column = -1
if by_transcript:
	match_column = 15
else:
	match_column = 14
with open(outfile, 'w') as open_outfile:
	with open(locus_characteristics, 'r') as open_locus_data:
		curr_lines = ''
		counter = 0
		for l in open_locus_data.readlines():
			any_ad = -1
			l_split = l.rstrip('\n').split('\t')
			if l_split[match_column] != '.':
				curr_connections = l_split[match_column].split(',')
				curr_loeuf_scores = [loeuf_dict[x] for x in curr_connections if x in loeuf_dict]
				if not curr_loeuf_scores:
					curr_lowest_loeuf = '-1'
				else:
					curr_lowest_loeuf = str(min(curr_loeuf_scores))
				any_ad = int(any([x in ad_genes for x in curr_connections]))
			else:
				curr_lowest_loeuf = '-1'
			curr_lines += '\t'.join(l_split + [curr_lowest_loeuf, str(any_ad)]) + '\n'
			counter += 1
			if counter == 1000:
				open_outfile.write(curr_lines)
				curr_lines = ''
				counter = 0
		open_outfile.write(curr_lines)
