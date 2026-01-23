#!/usr/bin/env python

import pdb

loeuf_file = snakemake.input[0]
ad_gene_list = snakemake.input[1]
locus_characteristics = snakemake.input[2]

loeuf_dict = {}

with open(loeuf_file, 'r') as open_loeuf:
	for l in open_loeuf.readlines():
		l_split = l.rstrip('\n').split('\t')
		curr_gene = l_split[3]
		curr_loeuf = float(l_split[5])
		if curr_gene in loeuf_dict:
			if curr_loeuf < loeuf_dict[curr_gene]:
				loeuf_dict[curr_gene] = curr_loeuf # finding the lowest scored transcript
		else:
			loeuf_dict[curr_gene] = curr_loeuf

ad_genes = []
with open(ad_gene_list, 'r') as open_ad_genes:
	ad_genes = open_ad_genes.read().split('\n')
	ad_genes.pop(-1) # empty line

outfile = snakemake.output[0]
with open(outfile, 'w') as open_outfile:
	with open(locus_characteristics, 'r') as open_locus_data:
		curr_lines = ''
		counter = 0
		for l in open_locus_data.readlines():
			any_ad = -1
			l_split = l.rstrip('\n').split('\t')
			is_elite = l_split[11].split(',')
			if '1' in is_elite:
				curr_connections = [y for x in l_split[12].split(',') for y in x.split(';')]
				curr_elite_connections = [x.split('|')[0] for x in curr_connections if x.split('|')[2] == '1']
				curr_loeuf_scores = [loeuf_dict[x] for x in curr_elite_connections if x in loeuf_dict]
				any_ad = int(any([x in ad_genes for x in curr_elite_connections]))
				if not curr_loeuf_scores:
					curr_lowest_loeuf = '-1'
				else:
					curr_lowest_loeuf = str(min(curr_loeuf_scores))
			else:
				curr_lowest_loeuf = '-1'
			# cutting down on size by removing the bulky loeuf line
			curr_lines += '\t'.join(l_split[:12] + [curr_lowest_loeuf, str(any_ad)]) + '\n'
			counter += 1
			if counter == 1000:
				open_outfile.write(curr_lines)
				curr_lines = ''
				counter = 0
		open_outfile.write(curr_lines)
