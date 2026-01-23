#!/usr/bin/env python

import pdb

locus_characteristics = snakemake.input[0]
tf_chip_table = snakemake.input[1]

output_file = snakemake.output[0]

tf_dict = {}

with open(tf_chip_table, 'r') as open_tf_table:
	header = open_tf_table.readline()
	for l in open_tf_table.readlines():
		l_split = l.rstrip('\n').split(',')
		tf_dict[l_split[2]] = l_split[0]

with open(locus_characteristics, 'r') as open_locus_characteristics:
	with open(output_file, 'w') as open_outfile:
		curr_allele = None
		curr_tf_list = []
		outlines = ''
		counter = 0
		last_line = None
		for l in open_locus_characteristics.readlines():
			l_split = l.rstrip('\n').split('\t')
# 			if l_split[0] == 'chr22' and l_split[1] == '50733151' and l_split[4] == 'TATTATTATTATTGATAATAATTATAATATTATTAT' and l_split[18] == '/scratch/ucgd/lustre-labs/quinlan/u6045601/interruptions/tf_chip_idr_peaks/ENCFF721NYV.bed.gz':
# 				pdb.set_trace()
			if not curr_allele:
				curr_allele = (l_split[0], l_split[1], l_split[4])
			if (l_split[0], l_split[1], l_split[4]) != curr_allele:
				curr_allele = (l_split[0], l_split[1], l_split[4])
				tf_to_add = '.'
				if curr_tf_list:
					tf_to_add = ','.join(curr_tf_list)
					curr_tf_list = []
				outlines += '\t'.join(last_line[:-2] + [tf_to_add]) + '\n'
				counter += 1
				if counter == 1000:
					open_outfile.write(outlines)
					outlines = ''
					counter = 0
				continue
				# process here
			if l_split[19] != '0':
				curr_tf_list.append(tf_dict[l_split[18].split('/')[-1].split('.')[0]])
			last_line = l_split
		tf_to_add = '.'
		if curr_tf_list:
			tf_to_add = ','.join(tf_list)
		outlines += '\t'.join(last_line[:-2] + [tf_to_add]) + '\n'
		open_outfile.write(outlines)

