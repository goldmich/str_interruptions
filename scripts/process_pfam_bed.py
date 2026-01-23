#!/usr/bin/env python

import pdb

infile = snakemake.input[0]
outfile = snakemake.output[0]

with open(outfile, 'w') as open_outfile:
	with open(infile, 'r') as open_infile:
		outlines = ''
		counter = 0
		for l in open_infile.readlines():
			l_split = l.rstrip('\n').split('\t')
			curr_start = int(l_split[1])
			curr_segs = zip(l_split[10].rstrip(',').split(','), l_split[11].rstrip(',').split(','))
			outlines += '\n'.join(['\t'.join([l_split[0], 
												str(curr_start + int(seg_offset)), 
												str(curr_start + int(seg_offset) + int(seg_length)),
												l_split[3]]) 
							for seg_length, seg_offset in curr_segs])
			outlines += '\n'
			counter += 1
			if counter > 1000:
				counter = 0
				open_outfile.write(outlines)
				outlines = ''
		open_outfile.write(outlines)

	