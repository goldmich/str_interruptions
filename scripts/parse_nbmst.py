#!/usr/bin/env python

import common
import pdb
from pathlib import Path
import numpy as np

# one list where index is floor(order in fasta / 2) and element is title (probably only need chrom:pos)
# then once I know the number of elements (* 2) make np array of 0s with 12 columns and len(list) rows
# iterate through each motif type in order of a list defined here
# every time the element number is even, add to even column
# if element number is odd, add to odd column
# then once processing is done, subtract odd columns from even columns to get delta for odd columns
# 	i want -3 to represent that the interrupted allele lost three bases of non-B forming potential and 3 to represent that it gained

fasta_file = snakemake.input['fasta_file']
# fasta_file = '/scratch/ucgd/lustre-labs/quinlan/u6045601/interruptions/tmp/allele_perfect_seq.chr20.fa'

motif_file_list = snakemake.input['motif_files']
# tmp_path = Path('/scratch/ucgd/lustre-labs/quinlan/u6045601/interruptions/tmp')
# motif_file_list = list(tmp_path.glob('*allele_perfect_seq.chr20*.tsv'))

seq_list = []
with open(fasta_file, 'r') as open_fasta:
	for l in open_fasta.readlines():
		if l.startswith(">") and l.endswith("obs\n"):
			seq_list.append(l.rstrip('\n'))

motif_table = np.zeros((len(seq_list), 12))
motif_list = ['DR', 'MR', 'IR', 'APR', 'GQ', 'Z']
motif_dict = dict(zip(motif_list, range(6)))

# sequences are 1-indexed
curr_motif = ''
curr_index = -1
for motif_file in motif_file_list:
	curr_motif = motif_file.split('motif_')[1].split('.')[0]
	curr_index = motif_dict[curr_motif]
	with open(motif_file, 'r') as open_motif_file:
		header = open_motif_file.readline()
		seq_index = 0
		seq_len = 0
		seq_shape_pos = set([])
		curr_seq_start = 151
		curr_seq_end = -1 # will need to initialize this
		for l in open_motif_file.readlines():
			l_split = l.rstrip('\n').split('\t')
			curr_seq_index = int(l_split[0].split('_')[1])
			if seq_index == 0:
				seq_index = curr_seq_index
				seq_len = len(seq_list[int((seq_index - 1) / 2)].split('|')[2])
				curr_seq_start = 151
				curr_seq_end = 150 + seq_len # inclusive, 1-indexed end
			if seq_index != curr_seq_index:
				if len(seq_shape_pos) > seq_len:
					pdb.set_trace()
				motif_table[int(np.floor((seq_index - 1) / 2)), (curr_index * 2) + ((seq_index - 1) % 2)] += len(seq_shape_pos)
				# iterate
				seq_index = curr_seq_index
				seq_len = len(seq_list[int(np.floor((seq_index - 1) / 2))].split('|')[2])
				curr_seq_start = 151
				curr_seq_end = 150 + seq_len # inclusive, 1-indexed end
				seq_shape_pos = set([])
			curr_shape_start = int(l_split[3])
			curr_shape_end = int(l_split[4])
			if curr_shape_start > curr_seq_end or curr_shape_end < curr_seq_start:
				continue
			curr_overlap_start = max(curr_seq_start, curr_shape_start)
			curr_overlap_end = min(curr_seq_end, curr_shape_end)
			seq_shape_pos = seq_shape_pos.union(set(range(curr_overlap_start, curr_overlap_end + 1)))
			
		# sequence positions are 1-indexed and end is inclusive i.e. 151 - 172 is 22bp long segment

# need to subtrack each odd indexed column from even
for x in range(len(motif_dict)):
	motif_table[:, (x * 2) + 1] = motif_table[:, (x * 2)] - motif_table[:, (x * 2) + 1]

with open(snakemake.output[0], 'w') as open_output:
	header = '#chrom\tpos\t' + '\t'.join([z for x in sorted(motif_dict.items(), key = lambda y: y[1]) for z in [x[0] + '_obs', x[0] + '_delta_perfect']])
	open_output.write(header + '\n')
	for i, x in enumerate(seq_list):
		x_split = x.split('|')
		curr_chrom = x_split[0][1:]
		curr_pos = x_split[1]
		open_output.write(curr_chrom + '\t' + curr_pos + '\t' + '\t'.join([str(int(y)) for y in motif_table[i, :]]) + '\n')

