#!/usr/bin/env python

import pdb
import random
import common
from Levenshtein import distance


coding_annotations = snakemake.input[0]
locus_characteristics = snakemake.input[1]

vep_dict = {}

with open(coding_annotations, 'r') as open_coding_annotations:
	for l in open_coding_annotations.readlines():
		l_split = l.rstrip('\n').split('\t')
		if l_split[0] not in vep_dict:
			vep_dict[l_split[0]] = {}
		if l_split[1] not in vep_dict[l_split[0]]:
			vep_dict[l_split[0]][l_split[1]] = {}
		vep_dict[l_split[0]][l_split[1]][l_split[2]] = l_split[3]

outfile = snakemake.output[0]
curr_lines = ''
counter = 0
with open(outfile, 'w') as open_outfile:
	with open(locus_characteristics, 'r') as open_locus_data:
		for l in open_locus_data.readlines():
			l_split = l.rstrip('\n').split('\t')
			to_add = '.'
			# in a genic region, but has a LOEUF score and is not perfect
			if l_split[14] != '.' and l_split[16] != '-1' and l_split[7] != '0':
				curr_pos = l_split[1]
				curr_allele = l_split[4]
				if curr_pos not in vep_dict[l_split[0]]: # there are some off-by-one errors that vep introduces...
					if str(int(l_split[1]) + 1) in vep_dict[l_split[0]]:
						curr_pos = str(int(l_split[1]) + 1)
						curr_allele = curr_allele[1:]
					elif str(int(l_split[1]) - 1) in vep_dict[l_split[0]]:
						curr_pos = str(int(l_split[1]) - 1)
						pdb.set_trace()
					else:
						pass
				# catching the missing loci and missing alleles
				if curr_pos not in vep_dict[l_split[0]] or curr_allele not in vep_dict[l_split[0]][curr_pos]:
# 					pdb.set_trace() # there are more issues here but they are rare and tend to involve multiple individual "loci" at the same position... i am choosing to ignore them for now
					pass
				elif vep_dict[l_split[0]][curr_pos][curr_allele] == 'ref_allele':
					list_of_possible_perfect_alleles = [(ru_to_repeat * (len(curr_allele) + 1))[x : (len(curr_allele) + x)] for x in range(len(l_split[3])) for ru_to_repeat in [l_split[3], common.revcomp(l_split[3])]]
					min_distance_combo = min([(distance(curr_allele, x), x) for x in list_of_possible_perfect_alleles], key = lambda y: y[0])
					if min_distance_combo[1] not in vep_dict[l_split[0]][curr_pos]:
						if curr_pos != l_split[1]: # this is for the vep off-by-one errors
							list_of_possible_perfect_alleles_obo = [(ru_to_repeat * (len(l_split[4]) + 1))[x : (len(l_split[4]) + x)] for x in range(len(l_split[3])) for ru_to_repeat in [l_split[3], common.revcomp(l_split[3])]]
							min_distance_combo_obo = min([(distance(l_split[4], x), x) for x in list_of_possible_perfect_alleles_obo], key = lambda y: y[0])
							if min_distance_combo_obo[1][1:] not in vep_dict[l_split[0]][curr_pos]:
								pdb.set_trace()
							else:
								to_add = vep_dict[l_split[0]][curr_pos][min_distance_combo_obo[1][1:]]
						else:
							pdb.set_trace() # the perfect allele should be here
					else:
						to_add = vep_dict[l_split[0]][curr_pos][min_distance_combo[1]]
				else:
					to_add = vep_dict[l_split[0]][curr_pos][curr_allele]
			curr_lines += '\t'.join(l_split + [to_add]) + '\n' # i still need to keep the allele
			counter += 1
			if counter == 1000:
				open_outfile.write(curr_lines)
				curr_lines = ''
				counter = 0
		open_outfile.write(curr_lines)
			
