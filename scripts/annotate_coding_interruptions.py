#!/usr/bin/env python

import pdb
import random
import common
import pysam
import math
from Levenshtein import distance

locus_characteristics = snakemake.input[0]
gtf_bed = snakemake.input[1]
ref_genome = snakemake.input[2]

output_file = snakemake.output[0]

exon_dict = {}

with open(gtf_bed, 'r') as open_gtf_bed:
	for line in open_gtf_bed.readlines():
		l_split = line.rstrip('\n').split('\t')
		if l_split[0] not in exon_dict:
			exon_dict[l_split[0]] = {} # by transcript
		if l_split[4] not in exon_dict[l_split[0]]:
			exon_dict[l_split[0]][l_split[4]] = []
		exon_dict[l_split[0]][l_split[4]].append([int(l_split[1]), int(l_split[2]), l_split[5], int(l_split[6])])

ref_complete = pysam.FastaFile(ref_genome)

with open(output_file, 'w') as open_output_file:
	output_list = []
	n_multiple = 0
	n_stop_gains = 0
	with open(locus_characteristics, 'r') as open_str_file:
		curr_chr = ''
		ref = ''
		fields_to_add = ['.', '.', '.']
		for line in open_str_file.readlines():
			# dump
			if len(output_list) > 1000:
				open_output_file.write('\n'.join(['\t'.join(x) for x in output_list]) + '\n')
				output_list = []
			fields_to_add = ['.', '.', '.', '.', '.']
			l_split = line.rstrip('\n').split('\t')
			# check modulo 3 here too
			if l_split[15] == '.' or len(l_split[3]) % 3 != 0: # now just looking at coding loci only; these have to be modulo 3
				continue
			if l_split[0] != curr_chr:
				curr_chr = l_split[0]
				ref = ref_complete[curr_chr]
			curr_start = int(l_split[1])
			curr_end = int(l_split[2])
			if ((curr_end - curr_start) % 3) != (len(l_split[4]) % 3):
				# there has been a frameshift relative to reference, i think i need to get rid of these
				output_list.append(l_split + fields_to_add)
				continue
			# get perfect STR here
			curr_ru = l_split[3]
			curr_allele = l_split[4]
			list_of_possible_perfect_alleles = [(ru_to_repeat * (len(curr_allele) + 1))[x : (len(curr_allele) + x)]
													 for x in range(len(curr_ru)) 
													 for ru_to_repeat in [curr_ru, common.revcomp(curr_ru)]]
# 			min_distance_combo = min([(distance(curr_allele, x), x) for x in list_of_possible_perfect_alleles], key = lambda y: y[0])
# 			min_hamming_distance_combo = min([(common.hamming_distance(curr_allele, x), x) for x in list_of_possible_perfect_alleles], key = lambda y: y[0])
# 			curr_perfect_allele = min_distance_combo[1]
			# which exon(s)?
			# i think i should be using curr start and curr end rather than the allele length given that it's ultimately mapped to the same coordinates as the gtf
			if ',' in l_split[15]:
				output_list.append(l_split + fields_to_add)
				n_multiple += 1
				continue
# 				pdb.set_trace() # handle this later
			for curr_tran in l_split[15].split(','):
				curr_exon_list = exon_dict[l_split[0]][curr_tran]
				exon_counter = 0
				while curr_exon_list[exon_counter][1] < curr_start:
					exon_counter += 1
				to_add_5prime = 0
				to_add_3prime = 0
				add_seq_5prime = ''
				add_seq_3prime = ''
				# remember: phase here is funny! 1 means that there is one extra 5' base before the first full codon begins; 2 means that there are 2 extra bases
				# i may run into some issues with STRs close to stop codons as the stop codon for the given transcript is not included in the coordinates of the 3'-most exon
				if curr_exon_list[exon_counter][2] == '+':
					to_add_5prime = (curr_start - curr_exon_list[exon_counter][0] + (abs(curr_exon_list[exon_counter][3] - 3) % 3)) % 3
					to_add_3prime = abs(((to_add_5prime + (curr_end - curr_start)) % 3) - 3) % 3
					if to_add_5prime > (curr_start - curr_exon_list[exon_counter][0]):
						add_seq_5prime_prior_exon = to_add_5prime - (curr_start - curr_exon_list[exon_counter][0])
						add_seq_5prime = ref[
												(curr_exon_list[exon_counter - 1][1] - add_seq_5prime_prior_exon): # this is really + 1 - 1 accounting for both gtf and python
												(curr_exon_list[exon_counter - 1][1])
												]
						to_add_5prime -= len(add_seq_5prime)
					add_seq_5prime += ref[(curr_start - to_add_5prime - 1):(curr_start - 1)]
					if to_add_3prime > ((curr_exon_list[exon_counter][1] + 1) - curr_end):
						add_seq_3prime_next_exon = to_add_3prime - ((curr_exon_list[exon_counter][1] + 1) - curr_end)
						add_seq_3prime = ref[
												(curr_exon_list[exon_counter + 1][0] - 1):
												(curr_exon_list[exon_counter + 1][0] - 1 + add_seq_3prime_next_exon)
												]
						to_add_3prime -= len(add_seq_3prime)
						# need to get next exon. remember the end of the gtf file segments are inclusive but not the case for the vcf
					add_seq_3prime = ref[(curr_end - 1):(curr_end - 1 + to_add_3prime)] + add_seq_3prime
					# make perfect version and actual version
				else:
					to_add_5prime = (curr_exon_list[exon_counter][1] + 1 - curr_end + (abs(curr_exon_list[exon_counter][3] - 3) % 3)) % 3 # length of codon before locus, then accounting for phase
					to_add_3prime = abs(((to_add_5prime + (curr_end - curr_start)) % 3) - 3) % 3
					if to_add_5prime > ((curr_exon_list[exon_counter][1] + 1) - curr_end):
						add_seq_5prime_prior_exon = to_add_5prime - ((curr_exon_list[exon_counter][1] + 1) - curr_end) # how much do i need from the prior exon
						add_seq_5prime = ref[
												(curr_exon_list[exon_counter + 1][0] - 1):
												(curr_exon_list[exon_counter + 1][0] - 1 + add_seq_5prime_prior_exon) # maybe add 1 more here?
												]
						to_add_5prime -= len(add_seq_5prime)
					add_seq_5prime = ref[(curr_end - 1):(curr_end - 1 + to_add_5prime)] + add_seq_5prime
					if to_add_3prime > (curr_start - curr_exon_list[exon_counter][0]):
						add_seq_3prime_next_exon = to_add_3prime - (curr_start - curr_exon_list[exon_counter][0]) # how much do i need from the next exon
						add_seq_3prime = ref[
												(curr_exon_list[exon_counter - 1][1] - add_seq_3prime_next_exon):
												(curr_exon_list[exon_counter - 1][1])
												]
						to_add_3prime -= len(add_seq_3prime)
						if to_add_3prime > 0:
							pdb.set_trace()
					add_seq_3prime += ref[(curr_start - to_add_3prime - 1):(curr_start - 1)]
					add_seq_5prime = common.revcomp(add_seq_5prime)
					add_seq_3prime = common.revcomp(add_seq_3prime)
					curr_allele = common.revcomp(curr_allele)
					list_of_possible_perfect_alleles = [common.revcomp(x) for x in list_of_possible_perfect_alleles]
# 					curr_perfect_allele = common.revcomp(curr_perfect_allele)
					
				curr_allele_in_frame = add_seq_5prime + curr_allele + add_seq_3prime
# 				curr_perfect_allele_in_frame = add_seq_5prime + curr_perfect_allele + add_seq_3prime
				list_of_possible_perfect_alleles = [add_seq_5prime + x + add_seq_3prime for x in list_of_possible_perfect_alleles]
				n_codons = int(len(curr_allele_in_frame) / 3)
				# should i handle stop gains?
				annotated_list_of_possible_perfect_alleles = []
				for poss_allele in list_of_possible_perfect_alleles:
					n_missense = 0
					n_silent = 0
					has_stop_gain = False
					has_stop_loss = False
					for i in range(n_codons):
						curr_actual_codon = curr_allele_in_frame[i * 3: (i + 1) * 3]
						curr_perfect_codon = poss_allele[i * 3: (i + 1) * 3]
						if curr_actual_codon != curr_perfect_codon:
							if common.codontab[curr_actual_codon] != common.codontab[curr_perfect_codon]:
								n_missense += 1
							else:
								n_silent += 1
						if common.codontab[curr_actual_codon] == "*":
							has_stop_gain = True
						if common.codontab[curr_perfect_codon] == "*":
							has_stop_loss = True
					if not has_stop_gain and not has_stop_loss:
						annotated_list_of_possible_perfect_alleles.append([poss_allele, n_missense, n_silent])
				if not annotated_list_of_possible_perfect_alleles:
					output_list.append(l_split + fields_to_add)
					continue
				min_aa_distance = min([x[1] for x in annotated_list_of_possible_perfect_alleles])
# 				sorted_by_nuc_distance = sorted([x + [distance(x[0], curr_allele_in_frame)] for x in annotated_list_of_possible_perfect_alleles], key = lambda x: x[3])
# 				if sorted_by_nuc_distance[0][1] != min_aa_distance and sorted_by_nuc_distance[0][3] != sorted_by_nuc_distance[1][3]:
# 					pdb.set_trace()
				best_possible_perfect_alleles = list(set([tuple(x) for x in annotated_list_of_possible_perfect_alleles if x[1] == min_aa_distance]))
# 				if curr_chr == "chr3" and curr_start == 138945766 and float(l_split[8]) > 0.9: # complex special case with multiple possible ancestral alleles with equal number of missense but unequal silent
# 					pdb.set_trace()
				curr_best_allele = sorted([(x, distance(x[0], curr_allele_in_frame)) for x in best_possible_perfect_alleles], key = lambda x: (x[1], x[0][2]))[0][0]
				n_missense_space = 0
				n_silent_space = 0
				for pos in range(len(add_seq_5prime), len(curr_best_allele[0]) - len(add_seq_3prime)):
					curr_codon = math.floor(pos / 3)
					curr_perfect_codon = curr_best_allele[0][curr_codon * 3: (curr_codon + 1) * 3]
					curr_nuc = curr_best_allele[0][pos]
					for mut in [x for x in ['A', 'C', 'G', 'T'] if x != curr_nuc]:
						curr_mut_codon = curr_perfect_codon[: pos % 3] + mut + curr_perfect_codon[(pos % 3) + 1:]
						if common.codontab[curr_mut_codon] == '*' or common.codontab[curr_perfect_codon] == "*":
							continue
						if common.codontab[curr_mut_codon] != common.codontab[curr_perfect_codon]: # maybe also count by the fraction of times a codon is missense vs silent? then take average for ratio?
							n_missense_space += 1
						else:
							n_silent_space += 1 # we can think about weighting these counts by ts/tv ratios eventually?
				# note: there are a very small number of alleles (4ish?) for which multiple different possible ancestral alleles require an equal number of missense and silent mutations 
				#      and have equivalent Levenshtein distance, but different possible missense and silent. I'm purposely leaving the choice of ancestral sequence as arbitrary between 
				#      the two as there does not seem to be a good logical explanation for one way or another
# 				if curr_chr == 'chr2' and curr_start == 184938704: # classic locus where observed major allele where AA parsimony is preferable to nucleotide
# 					pdb.set_trace()
				output_list.append(l_split + [str(n_codons), str(curr_best_allele[1]), str(curr_best_allele[2]), str(n_missense_space), str(n_silent_space)])
# 				n_stop_gains += has_stop_gain
				# save it here; remember to figure out what to do w the ones that overlap multiple transcripts
		print("# stop codons = %i" % (n_stop_gains))
		print("# multiple transcripts = %i" % (n_multiple))		
		open_output_file.write('\n'.join(['\t'.join(x) for x in output_list]) + '\n')


