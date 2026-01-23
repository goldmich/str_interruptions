#!/usr/bin/env python

import pysam
import pdb
import re
import random
import os
import common
import Levenshtein

# vcf_in = '/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/str_1kgp_h3a/ensemble_chr15_filtered.vcf.gz'
vcf_in = snakemake.input[0]
curr_chr = os.path.basename(vcf_in).split('_')[1]

output_list = []
output_allele_seq_list = []
with pysam.VariantFile(vcf_in) as open_vcf_in:
# 	counter = 0
	sample_list = open_vcf_in.header.samples
	for var in open_vcf_in.fetch():
		if var.info['METHODS'].split('|')[2] != "1": # need a hipSTR call
			continue
		if len(var.info['RU']) == 1:
			continue
		if var.filter:
			continue
		# These vcfs have an updated sample list - AC field is recalculated but REFAC is not, but i'm being safe and recalculating all ANs
		all_gts = [var.samples[x]['GT'] for x in sample_list if var.samples[x]['FILTER'] == 'PASS']
		all_haplotypes = [x for y in all_gts for x in y]
		curr_an_list = [sum([x == y for x in all_haplotypes]) for y in range(len(var.alleles))]
		for n, allele in enumerate(var.alleles):
			if curr_an_list[n] < 2:
				continue
			curr_matches = common.find_longest_homologous_seq(allele, var.info['RU'])
			curr_revcomp_matches = common.find_longest_homologous_seq(allele, common.revcomp(var.info['RU']))
			best_matches = None
			is_revcomp = False
			if curr_matches:
				if curr_revcomp_matches:
					if curr_matches[0][0] > curr_revcomp_matches[0][0]:
						best_matches = curr_matches
					else:
						best_matches = curr_revcomp_matches
						is_revcomp = True
				else:
					best_matches = curr_matches
			else:
				if curr_revcomp_matches:
					best_matches = curr_revcomp_matches
					is_revcomp = True
				else:
					best_matches = [[0, -1]]
			ru_to_repeat = var.info['RU']
			if is_revcomp:
				ru_to_repeat = common.revcomp(var.info['RU'])
			min_distance = len(allele)
			perfect_version = ''
			for x in range(len(var.info['RU'])):
				curr_perfect = (ru_to_repeat * (len(allele) + 1))[x : (len(allele) + x)]
				curr_offset_start = x
				curr_offset_end = (len(allele) + x - 1) % len(ru_to_repeat)
				curr_distance_tuple = common.levenshtein_distance_semilocal(allele, curr_perfect, ru_to_repeat, curr_offset_start, curr_offset_end)
				if curr_distance_tuple[1] < min_distance:
					min_distance = curr_distance_tuple[1]
					perfect_version = curr_distance_tuple[0]
# 			min_distance = min([Levenshtein.distance(allele, (ru_to_repeat * (len(allele) + 1))[x : (len(allele) + x)]) for x in range(len(var.info['RU']))])
# 			curr_homology = float(sum([x[0] for x in best_matches]) * len(var.info['RU'])) / float(len(allele))
# 			if curr_homology == 1 and (best_matches[0][0] * len(var.info['RU'])) < 5:
# 				pdb.set_trace()
			curr_af = curr_an_list[n] / sum(curr_an_list)
			output_list.append([var.pos, var.pos + len(var.ref), var.info['RU'], allele, len(allele), best_matches[0][0], min_distance, curr_af])
			output_allele_seq_list.append([var.pos, var.pos + len(var.ref), allele, perfect_version, curr_af])
# 		counter += 1
# 		if counter == 1000:
# 			break

af_homology_stats_output_file = snakemake.output['af_homology_stats']
allele_perfect_seq_output_file = snakemake.output['allele_perfect_seq']

with open(af_homology_stats_output_file, 'w') as open_af_outfile:
	open_af_outfile.write('#chr\tstart\tend\tru\tallele_len\tallele\tmax_homo_ru\tdistance_from_perfect\taf\n')
	for i in output_list:
		open_af_outfile.write('\t'.join([curr_chr] + [str(x) for x in i]) + '\n')

with open(allele_perfect_seq_output_file, 'w') as open_seq_outfile:
	open_seq_outfile.write('#chr\tstart\tend\tallele\tperfect\n')
	for i in output_allele_seq_list:
		open_seq_outfile.write('\t'.join([curr_chr] + [str(x) for x in i]) + '\n')
