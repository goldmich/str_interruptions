#!/usr/bin/env python

import pdb
import random
import pysam
import common
from Levenshtein import distance

bcf_in = snakemake.input[0]

# maybe consider making each allele a separate entry then calculating the "perfect" version for each one?
with pysam.VariantFile(bcf_in) as open_bcf_in:
	bcf_out = pysam.VariantFile(snakemake.output[0], 'w', header = open_bcf_in.header)
	sample_list = open_bcf_in.header.samples
	for var in open_bcf_in:
		# These vcfs have an updated sample list - AC field is recalculated but REFAC is not, but i'm being safe and recalculating all ANs
		# i already did locus-level filtering when i grabbed the locus characteristics
		all_gts = [var.samples[x]['GT'] for x in sample_list if var.samples[x]['FILTER'] == 'PASS']
		all_haplotypes = [x for y in all_gts for x in y]
		curr_an_list = [sum([x == y for x in all_haplotypes]) for y in range(len(var.alleles))]
		# filter singleton alleles
		# take reference allele: is it interrupted? if not then just add the nonsingletons to the bcf
		# if it's interrupted, is there a perfect allele that matches the var_to_ad? if not then add the var_to_add
		# annoying but there are reference alleles with < 2 allele counts... i think i need to keep them for now
		nonsingletons = [allele for n, allele in enumerate(var.alleles) if curr_an_list[n] > 1 or allele == var.ref]
		curr_ref_allele = var.ref
		list_of_possible_perfect_alleles = [(ru_to_repeat * (len(curr_ref_allele) + 1))[x : (len(curr_ref_allele) + x)] for x in range(len(var.info['RU'])) for ru_to_repeat in [var.info['RU'], common.revcomp(var.info['RU'])]]
		# here you need to change
		min_distance_combo = min([(distance(curr_ref_allele, x), x) for x in list_of_possible_perfect_alleles], key = lambda y: y[0])
		alleles_for_var = nonsingletons
		# if the ref allele is perfect or there is a perfect nonsingleton no need to add a hypothetical perfect allele
		if min_distance_combo[0] != 0 and min_distance_combo[1] not in nonsingletons:
			alleles_for_var.append(min_distance_combo[1])
		if len(alleles_for_var) == 1:
			continue # at the moment i don't think i need to worry about monoallelic pure loci
		var_to_add = bcf_out.new_record(contig = var.contig, start = var.start, stop = var.stop, alleles = tuple(alleles_for_var))
		var_to_add.info['RU'] = var.info['RU']
		bcf_out.write(var_to_add)
		
# 		for n, allele in enumerate(var.alleles):
# 			# no singletons
# 			if curr_an_list[n] < 2:
# 				continue
# 			list_of_possible_perfect_alleles = [(ru_to_repeat * (len(allele) + 1))[x : (len(allele) + x)] for x in range(len(var.info['RU'])) for ru_to_repeat in [var.info['RU'], common.revcomp(var.info['RU'])]]
# 			min_distance_combo = min([(distance(allele, x), x) for x in list_of_possible_perfect_alleles])
# 			if min_distance_combo[0] == 0:
# 				continue # i don't think i need to interpret the interruptions at an allele that is not interrupted
# 			var_to_add = bcf_out.new_record(contig = var.contig, start = var.start, stop = var.stop, alleles = (allele, min_distance_combo[1]))
# 			var_to_add.info['RU'] = var.info['RU']
# 			bcf_out.write(var_to_add)
			


