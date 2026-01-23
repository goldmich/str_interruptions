#!/usr/bin/env python

import pdb
import random
import common
from Levenshtein import distance
from os.path import basename

vep_file = snakemake.input[0]

# check that the new pipeline succesfully got the loeuf scores for a few coding loci
# make a dict site: allele: max effect
# need to read all lines for each site into temporary dictionary
# need to save: is canonical transcript? (probably in reality need to use transcript i used for loeuf but whatever), allele, effect of allele
# for each (interrupted) allele: is there a coding variant relative to the reference? if so this is the most important; is it synonymous or nonsynonymous
# for the reference allele: is the perfect allele a coding variant relative to it?

vep_dict = {}

filter_strings_from_effect = ['intron', 'upstream', 'downstream', 'non_coding', 'UTR', 'splice', 'NMD']

with open(vep_file, 'r') as open_vep:
	curr_pos = -1
	curr_ref = ''
	curr_lines = {}
	other_refs = []
	for l in open_vep.readlines():
		if l.startswith("#"):
			continue
		l_split = l.rstrip('\n').split('\t')
		pos = int(l_split[0].split('_')[1])
		if curr_pos == -1:
			curr_pos = pos
			curr_ref = l_split[0].split('_')[2].split('/')[0]
		if pos != curr_pos:
			vep_dict[curr_pos] = {}
			for allele in curr_lines:
				transcript_variants = [x[1] for x in curr_lines[allele] if x[0] == 'Transcript' and x[3] == 'YES']
				if transcript_variants:
					transcript_variants = [x for x in transcript_variants if all([y not in x for y in filter_strings_from_effect])]
					if all([x == 'synonymous_variant' for x in transcript_variants]):
						vep_dict[curr_pos][allele] = 'synonymous_variant'
					elif 'start_lost' in transcript_variants:
						vep_dict[curr_pos][allele] = 'start_lost'
					elif 'stop_gained' in transcript_variants:
						vep_dict[curr_pos][allele] = 'stop_gained'
					elif 'frameshift_variant' in transcript_variants:
						vep_dict[curr_pos][allele] = 'frameshift_variant'
					elif 'stop_lost' in transcript_variants:
						vep_dict[curr_pos][allele] = 'stop_lost'
					elif 'missense_variant' in transcript_variants: # yes this rule is more flexible
						vep_dict[curr_pos][allele] = 'missense_variant'
					elif all([x == 'inframe_deletion' for x in transcript_variants]):
						vep_dict[curr_pos][allele] = 'inframe_deletion'
					elif all([x == 'inframe_insertion' for x in transcript_variants]):
						vep_dict[curr_pos][allele] = 'inframe_insertion'
					elif all([x == 'protein_altering_variant' or x == 'coding_sequence_variant' for x in transcript_variants]):
						vep_dict[curr_pos][allele] = 'other'
					else:
						pdb.set_trace()
						vep_dict[curr_pos][allele] = 'other' # i'll try to catch more of these later
					continue # next allele; i care mostly about coding variants
				else:
					vep_dict[curr_pos][allele] = 'regulatory' # deal w this later
			# it's easier to process the ref allele later - need repeat unit
			vep_dict[curr_pos][curr_ref] = 'ref_allele'
			for x in other_refs:
				vep_dict[curr_pos][x] = 'ref_allele'
				pdb.set_trace()
			curr_lines = {}
			curr_pos = pos
			curr_ref = l_split[0].split('_')[2].split('/')[0]
			other_refs = []
		ref = l_split[0].split('_')[2].split('/')[0]
		if ref != curr_ref: # at least one of these loci is super complex and has multiple possible ref alleles for the same start position... i don't want to match on
			curr_ref = ref
		if l_split[2] not in curr_lines:
			curr_lines[l_split[2]] = []
		curr_info = {x.split('=')[0]: x.split('=')[1] for x in l_split[13].split(';')}
		if l_split[5] == 'Transcript':
			curr_effects = l_split[6].split(',') # sometimes there are multiple effects separated by comma
			if 'CANONICAL' in curr_info:
				curr_lines[l_split[2]] += [[l_split[5], x, curr_info['IMPACT'], curr_info['CANONICAL']] for x in curr_effects]
			else:
				curr_lines[l_split[2]] += [[l_split[5], x, curr_info['IMPACT'], 'N/A'] for x in curr_effects]
		else:
			# skip for now
			continue
	# one more locus
	vep_dict[curr_pos] = {}
	for allele in curr_lines:
		transcript_variants = [x[1] for x in curr_lines[allele] if x[0] == 'Transcript' and x[3] == 'YES']
		if transcript_variants:
			transcript_variants = [x for x in transcript_variants if all([y not in x for y in filter_strings_from_effect])]
			if all([x == 'synonymous_variant' for x in transcript_variants]):
				vep_dict[curr_pos][allele] = 'synonymous_variant'
			elif 'start_lost' in transcript_variants:
				vep_dict[curr_pos][allele] = 'start_lost'
			elif 'stop_gained' in transcript_variants:
				vep_dict[curr_pos][allele] = 'stop_gained'
			elif 'frameshift_variant' in transcript_variants:
				vep_dict[curr_pos][allele] = 'frameshift_variant'
			elif 'stop_lost' in transcript_variants:
				vep_dict[curr_pos][allele] = 'stop_lost'
			elif 'missense_variant' in transcript_variants: # yes this rule is more flexible
				vep_dict[curr_pos][allele] = 'missense_variant'
			elif all([x == 'inframe_deletion' for x in transcript_variants]):
				vep_dict[curr_pos][allele] = 'inframe_deletion'
			elif all([x == 'inframe_insertion' for x in transcript_variants]):
				vep_dict[curr_pos][allele] = 'inframe_insertion'
			elif all([x == 'protein_altering_variant' or x == 'coding_sequence_variant' for x in transcript_variants]):
				vep_dict[curr_pos][allele] = 'other'
			else:
				pdb.set_trace()
				vep_dict[curr_pos][allele] = 'other' # i'll try to catch more of these later
			continue # next allele; i care mostly about coding variants
		else:
			vep_dict[curr_pos][allele] = 'regulatory' # deal w this later
	vep_dict[curr_pos][curr_ref] = 'ref_allele'
	for x in other_refs:
		vep_dict[curr_pos][x] = 'ref_allele'


outfile = snakemake.output[0]
curr_chrom = basename(outfile).split('.')[1]
with open(outfile, 'w') as open_outfile:
	for locus in sorted(vep_dict.keys()):
		open_outfile.write('\n'.join(['\t'.join([curr_chrom, str(locus), allele, vep_dict[locus][allele]]) for allele in vep_dict[locus]]))
		open_outfile.write('\n')



