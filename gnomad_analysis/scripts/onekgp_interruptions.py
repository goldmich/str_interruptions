import pysam
import pdb
import gzip
import math
import os
import numpy as np

# I'm currently just using an SSC vcf to get the repeat coordinates and reference lengths... i could figure out something better eventually
# vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_1.sorted.filtered.vcf.gz'
# gnomad_vcf = '/net/harris/vol1/data/gnomad_v3/gnomad.1.vcf.bgz'
vcf_filename = snakemake.input[0]
onekgp_vcf= snakemake.input[1]

curr_chr = os.path.basename(onekgp_vcf).split('.')[1][3:]

interruption_dict = {}
# var_qual_list = []
with pysam.VariantFile(vcf_filename) as vcf_in:
	with pysam.VariantFile(onekgp_vcf) as onekgp_vcf_in:
		for var in vcf_in.fetch():
			if var.info['PERIOD'] == 1:
				continue
			n_var = np.zeros((5,4), dtype = int)
			curr_region = (var.pos, var.pos + int(var.info['REF'] * var.info['PERIOD']))
			for onekgp_var in onekgp_vcf_in.fetch(curr_chr, curr_region[0], curr_region[1]):
				if 'PASS' in onekgp_var.filter and onekgp_var.ref in ['A', 'C', 'G', 'T'] and all([x in ['A', 'C', 'G', 'T'] for x in onekgp_var.alts]):
					max_af = max(onekgp_var.info['AF'])
					pdb.set_trace()
					pos_within_repeat = (onekgp_var.pos - curr_region[0]) / (curr_region[1] - curr_region[0])
					pos_quantile = int(pos_within_repeat * 5)
					# allowing the final bin to be inclusive
					if pos_quantile == 5:
						pos_quantile = 4
					if max_af > 0.1:
						n_var[pos_quantile, 3] += 1
					elif max_af > 0.001:
						n_var[pos_quantile, 2] += 1
					elif max(onekgp_var.info['AC']) > 1: # rare but not singletons
						n_var[pos_quantile, 1] += 1
					else: # singletons
						n_var[pos_quantile, 0] += 1
						
# 					var_qual_list.append(onekgp_var.qual)
# 				elif 'PASS' in onekgp_var.filter:
# 					pdb.set_trace()
			interruption_dict[var.pos] = [str(int(var.info['REF'])), var.info['RU'].upper(), n_var]


output_file = snakemake.output[0]

with open(output_file, 'w') as open_outfile:
	open_outfile.write('#chr\tpos\tref_ru\tru\tpos_quantile\tallele_frequency_bin\tn\n')
	for pos in sorted(interruption_dict.keys()):
		for pos_quantile in range(5):
			for allele_freq_bin in range(4):
				open_outfile.write('\t'.join([curr_chr, str(pos), 
											  interruption_dict[pos][0], 
											  interruption_dict[pos][1], 
											  str(pos_quantile), 
											  str(allele_freq_bin), 
											  str(interruption_dict[pos][2][pos_quantile, allele_freq_bin])
											  ])
											  )
				open_outfile.write('\n')

# pdb.set_trace()

