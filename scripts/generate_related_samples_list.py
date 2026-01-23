#!/usr/bin/env python

import pdb

unrelated_1kgp_samples_file = snakemake.input[0]
all_1kgp_ped_file = snakemake.input[1]

unrelated_samples = []
with open(unrelated_1kgp_samples_file, 'r') as open_unrelated_file:
	for l in open_unrelated_file.readlines():
		if l.startswith("#"):
			continue
		unrelated_samples.append(l.split('\t')[9])

all_ped_samples = []
with open(all_1kgp_ped_file, 'r') as open_ped_file:
	open_ped_file.readline() # header
	for l in open_ped_file.readlines():
		all_ped_samples.append(l.split('\t')[1])

related_samples = set(all_ped_samples) - set(unrelated_samples)

output_file = snakemake.output[0]
with open(output_file, 'w') as open_output_file:
	open_output_file.write('\n'.join(related_samples))
