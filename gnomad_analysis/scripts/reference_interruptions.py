import pysam
import pdb
import gzip
import math

vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_1.sorted.filtered.vcf.gz'
curr_chr = 'chr1'
curr_ref = '/net/eichler/vol26/home/mdnoyes/nobackups/reference/ucsc.hg38.no_alts.fasta'

def get_ref(curr_ref):
	with open(curr_ref, 'r') as open_file:
		ref_lines = []
		save_lines = False
		for l in open_file.readlines():
			if not save_lines and l.startswith('>%s' % curr_chr):
				save_lines = True
			elif save_lines:
				if l.startswith('>'):
					return(''.join([x.rstrip('\n') for x in ref_lines]))
				ref_lines.append(l)
	return None

curr_ref = get_ref(curr_ref)
interruption_table = []
with pysam.VariantFile(vcf_filename) as vcf_in:
	for var in vcf_in.fetch():
		matches = [curr_ref[var.pos - 1 + (var.info['PERIOD'] * x) : 
							var.pos - 1 + (var.info['PERIOD'] * (x + 1))].upper() == var.info['RU'].upper()
					for x in range(math.floor(var.info['REF']))]
		if sum(matches) != var.info['REF']:
			pdb.set_trace()
		interruption_table.append([curr_chr, str(var.pos), str(var.info['REF']), var.info['RU'].upper(), str(sum(matches))])

with open('/net/harris/vol1/project/simons_simplex/interruptions/tmp/%s.interruptions.txt' % curr_chr, 'w') as open_outfile:
	open_outfile.write('#chr\tpos\tref_ru\tru\tmatches\n')
	open_outfile.write('\n'.join(['\t'.join(x) for x in interruption_table]))
	open_outfile.write('\n')

