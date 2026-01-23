#!/usr/bin/env python

import Levenshtein
import pdb

revcomp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def revcomp(seq):
	return ''.join([revcomp_dict[x] for x in list(seq[::-1])])

def find_longest_homologous_seq(allele, ru):
	matches = []
	counter = 0
	curr_match = None
	while counter <= len(allele) - len(ru):
		if allele[counter:].startswith(ru):
			if not curr_match:
				curr_match = [0, counter]
			curr_match[0] += 1
			counter += len(ru)
		else:
			if curr_match:
				matches.append(curr_match)
				curr_match = None
			counter += 1
	if curr_match:
		matches.append(curr_match)
# 	if not matches:
# 		pdb.set_trace()
	return sorted(matches, reverse = True, key = lambda x: x[0])

def hamming_distance(seq1, seq2):
	return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def levenshtein_distance_semilocal(seq1, seq2, motif, offset_start, offset_end):
	curr_terminal_indels = [x for x in Levenshtein.editops(seq1, seq2) if x[0] in ['insert', 'delete'] and (x[2] == len(seq2) or x[1] == len(seq1) or x[2] == 0)]
	curr_seq_distance_tuple = (seq2, Levenshtein.distance(seq1, seq2))
	if not curr_terminal_indels:
		return curr_seq_distance_tuple
	next_seq = ''
# 	pdb.set_trace()
	while curr_terminal_indels:
		if curr_terminal_indels[0][0] == 'insert':
			if curr_terminal_indels[0][2] == 0:
				next_seq = seq2[1:]
			else:
				next_seq = seq2[:-1]
		else:
			motif_string = motif * 2
			if curr_terminal_indels[0][2] == 0:
				to_add = motif_string[offset_start + len(motif) - 1]
				next_seq = to_add + seq2
				offset_start = (offset_start + len(motif) - 1) % len(motif)
			else:
				to_add = motif_string[offset_end + 1]
				next_seq = seq2 + to_add
				offset_end = (offset_end + 1) % len(motif)
		next_seq_distance_tuple = (next_seq, Levenshtein.distance(seq1, next_seq))
		if next_seq_distance_tuple[1] >= curr_seq_distance_tuple[1]:
			break
		curr_seq_distance_tuple = next_seq_distance_tuple
		curr_terminal_indels = [x for x in Levenshtein.editops(seq1, next_seq) 
									if x[0] in ['insert', 'delete'] and (x[2] == len(next_seq) or x[1] == len(seq1) or x[2] == 0)]
		seq2 = next_seq
	return curr_seq_distance_tuple

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

if __name__ == '__main__':
	print(levenshtein_distance_semilocal('ACCAC', 'ACACA', 'AC', 0, 0))
	print(levenshtein_distance_semilocal('ACCAC', 'CACAC', 'AC', 1, 1))
	print(levenshtein_distance_semilocal('ACACCACAACAC', 'ACACACACACAC', 'AC', 0, 1))
	print(levenshtein_distance_semilocal('ACGACACGACGGACGAG', 'ACGACGACGACGACGAC', 'ACG', 0, 1))
	print(levenshtein_distance_semilocal('ACGACACGACGGACGAG', 'GACGACGACGACGACGA', 'ACG', 2, 0))
	print(levenshtein_distance_semilocal('ACGACACGACGGACGAG', 'CGACGACGACGACGACG', 'ACG', 1, 2))
	print('bingo bango')
	
	
	
