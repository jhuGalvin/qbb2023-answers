#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys
from fasta import readFASTA

input_sequences = readFASTA(open(sys.argv[1]))
scoring_matrix_file = sys.argv[2]
gap_penalty = float(sys.argv[3])
align_file = sys.argv[4]

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]

# nt gap penalty = -300, aa gap penalty = -10

# ./w3Hw.py CTCF_38_M27_AA.faa BLOSUM62.txt -10

scoring_matrix = pd.DataFrame(pd.read_csv(scoring_matrix_file, sep = '\s+'))
print(scoring_matrix)

# test code on aa file, not nt; it will be faster

F_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
traceback_matrix = np.zeros((len(sequence1) + 1, len(sequence2) + 1), dtype = str)

for i in range(len(sequence1) + 1):
	F_matrix[i,0] = i * gap_penalty

for j in range(len(sequence2) + 1):
	F_matrix[0,j] = j * gap_penalty

for i in range(1, F_matrix.shape[0]):
	traceback_matrix[i, 0] = 'v'

	for j in range(1, F_matrix.shape[1]):
		traceback_matrix[0, j] = 'h'

		match_score = float(scoring_matrix.loc[sequence1[i - 1], sequence2[j - 1]])
		d = F_matrix[i - 1, j - 1] + match_score
		h = F_matrix[i, j - 1] + gap_penalty
		v = F_matrix[i - 1, j] + gap_penalty

		F_matrix[i,j] = max(d, h, v)

		if F_matrix[i,j] == d:
			traceback_matrix[i, j] = 'd'
		elif F_matrix[i,j] == h:
			traceback_matrix[i, j] = 'h'
		else:
			traceback_matrix[i, j] = 'v'

print(F_matrix)
print(traceback_matrix)

seq1Al = ''
seq2Al = ''


while i != 0 or j != 0:
	if traceback_matrix[i, j] == 'd':
		seq1Al = sequence1[i - 1] + seq1Al
		seq2Al = sequence2[j - 1] + seq2Al
		i -= 1
		j -= 1

	elif traceback_matrix[i, j] == 'h':
		seq1Al = '-' + seq1Al
		seq2Al = sequence2[j - 1] + seq2Al
		j -= 1

	else:
		seq1Al = sequence1[i] + seq1Al
		seq2Al = '-' + seq2Al
		i -= 1

print(seq1Al)
print("Gaps in sequence 1:", seq1Al.count('-'))
print(seq2Al)
print("Gaps in sequence 2:", seq2Al.count('-'))

seq2Al.count('-')

file = open(align_file, 'w')
lines = ["Sequence 1:", seq1Al, "\nSequence 2:", seq2Al]
file.writelines(lines)
file.close()

# for protein alignment
# ./w3Hw.py CTCF_38_M27_AA.faa BLOSUM62.txt -300 alignmentAA.txt

# for nucleotide
# ./w3Hw.py CTCF_38_M27_DNA.fna HOXD70.txt -10 alignmentNT.txt

