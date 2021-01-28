#!/usr/bin/python3

"""
Normalise scores for GSEA
@Author: Juami van Gils
@Date: October 2020
"""

import pandas as pd
import os
from sys import argv
import numpy as np

os.chdir(argv[0])

### Files and folders
# Files with non-normalized scores
f_in1 = 'LHPSA_GSEA.rnk'
f_in2 = 'THSA_GSEA.rnk'
f_in3 = 'RHSA_GSEA.rnk'

# Output files
f_out1 = 'LHPSA_GSEA_cs.rnk'
f_out2 = 'THSA_GSEA_cs.rnk'
f_out3 = 'RHSA_GSEA_cs.rnk'


### Calculate normalized scores
# Largest patch
scores = []
with open(f_in1, 'r') as f:
	for line in f:
		scores.append(float(line.strip().split('\t')[1]))
median = np.median(scores)
sd = np.std(scores)
print('LHPSA', median, sd)

with open(f_in1, 'r') as infile:
	with open(f_out1, 'w') as outfile:
		for line in infile:
			splitline = line.strip().split('\t')
			name = splitline[0]
			score = float(splitline[1])
			newscore = (score - median - 200) / sd
			
			outfile.write('%s\t%f\n' % (name, newscore))

			
# THSA
scores = []
with open(f_in2, 'r') as f:
	for line in f:
		scores.append(float(line.strip().split('\t')[1]))
median = np.median(scores)
sd = np.std(scores)
print('THSA', median, sd)

with open(f_in2, 'r') as infile:
	with open(f_out2, 'w') as outfile:
		for line in infile:
			splitline = line.strip().split('\t')
			name = splitline[0]
			score = float(splitline[1])
			newscore = (score - median - 1900) / sd
			
			outfile.write('%s\t%f\n' % (name, newscore))
			
# RHSA
scores = []
with open(f_in3, 'r') as f:
	for line in f:
		scores.append(float(line.strip().split('\t')[1]))
median = np.median(scores)
sd = np.std(scores)
print('RHSA', median, sd)

with open(f_in3, 'r') as infile:
	with open(f_out3, 'w') as outfile:
		for line in infile:
			splitline = line.strip().split('\t')
			name = splitline[0]
			score = float(splitline[1])
			newscore = (score - median - 0.08) / sd
			
			outfile.write('%s\t%f\n' % (name, newscore))
