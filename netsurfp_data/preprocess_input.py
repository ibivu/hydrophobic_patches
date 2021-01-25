#!/usr/bin/python3

"""
Convert NetsurfP2.0 output to feature format used by the model
@Author: Juami van Gils
@Date: October 2020
"""

import pandas as pd
import os
from sys import argv

### Files and folders
# File with Uniprot info from human proteins (tab-separated)
f_in = '/home/juami/hydrophobic_patches_jan/netsurfp_dea/uniprot_9606.tab'
# Folder with NetsurfP2.0 output files
#d_in = '/mnt/c/Users/jhmva/surfdrive/PhD/hydrophobic_patches/NetsurfP2.0-predictions'
d_in = '/home/juami/hydrophobic_patches_jan/netsurfp_dea/all_csv_files'
netsurfp_files = [os.path.join(d_in, x) for x in os.listdir(d_in)]
# Output file
f_out = '/home/juami/hydrophobic_patches_jan/netsurfp_dea/uniprot_9606_formatted.tab'

with open(f_out, 'w') as f:
	f.write('\t'.join(['id', 'Gene name', 'PDB', 'length', 'hydr_count', 'polar_count', 'molecular_weight', 'helix', 'turn', 'sheet', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'fasta_sequence', 'thsa_netsurfp2', 'tasa_netsurfp2', 'q3_H', 'q3_E', 'q3_C', 'rhsa_netsurfp2', 'on_surface', 'disorder']) + '\n')

### Residues
alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
hydrophobic = ['A', 'C', 'F', 'L', 'I', 'M', 'V', 'W', 'Y']
polar = ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'T']
molecular_weight = {'A': 89.09, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2, 'G': 75.07, 'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2, 'M': 149.2, 'N': 132.1, 'P': 115.1, 'Q': 146.1, 'R': 174.2, 'S': 105.09, 'T': 119.1, 'V': 117.1, 'W': 204.2, 'Y': 181.2}
disorder_threshold_residue = 0.4 #If disorder prediction for a residue is larger than this, it is considered a disordered residue
#disorder_threshold_protein = 0.5 #If a larger fraction of the residues is predicted to be disordered, consider the protein disordered

### Create dataframes
uniprot_df = pd.read_csv(f_in, sep='\t')
netsurfp_df = pd.read_csv(netsurfp_files[0], sep=',')

for i in range(1, len(netsurfp_files)):
	tmp = pd.read_csv(netsurfp_files[i], sep=',')
	netsurfp_df = pd.concat([netsurfp_df, tmp])

rownames = [x.split('_')[0] for x in netsurfp_df.loc[:,'id']]
#rownames_dict = {}
#for i in netsurfp_df.index:
#	rownames_dict[i] = rownames[i]

#netsurfp_df = netsurfp_df.rename(index=rownames_dict)

#print(netsurfp_df.shape)
#print(sum([x == 'Q8N7X0' for x in rownames]))
#print(sum([x == 'Q8N7X0' for x in netsurfp_df.index]))
#print(sum(['Q8N7X0' in x for x in netsurfp_df.loc[:, 'id']]))
#print(sum([x.split('_')[1] == 'Q8N7X0' for x in netsurfp_df.loc[:, 'id']]))
#print(netsurfp_df.iloc[:5,:5])
#print(netsurfp_df.loc['Q8N7X0', ['id', 'seq', 'n']])
#print(netsurfp_df.loc['Q8N7X0', ['id', 'seq', 'n']].shape)

#print(netsurfp_df.index[:5])

#print(netsurfp_df.iloc[:5,:5])
#print(netsurfp_df.shape)

### For every Uniprot entry, get relevant columns, merge with NetsurfP info if available
for i in range(len(uniprot_df.index)):
	#print(uniprot_df.iloc[i,:5])
	# Uniprot info
	id          = uniprot_df.loc[i,'Entry']
	#print(i, id)
	gene_name   = uniprot_df.loc[i, 'Gene names  (primary )']
	pdb         = uniprot_df.loc[i,'Cross-reference (PDBsum)']
	length      = uniprot_df.loc[i,'Length']
	sequence    = uniprot_df.loc[i,'Sequence']
	hydr_count  = sum([x in hydrophobic for x in sequence])
	#print(hydr_count)
	polar_count = sum([x in polar for x in sequence])
	aa_counts = {}
	mw = 0
	for aa in alphabet:
		aa_counts[aa] = sum([x == aa for x in sequence])
		mw += molecular_weight[aa] * sum([x == aa for x in sequence])
	
	asa = 0
	hsa = 0
	rsa = 0
	on_surface = 0
	q3 = {'H': 0, 'E': 0, 'C': 0}
	disorder_frac = 0
	# NetsurfP2 info
	if id in rownames:
		#print(id)
		tmp = netsurfp_df.loc[[x == id for x in rownames], :]
		#print(netsurfp_df.shape, tmp.shape)
		#print(tmp.index)
		
		#print(length)
		q3_col = tmp.loc[:, 'q3']
		asa_col = tmp.loc[:, 'asa']
		seq_col = tmp.loc[:, 'seq']
		dio_col = tmp.loc[:, 'disorder']
		
		for j in tmp.index:
			#print(id, j)
			#print(q3_col[j])
			q3[q3_col[j]] += 1
			asa += float(asa_col[j])
			if seq_col[j] in hydrophobic:
				hsa += float(asa_col[j])
				on_surface += 1
			if dio_col[j] > disorder_threshold_residue:
				disorder_frac += 1
				
		rsa = float(hsa)/asa
	disorder_frac = disorder_frac/float(length)
	helix = q3['H'] / float(length)
	turn = q3['C'] / float(length)
	sheet = q3['E'] / float(length)
	
	
	with open(f_out, 'a') as f:
		f.write('\t'.join([str(value) for value in [id, gene_name, pdb, length, hydr_count, polar_count, mw, helix, turn, sheet, '\t'.join([str(g) for g in aa_counts.values()]), sequence, hsa, asa, q3['H'], q3['E'], q3['C'], rsa, on_surface, disorder_frac]]) + '\n')
		#['id', 'length', 'hydr_count', 'polar_count', 'molecular_weight', 'helix', 'turn', 'sheet', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'fasta_sequence', 'thsa_netsurfp2', 'tasa_netsurfp2', 'q3_H', 'q3_E', 'q3_C', 'rhsa_netsurfp2', 'on_surface']) + '\n')
	
	
	#break
	
