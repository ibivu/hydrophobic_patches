#!/usr/bin/python3

"""
Make hydrophobic patch predictions using the NetsurfP2.0-based model
@Author: Juami van Gils
@Date: October 2020
"""

import pickle
import pandas as pd
import os
from sys import argv


### Files
# Output of preprocess_input.py
#f_in = '/mnt/c/Users/jhmva/surfdrive/PhD/hydrophobic_patches/uniprot_9606_formatted.tab'
f_in = argv[1]
# LHPSA model as a pickle
#f_model = '/mnt/c/Users/jhmva/surfdrive/PhD/Internship_Jan/stage1/model/lhpsa_nsp2.model'
f_model = argv[2]
# Output file
#f_out = '/mnt/c/Users/jhmva/surfdrive/PhD/hydrophobic_patches/model_predictions.tsv'
f_out = argv[3]

### Load data
input = pd.read_csv(f_in, sep='\t')
model = pickle.load(open(f_model, 'rb'))

### Run predictions
input_cols = input.loc[:, ['rhsa_netsurfp2','thsa_netsurfp2','tasa_netsurfp2']]
id = input.loc[:, 'id']

predictions = model.predict(input_cols)

### Write output file
pd.DataFrame({'id':id, 'prediction':predictions}).to_csv(f_out, sep='\t', index=False)
