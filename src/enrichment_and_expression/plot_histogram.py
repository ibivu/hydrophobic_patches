#!/usr/bin/python3

"""
Normalise scores for GSEA
@Author: Juami van Gils
@Date: October 2020
"""

import os
from sys import argv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import pandas as pd
#import matplotlib.pyplot as plt

os.chdir(argv[0])

# Input files
f_in1 = 'LHPSA_GSEA_cs.rnk'
f_in2 = 'THSA_GSEA_cs.rnk'
f_in3 = 'RHSA_GSEA_cs.rnk'

# Output files
f_out1 = 'LHPSA_histogram.png'
f_out2 = 'THSA_histogram.png'
f_out3 = 'RHSA_histogram.png'

### Create plots
plt.rc('font', size=18)
# LHPSA
df = pd.read_csv(f_in1, sep='\t')
ax = df.plot.density().get_figure()
plt.vlines(0, 0, 1.2, linestyles='dashed')
ax.savefig(f_out1)

# THSA
df = pd.read_csv(f_in2, sep='\t')
ax = df.plot.density().get_figure()
plt.vlines(0, 0, 1.2, linestyles='dashed')
ax.savefig(f_out2)

# RHSA
df = pd.read_csv(f_in3, sep='\t')
ax = df.plot.density().get_figure()
plt.vlines(0, 0, 0.8, linestyles='dashed')
ax.savefig(f_out3)