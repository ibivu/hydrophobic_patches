#!/bin/bash

#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -e hydrophobic_patches.err
#SBATCH -o hydrophobic_patches.out

module load 2019
#module load biopython
module load Miniconda3
source activate biopython

echo "making predictions on the human genome (LHPSA netsurfp2 model)"
cd ~/hydrophobic_patches/netsurfp_data
#python3 preprocess_input.py uniprot_9606.tab all_csv_files uniprot_9606_formatted.tab
python3 lhpsa_predictions.py NSP2_complete.tab ../model/lhpsa_nsp2.model lhpsa_netsurfp_model_predictions.tsv

