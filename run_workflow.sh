#!/bin/bash

#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -e hydrophobic_patches.err
#SBATCH -o hydrophobic_patches.out

#module load 2019
#module load biopython
#module load Miniconda3
#source activate biopython

echo "preparing files"
cd <...working_directory...>
# Download PDB files
python3 download_pdb.py
# Select chains
python3 create_chain_pdb.py
# Create DSSP files
python3 create_dssp.py

echo "processing data files"
cd ../processing
#echo "DSSP"
python3 process_dssp_data.py
#echo "fasta"
python3 process_fasta_data.py
#echo "netsurfp2"
python3 process_netsurfp2.py
#echo "global"
python3 process_global_sequence_data.py
#echo "monomer"
python3 process_monomer_data.py
#echo "TMHMM"
python3 process_TMHMM.py

#echo "final preprocessing for predictions"
python3 process_all_data_for_prediction.py

echo "building models"
cd ../modeling
python3 model.py

echo "making predictions on the human genome (LHPSA netsurfp2 model)"
cd ../../netsurfp_data
python3 preprocess_input.py uniprot_9606.tab all_csv_files uniprot_9606_formatted.tab
python3 lhpsa_predictions.py uniprot_9606_formatted.tab ../model/lhpsa_nsp2.model lhpsa_netsurfp_model_predictions.tsv
