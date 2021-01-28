# SRC

This directory has three main subdirectories.

* /enrichment_and_expression
  * normalize_scores.py --> plot_histogram.py; normalize the predicted values for the measures for GSEA, and plot the new distributions
  * Human_Proteome_Mapping.ipynb; code for the functional analyses in this work
* /modeling
  * src for the creation of the ML models
* /molpatch
  * src code for MolPatch
* /preparation
  * src for downloading and generating all needed files for global feature processing (PDB, DSSP)
* /processing
  * code to generate global feature data

## How to use

First fill in the root directory in the config.yml file

#### preparation

* Run download_pdb.py to get all the full pdbs in /data/pdb/proteins
* Run create_chain.py to create pdb files of a selected pdb in /data/pdb/chains
* Run create_dssp.py to create dssp and dssp_csv files

#### processing

* All processed data files are in /data/processed_data already
* to run process_all_data_for_prediction.py you should first run all other preprocessing files

#### modeling

* TFM models are created in TFM.R since Cubist regression is an R methods
* all other models are created in model.py

#### molpatch

* calculate_patches_pisite_atom.py and calculate_patches_pisite_res.py are examples of how the ProteinPatch class can be used to find the largest patch
* PiSITE parser can be used to check what part of the patch are protein interacting residues
