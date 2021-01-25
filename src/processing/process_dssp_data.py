from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from os import listdir
from os.path import isfile, join
import pandas as pd
import yaml

config = yaml.safe_load(open("../config.yml"))

hydr = config['hydrophobic']

def correct_sequence(sequence):
    """Correct ambigious aa's or cystine bridges

    Parameters
    ----------
    sequence : str
        protien sequence

    Returns
    -------
    str
        corrected protein sequence
    """
    corrected_sequence = ''

    for aa_type in sequence:
        if aa_type.islower() or aa_type in ['B', 'Z']:
            print(aa_type)
        # cysteines forming a cysteine bond
        elif aa_type.islower():
            aa_type = 'C'
        # aspartic acid or asparagine
        elif aa_type == 'B':
            aa_type = 'D'
        # Z = glutamic acid or glutamate
        elif aa_type == 'Z':
            aa_type = 'E'

        corrected_sequence += aa_type

    return corrected_sequence

def process_dssp_data(DSSP_CSV_PATH, PROCESSED_DATA_PATH):
    """calculate the thsa, rhsa and tasa and write this to a csv

    Parameters
    ----------
    DSSP_CSV_PATH : str
        path to the dssp_csv files
    PROCESSED_DATA_PATH : str
        output file

    """
    dssp_csv_files = [f for f in listdir(DSSP_CSV_PATH) if isfile(join(DSSP_CSV_PATH, f))]
    results = {
            'id':[],
            'tasa':[],
            'thsa':[],
            'dssp_sequence':[]
            }

    # add the TASA, THSA and RHSA to the dictionary for every protein
    for dssp_csv_file in dssp_csv_files:
        df = pd.read_csv(DSSP_CSV_PATH+dssp_csv_file, index_col=0)
        results['id'].append(dssp_csv_file[:-4])
        corrected_sequence = correct_sequence(df['2'].sum())
        results['dssp_sequence'].append(corrected_sequence)
        results['tasa'].append(df['4'].sum())
        results['thsa'].append(df[df['2'].isin(hydr)]['4'].sum())

    df = pd.DataFrame(results)
    df['rhsa'] = df['thsa']/df['tasa']
    df.to_csv(PROCESSED_DATA_PATH, index=False)

if __name__ == '__main__':

    config = yaml.safe_load(open("../config.yml"))

    process_dssp_data(config['path']['dssp_csv'], config['path']['processed_data']+'dssp_data.csv')
