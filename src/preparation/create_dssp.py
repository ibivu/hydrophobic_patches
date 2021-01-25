import os
from os import listdir
from os.path import isfile, join
import pandas as pd
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

def create_dssp_dssp(pdb_chain_file, dssp_dssp_file):
    """create a dssp file

    Parameters
    ----------
    pdb_chain_file : str
        The file location of the pdb chain
    dssp_dssp_file : str
        The file location of the output dssp
    """

    os.system("mkdssp -i "+pdb_chain_file+" -o "+dssp_dssp_file)
    print('Created', dssp_dssp_file)

def create_dssp_csv(pdb_chain_file, dssp_csv_file):
    """create a dssp csv file

    Parameters
    ----------
    pdb_chain_file : str
        The file location of the pdb chain
    dssp_csv_file : str
        The file location of the output dssp_csv
    """
    print(pdb_chain_file)
    values, keys = dssp_dict_from_pdb_file(pdb_chain_file)
    data = [x+y for x, y in zip(keys, values.values())]
    pd.DataFrame(data).to_csv(dssp_csv_file, index=False)

def create_dssp(PDB_CHAIN_PATH, DSSP_DSSP_PATH, DSSP_CSV_PATH):
    """create a dssp csv file

    Parameters
    ----------
    PDB_CHAIN_PATH : str
        The directory location of chain pdb
    DSSP_DSSP_PATH : str
        The directory location of the output dssp
    DSSP_CSV_PATH : str
        The directory location of the output dssp_csv
    """

    #get all pdb chains in directory
    chain_files = [f for f in listdir(PDB_CHAIN_PATH) if isfile(join(PDB_CHAIN_PATH, f)) and not '.git' in f]

    #create dssp and dssp_csv for every chain pdb
    for chain_file in chain_files:
        chain_id = chain_file[:-4]
        pdb_chain_file = PDB_CHAIN_PATH+chain_file
        dssp_dssp_file = DSSP_DSSP_PATH+chain_id+".dssp"
        dssp_csv_file = DSSP_CSV_PATH+chain_id+".csv"

        if os.path.isfile(dssp_dssp_file) and os.path.isfile(dssp_csv_file):
            print ("Files exist", chain_id)
        elif os.path.isfile(dssp_csv_file):
            create_dssp_dssp(pdb_chain_file, dssp_dssp_file)
            print("DSSP created", chain_id)
        elif os.path.isfile(dssp_dssp_file):
            create_dssp_csv(pdb_chain_file, dssp_csv_file)
            print("DSSP CSV created", chain_id)
        else:
            create_dssp_dssp(pdb_chain_file, dssp_dssp_file)
            create_dssp_csv(pdb_chain_file, dssp_csv_file)
            print("DSSP CSV and DSSP created", chain_id)

if __name__ == '__main__':
    import yaml
    import sys

    config = yaml.safe_load(open("../config.yml"))

    create_dssp(config['path']['chain'], config['path']['dssp_dssp'], config['path']['dssp_csv'])
