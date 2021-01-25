from Bio import SeqIO
import pandas as pd
from os import listdir
from os.path import isfile, join
from Bio.PDB.MMCIFParser import MMCIFParser

def process_monomer_data(PDB_PROTEIN_PATH, PROCESSED_DATA_PATH):
    """check if pdb protein is a monomer

    Parameters
    ----------
    PDB_PROTEIN_PATH : str
        location of the directory with full pdb files
    PROCESSED_DATA_PATH : str
        output file

    """

    protein_files = [f for f in listdir(PDB_PROTEIN_PATH) if isfile(join(PDB_PROTEIN_PATH, f))]
    results = {'id':[], 'monomer':[]}

    #for all files in directory
    for protein_file in protein_files:
        protein_id = protein_file[:-4]
        pdb_protein_file = PDB_PROTEIN_PATH+protein_file

        # if only 1 chain, pdb = monomer
        try:
            p = MMCIFParser(QUIET=1)
            structure = p.get_structure(protein_id, pdb_protein_file)
            results['id'].append(protein_id)
            if len(structure[0]) == 1:
                results['monomer'].append(True)
            else:
                results['monomer'].append(False)
        except:
            continue

    df = pd.DataFrame(results)
    df.to_csv(PROCESSED_DATA_PATH, index=False)

if __name__ == '__main__':
    import yaml
    config = yaml.safe_load(open("../config.yml"))
    process_monomer_data(config['path']['protein'], config['path']['processed_data']+'monomer_data.csv')
