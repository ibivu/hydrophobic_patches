from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from contextlib import closing
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
import pandas as pd
import os.path

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0

def create_chain_pdb(FASTA_FILE, PDB_PROTEIN_PATH, PDB_CHAIN_PATH):
    """Create a new pdb file of one chain

    Parameters
    ----------
    FASTA_FILE : str
        The file location of the pisces file
    PDB_PROTEIN_PATH : str
        The directory location of full pdb
    PDB_CHAIN_PATH : str
        The directory location of chain pdb
    """
    # get fasta records
    records = list(SeqIO.parse(FASTA_FILE, "fasta"))

    for record in records:
        chain_char = record.id[-1]
        protein_id = record.id[:-1]
        chain_id = record.id
        pdb_protein_file = (PDB_PROTEIN_PATH+protein_id+'.cif').lower()
        pdb_chain_file = (PDB_CHAIN_PATH+record.id+'.ent')

        if os.path.isfile(pdb_chain_file):
            print ("File exist", pdb_chain_file)
            continue

        try:
            #Filter for one chain and create a new pdb file
            p = MMCIFParser(QUIET=1)
            structure = p.get_structure(protein_id, pdb_protein_file)
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save(pdb_chain_file, ChainSelect(chain_char))
            print('Saved', chain_id)
        except:
            print('Error', pdb_protein_file)

if __name__ == '__main__':
    import yaml
    import sys

    config = yaml.safe_load(open("../config.yml"))

    create_chain_pdb(config['path']['fasta'], config['path']['protein'], config['path']['chain'])
