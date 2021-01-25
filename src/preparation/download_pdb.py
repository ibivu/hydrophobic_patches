from Bio.PDB import PDBList
from Bio import SeqIO

def download_pdb_files(FASTA_FILE, PDB_OUTPUT):
    """Downloads all pdb files from a given fasta file with pdb ids

    Parameters
    ----------
    FASTA_FILE : str
        The file location of the .fasta
    PDB_OUTPUT : str
        The directory for pdb outputs

    """

    # get protein ids
    records = list(SeqIO.parse(FASTA_FILE, "fasta"))
    protein_ids = set([x.id[:-1] for x in records])

    # download pdbs
    pdbl = PDBList()
    pdbl.download_pdb_files(protein_ids, pdir=PDB_OUTPUT)

if __name__ == '__main__':
    import yaml
    import sys

    # retrieve config values
    config = yaml.safe_load(open("../config.yml"))
    
    download_pdb_files(config['path']['fasta'], config['path']['protein'])
