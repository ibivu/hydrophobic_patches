from Bio import SeqIO
import pandas as pd

def process_fasta_data(FASTA_FILE, PROCESSED_DATA_PATH):
    """create csv with full sequence

    Parameters
    ----------
    FASTA_FILE : str
        location of the fasta file
    PROCESSED_DATA_PATH : str
        output file

    """
    records = list(SeqIO.parse(FASTA_FILE, "fasta"))
    results = {
            'id':[x.id for x in records],
            'fasta_sequence':[str(x.seq) for x in records],
            }

    df = pd.DataFrame(results)
    df.to_csv(PROCESSED_DATA_PATH, index=False)

if __name__ == '__main__':
    import yaml
    config = yaml.safe_load(open("../config.yml"))

    process_fasta_data(config['path']['fasta'], config['path']['processed_data']+'fasta_data.csv')
