import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import math
import yaml
config = yaml.safe_load(open("../config.yml"))

hydrophobic_proteins = config['hydrophobic']
polar_proteins = config['polar']

buried = {'A':86.600, 'R':162.200, 'N':103.300, "D":97.800, "C":132.300,
            "Q":119.200, "E":113.900, "G":62.900, "H":155.800, "I":158.000,
            "L":164.100, "K":115.500, "M":172.900, "F":194.100, "P":92.900,
            "S":85.600, "T":106.500, "W":224.600, "Y":177.700, "V":141.000}

def entropy(string):
    # "Calculates the Shannon entropy of a string"
    # get probability of chars in string
    prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string)) ]
    # calculate the entropy
    entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])
    return entropy

def entropy_ideal(length):
    # "Calculates the ideal Shannon entropy of a string with given length"
    prob = 1.0 / length
    return -1.0 * length * prob * math.log(prob) / math.log(2.0)

def get_features(seq):
    """get global features from a protein sequence

    Parameters
    ----------
    seq : str
        protein sequence

    Return
    ----------
    dictionary:
        global features of the protein sequence

    """

    features = {}
    features['undefined_count'] = len([x for x in seq if x in ['X','B','Z',"'",'O','U']])
    features['length'] = len(seq)
    features['perc_undefined_count'] = features['undefined_count']/features['length']
    features['entropy'] = entropy(seq)
    features['ideal_entropy'] = entropy_ideal(len(seq))
    features['perc_entropy'] = features['entropy']/features['ideal_entropy']
    features['hydr_count'] = sum(1 for x in seq if x in hydrophobic_proteins)
    features['polar_count'] = sum(1 for x in seq if x in polar_proteins)
    features['buried'] = sum(buried[x] for x in seq if x in hydrophobic_proteins)

    seq = ''.join([x for x in seq if x not in ['X','B','Z',"'",'O','U']])

    protein = ProteinAnalysis(seq)
    features['gravy'] = protein.gravy()
    features['molecular_weight'] = protein.molecular_weight()
    features['aromaticity'] = protein.aromaticity()
    features['instability_index'] = protein.instability_index()
    features['isoelectric_point'] = protein.isoelectric_point()
    features['helix'], features['turn'], features['sheet'] = protein.secondary_structure_fraction()

    features.update(protein.count_amino_acids())
    # features.update(protein.get_amino_acids_percent())
    return features

def process_global_sequence_data(PROCESSED_DSSP_DATA, PROCESSED_GLOBAL_DATA):
    """create csv with full sequence

    Parameters
    ----------
    PROCESSED_DSSP_DATA : str
        location of the processed dssp file
    PROCESSED_GLOBAL_DATA : str
        output file

    """
    df_dssp = pd.read_csv(PROCESSED_DSSP_DATA)
    features = lambda x: get_features(x)
    protein_features = df_dssp['dssp_sequence'].apply(features).apply(pd.Series)
    protein_features['id'] = df_dssp['id']

    protein_features.to_csv(PROCESSED_GLOBAL_DATA,  index=False)

if __name__ == '__main__':
    config = yaml.safe_load(open("../config.yml"))
    process_global_sequence_data(config['path']['processed_data']+'dssp_data.csv', config['path']['processed_data']+'global_seq_features_data.csv')
