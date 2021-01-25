import pandas as pd
import yaml
from functools import reduce
import os
import glob

config = yaml.safe_load(open("../config.yml"))
hydrophobic_aa = config['hydrophobic']

def process_all_data_for_predictions():
    """combine preprocessed features in one DataFrame

    Returns
    -------
    pandas.DataFrame
        all global features combined
    """
    df_DSSP = pd.read_csv(config['path']['processed_data']+'dssp_data.csv')
    df_TMHMM = pd.read_csv(config['path']['processed_data']+'tmhmm_data.csv')
    df_fasta = pd.read_csv(config['path']['processed_data']+'fasta_data.csv')
    df_nsp2 = pd.read_csv(config['path']['processed_data']+'netsurp2_data.csv')
    df_monomer = pd.read_csv(config['path']['processed_data']+'monomer_data.csv')
    df_global_seq_features = pd.read_csv(config['path']['processed_data']+'global_seq_features_data.csv')

    df_patches = pd.read_csv(config['path']['patch']+'lp_residue_all.csv')
    df_patches = df_patches[df_patches['rank'] == 1]

    data_frames = [df_DSSP, df_TMHMM, df_global_seq_features, df_fasta, df_nsp2, df_patches]

    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['id'],
                                                how='inner'), data_frames)

    df_no_tmp = df_merged[~df_merged['tmp']]

    df_no_tmp.to_csv(config['path']['processed_data']+'ready_to_use_data.csv', index=False)

    return df_no_tmp

def process_all_data_for_predictions_test_train():
    """split combined features in train and test set

    """

    train_file = config['path']['data'] + 'csv/train/train.csv'
    test_file = config['path']['data'] + 'csv/test/test.csv'

    train_old = pd.read_csv(train_file)
    test_old = pd.read_csv(test_file)
    df = process_all_data_for_predictions()

    train = df[df['id'].isin(train_old['id'])]
    test = df[df['id'].isin(test_old['id'])]

    train.to_csv(config['path']['processed_data']+'ready_to_use_data_train.csv', index=False)
    test.to_csv(config['path']['processed_data']+'ready_to_use_data_test.csv', index=False)

    print(len(train), len(test))

if __name__ == '__main__':
    import yaml
    config = yaml.safe_load(open("../config.yml"))
    process_all_data_for_predictions_test_train()
