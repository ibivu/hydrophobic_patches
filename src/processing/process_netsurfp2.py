import pandas as pd
import yaml
import os
import glob

config = yaml.safe_load(open("../config.yml"))
hydrophobic_aa = config['hydrophobic']

def process_netsurfp2(NSP2_PATH, PROCESSED_DATA_PATH):
    """retrieve the predicted TASA, THSA and RHSA from NetSurfP2 output files

    Parameters
    ----------
    NSP2_PATH : str
        location of the directory with netsurfp2 files
    PROCESSED_DATA_PATH : str
        output file

    """

    extension = 'csv'
    os.chdir(NSP2_PATH)
    result = glob.glob('*.{}'.format(extension))
    print(result)
    df = pd.concat([pd.read_csv(x) for x in result])
    df['seq_hydr'] = df['seq'].isin(hydrophobic_aa)
    df2 = df[df['seq_hydr'] == True].groupby(['id']).sum()['asa'].reset_index()
    df3 = df2.merge(df.groupby(['id']).sum()['asa'].reset_index(), on='id', how='inner')
    df4 = df3.merge(df.groupby(['id']).sum()[['p[q3_H]','p[q3_E]','p[q3_C]']].reset_index(), on='id', how='inner')
    df4 = df4.rename(columns={
                            'asa_x':'thsa_netsurfp2',
                            'asa_y':'tasa_netsurfp2',
                            'p[q3_H]':'q3_H',
                            'p[q3_E]':'q3_E',
                            'p[q3_C]':'q3_C'
                            })
    df4['rhsa_netsurfp2'] = df4['thsa_netsurfp2']/df4['tasa_netsurfp2']

    df4.to_csv(PROCESSED_DATA_PATH, index=False)

if __name__ == '__main__':
    import yaml
    config = yaml.safe_load(open("../config.yml"))

    process_netsurfp2(config['path']['netsurfp2'], config['path']['processed_data']+'netsurp2_data.csv')
