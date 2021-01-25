from AtomBased.ProteinPatch import ProteinPatch
from PisiteParser import PisiteParser
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from os import listdir
from os.path import isfile, join
import pandas as pd

config = yaml.safe_load(open("../config.yml"))

atoms = ['C', 'S']

path = '/home/jan/Documents/BioInformatics/final_project_patch/data/pdb/hoh/'
pisite_path = '/home/jan/Documents/BioInformatics/final_project_patch/data/pisite/nr/'
csv_file = '../../data/patches/lp_pisite_atom.csv'

files = [f for f in listdir(path) if isfile(join(path, f))]

df = pd.read_csv(csv_file)
for file in files:
    try:
        pisite = file[:-5].lower()
        pisite += file[4].upper()+'.pisite'

        pisiteParser = PisiteParser(pisite_path+pisite)
        print(file, 'interaction_sites', pisiteParser.interaction_sites())
        if not pisiteParser.interaction_sites():
            print(file, 'no interaction sites')
            continue

        pisite_data = pisiteParser.get_data()
        print(file)
        proteinPatches = ProteinPatch(id,path+file,atoms)
        patches = proteinPatches.patches()
        patches = sorted(patches, key=(lambda x:x.size()), reverse=True)

        for i,patch in enumerate(patches[:20]):

            result_dict = patch.check_pisite(pisite_data)
            result_dict['id'] = file[:-4]
            result_dict['rank'] = int(i+1)
            result_dict['total_ppis'] = pisiteParser.interaction_sites()
            result_dict['random'] = patch.random_patch(pisite_data)
            df = df.append(pd.Series(result_dict), ignore_index=True)

        df.to_csv(csv_file, index=False)
    except:
        print(file, 'failed')
