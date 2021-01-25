from ResidueBased.ProteinPatch import ProteinPatch
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from os import listdir
from os.path import isfile, join
from PisiteParser import PisiteParser
import pandas as pd

config = yaml.safe_load(open("../config.yml"))

hydr_residues = config['hydrophobic']

path = '/home/jan/Documents/BioInformatics/final_project_patch/data/pdb/chain/'
pisite_path = '/home/jan/Documents/BioInformatics/final_project_patch/data/pisite/nr/'
csv_file = '../../data/patches/lp_pisite_residue_with_aa.csv'

files = [f for f in listdir(path) if isfile(join(path, f))]

df = pd.read_csv(csv_file)

for file in files[0:]:
    try:
        pisite = file[:-5].lower()
        pisite += file[4].upper()+'.pisite'

        pisiteParser = PisiteParser(pisite_path+pisite)
        print(file, 'interaction_sites', pisiteParser.interaction_sites())
        if not pisiteParser.interaction_sites():
            print(file, 'no interaction sites')
            continue

        pisite_data = pisiteParser.get_data()
        proteinPatches = ProteinPatch(id,path+file,hydr_residues)
        patches = proteinPatches.patches
        for i,patch in enumerate(patches[:20]):
            result_dict = {}
            result_dict['ppi_in_patch'] = patch.check_pisite(pisite_data)
            result_dict['patch_size'] = patch.patch_length()
            result_dict['on_surface'] = len(patch.residue_on_surface())
            result_dict['id'] = file[:-4]
            result_dict['rank'] = int(i+1)
            result_dict['residues'] = ''.join(patch.residues())
            result_dict['total_ppis'] = pisiteParser.interaction_sites()
            result_dict['random_patch_size'] = patch.random_patch_ppis(pisite_data)
            result_dict['size'] = patch.size()
            df = df.append(pd.Series(result_dict), ignore_index=True)

        df.to_csv(csv_file, index=False)
    except:
        print(file, 'failed')
