from Bio.SeqUtils import seq1
import random

class ResiduePatch():
    """
    Parse pisite files to dict

    Attributes
    ----------
    ids : list
        residue ids
    file: str
        pdb file
    dssp: dict
        dssp dict
    dssp_dict_keys: dict
        dssp dict keys
    """

    def __init__(self, ids, dssp, dssp_dict_keys):
        self.ids = ids
        self.dssp = dssp
        self.dssp_dict_keys = dssp_dict_keys

    def get_ids(self):
        return self.ids

    def size(self):
        return sum(self.dssp[i[-2:]][2] for i in self.ids)

    def residues(self):
        return [self.dssp[i[-2:]][0] for i in self.ids]

    def patch_length(self):
        return len(self.ids)

    def residue_on_surface(self):
        """
        Get all residues accessible from the surface

        Return
        ------
        list
            residue ids
        """
        residue_list = []
        for x in self.dssp_dict_keys:
            if x[1][0] != ' ' or  x[1][2] != ' ':
                continue
            if self.dssp[x][2] > 0:
                residue_list.append(x)
        return residue_list

    def check_pisite(self, pisite_dict):
        """
        Get all interaction sites

        Attributes
        ----------
        pisite_dict: pisite dict

        Return
        ------
        int
            interaction sites
        """
        interaction_sites = 0
        for i in self.ids:
            id = i[-1][1]
            if self.dssp[i[-2:]][0] == pisite_dict[id][1]:
                if int(pisite_dict[id][2]) >= 1:
                    interaction_sites += 1

        return interaction_sites

    def random_patch_ppis(self, pisite_dict):
        """
        Get all interaction sites

        Attributes
        ----------
        pisite_dict: pisite dict

        Return
        ------
        int
            fraction of interaction site in a random patch
        """
        residues_in_patch = self.patch_length()
        random_selected = [0]*(len(self.residue_on_surface())-residues_in_patch) + [1]*residues_in_patch

        random.shuffle(random_selected)

        random_patch_size = 0

        for i,j in enumerate(self.residue_on_surface()):
            id = j[-1][1]
            if self.dssp[j][0] == pisite_dict[id][1]:
                if random_selected[i] == 1 and int(pisite_dict[id][2]) == 1:
                    random_patch_size += 1
        return random_patch_size
