from Bio.SeqUtils import seq1
import random

class AtomPatch():
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

    def __init__(self, ids, structure, result):
        self.ids = ids
        self.structure = structure
        self.result = result
        self.residue_dict = self.create_residue_dict()

    def get_ids(self):
        return self.ids

    def size(self):
        return sum(self.result.atomArea(i) for i in self.ids)

    def residue_on_surface(self):
        """
        Get all residues accessible from the surface

        Return
        ------
        list
            residue ids
        """
        s = self.structure
        n_atoms = s.nAtoms()
        residue_dict = []
        for i in range(n_atoms):
            if self.result.atomArea(i) > 0:
                residue_dict.append(s.residueNumber(i))
        return len(set(residue_dict))

    def create_residue_dict(self):
        """
        Get the residues corresponding to the atom

        Return
        ------
        dict
            {atom id: residue attributes}
        """
        s = self.structure
        n_atoms = s.nAtoms()
        residue_dict = {}
        for i in range(n_atoms):
            number = int(s.residueNumber(i))
            if number not in residue_dict.keys():
                residue_dict[number] = {'name':s.residueName(i),
                                        'char':seq1(s.residueName(i)),
                                        'selected':False}
            if i in self.ids:
                residue_dict[number]['selected'] = True

        return residue_dict

    def print_residues(self):
        """
        Get all residues accessible from the surface

        Return
        ------
        list
            residue ids
        """
        residues = ''
        rd = self.residue_dict

        for i in rd:
            if rd[i]['selected']:
                residues += rd[i]['char']
            else:
                residues += rd[i]['char'].lower()

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
        interaction_site = 0
        in_patch = 0
        total_interaction_sites = 0

        for i in self.residue_dict:
            if self.residue_dict[i]['char'] == pisite_dict[i][1]:
                if int(self.residue_dict[i]['selected']) == 1 and int(pisite_dict[i][2]) == 1:
                    interaction_site += 1
                if self.residue_dict[i]['selected']:
                    in_patch += 1
                if int(pisite_dict[i][2]) == 1:
                    total_interaction_sites += 1

        return {'len':len(self.residue_dict),
                'residues_in_patch':in_patch,
                'residues_on_surface':self.residue_on_surface(),
                'interaction_site':interaction_site,
                'total_interaction_sites':total_interaction_sites,
                'size':self.size()}

    def random_patch(self, pisite_dict):
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
        residues_in_patch = sum(1 for x in self.residue_dict if self.residue_dict[x]['selected'])
        my_list = [0]*(len(self.residue_dict)-residues_in_patch) + [1]*residues_in_patch
        random.shuffle(my_list)

        random_patch_size = 0

        for j,i in enumerate(self.residue_dict):
            if self.residue_dict[i]['char'] == pisite_dict[i][1]:
                if my_list[j] == 1 and int(pisite_dict[i][2]) == 1:
                    random_patch_size += 1
        return random_patch_size
