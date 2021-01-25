import Bio.PDB.ResidueDepth as ResidueDepth
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB import MMCIFParser
from scipy.spatial import KDTree
import networkx as nx
from Bio.SeqUtils import seq1, seq3
import numpy as np
from mayavi import mlab
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from AtomBased.Patch import AtomPatch
import freesasa
import tempfile
import os

class ProteinPatch():
    """
    Parse pisite files to dict

    ...

    Attributes
    ----------
    id : str
        protein id
    file: str
        pdb file
    atom_selection: list
        list of allowed atoms in patch
    r: float
        search radius
    msms: str
        msms command
    """
    def __init__(self, id, file, atom_selection, msms = 'msms -density 1.5', r=1.25):

        if file.endswith('.cif'):
            parser = MMCIFParser()
        else:
            parser = PDBParser()

        structure = parser.get_structure(id, file)
        self.model = structure[0]
        self.structure = freesasa.structureFromBioPDB(self.model)
        self.result = freesasa.calc(self.structure)
        self.atom_selection = atom_selection
        self.r = r
        self.msms = msms
        self.G = self.labeled_dot_cloud_graph(self.model, self.msms)
        self.G = self.patch_network(self.G)

    def labeled_dot_cloud_graph(self, model, msms):
        """
        Create a dotted cloud of the proteins surface and label the dots
        hydrophobic or hydrophilic
        ...

        Attributes
        ----------
        model : BIO.pdb
            PDB model
        msms: str
            msms command

        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic or not
        """

        self.surface_vertices = get_surface(model, MSMS=msms)
        atom_list = Selection.unfold_entities(model, "A")
        atom_vectors = [list(v.get_vector()) for v in atom_list]

        # create KDtree to search for the closest atoms to a surface vector
        T = KDTree(atom_vectors)
        closest_atoms = T.query(self.surface_vertices, k=2)

        G = nx.Graph()
        # for every surface dot
        for i in range(len(self.surface_vertices)):
            score = 0
            G.add_node(i)
            close_atoms = closest_atoms[1][i]
            # add node attributes
            G.node[i]['surface_vector_pos'] = self.surface_vertices[i]
            G.node[i]['closest_atom_id'] = close_atoms[0]
            G.node[i]['residue_id'] = atom_list[close_atoms[0]].get_full_id()[3][1]
            G.node[i]['residue'] = atom_list[close_atoms[0]].get_parent().get_resname()

            #check if the 2 closest atoms are hydrophobic
            score = sum(1 for x in close_atoms if atom_list[x].get_name()[0] in self.atom_selection)

            if score >= 2:
                G.node[i]['selected'] = 1
            else:
                G.node[i]['selected'] = 0

        return G

    def patch_network(self, G):
        """
        Create a edges between hydrophobic nodes if they are within r
        hydrophobic or hydrophilic
        ...

        Attributes
        ----------
        G : networkx.Graph
            Graph with nodes labeled as hydrophobic or not

        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic and edges between nodes within range r
        """
        node_list = [i for i in G.nodes if G.node[i]['selected']]
        x = [G.nodes[i]['surface_vector_pos'] for i in G.nodes if G.node[i]['selected']]
        T = KDTree(x)
        pairs = T.query_pairs(self.r)
        G.add_edges_from([(node_list[x[0]],node_list[x[1]]) for x in pairs])
        return G

    def largest_patch(self):
        """
        get largest patch

        Return
        ------
        networkx.Graph
            Graph of the largest patch
        """
        patches = self.patches()
        largest_patch = max(patches, key=(lambda x: x.size()))
        return largest_patch

    def patches(self):
        """
        Get the components of the graph

        Return
        ------
        list
            list of Graph components where an item is a patch
        """
        patched_G = self.patch_network(self.G)
        components = nx.connected_component_subgraphs(patched_G)
        patch_dict = []
        for component in components:
            if len(component.nodes) <= 1:
                continue
            atom_ids_in_patch = list(set([component.nodes[i]['closest_atom_id'] for i in component.nodes]))
            patch_dict.append(AtomPatch(atom_ids_in_patch, self.structure, self.result))
        return patch_dict

    def plot_largest_patches(self):
        """
        Plot the largest patch
        """
        largest_patch = self.largest_patch()
        print('largest_patch',largest_patch.size())
        #change color of largest patch
        for node in self.G.nodes:
            if self.G.node[node]['closest_atom_id'] in largest_patch.get_ids():
                self.G.node[node]['selected'] = 2

        #plot graph in three colors
        G = self.G
        xyz = np.array([G.node[v]['surface_vector_pos'] for v in sorted(G)])
        # scalar colors
        scalars = np.array([int(G.node[node]['selected']) for node in G.nodes]) + 2
        mlab.figure(1, bgcolor=(0.5, 0.5, 0.5))
        mlab.clf()
        pts = mlab.points3d(xyz[:, 0], xyz[:, 1], xyz[:, 2],
                            scalars,
                            scale_factor=0.25,
                            scale_mode='none',
                            resolution=20,
                            colormap='coolwarm')

        pts.mlab_source.dataset.lines = np.array(list(G.edges()))
        tube = mlab.pipeline.tube(pts, tube_radius=0.05)
        mlab.pipeline.surface(tube,colormap='Reds')
        mlab.show()
