from os.path import join

import networkx as nx
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import distance_array, self_distance_array

from MoleculeSystem import MoleculeSystem
object

class SystemGraph(MoleculeSystem):

    def __init__(self, top_path: str, trj_path=None):
        self.ag = None
        if not trj_path:
            self.u = mda.Universe(top_path)
        else:
            self.u = mda.Universe(top_path, trj_path)

    def select_atoms(self, selection: str) -> None:
        self.ag = self.u.select_atoms(selection)

    def graph(self):
        matrix = distance_array(self.ag.positions, self.ag.positions)
        processed = self._process_array(matrix, cutoff=4.7)
        np.fill_diagonal(processed, 0)

        n_mol = self.ag.n_residues
        n_atoms = self.ag.n_atoms // n_mol

        temp = np.zeros(shape=(n_mol, n_mol))

        for i in range(0, n_mol):
            for j in range(0, n_mol):
                if i == j:
                    continue
                temp[i, j] = processed[i * n_atoms: (i + 1) * n_atoms, j * n_atoms: (j + 1) * n_atoms].sum()

        # print("Median of temp:", np.median(temp))
        # print("Max of temp:", temp.max())
        # print("Median / Max perc:", np.median(temp) / temp.max() * 100)  # 3%
        # print("25th percentile:", np.percentile(temp, 50))

        # c_off = np.ceil(np.mean(temp) / temp.max() * 100) / 100
        # print(c_off)
        np.savetxt("../temp.txt", temp, fmt='%.i')

        # adj_matrix = self._process_array(temp, cutoff=c_off)
        # np.savetxt("adj_matrix_compact.txt", adj_matrix, fmt='%.i')
        # np.fill_diagonal(adj_matrix, 0)
        G = nx.Graph(temp)
        clusters = list(nx.connected_components(G))
        for i, cluster in enumerate(clusters):
            print(f'Cluster {i}: ', *[c + 1 for c in cluster])

        pos = nx.spring_layout(G, k=0.15, iterations=20)
        nx.draw(G, pos=pos, with_labels=True)
        plt.show()

    @staticmethod
    @np.vectorize
    def _process_array(a, if_value=1, else_value=0, cutoff=5):
        return if_value if a <= cutoff else else_value

str
if __name__ == "__main__":
    BASE_PATH = r"C:\Users\hrach\Documents\Simulations\naol_micelles"
    SYSTEM_NAME = "200-ctab-micelle"
    # atom_names = "name O8 O18 O28 O38 O48 O58 O68 O3 O13 O23 O33 O43 O53 O63 C1 C22 C39 C56 C73 C90 C107"
    # atom_names = "name C1 C22 C39 C56 C73 C90 C107"  # core only
    #    atom_names = """name C51 C52 C53 C54 C55 C56 C57 C58 C59 C60 C61 C62 C63 C64 C65 C103 C104 C105 C106 C107 C108 C109
    # C110 C111 C112 C113 C114 C115 C116 C69 C70 C71 C72 C73 C74 C75 C76 C77 C78 C79 C80 C81 C82 C68 C1 C2 C3 C4 C5 C6
    # C7 C8 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C34 C35 C36 C37 C38 C39 C40
    #  C41 C42 C43 C44 C86 C87 C88 C89 C90 C91 C92 C93 C94 C95 C96 C97 C98 C99 C47 C46 C45 O1 C48 C85 C102"""
    atom_names = """resname CTB16"""

    atoms_per_mol = len(atom_names)
    top = join(BASE_PATH, SYSTEM_NAME, "clustered.gro")
    trj = join(BASE_PATH, SYSTEM_NAME, "clustered_skip5.xtc")

    m = SystemGraph(top)  # , trj)
    m.select_atoms(atom_names)
    m.graph()

# Cluster 0:  1 4 5 12 13 16 18 19 21
# Cluster 1:  2 3 6 7 8 20
# Cluster 2:  9 10 11 14
# Cluster 3:  15
# Cluster 4:  17


# Cluster 0:  1 2 3 6 9 10 14 18
# Cluster 1:  4 5 12 13 16 19 21
# Cluster 2:  7
# Cluster 3:  8
# Cluster 4:  11
# Cluster 5:  15
# Cluster 6:  17
# Cluster 7:  20