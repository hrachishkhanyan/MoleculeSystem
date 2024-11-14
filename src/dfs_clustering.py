import pickle

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array

from MoleculeSystem import MoleculeSystem

cutoff = 4.7


class ClusterSearch(MoleculeSystem):

    def __init__(self, top_path: str, trj_path=None):
        self.ag = None
        if not trj_path:
            self.u = mda.Universe(top_path)
        else:
            self.u = mda.Universe(top_path, trj_path)

    def select_atoms(self, selection: str) -> None:
        self.ag = self.u.select_atoms(selection)

    def _dfs(self, index: int, positions: np.ndarray, cluster: list, checked: list) -> list:
        checked.append(index)

        for i in range(len(positions)):

            if i != index and i not in checked:
                distances = distance_array(positions[index], positions[i])

                if (distances < cutoff).any():
                    cluster.append(i)
                    self._dfs(i, positions, cluster, checked)

        return cluster

    def find_clusters(self):
        # u = mda.Universe(top, trj)
        # ag = u.select_atoms(selection)
        n_res = self.ag.n_residues
        resids = np.unique(self.ag.resids)
        n_atoms = self.ag.n_atoms // n_res
        clusters = []

        positions = np.zeros((len(self.u.trajectory), n_res, n_atoms, 3))  # positions of carbons in rings
        for ts in self.u.trajectory:
            current_frame_clusters = []
            checked_residues = []

            for i in range(n_res):
                for j in range(n_atoms):
                    positions[ts.frame, i, j] = self.ag[i * n_atoms + j].position  # WRONG!!!  ????

            for i, position in enumerate(positions[ts.frame]):
                if i not in checked_residues:
                    current_frame_clusters.append(resids[self._dfs(i, positions[ts.frame], [i], checked_residues)])
            clusters.append(current_frame_clusters)
        return clusters


if __name__ == '__main__':
    from os.path import join

    BASE_PATH = r"C:\Users\hrach\Documents\Simulations\naol_micelles"
    SYSTEM_NAME = "200-ctab-micelle"

    atom_names = """resname CTB16 and not type H"""

    atoms_per_mol = len(atom_names)
    top = join(BASE_PATH, SYSTEM_NAME, "clustered.gro")
    trj = join(BASE_PATH, SYSTEM_NAME, "clustered_skip5.xtc")
    s = ClusterSearch(top, trj)
    s.select_atoms(atom_names)
    clusters = s.find_clusters()

    with open('../ctab_micelle_clusters.pkl', 'wb') as f:
        pickle.dump(clusters, f)

    # [print(f"Cluster {i}: ", *c) for i, c in enumerate(clusters)]
