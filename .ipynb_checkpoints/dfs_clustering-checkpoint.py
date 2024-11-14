import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array

from src.MoleculeSystem import MoleculeSystem

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
        n_atoms = self.ag.n_atoms // n_res
        clusters = []

        positions = np.zeros((len(self.u.trajectory), n_res, n_atoms, 3))  # positions of carbons in rings
        for ts in self.u.trajectory:
            checked_residues = []

            for i in range(n_res):
                for j in range(n_atoms):
                    positions[ts.frame, i, j] = self.ag[i * n_atoms + j].position  # WRONG!!!

            for i, position in enumerate(positions[ts.frame]):
                if i not in checked_residues:
                    clusters.append(self._dfs(i, positions[ts.frame], [i], checked_residues))
            break
        return clusters


if __name__ == '__main__':
    from os.path import join

    BASE_PATH = r"E:\C_backup\Documents-2\Documents\Simulations\tyloxapol_tx\tyl_7"
    SYSTEM_NAME = "100pc_tyl"

    atom_names = """name 
    C0 C1 C2 C3 C4 C5 C6 C7 C8
    C17 C18 C19 C20 C21 C22 C23
    C34 C35 C36 C37 C38 C39 C40
    C51 C52 C53 C54 C55 C56 C57
    C68 C69 C70 C71 C72 C73 C74
    C85 C86 C87 C88 C89 C90 C91 
    C102 C103 C104 C105 C106 C107 C108"""

    atoms_per_mol = len(atom_names)
    top = join(BASE_PATH, SYSTEM_NAME, "centered.gro")
    trj = join(BASE_PATH, SYSTEM_NAME, "centered.xtc")
    s = ClusterSearch(top, trj)
    s.select_atoms(atom_names)
    clusters = s.find_clusters()
    print(clusters)
    for i, c in enumerate(clusters):
        clusters[i] = [e + 1 for e in c]
        clusters[i].sort()

    [print(f"Cluster {i}: ", *c) for i, c in enumerate(clusters)]
