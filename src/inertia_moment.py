from typing import Union
from numpy import linalg as lng
import numpy as np
import MDAnalysis as mda

from src.MoleculeSystem import MoleculeSystem


class Inertia(MoleculeSystem):
    shape_mapping = {
        0: "Spherical",
        1: "Oblate",
        2: "Prolate",
        3: "Triaxial"
    }

    def __init__(self, top_path: str, trj_path=None):
        self.ag = None
        if not trj_path:
            self.u = mda.Universe(top_path)
        else:
            self.u = mda.Universe(top_path, trj_path)
        self.moments: Union[np.ndarray, None] = None

    def select_atoms(self, selection: str) -> None:
        self.ag = self.u.select_atoms(selection)

    def _calculate_moment(self) -> np.ndarray:
        moment_of_inertia = self.ag.moment_of_inertia()
        eigenvalues = lng.eig(moment_of_inertia)[0]
        return eigenvalues

    def _shape(self, tolerance=None) -> int:

        i_1, i_2, i_3 = self.moments.mean(axis=0)
        print(i_1, i_2, i_3)
        if not tolerance:
            tolerance = 0.05 * np.mean((i_1, i_2, i_3))

        if abs(i_1 - i_2) < tolerance and abs(i_1 - i_3) < tolerance and abs(i_2 - i_3) < tolerance:
            return 0
        if ((abs(i_1 - i_2) < tolerance < i_3 - i_1 and tolerance < i_3 - i_2 and i_3 > i_1 and i_3 > i_2) or
                (abs(i_3 - i_2) < tolerance < i_1 - i_3 and tolerance < i_1 - i_2 and i_1 > i_3 and i_1 > i_2) or
                (abs(i_3 - i_1) < tolerance < i_2 - i_1 and tolerance < i_2 - i_3 and i_2 > i_3 and i_2 > i_1)):
            return 1
        if ((abs(i_1 - i_2) < tolerance < i_1 - i_3 and tolerance < i_2 - i_3 and i_3 < i_1 and i_3 < i_2) or
                (abs(i_3 - i_2) < tolerance < i_2 - i_1 and tolerance < i_3 - i_1 and i_1 < i_3 and i_1 < i_2) or
                (abs(i_3 - i_1) < tolerance < i_1 - i_2 and tolerance < i_3 - i_2 and i_2 < i_3 and i_2 < i_1)):
            return 2

        return 3

    def determine_shape(self, tolerance=None) -> str:
        return self.shape_mapping[self._shape(tolerance)]

    def run(self, start=0, end=None, skip=1) -> None:
        print("Calculating moment of inertia")
        if not end:
            end = self.u.trajectory.n_frames

        moments: list = []

        for _ in self.u.trajectory[slice(start, end, skip)]:
            moments.append(self._calculate_moment())

        self.moments = np.array(moments)

    def eccentricity(self):

        if self.moments is None:
            raise ValueError("Use the run() method first!")
        eccentricities = []

        for moment in self.moments:
            ecc = 1 - moment.min() / moment.mean()
            eccentricities.append(ecc)

        return eccentricities
