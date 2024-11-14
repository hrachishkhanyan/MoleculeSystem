import multiprocessing
import time
from functools import partial

import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.base import AnalysisFromFunction
import numpy as np
from MDAnalysis.lib.distances import distance_array
from freesasa import calcCoord
from MDAnalysis.topology.tables import vdwradii
import os
from scipy.signal import savgol_filter
from tqdm.contrib.concurrent import process_map

from src.MoleculeSystem import MoleculeSystem

TQDM_BAR_FORMAT = '\033[37m{percentage:3.0f}%|{bar:30}\033[37m|[Estimated time remaining: {remaining}]\033[0m'
CPU_COUNT = multiprocessing.cpu_count()


class StandardParams(MoleculeSystem):
    clustering_cutoff = 4.7

    def __init__(self, top_path: str, trj_path=None):
        self.clusters = None
        self.cluster_count = 1
        self.file_name = os.path.splitext(os.path.split(trj_path)[1])[0]
        self.u = mda.Universe(top_path, trj_path)
        self.ag = None
        self.rms_data = None
        self.selection = None
        self.sasa_data = []  # np.zeros((len(self.u.trajectory, )))
        self.gyr_data = None

    def get_sasa_data(self):
        return self.sasa_data

    def get_gyr_data(self):
        return self.gyr_data.timeseries[:, 0]

    def get_rms_data(self):
        return self.rms_data

    def select_atoms(self, sel):
        self.ag = self.u.select_atoms(sel)
        self.selection = sel

    def rms(self, resids, start=0, skip=1, end=None, cpu_count=CPU_COUNT):
        reference_pos = self.ag.select_atoms(f'resid  {" ".join(map(str, resids))}').positions

        n_frames = self.u.trajectory.n_frames if end is None else end
        frame_range = range(start, n_frames, skip)

        rms_func = partial(self._rms, reference_pos=reference_pos)
        res = process_map(self._rms, frame_range,
                          max_workers=cpu_count,
                          bar_format=TQDM_BAR_FORMAT
                          )
        return res

    def _rms(self, frame, reference_pos):
        self.u.trajectory[frame]

        rms_result = {i: [] for i in range(len(self.clusters[frame]))}

        for i, cluster in enumerate(rms_result):
            pos = self.ag.select_atoms(f'resid  {" ".join(map(str, self.clusters[frame][i]))}').positions
            print(rmsd(pos, pos))
            rms_result[cluster].append(rmsd(pos, reference_pos))

        return rms_result

    def _radgyr(self, atomgroup, masses, total_mass=None):
        # coordinates change for each frame
        coordinates = atomgroup.positions
        center_of_mass = atomgroup.center_of_mass()

        # get squared distance from center
        ri_sq = (coordinates - center_of_mass) ** 2

        # sum the unweighted positions
        sq = np.sum(ri_sq, axis=1)
        sq_x = np.sum(ri_sq[:, [1, 2]], axis=1)  # sum over y and z
        sq_y = np.sum(ri_sq[:, [0, 2]], axis=1)  # sum over x and z
        sq_z = np.sum(ri_sq[:, [0, 1]], axis=1)  # sum over x and y

        # make into array
        sq_rs = np.array([sq, sq_x, sq_y, sq_z])

        # weight positions
        rog_sq = np.sum(masses * sq_rs, axis=1) / total_mass
        # square root and return
        return np.sqrt(rog_sq)

    def gyr(self):
        print('Calculating radius of gyration...')
        gyradius = AnalysisFromFunction(self._radgyr, self.ag, self.ag.masses, total_mass=np.sum(self.ag.masses))
        gyradius.run()
        print('Done!')
        self.gyr_data = gyradius.results

        return self

    def sasa(self, start=0, skip=1, end=None, cpu_count=CPU_COUNT):
        print('Calculating SASA...')
        radii = np.array([vdwradii[ag.type] for ag in self.ag])

        n_frames = self.u.trajectory.n_frames if end is None else end
        frame_range = range(start, n_frames, skip)
        partial_sasa = partial(self._sasa, radii=radii)

        self.sasa_data = self._mp_calc(partial_sasa, frame_range, cpu_count)
        print('Done!')

    def _sasa(self, frame, radii):
        self.u.trajectory[frame]

        coords = self.ag.positions.flatten()
        sasa_res = calcCoord(coords, radii)
        return sasa_res.totalArea()

    def _smooth(self, data, windows=5, order=2):
        return savgol_filter(data, windows, order)

    def _cluster(self, frame):
        self.u.trajectory[frame]

        n_res = self.ag.n_residues
        resids = np.unique(self.ag.resids)
        n_atoms = self.ag.n_atoms // n_res
        clusters = []

        positions = np.zeros((n_res, n_atoms, 3))  # positions of carbons in rings
        checked_residues = []

        for i in range(n_res):
            for j in range(n_atoms):
                positions[i, j] = self.ag[i * n_atoms + j].position  # WRONG!!!  ????

        for i, position in enumerate(positions):
            if i not in checked_residues:
                clusters.append(resids[self._dfs(i, positions, [i], checked_residues)])

        if len(clusters) > self.cluster_count:
            self.cluster_count = len(clusters)

        return clusters

    def _dfs(self, index, positions: np.ndarray, cluster, checked) -> list:
        checked.append(index)

        for i in range(len(positions)):

            if i != index and i not in checked:
                distances = distance_array(positions[index], positions[i])

                if (distances < self.clustering_cutoff).any():
                    cluster.append(i)
                    self._dfs(i, positions, cluster, checked)

        return cluster

    def cluster_molecules(self, start=0, skip=1, end=None, cpu_count=multiprocessing.cpu_count()):
        n_frames = self.u.trajectory.n_frames if end is None else end
        frame_range = range(start, n_frames, skip)

        self.clusters = self._mp_calc(self._cluster, frame_range, cpu_count)

    @staticmethod
    def _mp_calc(func, frame_range, cpu_count, **kwargs) -> np.ndarray:
        """
        This method handles multiprocessing
        :param func:
        :param frame_range:
        :param cpu_count:
        :param kwargs:
        :return:
        """
        start = time.perf_counter()
        per_frame_func = partial(func, **kwargs)
        res = process_map(per_frame_func, frame_range,
                          max_workers=cpu_count,
                          bar_format=TQDM_BAR_FORMAT
                          )
        print(f'Execution time for {len(frame_range)} frames:', time.perf_counter() - start)
        return np.array(res)
