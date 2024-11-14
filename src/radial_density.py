import numpy as np
from matplotlib import pyplot as plt

from src.MoleculeSystem import MoleculeSystem
import MDAnalysis as mda


class RadialDensity(MoleculeSystem):

    def __init__(self, top_path: str, trj_path=None):
        self.ag = None
        if not trj_path:
            self.u = mda.Universe(top_path)
        else:
            self.u = mda.Universe(top_path, trj_path)

    def select_atoms(self, selection: str) -> None:
        self.ag = self.u.select_atoms(selection)

    def calculate_radial_density(self, selection, ref, layer_count=15):
        pbc = self.u.dimensions[0]
        density_ag = self.u.select_atoms(selection)

        densities = np.zeros(layer_count)
        distances = np.zeros(layer_count)

        ls, step = np.linspace(0, pbc, layer_count, retstep=True)

        for index, i in enumerate(ls):
            layer_volume = 3 / 4 * np.pi * ((i + step) ** 3 - i ** 3)
            n_layer_atoms = density_ag.select_atoms(f'sphlayer {i} {i + step} ( {ref} )').n_atoms
            densities[index] = (n_layer_atoms / layer_volume)
            distances[index] = i

        return distances, densities


if __name__ == '__main__':
    import matplotlib
    matplotlib.rcParams.update({'font.size': 14})
    layer_count = 20
    rd = RadialDensity(
        r'C:\Users\hrach\PycharmProjects\md_grid_project\PUCHIK\test\test_structures\InP_sphere_r_29.pdb')
    distances, densities = rd.calculate_radial_density('resname UNL and type P', ref='resname UNL',
                                                       layer_count=layer_count)
    rd_1 = RadialDensity(r'C:\Users\hrach\PycharmProjects\md_grid_project\PUCHIK\test\test_structures\InP_cylinder.pdb')
    distances_1, densities_1 = rd_1.calculate_radial_density('resname UNL and type P', ref='resname UNL',
                                                             layer_count=layer_count)

    plt.plot(distances - 27, densities, label='Sphere')
    plt.plot(distances_1 - 27, densities_1, label='Cylinder')
    plt.legend()
    plt.ylabel('Number density (#/$\AA^3$)')
    plt.xlabel('Distance ($\AA$)')
    plt.xlim(-20, 20)
    plt.ylim(0, 0.04)

    plt.savefig(r'C:\Users\hrach\OneDrive\Pictures\Graphs\PUCHIK\test_rad_density.png', bbox_inches='tight')
    plt.show()

