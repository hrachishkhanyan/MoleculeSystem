import pickle

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
from os.path import join

from src.standard_params import StandardParams

mpl.rcParams['font.size'] = 14
TQDM_BAR_FORMAT = '\033[37m{percentage:3.0f}%|{bar:30}\033[37m|[Estimated time remaining: {remaining}]\033[0m'

def meh():
    font = {
        'size': 14}
    matplotlib.rc('font', **font)
    surf = 'tyl_3'
    BASE_PATH = rf"E:\Simulation\tyloxapol_tx\{surf}"
    SYSTEM_NAME = "75_25"

    top = join(BASE_PATH, SYSTEM_NAME, "centered.gro")
    trj = join(BASE_PATH, SYSTEM_NAME, "centered_whole_skip_100.xtc")

    m = StandardParams(top, trj)
    m.select_atoms('(resname TY39 or resname TX0) and not type H')
    m.sasa(420, 600)
    np.save(f'sasa_{surf}_{SYSTEM_NAME}.npy', m.get_sasa_data())
    print(m.get_sasa_data().mean(), m.get_sasa_data().std())
    # np.save(f'data/{surf}_{SYSTEM_NAME}_ellipticity.npy', m.moments)
    # m = StandardParams(top, trj)
    # m.select_atoms("resname OLE or resname CTB16")
    # m.sasa()
    # np.save('sasa_naol_ctab.npy', m.get_sasa_data())
    # systems = [
    #     '100',
    #     '75_25',
    #     '50_50',
    #     '25_75'
    # ]
    # data = []
    # for s in systems:
    #     data.append(np.array([max(d) / min(d) for d in np.load(f'data/{surf}_{s}_ellipticity.npy')]))
    #
    # ellipt = [max(d) / min(d) for d in data]
    # plt.hist(ellipt)
    # plt.figure(figsize=(10, 6))
    # for i, d in enumerate(data):
    #     if i == 0:
    #         plt.plot(d, label=f'{systems[i]}% TYL7', linewidth=2)
    #     else:
    #         tyl_pc, tx_pc = systems[i].split('_')
    #         plt.plot(d, label=f'{tyl_pc}% TYL7 - {tx_pc}% TX100', linewidth=2)
    # #     plt.hist(d[420:], rwidth=.9, alpha=.9)
    # plt.legend()
    # plt.xlim(0, 600)
    # plt.ylim(1, 5.5)
    # plt.legend()
    # plt.xlabel('Time (ns)')
    # plt.ylabel('Ellipticity')
    # plt.savefig(f'{BASE_PATH}/TYL7_ellipt.png')
    # plt.show()

str
def main():
    BASE_PATH = r"C:\Users\hrach\Documents\Simulations\naol_micelles"
    SYSTEM_NAME = "200-ctab-micelle"
    top = join(BASE_PATH, SYSTEM_NAME, "clustered.gro")
    trj = join(BASE_PATH, SYSTEM_NAME, "clustered_skip5.xtc")
    start, skip, end = 0, 1, 1000
    s = StandardParams(top, trj)
    s.select_atoms('resname CTB16 and not type H')
    s.sasa(start, skip, end)
    plt.plot(s.get_sasa_data())
    plt.show()
    # with open('ctab_micelle_clusters_5000frames.pkl', 'wb') as file:
    #     pickle.dump(s.clusters, file)
    # r = s.rms(start, skip, end)
    # print(r)


if __name__ == '__main__':
    main()
