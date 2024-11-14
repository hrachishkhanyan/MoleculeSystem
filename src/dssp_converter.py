import sys

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mode

matplotlib.rcParams.update({'font.size': 16})


def _convert_dat(path: str):
    s_s_letters: list[str] = ['=', 'P', 'S', 'T', 'G', 'H', 'I', 'B', 'E', '~']

    s_s_values: list[float] = [0.0, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    raw_data = np.loadtxt(path, dtype=str)

    for i, row in enumerate(raw_data):
        raw_data[i] = row.replace('=', '')

    shape = raw_data.shape + (len(raw_data[0]),)
    clean_data = np.zeros(shape=shape)
    s_s_existence: list[bool] = [False] * len(s_s_letters)

    for i, row in enumerate(raw_data):
        for j, symbol in enumerate(row):
            if symbol != '=':
                clean_data[i, j] = s_s_values[s_s_letters.index(symbol)]
                if not s_s_existence[s_s_letters.index(symbol)]:
                    s_s_existence[s_s_letters.index(symbol)] = True

    clean_data = clean_data[::10].T
    # Hardcoding section
    clean_data.shape = (116, 6, -1)
    reduced_array = mode(clean_data, axis=1).mode.squeeze()
    print(reduced_array.shape)
    return reduced_array, s_s_existence


def dat2png(path: str, output_file: str):
    s_s_descriptions: list[str] = ['Break', 'κ-Helix', 'Bend', 'Turn', '3₁₀-Helix', 'α-Helix', 'π-Helix', 'β-Bridge',
                                   'β-Strand', 'Loop']
    s_s_values: list[float] = [0.0, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    clean_data, s_s_existence = _convert_dat(path)

    fig, ax0 = plt.subplots(1, 1, figsize=(10, 6))
    c = ax0.pcolor(clean_data, linewidths=0, edgecolors='k', cmap='nipy_spectral', vmin=0.0, vmax=1.0)
    ax0.set_xlabel('Time (ns)')
    ax0.set_ylabel('Protein residue number')
    colors: list = [c.cmap(c.norm(value)) for value in s_s_values]
    patches = [mpatches.Patch(color=colors[i], label=f'{s_s_descriptions[i]}') for i in range(len(colors)) if
               s_s_existence[i]]
    ax0.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    salt = 'urea_kh2po4'
    dat2png(fr'E:\Simulation\spider_silk\oligomer\{salt}\dssp.dat',
            fr'C:\Users\hrach\OneDrive\Pictures\Graphs\Spider silk\oligo_{salt}_dssp.png')
