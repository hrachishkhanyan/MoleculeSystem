# MoleculeSystem

A Python toolkit for molecular dynamics (MD) simulation analysis, focusing on conformational clustering using dimensionality reduction techniques.

## Overview

This project provides tools for analyzing molecular dynamics trajectories, with functionality for:

- **Conformational clustering** using UMAP dimensionality reduction and HDBSCAN clustering
- **Distance-based analysis** of molecular conformations
- **Micelle processing** utilities

## Project Structure

```
MoleculeSystem-main/
├── UMAP-and_HDBSCAB_dist.ipynb    # Main clustering analysis notebook
├── src/
│   └── make_micelle_whole.ipynb   # Micelle PBC unwrapping utilities
└── umap_files/                     # Output directory for UMAP results
```

## Dependencies

- [MDAnalysis](https://www.mdanalysis.org/) - MD trajectory analysis
- [UMAP](https://umap-learn.readthedocs.io/) - Dimensionality reduction
- [HDBSCAN](https://hdbscan.readthedocs.io/) - Density-based clustering
- [NumPy](https://numpy.org/) - Numerical computations
- [Matplotlib](https://matplotlib.org/) - Visualization
- [scikit-learn](https://scikit-learn.org/) - t-SNE and ML utilities
- [joblib](https://joblib.readthedocs.io/) - Model persistence

## Installation

```bash
pip install MDAnalysis umap-learn hdbscan numpy matplotlib scikit-learn joblib
```

## Usage

### Load Trajectory Data

```python
import MDAnalysis as mda

structure = 'path/to/structure.gro'
trajectory = 'path/to/trajectory.xtc'
u = mda.Universe(structure, trajectory)
```

### Extract Distances and Run UMAP

```python
from umap import UMAP

# Select reference atoms
atom_names = """name C7 C36 C70 C102"""
ag = u.select_atoms(atom_names)

# Compute distances and run UMAP
reducer = UMAP(n_components=2, n_neighbors=24, random_state=42)
embedded = reducer.fit_transform(distances)
```

### Cluster with HDBSCAN

```python
from hdbscan import HDBSCAN

clusterer = HDBSCAN(min_cluster_size=50)
labels = clusterer.fit_predict(embedded)
```

## Output Files

- `embedded_distances.npy` - 2D UMAP embeddings
- `cluster_labels.npy` - Cluster assignments per frame
- `2d_embedding.joblib` - Saved UMAP model
- `clusters.joblib` - Saved HDBSCAN model
- `clusters.png` - Cluster visualization

## License

MIT License
