# Source Module Documentation

This directory contains Python modules for molecular dynamics simulation analysis. All analysis classes inherit from the abstract `MoleculeSystem` base class.

---

## Table of Contents

1. [MoleculeSystem.py](#moleculesystempy) - Abstract base class
2. [standard_params.py](#standard_paramspy) - SASA, RMSD, radius of gyration
3. [dfs_clustering.py](#dfs_clusteringpy) - DFS-based molecular clustering
4. [graph_clustering.py](#graph_clusteringpy) - Graph-based clustering with NetworkX
5. [inertia_moment.py](#inertia_momentpy) - Moment of inertia and shape analysis
6. [radial_density.py](#radial_densitypy) - Radial density profiles
7. [dssp_converter.py](#dssp_converterpy) - DSSP secondary structure visualization
8. [strip_water.py](#strip_waterpy) - Solvent removal utility

---

## MoleculeSystem.py

### Description
Abstract base class that defines the interface for all analysis modules. All other analysis classes must implement `__init__` and `select_atoms` methods.

### Class: `MoleculeSystem`

```python
from abc import ABC, abstractmethod

class MoleculeSystem(ABC):
    @abstractmethod
    def __init__(self, top_path: str, trj_path=None):
        ...

    @abstractmethod
    def select_atoms(self, selection: str) -> None:
        ...
```

### Usage
This class should not be instantiated directly. Inherit from it to create new analysis modules.

---

## standard_params.py

### Description
Calculates standard structural parameters including Solvent Accessible Surface Area (SASA), RMSD, and radius of gyration. Supports multiprocessing for efficient trajectory analysis.

### Class: `StandardParams`

#### Methods

| Method | Description |
|--------|-------------|
| `select_atoms(sel)` | Select atoms using MDAnalysis selection syntax |
| `sasa(start, skip, end, cpu_count)` | Calculate SASA over trajectory frames |
| `gyr()` | Calculate radius of gyration over trajectory |
| `rms(resids, start, skip, end, cpu_count)` | Calculate RMSD for specified residues |
| `cluster_molecules(start, skip, end, cpu_count)` | Identify molecular clusters using DFS |
| `get_sasa_data()` | Returns SASA data array |
| `get_gyr_data()` | Returns gyration radius timeseries |
| `get_rms_data()` | Returns RMSD data |

#### Example

```python
from src.standard_params import StandardParams
from os.path import join

BASE_PATH = r"path/to/simulations"
SYSTEM_NAME = "my_system"

top = join(BASE_PATH, SYSTEM_NAME, "structure.gro")
trj = join(BASE_PATH, SYSTEM_NAME, "trajectory.xtc")

# Initialize and select atoms
params = StandardParams(top, trj)
params.select_atoms('resname MOL and not type H')

# Calculate SASA (frames 0-1000, skip every frame)
params.sasa(start=0, skip=1, end=1000)
sasa_values = params.get_sasa_data()

# Calculate radius of gyration
params.gyr()
gyr_values = params.get_gyr_data()

# Find molecular clusters
params.cluster_molecules(start=0, skip=10, end=1000)
```

#### Dependencies
- MDAnalysis
- FreeSASA
- NumPy
- tqdm

---

## dfs_clustering.py

### Description
Identifies molecular clusters in a trajectory using depth-first search (DFS) algorithm. Molecules within a cutoff distance (default 4.7 Å) are considered part of the same cluster.

### Class: `ClusterSearch`

#### Attributes
- `cutoff` (module level): Distance cutoff for clustering (default: 4.7 Å)

#### Methods

| Method | Description |
|--------|-------------|
| `select_atoms(selection)` | Select atoms using MDAnalysis selection syntax |
| `find_clusters()` | Run DFS clustering over all trajectory frames |

#### Example

```python
from src.dfs_clustering import ClusterSearch
from os.path import join
import pickle

BASE_PATH = r"path/to/simulations"
SYSTEM_NAME = "micelle_system"

top = join(BASE_PATH, SYSTEM_NAME, "structure.gro")
trj = join(BASE_PATH, SYSTEM_NAME, "trajectory.xtc")

# Initialize and select atoms
searcher = ClusterSearch(top, trj)
searcher.select_atoms('resname CTB16 and not type H')

# Find clusters for each frame
clusters = searcher.find_clusters()

# Save results
with open('clusters.pkl', 'wb') as f:
    pickle.dump(clusters, f)

# Output format: list of lists
# clusters[frame] = [cluster1_resids, cluster2_resids, ...]
```

---

## graph_clustering.py

### Description
Graph-based molecular clustering using NetworkX. Creates an adjacency matrix based on inter-molecular distances and identifies connected components as clusters.

### Class: `SystemGraph`

#### Methods

| Method | Description |
|--------|-------------|
| `select_atoms(selection)` | Select atoms using MDAnalysis selection syntax |
| `graph()` | Build molecular graph and identify clusters |

#### Example

```python
from src.graph_clustering import SystemGraph
from os.path import join

BASE_PATH = r"path/to/simulations"
SYSTEM_NAME = "micelle_system"

top = join(BASE_PATH, SYSTEM_NAME, "structure.gro")

# Initialize (single frame analysis)
graph = SystemGraph(top)
graph.select_atoms('resname CTB16')

# Build graph and visualize clusters
graph.graph()
# Outputs: Cluster assignments and NetworkX visualization
```

#### Output
- Prints cluster assignments to console
- Displays NetworkX spring layout visualization
- Saves adjacency matrix to `temp.txt`

---

## inertia_moment.py

### Description
Calculates the moment of inertia tensor eigenvalues to determine molecular/aggregate shape and eccentricity over trajectory frames.

### Class: `Inertia`

#### Shape Classifications
| Code | Shape |
|------|-------|
| 0 | Spherical |
| 1 | Oblate (disk-like) |
| 2 | Prolate (rod-like) |
| 3 | Triaxial |

#### Methods

| Method | Description |
|--------|-------------|
| `select_atoms(selection)` | Select atoms using MDAnalysis selection syntax |
| `run(start, end, skip)` | Calculate moments of inertia over frames |
| `determine_shape(tolerance)` | Classify shape based on eigenvalue ratios |
| `eccentricity()` | Calculate eccentricity for each frame |

#### Example

```python
from src.inertia_moment import Inertia
from os.path import join

BASE_PATH = r"path/to/simulations"
SYSTEM_NAME = "micelle_system"

top = join(BASE_PATH, SYSTEM_NAME, "structure.gro")
trj = join(BASE_PATH, SYSTEM_NAME, "trajectory.xtc")

# Initialize
inertia = Inertia(top, trj)
inertia.select_atoms('resname MOL')

# Run analysis over frames 0-1000
inertia.run(start=0, end=1000, skip=1)

# Get shape classification
shape = inertia.determine_shape(tolerance=None)  # Auto-calculates tolerance
print(f"Aggregate shape: {shape}")

# Calculate eccentricity time series
ecc = inertia.eccentricity()
print(f"Mean eccentricity: {np.mean(ecc):.3f}")
```

#### Eccentricity Formula
$$e = 1 - \frac{I_{min}}{I_{mean}}$$

Where $I_{min}$ is the smallest eigenvalue and $I_{mean}$ is the mean of all eigenvalues.

---

## radial_density.py

### Description
Calculates radial number density profiles from a reference point or selection. Useful for analyzing density distributions in spherical or micellar systems.

### Class: `RadialDensity`

#### Methods

| Method | Description |
|--------|-------------|
| `select_atoms(selection)` | Select atoms using MDAnalysis selection syntax |
| `calculate_radial_density(selection, ref, layer_count)` | Calculate density in spherical shells |

#### Parameters
- `selection`: Atoms to count in each shell
- `ref`: Reference selection for center
- `layer_count`: Number of spherical shells (default: 15)

#### Example

```python
from src.radial_density import RadialDensity
import matplotlib.pyplot as plt

# Initialize
rd = RadialDensity('structure.gro', 'trajectory.xtc')

# Calculate radial density of phosphorus atoms
distances, densities = rd.calculate_radial_density(
    selection='resname LIPID and type P',
    ref='resname CORE',
    layer_count=20
)

# Plot results
plt.plot(distances, densities)
plt.xlabel('Distance (Å)')
plt.ylabel('Number density (#/Å³)')
plt.savefig('radial_density.png')
plt.show()
```

#### Shell Volume Calculation
$$V_{shell} = \frac{4}{3}\pi\left((r + \Delta r)^3 - r^3\right)$$

---

## dssp_converter.py

### Description
Converts DSSP secondary structure output files (`.dat`) to visualization plots. Maps secondary structure symbols to colored representations.

### Secondary Structure Codes
| Symbol | Value | Description |
|--------|-------|-------------|
| `=` | 0.0 | Break |
| `P` | 0.1 | κ-Helix |
| `S` | 0.2 | Bend |
| `T` | 0.3 | Turn |
| `G` | 0.5 | 3₁₀-Helix |
| `H` | 0.6 | α-Helix |
| `I` | 0.7 | π-Helix |
| `B` | 0.8 | β-Bridge |
| `E` | 0.9 | β-Strand |
| `~` | 1.0 | Loop |

### Functions

| Function | Description |
|----------|-------------|
| `dat2png(path, output_file)` | Convert DSSP .dat file to PNG visualization |

#### Example

```python
from src.dssp_converter import dat2png

# Convert DSSP output to image
dat2png(
    path='path/to/dssp.dat',
    output_file='secondary_structure.png'
)
```

#### Output
Creates a 2D heatmap showing:
- X-axis: Time (frames/ns)
- Y-axis: Residue number
- Color: Secondary structure type

---

## strip_water.py

### Description
Utility function to remove water (or other residues) from trajectories and save the stripped system.

### Function: `strip_water`

#### Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `top` | str | required | Path to topology file |
| `trj` | str | None | Path to trajectory file |
| `file_name` | str | "no_water" | Output filename (without extension) |
| `save_path` | str | "." | Directory for output files |
| `water_resname` | str | "TIP3" | Residue name to remove |

#### Example

```python
from src.strip_water import strip_water

# Remove TIP3 water from trajectory
strip_water(
    top='solvated.gro',
    trj='trajectory.xtc',
    file_name='no_solvent',
    save_path='./output',
    water_resname='TIP3'
)
# Creates: ./output/no_solvent.gro and ./output/no_solvent.xtc

# Remove ions (e.g., sodium)
strip_water(
    top='system.gro',
    trj='trajectory.xtc',
    file_name='no_ions',
    water_resname='NA'
)
```

#### Output Files
- `{file_name}.gro` - Topology without specified residue
- `{file_name}.xtc` - Trajectory without specified residue (if `trj` provided)

---

## Dependencies Summary

| Module | Required Packages |
|--------|-------------------|
| MoleculeSystem | abc (built-in) |
| standard_params | MDAnalysis, FreeSASA, NumPy, SciPy, tqdm |
| dfs_clustering | MDAnalysis, NumPy |
| graph_clustering | MDAnalysis, NumPy, NetworkX, Matplotlib |
| inertia_moment | MDAnalysis, NumPy |
| radial_density | MDAnalysis, NumPy, Matplotlib |
| dssp_converter | NumPy, Matplotlib, SciPy |
| strip_water | MDAnalysis |

## Installation

```bash
pip install MDAnalysis numpy matplotlib scipy networkx freesasa tqdm
```
