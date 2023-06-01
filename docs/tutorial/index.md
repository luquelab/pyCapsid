---
layout: default
title: Tutorial
nav_order: 3
---
* 
{:toc}

# Interactive Usage
The simplest way to follow this tutorial is using this [colab notebook](https://colab.research.google.com/drive/1p4E1ne8t47yGaiKC6NXpDo4TTnPLOrZ5?usp=sharing).
This tutorial covers the step by step use pyCapsid to identify the quasi-rigid subunits of an example capsid. This tutorial also comes
in the form of a jupyter notebook for those who wish to run it locally.
The example notebook is also provided for local use in the [notebooks folder](https://github.com/luquelab/pyCapsid/tree/main/notebooks).

Once the package and other dependencies are [installed](https://luquelab.github.io/pyCapsid/installation/), download the 
notebook and run the following command in its directory to launch jupyter notebook.

~~~~
jupyter notebook
~~~~

## Fetch and load PDB
This code acquires the pdb file from the RCSB databank, loads the necessary information, and saves copies for possible use in visualization in other software.

```python
from pyCapsid.PDB import getCapsid
pdb = '4oq8'
capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(pdb)
```

## Build ENM Hessian
This code builds a hessian matrix using an elastic network model defined by the given parameters. The types of model and the meaning of the parameters are provided in the documentation.

```python
from pyCapsid.CG import buildENMPreset
kirch, hessian = buildENMPreset(coords, preset='U-ENM')
```

## Perform NMA
This code calculates the n lowest frequency modes of the system by calculating the eigenvalues and eigenvectors of the hessian matrix.

```python
from pyCapsid.NMA import modeCalc
evals, evecs = modeCalc(hessian)
```

## Predict, fit, and compare b-factors
This code uses the resulting normal modes and frequencies to predict the b-factors of each alpha carbon, fits these results to experimental values from the pdb entry, and plots the results for comparison.

```python
from pyCapsid.NMA import fitCompareBfactors
evals_scaled, evecs_scaled = fitCompareBfactors(evals, evecs, bfactors, pdb)
```

![capsid_chx](4oq8_bfactorplot.png){: width="500"}

## Perform quasi-rigid cluster identification (QRC)

```python
from pyCapsid.NMA import calcDistFlucts
from pyCapsid.QRC import findQuasiRigidClusters

dist_flucts = calcDistFlucts(evals, evecs, coords)

cluster_start = 4
cluster_stop = 130
cluster_step = 2
labels, score, residue_scores  = findQuasiRigidClusters(pdb, dist_flucts, cluster_start=cluster_start, cluster_stop=cluster_stop, cluster_step=cluster_step)
```

![capsid_chx](4oq8_score_profile.png){: width="500"}

## Visualize in ChimeraX
If ChimeraX (https://www.cgl.ucsf.edu/chimerax/download.html) is installed you may provide a path to the chimerax 
executable file to automatically visualize the results in chimerax. This is done using the runscript command in chimerax 
and this python script: (https://github.com/luquelab/pyCapsid/blob/main/src/pyCapsid/scripts/chimerax_script.py).

```python
from pyCapsid.VIS import chimeraxViz
chimeraxViz(labels, pdb, chimerax_path='C:\\Program Files\\ChimeraX\\bin')
```

![capsid_chx](4oq8_chimerax.png){: width="500"}

Running the same code but replacing labels with residue_scores and adding rwb_scale=True visualizes the quality score of 
each residue. This is a measure of how rigid each residue is with respect to its cluster. Blue residues make up the 
cores of rigid clusters, and red residues represent borders between clusters. 

```python
from pyCapsid.VIS import chimeraxViz
chimeraxViz(residue_scores, pdb, chimerax_path='C:\\Program Files\\ChimeraX\\bin', rwb_scale=True)
```

![capsid_score_chx](4oq8_score_cx.png){: width="500"}

## Visualize in jupyter notebook with nglview
You can visualize the results in a jupyter notebook with nglview. The following function returns an nglview view with the 
results colored based on cluster. See the nglview documentation for further info 
(http://nglviewer.org/nglview/release/v2.7.7/index.html)
This cell will create an empty view, which the next cell will modify to create the final result.
```python
from pyCapsid.VIS import createCapsidView
view_clusters = createCapsidView(pdb, capsid)
view_clusters
```

![capsid_ngl](4oq8_nglview.png){: width="500"}

If the above view doesn't change coloration, run this cell again. In general do not run this cell until the above cell 
has finished rendering.

```python
from pyCapsid.VIS import createClusterRepresentation
createClusterRepresentation(pdb, labels, view_clusters)
```

Once you've done this use this code to download the results

```python
view_clusters.download_image()
```

![capsid_ngl](4oq8_score.png){: width="500"}


# Running pyCapsid using a simple config.toml file
This tutorial also has a corresponding [colab notebook](https://colab.research.google.com/drive/1Ct9Lh6w5qpO_9vGRXt7_9YHcUB9aJL8Z?usp=sharing).
This is a simpler and faster way to run the entire pyCapsid pipeline and save the results by setting the parameters ahead
of time in a text file. To do this download [this example](https://github.com/luquelab/pyCapsid/blob/main/docs/tutorial/config_simple.toml) 
from our github or copy and paste the following into a text editor and save the output as 'config.toml'

### A simple config.toml example

```toml
[PDB]
pdb = '4oq8' # PDB ID of structure
save_all_path = './4oq8/' # where to save the results

[CG]
preset = 'U-ENM' # Model Preset To Use
save_hessian = true # Whether to save the hessian matrix

[NMA]
n_modes = 200 # Number of low frequency modes to calculate
eigen_method = 'eigsh' # eigen method to use

[b_factors]
fit_modes = true # Whether to select the number of modes used to maximize correlation

[QRC]

[VIS]
chimerax_path = 'C:\Program Files\ChimeraX\bin\ChimeraX.exe'

[plotting]

```

Once you've created the 'config.toml' file, the following python code will run the entire pyCapsid pipeline using the
specified settings. Make sure to either run python in the same directory as the config file or properly include the
path to the file in the python code.

```python
from pyCapsid import run_capsid
run_capsid('config.toml')
```

### A more complex config.toml

```toml
[PDB]
pdb = '4oq8' # PDB ID of structure
pdbx = false
local = false
save_full_pdb = true
save_all_path = './4oq8/' # will be prepended to all other save_paths
save_pdb_path = 'pdb/'

[CG]
preset = 'U-ENM'
save_hessian = true
save_kirchhoff = true
save_cg_path = 'matrices/'

[NMA]
n_modes = 200
eigen_method = 'eigsh'
shift_invert = true
save_modes = true
save_mode_path = 'modes/'

[b_factors]
fit_modes = true
plot_modes = false
force_ico = true
ico_tol = 0.002
save_bfactors = true
save_bfactors_path = 'bfactors/'


[QRC]
cluster_start = 10
cluster_stop = 100
cluster_step = 2
cluster_method = 'discretize'
score_method = 'median'
return_type = 'final'
save_results =  true
save_results_path =  'results/'

[VIS]
chimerax_path = 'C:\Program Files\ChimeraX\bin\ChimeraX.exe'

[plotting]
suppress_plots = true

```

# ProDy Integration
One can make use of pyCapsids faster ENM and NMA module while still being able to use ProDy's other features by performing
the calculations using pyCapsid and passing the results to ProDy.

```python
from prody import ANM, parsePDB
capsid = parsePDB('7kq5', biomol=True)
calphas = capsid.select('protein and name CA')
anm = ANM('T. maritima ANM analysis')

from pyCapsid.CG import buildENM
coords = calphas.getCoords()
kirch, hessian = buildENM(coords, cutoff=10)

from pyCapsid.NMA import modeCalc
evals, evecs = modeCalc(hessian)

anm._hessian = hessian
anm._array = evecs
anm._eigvals = evals
```