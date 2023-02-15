---
layout: default
title: Tutorial
nav_order: 3
---

# Tutorial: Quasi-rigid subunits of Satellite Tobacco Mosaic Virus

This tutorial covers the use pyCapsid to identify the quasi-rigid subunits of an example capsid. An example 
notebook is provided with all of the example code in the [notebooks folder](https://github.com/luquelab/pyCapsid/tree/main/notebooks). 

Once the package and other dependecies are [installed](https://luquelab.github.io/pyCapsid/installation/), download the 
notebook and run the following command in its directory to launch jupyter lab.

~~~~
jupyter lab
~~~~

### Fetch and load PDB
This code acquires the pdb file from the RCSB databank, loads the necessary information, and saves copies for possible use in visualization in other software.

```python
from pyCapsid.PDB import getCapsid
pdb = '4oq8'
capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(pdb, save=True)
```

### Build ENM Hessian
This code builds a hessian matrix using an elastic network model defined by the given parameters. The types of model and the meaning of the parameters are provided in the documentation.

```python
from pyCapsid.ENM import buildENMPreset
kirch, hessian = buildENMPreset(coords, preset='U-ENM')
```

### Perform NMA
This code calculates the n lowest frequency modes of the system by calculating the eigenvalues and eigenvectors of the hessian matrix.

```python
from pyCapsid.NMA import modeCalc
evals, evecs = modeCalc(hessian)
```

### Predict, fit, and compare b-factors
This code uses the resulting normal modes and frequencies to predict the b-factors of each alpha carbon, fits these results to experimental values from the pdb entry, and plots the results for comparison.

```python
from pyCapsid.NMA import fitCompareBfactors
evals_scaled, evecs_scaled = fitCompareBfactors(evals, evecs, bfactors, pdb, fitModes=False)
```

### Perform quasi-rigid cluster identification (QRC)

```python
from pyCapsid.NMA import calcDistFlucts
from pyCapsid.QRC import findQuasiRigidClusters
import numpy as np

dist_flucts = calcDistFlucts(evals_scaled, evecs, coords)

n_cluster_max = 62
n_range = np.arange(4, n_cluster_max, 2)
labels, score  = findQuasiRigidClusters(pdb, dist_flucts, n_range)
```

## Visualize in ChimeraX
If ChimeraX (https://www.cgl.ucsf.edu/chimerax/download.html) is installed you may provide a path to the chimerax 
executable file and the path to the script (https://github.com/luquelab/pyCapsid/blob/main/src/pyCapsid/scripts/chimerax_script.py) 
to automatically visualize the results in chimerax.

```python
from pyCapsid.viz_util import chimeraxViz
chimeraxViz(labels, pdb, chimerax_path='C:\\Program Files\\ChimeraX\\bin')
```

![capsid_chx](4oq8_chimerax.png){: width="500"}

## Visualize in jupyter notebook with nglview
You can visualize the results in a jupyter notebook with nglview. The following function returns an nglview view with the 
results colored based on cluster. See the nglview documentation for further info 
(http://nglviewer.org/nglview/release/v2.7.7/index.html)

```python
from pyCapsid.viz_util import view_pdb_ngl
view = view_pdb_ngl(pdb, capsid, labels)
view
```

![capsid_ngl](4oq8_nglview.png){: width="500"}