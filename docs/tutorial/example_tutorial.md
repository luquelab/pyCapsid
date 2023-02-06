---
layout: default
title: Tutorial
nav_order: 4
---

# Example Tutorial
This tutorial covers an example of the process to use PyCapid.


## Fetch and load PDB
This code acquires the pdb file from the RCSB databank, loads the necessary information, and saves copies for possible use in visualization in other software.

```python
from pyCapsid.PDB import getCapsid
pdb = '4oq8'
capsid, calphas, coords, bfactors, chain_starts, title = getCapsid(pdb, save=True)

```

## Build ENM Hessian
This code builds a hessian matrix using an elastic network model defined by the given parameters. The types of model and the meaning of the parameters are provided in the documentation.

```python
from pyCapsid.ENM import buildENMPreset
kirch, hessian = buildENMPreset(coords, preset='U-ENM')
```

## Perform NMA
### Calculate low frequency modes
This code calculates the n lowest frequency modes of the system by calculating the eigenvalues and eigenvectors of the hessian matrix.

```python
from pyCapsid.NMA import modeCalc
n_modes = 200
eigmethod = 'eigsh'

evals, evecs = modeCalc(hessian, n_modes, eigmethod=eigmethod)
```

## Predict, fit, and compare b-factors
This code uses the resulting normal modes and frequencies to predict the b-factors of each alpha carbon, fits these results to experimental values from the pdb entry, and plots the results for comparison.

```python
from pyCapsid.bfactorfit import plotBfactors
plotBfactors(evals, evecs, bfactors, pdb, is3d=True, fitModes=True, plotModes=True, forceIc
```

## Perform quasi-rigid cluster identification (QRC)
### Build weighted graph based on distance fluctuations
This code calculates the distance fluctuations between residues within a cutoff distance from each other and transforms those distance fluctuations into a similarity matrix representing a sparse weighted graph.

```python
from pyCapsid.NMA import calcDistFlucts

n_modes = 200
fluct_cutoff = 10

dist_flucts = calcDistFlucts(evals, evecs, n_modes, coords, fluct_cutoff, is3d=True)
print(dist_flucts.data)
```

## Calculate the spectral embedding of the graph
This code calculates the spectral embedding (eigenvectors) of the sparse weighted graph.

```python
from pyCapsid.QRC import calcEmbedding, cluster_embedding
from pyCapsid import clustering_util
import numpy as np

n_cluster_max = 130
embedding = calcEmbedding(sims, n_cluster_max)
```

## Cluster the embedded points
This code performs clustering for a set of # of clusters

```python
n_range = np.arange(4, n_cluster_max, 2)

labels, scores, variances, numtypes = cluster_embedding(n_range, embedding, method='discretize')
```


```python
from pyCapsid.clustering_util import plotScores
plotScores(pdb, n_range, scores, variances, numtypes)
clusters = labels[np.argmax(scores)]
```




