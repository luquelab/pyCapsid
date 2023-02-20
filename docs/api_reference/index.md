---
layout: default
title: API Reference
nav_order: 5
---

# pycapsid


* [pyCapsid package]()


    * [pyCapsid.ENM module](#module-pyCapsid.ENM)


        * [`buildBackbone()`](#pyCapsid.ENM.buildBackbone)


        * [`buildENM()`](#pyCapsid.ENM.buildENM)


        * [`buildENMPreset()`](#pyCapsid.ENM.buildENMPreset)


        * [`hessCalc()`](#pyCapsid.ENM.hessCalc)


        * [`kirchGamma()`](#pyCapsid.ENM.kirchGamma)


    * [pyCapsid.NMA module](#module-pyCapsid.NMA)


        * [`calcCovMat()`](#pyCapsid.NMA.calcCovMat)


        * [`calcDistFlucts()`](#pyCapsid.NMA.calcDistFlucts)


        * [`con_c()`](#pyCapsid.NMA.con_c)


        * [`cov()`](#pyCapsid.NMA.cov)


        * [`distFluctFromCov()`](#pyCapsid.NMA.distFluctFromCov)


        * [`fitCompareBfactors()`](#pyCapsid.NMA.fitCompareBfactors)


        * [`fluctPlot()`](#pyCapsid.NMA.fluctPlot)


        * [`gCon_c()`](#pyCapsid.NMA.gCon_c)


        * [`gCov()`](#pyCapsid.NMA.gCov)


        * [`modeCalc()`](#pyCapsid.NMA.modeCalc)


    * [pyCapsid.PDB module](#module-pyCapsid.PDB)


        * [`downloadPDB()`](#pyCapsid.PDB.downloadPDB)


        * [`getCapsid()`](#pyCapsid.PDB.getCapsid)


        * [`getProdyChainStarts()`](#pyCapsid.PDB.getProdyChainStarts)


        * [`loadPDB()`](#pyCapsid.PDB.loadPDB)


        * [`loadPDBx()`](#pyCapsid.PDB.loadPDBx)


    * [pyCapsid.QRC module](#module-pyCapsid.QRC)


        * [`calcEmbedding()`](#pyCapsid.QRC.calcEmbedding)


        * [`cluster_embedding()`](#pyCapsid.QRC.cluster_embedding)


        * [`discretize()`](#pyCapsid.QRC.discretize)


        * [`findQuasiRigidClusters()`](#pyCapsid.QRC.findQuasiRigidClusters)


        * [`fluctToSims()`](#pyCapsid.QRC.fluctToSims)


    * [pyCapsid.bfactorfit module](#module-pyCapsid.bfactorfit)


        * [`checkIcoFlucts()`](#pyCapsid.bfactorfit.checkIcoFlucts)


        * [`fastFlucts()`](#pyCapsid.bfactorfit.fastFlucts)


        * [`fitBfactors()`](#pyCapsid.bfactorfit.fitBfactors)


        * [`fitPlotBfactors()`](#pyCapsid.bfactorfit.fitPlotBfactors)


        * [`fluctModes()`](#pyCapsid.bfactorfit.fluctModes)


        * [`plotByMode()`](#pyCapsid.bfactorfit.plotByMode)


        * [`springFit()`](#pyCapsid.bfactorfit.springFit)


    * [pyCapsid.clustering_util module](#module-pyCapsid.clustering_util)


        * [`calcCentroids()`](#pyCapsid.clustering_util.calcCentroids)


        * [`calcCosCentroids()`](#pyCapsid.clustering_util.calcCosCentroids)


        * [`cluster_types()`](#pyCapsid.clustering_util.cluster_types)


        * [`median_score()`](#pyCapsid.clustering_util.median_score)


        * [`plotScores()`](#pyCapsid.clustering_util.plotScores)


## pyCapsid.PDB module

Module with functions for downloading and dealing with PDB/PDBx files.


### pyCapsid.PDB.downloadPDB(pdb, dir='.', pdbx=False)
Downloads pdb and returns the filename


* **Parameters**

    
    * **pdb** – PDB id of entry to download. Can also be the name of a local file


    * **dir** – Target directory where download will be placed


    * **pdbx** – Whether the target structure should be acquired in pdbx/mmcif format



### pyCapsid.PDB.getCapsid(pdb, dir='.', pdbx=False, local=False, save=False, chains='', chains_clust='')
Downloads and opens molecular data from a PDB entry or loads data from a local file.


* **Parameters**

    
    * **pdb** – PDB id of entry to download. Can also be the name of a local file


    * **dir** – Target directory where download will be placed


    * **pdbx** – Whether the target structure should be acquired in pdbx/mmcif format


    * **local** – Whether to instead load a local file


    * **save** – Whether to save a copy of the complete assembly as pdb/pdbx. Necessary if visualizing in external software.


    * **chains** – List of chains from the entry to include in the ENM model


    * **chains_clust** – List of chains that will be assigned to quasi-rigid clusters. Must be a subset of ‘chains’



### pyCapsid.PDB.getProdyChainStarts(calphas_asym, n_units=60)

* **Parameters**

    
    * **calphas_asym** – 


    * **n_units** – 



* **Returns**

    


### pyCapsid.PDB.loadPDB(filename, pdb, save)
Loads PDBx data from a file


* **Parameters**

    
    * **filename** – Name of local file


    * **pdb** – PDB id of entry


    * **save** – Whether to save a copy of the complete assembly as pdb/pdbx



### pyCapsid.PDB.loadPDBx(filename, pdb, save)
Loads PDBx data from a file


* **Parameters**

    
    * **filename** – Name of local file


    * **pdb** – PDB id of entry


    * **save** – Whether to save a copy of the complete assembly as pdb/pdbx.


## pyCapsid.ENM module

Module with functions for building hessian matrices of different elastic network models.


### pyCapsid.ENM.buildBackbone(bblen, bbgamma, kirch, chain_starts)
Modifies the kirchhoff matrix to include backbone terms


* **Parameters**

    
    * **bblen** – The number of backbone neighbors in either direction with strengthened spring constants.


    * **bbgamma** – The relative strength to use for the backbone interactions.


    * **kirch** – The sparse kirchhoff matrix


    * **chain_starts** – Indices where protein chains begin/end.



### pyCapsid.ENM.buildENM(coords, cutoff=10, gnm=False, fanm=1, wfunc='power', base_dist=1, d_power=0, backbone=False, k_backbone=1, l_backbone=1, chain_starts=None)
Builds a hessian matrix representing an ENM based on the provided parameters.


* **Parameters**

    
    * **coords** – Cartesian of alpha carbon (or choice of representation) atoms


    * **cutoff** – Cutoff distance for long range interactions in the ENM


    * **gnm** – Whether to use only the kirchhoff matrix (Isotropic GNM), otherwise use full hessian


    * **fanm** – Parameter representing degree of anisotropy for U-ENM


    * **wfunc** – Weight function to assign spring constant based on distance. Use ‘power’ or ‘exp’


    * **base_dist** – In wfunc, divide distance by the base_distance


    * **d_power** – In wfunc, use this power of the distance


    * **backbone** – Whether to use stronger interactions for residues connected along the backbone


    * **k_backbone** – Relative strength of backbone interaction


    * **l_backbone** – How many steps along the backbone to give stronger interactions


    * **chain_starts** – Used for defining backbone interactions



* **Returns**

    A tuple of sparse matrices. The kirchoff matrix and the hessian matrix



* **Return type**

    (scipy.sparse.csr_matrix, scipy.sparse.csr_matrix)



### pyCapsid.ENM.buildENMPreset(coords, preset='ANM', \*\*kwargs)
Builds a hessian matrix representing an ENM based on one of several presets.


* **Parameters**

    
    * **coords** (*ndarray*) – Cartesian coordinates of alpha carbon atoms(or choice of alternate representation).


    * **preset** (*str*) – The specific model preset to use. Only accepts the following values:
    - ‘ANM’: Anisotropic Network Model with a cutoff of 15Å and no distance weighting.
    - ‘GNM’: Gaussian Network Model with a cutoff of 7.5Å and no distance weighting.
    - ‘U-ENM’: Unified Elastic Network Model with a cutoff of 7.5Å and and f_anm parameter of 0.1. [ref]
    - ‘bbENM’: Backbone-enhanced Elastic Network Model with a cutoff of 7.5Å and no distance weighting.
    - ‘betaENM’: Not yet implemented



* **Returns**

    A tuple of sparse matrices. The kirchoff matrix and the hessian matrix



* **Return type**

    (scipy.sparse.csr_matrix, scipy.sparse.csr_matrix)



### pyCapsid.ENM.hessCalc(row, col, kGamma, coords)

* **Parameters**

    
    * **row** – 


    * **col** – 


    * **kGamma** – 


    * **coords** – 



* **Returns**

    


### pyCapsid.ENM.kirchGamma(dists, \*\*kwargs)
Calculates the kirchhoff matrix (spring constant matrix) from a sparse distance matrix.


* **Parameters**

    
    * **dists** – Sparse distance matrix


    * **kwargs** – 



* **Returns**

    The sparse kirchhoff matrix.


## pyCapsid.NMA module

Module with functions for calculating the normal modes and frequencies of a given hessian. Eigenvalue/vector functions.
Also contains functions for calculating mechanical properties from NMA results. I.e. squared fluctuations, distance
fluctuations, compressibility, collectivity etc.


### pyCapsid.NMA.calcCovMat(evals, evecs, n_modes, coords, fluct_cutoff, is3d=True)
Calculates a sparse covariance matrix from the low frequency modes.


* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **n_modes** – 


    * **coords** – 


    * **fluct_cutoff** – 


    * **is3d** – 



* **Returns**

    


### pyCapsid.NMA.calcDistFlucts(evals, evecs, coords, fluct_cutoff=7.5, is3d=True)
Calculates a sparse covariance matrix from the low frequency modes.


* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **n_modes** – 


    * **coords** – 


    * **fluct_cutoff** – 


    * **is3d** – 



* **Returns**

    


### pyCapsid.NMA.con_c(evals, evecs, row, col)

* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **row** – 


    * **col** – 



* **Returns**

    


### pyCapsid.NMA.cov(evals, evecs, i, j)

* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **i** – 


    * **j** – 



* **Returns**

    


### pyCapsid.NMA.distFluctFromCov(c_diag, c_data, row, col)

* **Parameters**

    
    * **c_diag** – 


    * **c_data** – 


    * **row** – 


    * **col** – 



* **Returns**

    


### pyCapsid.NMA.fitCompareBfactors(evals, evecs, bfactors, pdb, is3d=True, fitModes=True, plotModes=False, forceIco=True, icotol=0.002)

* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **bfactors** – 


    * **pdb** – 


    * **is3d** – 


    * **fitModes** – 


    * **plotModes** – 


    * **forceIco** – 


    * **icotol** – 



* **Returns**

    


### pyCapsid.NMA.fluctPlot(d, title, pdb)

* **Parameters**

    
    * **d** – 


    * **title** – 


    * **pdb** – 



### pyCapsid.NMA.gCon_c(evals, evecs, row, col)

* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **row** – 


    * **col** – 



* **Returns**

    


### pyCapsid.NMA.gCov(evals, evecs, i, j)

* **Parameters**

    
    * **evals** – 


    * **evecs** – 


    * **i** – 


    * **j** – 



* **Returns**

    


### pyCapsid.NMA.modeCalc(hess, n_modes=200, eigmethod='eigsh', is3d=True)
Calculate the ‘n_modes’ lowest frequency modes of the system by calculating the smallest eigenvalues and eigenvectors
of the hessian matrix.


* **Parameters**

    
    * **hess** – Sparse hessian matrix


    * **n_modes** – Integer number of low-frequency modes to calculate.


    * **eigmethod** – Choice of method for solving the eigenvalue problem.



* **Returns**

    

## pyCapsid.QRC module

Module with functions for identifying quasi-rigid clusters given a distance fluctuation matrix, either calculated using
ENM or provided.


### pyCapsid.QRC.calcEmbedding(sims, n_vecs)

* **Parameters**

    
    * **sims** – 


    * **n_vecs** – 



* **Returns**

    


### pyCapsid.QRC.cluster_embedding(n_range, maps, method='discretize', score_method='median')

* **Parameters**

    
    * **n_range** – 


    * **maps** – 


    * **method** – 


    * **score_method** – 



* **Returns**

    


### pyCapsid.QRC.discretize(vectors, \*, copy=False, max_svd_restarts=300, n_iter_max=2000, random_state=None)
Search for a partition matrix which is closest to the eigenvector embedding.
This implementation was proposed in .


* **Parameters**

    
    * **vectors** – array-like of shape (n_samples, n_clusters)
    The embedding space of the samples.


    * **copy** – bool, default=True
    Whether to copy vectors, or perform in-place normalization.


    * **max_svd_restarts** – int, default=30
    Maximum number of attempts to restart SVD if convergence fails


    * **n_iter_max** – int, default=30
    Maximum number of iterations to attempt in rotation and partition
    matrix search if machine precision convergence is not reached


    * **random_state** – int, RandomState instance, default=None
    Determines random number generation for rotation matrix initialization.
    Use an int to make the randomness deterministic.
    See Glossary.


### Returns

labels: array of integers, shape: n_samples

    The labels of the clusters.

### References

### Notes

The eigenvector embedding is used to iteratively search for the
closest discrete partition.  First, the eigenvector embedding is
normalized to the space of partition matrices. An optimal discrete
partition matrix closest to this normalized embedding multiplied by
an initial rotation is calculated.  Fixing this discrete partition
matrix, an optimal rotation matrix is calculated.  These two
calculations are performed until convergence.  The discrete partition
matrix is returned as the clustering solution.  Used in spectral
clustering, this method tends to be faster and more robust to random
initialization than k-means.


### pyCapsid.QRC.findQuasiRigidClusters(pdb, dist_flucts, n_range, cluster_method='discretize', return_type='final', score_method='median', save=False, dir='.')
Uses spectral clustering to split the residues into clusters with minimal internal distance fluctuations.


* **Parameters**

    
    * **pdb** (*str*) – 


    * **dist_flucts** – 


    * **n_range** (*list*) – 


    * **cluster_method** (*str*) – 


    * **score_method** (*str*) – 


    * **save** (*bool*) – 


    * **dir** (*str*) – 



* **Returns**

    


### pyCapsid.QRC.fluctToSims(d)
Transforms a distance fluctuation matrix into a similarity matrix


* **Parameters**

    **d** – Sparse distance fluctuation matrix



* **Returns**

    Sparse similarity matrix
