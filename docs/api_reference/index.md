---
layout: default
title: API Reference
nav_order: 4
---

# pycapsid


* [pyCapsid package]()


    * [pyCapsid.ENM module](#module-pyCapsid.ENM)


        * [`buildENM()`](#pyCapsid.ENM.buildENM)


    * [pyCapsid.NMA module](#module-pyCapsid.NMA)


        * [`modeCalc()`](#pyCapsid.NMA.modeCalc)


    * [pyCapsid.PDB module](#module-pyCapsid.PDB)


        * [`downloadPDB()`](#pyCapsid.PDB.downloadPDB)


        * [`getCapsid()`](#pyCapsid.PDB.getCapsid)


        * [`loadPDB()`](#pyCapsid.PDB.loadPDB)


        * [`loadPDBx()`](#pyCapsid.PDB.loadPDBx)


    * [pyCapsid.QRC module](#module-pyCapsid.QRC)


    * [pyCapsid.bfactorfit module](#module-pyCapsid.bfactorfit)


    * [pyCapsid.clustering_util module](#module-pyCapsid.clustering_util)


    * [Module contents](#module-pyCapsid)

# pyCapsid package

## pyCapsid.ENM module

Module with functions for building hessian matrices of different elastic network models.


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


## pyCapsid.NMA module

Module with functions for calculating the normal modes and frequencies of a given hessian. Eigenvalue/vector functions.
Also contains functions for calculating mechanical properties from NMA results. I.e. squared fluctuations, distance
fluctuations, compressibility, collectivity etc.


### pyCapsid.NMA.modeCalc(hess, kirch, n_modes, eigmethod='eigsh', gnm=False)

* **Parameters**

    
    * **hess** – 


    * **kirch** – 


    * **n_modes** – 


    * **eigmethod** – 


    * **gnm** – 



* **Returns**

    
    * x -



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


## pyCapsid.QRC module

Module with functions for identifying quasi-rigid clusters given a distance fluctuation matrix, either calculated using
ENM or provided.

## pyCapsid.bfactorfit module

## pyCapsid.clustering_util module

## Visualization: ChimeraX Bundle
The following bundle will be adapted to work for any user.
https://github.com/colintravisbrown/show-subdivision