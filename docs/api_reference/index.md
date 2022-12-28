---
layout: default
title: API Reference
nav_order: 4
---

# API Reference

Document the various modules and functions in the package.

## PDB Utility
Module with functions for downloading and dealing with PDB/PDBx files.

## Elastic Network Models
Module with functions for building hessian matrices of different elastic network models.

## Normal Mode Analysis
Module with functions for calculating the normal modes and frequencies of a given hessian. Eigenvalue/vector functions.
Also contains functions for calculating mechanical properties from NMA results. I.e. squared fluctuations, distance 
fluctuations, compressibility, collectivity etc.

## Quasi-rigid Clustering
Module with functions for identifying quasi-rigid clusters given a distance fluctuation matrix, either calculated using 
ENM or provided.

## Visualization: nglview
Functions for visualizing results in jupyter notebooks.

## Visualization: ChimeraX Bundle
The following bundle will be adapted to work for any user.
https://github.com/colintravisbrown/show-subdivision