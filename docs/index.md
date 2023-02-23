---
layout: default
title: Home
nav_order: 1
---

# pyCapsid Documentation

[pyCapsid](https://github.com/luquelab/pyCapsid) is a high-performance Python package for analyzing protein shells and other protein complexes. Given a molecular model in [PDB format](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction), it identifies the dominant dynamics and quasi-rigid mechanical units in the structure. The figure displayed below and [Gallery page](https://luquelab.github.io/pyCapsid/gallery/) illustrate the main steps of the analysis and standard outputs. The results from pyCapsid have been generated having in mind the visualization software [NGLView](http://nglviewer.org/nglview/latest/) and [ChimeraX](https://www.cgl.ucsf.edu/chimerax/), but they can be visualized in other platforms too. The package is available on [GitHub](https://github.com/luquelab/pyCapsid), [PIP](https://pypi.org/project/pyCapsid/), and [Conda](https://anaconda.org/luque_lab/pycapsid), and the package installation steps are described in the [Installation page](https://luquelab.github.io/pyCapsid/installation/).

generate comprehensive outputs to interpret the results. The development of pyCapsid is a result of our lab's need for
high performance computational tools when investigating macromolecular protein shells. Such protein shells include viral capsids, cellular 
protein compartments like encapsulins, and capsids of gene transfer agents (GTAs), gene vectors or antivirals.

pyCapsid uses coarse-grained elastic network models in conjunction with normal mode analysis to model the motions of large 
protein complexes. These methods are commonly available in existing packages such as [ProDy](http://prody.csb.pitt.edu/)
but pyCapsid is significantly more optimized for analyzing large protein complexes. pyCapsid can perform NMA in as low as
1/5th the time for even small protein shells.

![myimg](figure_process_overview_07_13_CB.png){: width="800"; style="display:block; margin-left:auto; margin-right:auto"}



## Tutorials
A tutorial notebook is provided in the [notebooks folder](https://github.com/luquelab/pyCapsid/tree/main/notebooks).
An accompanying tutorial is provided in the documentation.

## Project History
This is an evolving repository
Started: 2022-10-24
