---
layout: default
title: Installation
nav_order: 2
---

# Installation

Via pip:
~~~~
pip install pyCapsid
~~~~

Via conda:

First, create a new virtual environment using conda and then activate it.

~~~~
conda create -n pycapsid python=3.10 -y
conda activate pycapsid
~~~~
~~~~
conda install -c luque_lab -c conda-forge pycapsid
~~~~

## Visualization in jupyter lab
To install the necessary packages for visualizing the results in a jupyter notebook, 
Via pip:
~~~~
pip install jupyterlab ipywidgets==7.7.2 nglview
~~~~

Via conda:
~~~~
conda install -c conda-forge jupyterlab ipywidgets==7.7.2 nglview
~~~~

## Visualization in ChimeraX
For higher quality visualization of the results you need to download and install [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html).
We reccomend ChimeraX 1.5
