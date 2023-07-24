---
layout: default
title: Installation
nav_order: 2
---

# Running remotely with Google Colab
The simplest way to use pyCapsid is using this [colab notebook](https://colab.research.google.com/github/luquelab/pyCapsid/blob/main/notebooks/pyCapsid_colab_notebook.ipynb)
which runs pyCapsid on a free Google cloud-based platform in a Jupyter environment. The Colab notebook is self-documenting 
and is designed to be simple to use. The Colab notebook has built in methods for visualizing the results for small 
structures using NGLView, but we recommend installing molecular visualization software 
[UCSF ChimeraX](#visualization-in-chimerax) locally for high-quality visualizations of larger structures.


# Requirements
Requires python 3.8-3.10. If you have an older version of python you may download a new version of python from 
[here](https://www.python.org/downloads/). Alternatively, you may use 
[conda](https://docs.conda.io/projects/conda/en/stable/), a python package manager, to create a new virtual environment
in which to install pyCapsid. You can get conda [here](https://docs.conda.io/en/latest/miniconda.html).

# Installation

## Via pip:
If you have python 3.8-3.10 installed, simply use pip, which will install pyCapsid and all of it's dependencies.
~~~~
pip install pyCapsid
~~~~
If you have a version of pyCapsid already installed, use the '--upgrade' flag to update it.
~~~~
pip install --upgrade pyCapsid
~~~~

## Via conda:
If you have an incompatible python version and don't want to upgrade we recommend installing and using [conda](https://docs.conda.io/en/latest/miniconda.html) 
to create a virtual environment with its own python version and install pyCapsid. First, create a new virtual environment 
with its own Python version using conda and then activate it.
~~~~
conda create -n pycapsid python=3.10 -y
conda activate pycapsid
~~~~
~~~~
conda install -c luque_lab -c conda-forge pycapsid
~~~~

# Visualization in ChimeraX

For the highest quality visualization of the results we recommend downloading and installing [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html).

pyCapsid has been tested on ChimeraX versions 1.5 and 1.6.

# Visualization in Jupyter Notebook with NGLView
To install the necessary packages for visualizing the results in a Jupyter notebook, 
Via pip:
~~~~
pip install jupyterlab ipywidgets==7.7.2 nglview
~~~~

Via conda:
~~~~
conda install -c conda-forge jupyterlab ipywidgets==7.7.2 nglview
~~~~


