---
layout: default
title: Tutorial
nav_order: 3
---

# Tutorial Notebook

An example notebook is provided in the [notebooks folder](https://github.com/luquelab/pyCapsid/tree/main/notebooks). 
This tutorial covers what is necessary to run the notebook.

First, create a new virtual environment using conda and then activate it.

~~~~
conda create -n pycapsid python=3.10
conda activate pycapsid
~~~~

In the virtual environment, install pycapsid using pip.
~~~~
pip install pycapsid
~~~~

To use the notebook and visualize the results, you will need jupyter lab with widgets and the nglview package. Install
these using pip.
~~~~
pip install jupyterlab ipywidgets==7.7.2 nglview
~~~~

Once these are all installed, download the notebook and run the following command in its directory to launch jupyter lab.

~~~~
jupyter lab
~~~~

The remaining steps are documented in the notebook.
