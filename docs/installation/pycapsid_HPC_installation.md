---
layout: default
title: Installation
nav_order: 2
---

> ## Running remotely with Google Colab
> The simplest way to use pyCapsid is using this [colab notebook](https://colab.research.google.com/github/luquelab/pyCapsid/blob/main/notebooks/pyCapsid_colab_notebook.ipynb)
which runs pyCapsid on a free Google cloud-based platform in a Jupyter environment. The Colab notebook is self-documenting 
and is designed to be simple to use. The Colab notebook has built in methods for visualizing the results for small 
structures using NGLView, but we recommend installing molecular visualization software 
[UCSF ChimeraX](#visualization-in-chimerax) locally for high-quality visualizations of larger structures.

# Installation

## Requirements
Requires python 3.8-3.10. If you have an older version of python you may download a new version of python from [here](https://www.python.org/downloads/). Alternatively, you may use [conda](https://docs.conda.io/projects/conda/en/stable/), a python package manager, to create a new virtual environment in which to install pyCapsid. You can get conda [here](https://docs.conda.io/en/latest/miniconda.html).

### Via pip:
If you have python 3.8-3.10 installed, simply use pip, which will install pyCapsid and all of it's dependencies.
```bash
pip install pyCapsid
```
If you have a version of pyCapsid already installed, use the '--upgrade' flag to update it.
```bash
pip install --upgrade pyCapsid
```

### Via conda:
If you have an incompatible python version and don't want to upgrade we recommend installing and using [conda](https://docs.conda.io/en/latest/miniconda.html) to create a virtual environment with its own python version and install pyCapsid. First, create a new virtual environment 
with its own Python version using conda and then activate it.
```bash
conda create -n pycapsid -y
conda activate pycapsid
conda install -c luque_lab -c conda-forge pycapsid
```

### GPU acceleration with CuPy (experimental)
[Cupy](https://cupy.dev/) may be used to accelerate the calculation of low-frequency modes on GPUs using CUDA. To install
cupy alongside pyCapsid, use the following command to install them together:
```
conda install -c luque_lab -c conda-forge pycapsid cupy
```
This may take some time to install. For further information on installing cupy see cupy's [installation documentation](https://docs.cupy.dev/en/stable/install.html).
The speed improvement from GPU acceleration will depend heavily on your GPU, and will be limited by the memory available 
to your GPU.

## Visualization

### Visualization in ChimeraX

For the highest quality visualization of the results we recommend downloading and installing [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html). pyCapsid has been tested on ChimeraX versions 1.5 and 1.6.

### Visualization in Jupyter Notebook with NGLView
To install the necessary packages for visualizing the results in a Jupyter notebook. Via pip:
```bash
pip install jupyterlab ipywidgets==7.7.2 nglview
```

Via conda:
```bash
conda install -c conda-forge jupyterlab ipywidgets==7.7.2 nglview
```


## High Performance Computing (HPC) Installation

The most recommended way to install PyCapsid on an HPC is by creating a conda environment to manage its dependencies.

A good practice when using servers is to avoid redundancy in software installations. Before installing conda, it is advisable to check if it is already available on the server through environment modules. This can be done with the following command once logged into the HPC:

```bash
module avail anaconda
```

Una vez que comprobamos que anaconda3 está disponible en nuestro servicio HPC, podemos cargarlo como un módulo.

```bash
module load anaconda3
```
or if you know the path to the anaconda installation:

```bash
source /share/apps/anaconda/anaconda3_build/bin/activate
```
After this, we will know that anaconda has been loaded into our session when we see the text `(base)` preceding our prompt or by evaluating the output of a command like `conda --version` or `which python`.

With conda running in our session on the HPC, we proceed to create and activate an environment for pyCapsid with the recommended version of Python (3.10):

```bash
conda create -n pycapsid_HPC python=3.10
conda activate pycapsid_HPC
```
From there the simplest way to install pyCapsid is with the following command:

```bash
conda install -c luque_lab -c conda-forge pycapsid
```
Alternatively, to ensure compatibility between dependency versions, these can be explicitly declared as follows:

```bash
conda install -c luque_lab -c conda-forge \
pycapsid=0.5.9 \
biotite=0.39.0 \
numba=0.58.1 \
numpy=1.23.5 \
scipy=1.11.4 \
matplotlib=3.8.2 \
seaborn=0.13.1 \
scikit-learn=1.3.2 \
statsmodels=0.14.1 \
markdown=3.5.2 \
pandas=2.1.4
```
### Tested environment from UM Pegasus HPC

For convenience, the functional conda environment on the HPC (x86_64 CentOS 7) Pegasus at the University of Miami has been exported to the following YAML file:

```yaml
# Conda environment that include all dependencies for PyCapsid
# with explicit versions

name: pycapsid_HPC

channels:
  - luque_lab
  - conda-forge
  - defaults

dependencies:
  - python=3.10.0
  - biotite=0.39.0
  - numpy=1.21.4
  - pandas=1.3.4
  - matplotlib=3.4.3
  - scipy=1.7.1
  - scikit-learn=0.24.2
  - requests=2.26.0
  - pytest=6.2.5 
  - pycapsid=0.5.9

```

Assuming the content of this file is contained in **pycapsid_environment.yml**, executing the following command from the same folder will replicate the environment:

```bash
conda env create -f pycapsid_environment.yml
```
To verify the installation of pyCapsid, load the environment and check that it can be imported from Python:

```bash
conda activate pycapsid_HPC
conda list | grep "pycapsid"
```