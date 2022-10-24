"""Install script for setuptools."""

from setuptools import find_packages
from setuptools import setup

setup(
    name='pycapsid',
    version='0.0.1',
    description='A set of computational tools written in python for the analysis of viral capsids',
    author='Luque Lab, Colin Brown',
    author_email='colintravisbrown@gmail.com',
    license='Apache License, Version 2.0',
    url='https://github.com/luquelab/pycapsid',
    packages=find_packages(),
    install_requires=[ # this will be updated down the line
        'absl-py',
        'biopython',
        'chex',
        'dm-haiku',
        'dm-tree',
        'docker',
        'immutabledict',
        'jax',
        'ml-collections',
        'numpy',
        'pandas',
        'scipy',
        'tensorflow-cpu',
    ],
    tests_require=[
        'matplotlib',  # For notebook_utils_test.
        'mock',
    ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License', # incorrect license
        'Operating System :: POSIX :: Linux',  # needs to be changed
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
