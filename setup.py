"""Install script for setuptools."""

from setuptools import find_packages
from setuptools import setup

setup(
    name='alphafold',
    version='2.2.4',
    description='An implementation of the inference pipeline of AlphaFold v2.0.'
    'This is a completely new model that was entered as AlphaFold2 in CASP14 '
    'and published in Nature.',
    author='DeepMind',
    author_email='alphafold@deepmind.com',
    license='Apache License, Version 2.0',
    url='https://github.com/deepmind/alphafold',
    packages=find_packages(),
    install_requires=[
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
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
)
