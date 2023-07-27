"""Install script for setuptools."""

from setuptools import find_packages
from setuptools import setup
long_description = open("README.md").read()
setup(
    name='pyCapsid',
    version='0.5.3',
    description='A set of computational tools written in python for the analysis of viral capsids',
    long_description =long_description,
    long_description_content_type='text/markdown',
    author='Luque Lab, Colin Brown, Anuradha Agarwal',
    author_email='colintravisbrown@gmail.com',
    license='MIT License',
    url='https://github.com/luquelab/pycapsid',
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires = '>=3.7, <3.11',
    install_requires=[
        'biotite',
        'scikit-learn',
        'numpy',
        'scipy',
        'matplotlib',
        'pillow',
        'numba>=0.57',
        'statsmodels',
        'toml',
        'markdown'
    ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
