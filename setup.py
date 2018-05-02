# -*- coding: utf-8 -*-

"""
Prepare tvb Python package for setup.
"""

import shutil
from setuptools import setup, find_packages

requirements = (
    'numpy',
    'scipy',
    'scikit-learn',
    'matplotlib',
    'trimesh',
    'anytree',
    'Pegasus',
    'h5py',
    'pytest',
    'Cython',
    'gdist',
    'nibabel'
)

setup(
    name="tvb",
    description="Brain Network Models - Reconstruction tool from structural MR scans",
    packages=find_packages(),
    version="0.1",
    license="GPL",
    author="BNM Team",
    install_requires=requirements,
)

shutil.rmtree('tvb.egg-info', True)
